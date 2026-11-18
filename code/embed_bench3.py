"""Benchmark: one tiling pass per slide, three encoders (Phikon, UNI, UNI2-h) on the same tiles.
Measures per-encoder GPU seconds and a model-free decode baseline, so the cost of each
encoder subset can be projected. Usage: python embed_bench3.py slides.csv
"""
import os, sys, time, json, glob, csv, shutil, subprocess, threading, queue
import numpy as np, cv2, torch, h5py, timm
from torch.utils.data import Dataset, DataLoader
from wsidicom import WsiDicom
from transformers import ViTModel

OUT = 'out'; TMP = '/tmp/slides'; FEAT = '/data/features/bench3'
for p in (OUT, TMP): os.makedirs(p, exist_ok=True)
NCPU = int(os.environ.get('NCPU', '4')); BS = int(os.environ.get('BS', '64'))
MAXT = int(os.environ.get('MAX_TILES', '4000')); TARGET_MPP, TILE = 0.5, 224
MEAN = np.array([0.485, 0.456, 0.406], np.float32); STD = np.array([0.229, 0.224, 0.225], np.float32)
rows = list(csv.DictReader(open(sys.argv[1])))
print(f'{len(rows)} slides | batch {BS} | max_tiles {MAXT}', flush=True)

# ---- encoders (all three take ImageNet mean/std at 224 px; verified from their preprocessor configs) ----
t0 = time.time()
phikon = ViTModel.from_pretrained('/data/models/phikon', add_pooling_layer=False).to('cuda').half().eval()
uni = timm.create_model('vit_large_patch16_224', img_size=224, patch_size=16,
                        init_values=1e-5, num_classes=0, dynamic_img_size=True)
uni.load_state_dict(torch.load('/data/models/uni/pytorch_model.bin', map_location='cpu'), strict=True)
uni = uni.to('cuda').half().eval()
uni2_kwargs = dict(model_name='vit_giant_patch14_224', img_size=224, patch_size=14, depth=24,
                   num_heads=24, init_values=1e-5, embed_dim=1536, mlp_ratio=2.66667 * 2,
                   num_classes=0, no_embed_class=True, mlp_layer=timm.layers.SwiGLUPacked,
                   act_layer=torch.nn.SiLU, reg_tokens=8, dynamic_img_size=True)
uni2 = timm.create_model(pretrained=False, **uni2_kwargs)
uni2.load_state_dict(torch.load('/data/models/uni2h/pytorch_model.bin', map_location='cpu'), strict=True)
uni2 = uni2.to('cuda').half().eval()
print(f'models loaded in {time.time()-t0:.0f}s | params M: '
      f'phikon {sum(p.numel() for p in phikon.parameters())/1e6:.0f} '
      f'uni {sum(p.numel() for p in uni.parameters())/1e6:.0f} '
      f'uni2h {sum(p.numel() for p in uni2.parameters())/1e6:.0f}', flush=True)

ENC = {'phikon': lambda x: phikon(pixel_values=x).last_hidden_state[:, 0],
       'uni':    lambda x: uni(x),
       'uni2h':  lambda x: uni2(x)}
for k in ENC: os.makedirs(f'{FEAT}/{k}', exist_ok=True)

def fetch(r):
    b = r['slide_barcode']; dst = f'{TMP}/{b}'; os.makedirs(dst, exist_ok=True); t = time.time()
    for _ in range(3):
        p = subprocess.run(['s5cmd', '--no-sign-request', '--endpoint-url', 'https://s3.amazonaws.com',
                            'cp', '--show-progress=false', r['s3_url'].rstrip('/') + '/*', dst + '/'],
                           capture_output=True, text=True)
        if p.returncode == 0 and glob.glob(dst + '/*.dcm'): break
        time.sleep(5)
    return dst, time.time() - t

q = queue.Queue(maxsize=2)
def producer():
    for r in rows:
        try: q.put((r,) + fetch(r))
        except Exception as e: q.put((r, None, f'{type(e).__name__}: {e}'))
    q.put(None)
threading.Thread(target=producer, daemon=True).start()

def tissue_grid(w):
    base_mpp = w.levels[0].mpp.width; lo = w.levels[-1]
    thumb = np.asarray(w.read_region((0, 0), lo.level, (lo.size.width, lo.size.height)).convert('RGB'))
    hsv = cv2.cvtColor(thumb, cv2.COLOR_RGB2HSV)
    mask = ((hsv[..., 1] > 20) & (hsv[..., 2] > 30) & (hsv[..., 2] < 250)).astype(np.uint8)
    mask = cv2.morphologyEx(cv2.morphologyEx(mask, cv2.MORPH_OPEN, np.ones((3, 3), np.uint8)),
                            cv2.MORPH_CLOSE, np.ones((5, 5), np.uint8))
    read_px = int(round(TILE * TARGET_MPP / base_mpp))
    W0, H0 = w.levels[0].size.width, w.levels[0].size.height; sc = base_mpp / lo.mpp.width
    s = max(1, int(read_px * sc))
    coords = [(x, y) for y in range(0, H0 - read_px, read_px) for x in range(0, W0 - read_px, read_px)
              if mask[int(y * sc):int(y * sc) + s, int(x * sc):int(x * sc) + s].mean() >= 0.5]
    return np.array(coords, np.int64).reshape(-1, 2), read_px, base_mpp

class Tiles(Dataset):
    def __init__(s, path, coords, rp): s.path, s.coords, s.rp, s.w = path, coords, rp, None
    def __len__(s): return len(s.coords)
    def __getitem__(s, i):
        if s.w is None: s.w = WsiDicom.open(s.path)
        x, y = s.coords[i]
        im = np.asarray(s.w.read_region((int(x), int(y)), 0, (s.rp, s.rp)).convert('RGB'))
        im = cv2.resize(im, (TILE, TILE), interpolation=cv2.INTER_AREA).astype(np.float32) / 255.
        return torch.from_numpy(((im - MEAN) / STD).transpose(2, 0, 1))

def loader(path, coords, rp):
    return DataLoader(Tiles(path, coords, rp), batch_size=BS, num_workers=NCPU,
                      pin_memory=True, prefetch_factor=4)

timing = []; T0 = time.time()
while True:
    item = q.get()
    if item is None: break
    if len(item) == 3: r, path, dl = item
    else: r, path, dl = item[0], None, item[2]
    name = r['slide_barcode']; rec = {'slide': name, 'project': r.get('project_id', ''),
                                      'MB': r.get('MB', ''), 'dl_s': round(dl, 1) if path else 0}
    t_slide = time.time()
    try:
        if path is None: raise RuntimeError(dl)
        w = WsiDicom.open(path); coords, rp, mpp = tissue_grid(w); w.close(); n_all = len(coords)
        if n_all > MAXT:
            coords = coords[np.sort(np.random.default_rng(0).choice(n_all, MAXT, replace=False))]
        rec.update(n_tissue_tiles=int(n_all), n_tiles=int(len(coords)), read_px=rp, base_mpp=round(mpp, 4))
        if len(coords) == 0: raise RuntimeError('no tissue tiles')

        # pass 1 — decode only, no model: the data-side floor for this slide
        t = time.time()
        for xb in loader(path, coords, rp): xb.pin_memory()
        rec['decode_only_s'] = round(time.time() - t, 1)

        # pass 2 — same tiles through all three encoders, per-encoder GPU time
        feats = {k: [] for k in ENC}; gpu_s = {k: 0.0 for k in ENC}
        t_pass = time.time(); t_wait = 0.0; t_mark = time.time()
        with torch.no_grad():
            for xb in loader(path, coords, rp):
                t_wait += time.time() - t_mark
                xb = xb.to('cuda', non_blocking=True).half()
                for k, fn in ENC.items():
                    torch.cuda.synchronize(); t = time.time()
                    feats[k].append(fn(xb).float().cpu().numpy())
                    torch.cuda.synchronize(); gpu_s[k] += time.time() - t
                t_mark = time.time()
        rec['encode_pass_s'] = round(time.time() - t_pass, 1)
        rec['data_wait_s'] = round(t_wait, 1)
        for k in ENC:
            F = np.concatenate(feats[k]).astype(np.float16)
            rec[f'{k}_gpu_s'] = round(gpu_s[k], 1); rec[f'{k}_dim'] = int(F.shape[1])
            rec[f'{k}_nan'] = int(np.isnan(F.astype(np.float32)).sum())
            rec[f'{k}_std'] = round(float(F.astype(np.float32).std(0).mean()), 4)
            with h5py.File(f'{FEAT}/{k}/{name}.h5.tmp', 'w') as f:
                f['features'] = F; f['coords'] = coords
                f.attrs.update({'read_px_level0': rp, 'mpp': TARGET_MPP, 'model': k,
                                'n_tissue_tiles': n_all, 'max_tiles': MAXT})
            os.replace(f'{FEAT}/{k}/{name}.h5.tmp', f'{FEAT}/{k}/{name}.h5')
    except Exception as e:
        rec['error'] = f'{type(e).__name__}: {str(e)[:300]}'
    finally:
        if path: shutil.rmtree(path, ignore_errors=True)
    rec['slide_s'] = round(time.time() - t_slide, 1); rec['t_s'] = round(time.time() - T0, 1)
    timing.append(rec); print(json.dumps(rec), flush=True)
    json.dump(timing, open(f'{OUT}/bench3_timing.json', 'w'))
print(f'DONE in {time.time()-T0:.0f}s', flush=True)
