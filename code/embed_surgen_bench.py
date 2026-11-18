"""Benchmark: embed SurGen CZI slides with Phikon + UNI, same tiling rules as embed_stream2.py
(0.5 um/px, 224 px tiles, HSV tissue mask, >=50% tissue, random cap 4,000 tiles, seed 0).
Usage: python embed_surgen_bench.py bench_slides.csv   (csv: slide,rel_path,stage,dl,size_gb)"""
import os, sys, time, json, csv, shutil, subprocess
import numpy as np, cv2, torch, h5py, timm
from torch.utils.data import Dataset, DataLoader
from transformers import ViTModel
from pylibCZIrw import czi as pyczi

BASE = 'https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/285/S-BIAD1285/Files/'
OUT = 'out'; os.makedirs(OUT, exist_ok=True)
STAGE = {'ram': '/work/slides', 'vol': '/data/tmp_surgen'}
FEAT = {'phikon': '/data/features_surgen/phikon_max4k', 'uni': '/data/features_surgen/uni_max4k'}
for p in (*STAGE.values(), *FEAT.values()): os.makedirs(p, exist_ok=True)
NCPU, BS, MAXT, TARGET_MPP, TILE, OVER = 4, 64, 4000, 0.5, 224, 2
MEAN = np.array([0.485, 0.456, 0.406], np.float32); STD = np.array([0.229, 0.224, 0.225], np.float32)
timing = []; meta = {'gpu': torch.cuda.get_device_name(0), 'aria2c': bool(shutil.which('aria2c'))}
def flush(): json.dump({'meta': meta, 'slides': timing}, open(f'{OUT}/timing.json', 'w'), indent=1)

t = time.time()
phikon = ViTModel.from_pretrained('/data/models/phikon', add_pooling_layer=False).to('cuda').half().eval()
uni = timm.create_model('vit_large_patch16_224', img_size=224, patch_size=16, init_values=1e-5, num_classes=0, dynamic_img_size=True)
uni.load_state_dict(torch.load('/data/models/uni/pytorch_model.bin', map_location='cpu'), strict=True)
uni = uni.to('cuda').half().eval()
ENC = {'phikon': lambda x: phikon(pixel_values=x).last_hidden_state[:, 0], 'uni': lambda x: uni(x)}
meta['model_load_s'] = round(time.time() - t, 1); print('models ready', meta, flush=True); flush()

def download(url, dst, how):
    t = time.time()
    if how == 'aria' and shutil.which('aria2c'):
        cmd = ['aria2c', '-x', '8', '-s', '8', '-k', '16M', '-q', '--file-allocation=none', '-d', os.path.dirname(dst), '-o', os.path.basename(dst), url]
    else:
        cmd = ['curl', '-sS', '-L', '--fail', '-o', dst, url]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0: raise RuntimeError(f'download failed: {p.stderr[-300:]}')
    return time.time() - t, os.path.getsize(dst)

def geom(cz):
    bb = cz.total_bounding_rectangle
    dist = cz.metadata['ImageDocument']['Metadata']['Scaling']['Items']['Distance']
    mpp = [float(d['Value']) * 1e6 for d in dist if d.get('@Id') == 'X'][0]
    return (int(bb[0]), int(bb[1]), int(bb[2]), int(bb[3])), mpp

def tissue_grid(cz):
    (bx, by, bw, bh), mpp = geom(cz)
    zt = mpp / 8.0                                  # thumbnail at ~8 um/px
    thumb = np.ascontiguousarray(cz.read(roi=(bx, by, bw, bh), zoom=zt, plane={'C': 0})[..., ::-1])   # BGR -> RGB
    hsv = cv2.cvtColor(thumb, cv2.COLOR_RGB2HSV)
    mask = ((hsv[..., 1] > 20) & (hsv[..., 2] > 30) & (hsv[..., 2] < 250)).astype(np.uint8)
    mask = cv2.morphologyEx(cv2.morphologyEx(mask, cv2.MORPH_OPEN, np.ones((3, 3), np.uint8)), cv2.MORPH_CLOSE, np.ones((5, 5), np.uint8))
    rp = int(round(TILE * TARGET_MPP / mpp)); sc = thumb.shape[1] / bw; s = max(1, int(rp * sc))
    coords = [(x, y) for y in range(by, by + bh - rp, rp) for x in range(bx, bx + bw - rp, rp)
              if mask[int((y - by) * sc):int((y - by) * sc) + s, int((x - bx) * sc):int((x - bx) * sc) + s].mean() >= 0.5]
    return np.array(coords, np.int64).reshape(-1, 2), rp, mpp, thumb, (bw, bh)

class Tiles(Dataset):
    def __init__(s, path, coords, rp): s.path, s.coords, s.rp, s.cz = path, coords, rp, None
    def __len__(s): return len(s.coords)
    def raw(s, i):
        if s.cz is None: s.cz = pyczi.CziReader(s.path)
        x, y = s.coords[i]
        im = s.cz.read(roi=(int(x), int(y), s.rp, s.rp), zoom=OVER * TILE / s.rp, plane={'C': 0})[..., ::-1]   # ~448 px, BGR -> RGB
        return cv2.resize(np.ascontiguousarray(im), (TILE, TILE), interpolation=cv2.INTER_AREA)
    def __getitem__(s, i):
        im = s.raw(i).astype(np.float32) / 255.
        return torch.from_numpy(((im - MEAN) / STD).transpose(2, 0, 1))

pooled = {}
for k, r in enumerate(csv.DictReader(open(sys.argv[1]))):
    name = r['slide']; dst = f"{STAGE[r['stage']]}/{name}.czi"
    rec = {'slide': name, 'stage': r['stage'], 'dl': r['dl']}; t0 = time.time()
    try:
        dl_s, nbytes = download(BASE + r['rel_path'], dst, r['dl'])
        rec.update(dl_s=round(dl_s, 1), size_gb=round(nbytes / 1e9, 2), dl_MBps=round(nbytes / 1e6 / dl_s, 1))
        t1 = time.time(); cz = pyczi.CziReader(dst)
        coords, rp, mpp, thumb, wh = tissue_grid(cz); cz.close(); n_all = len(coords)
        if n_all > MAXT: coords = coords[np.sort(np.random.default_rng(0).choice(n_all, MAXT, replace=False))]
        rec.update(grid_s=round(time.time() - t1, 1), n_tissue_tiles=int(n_all), n_tiles=int(len(coords)), read_px=rp, base_mpp=round(mpp, 4),
                   level0_wh=list(wh), thumb_rgb_mean=[round(float(v), 1) for v in thumb.reshape(-1, 3)[thumb.reshape(-1, 3).sum(1) > 0].mean(0)])
        if len(coords) == 0: raise RuntimeError('no tissue tiles')
        ds = Tiles(dst, coords, rp)
        if k == 0:   # visual checks: thumbnail + 12 random tiles
            f = 900 / max(thumb.shape[:2]); cv2.imwrite(f'{OUT}/check_thumb.jpg', cv2.resize(thumb, None, fx=f, fy=f)[..., ::-1], [cv2.IMWRITE_JPEG_QUALITY, 80])
            ii = np.random.default_rng(1).choice(len(coords), 12, replace=False)
            mont = np.vstack([np.hstack([ds.raw(int(i)) for i in ii[a:a + 4]]) for a in (0, 4, 8)])
            cv2.imwrite(f'{OUT}/check_tiles.jpg', mont[..., ::-1], [cv2.IMWRITE_JPEG_QUALITY, 85]); ds.cz.close(); ds.cz = None
        t2 = time.time(); feats = {e: [] for e in ENC}
        with torch.no_grad():
            for xb in DataLoader(ds, batch_size=BS, num_workers=NCPU, pin_memory=True, prefetch_factor=4):
                xb = xb.to('cuda', non_blocking=True).half()
                for e, fn in ENC.items(): feats[e].append(fn(xb).float().cpu().numpy())
        rec['embed_s'] = round(time.time() - t2, 1)
        for e in ENC:
            F = np.concatenate(feats[e]); pooled[f'{name}__{e}'] = F.mean(0).astype(np.float32)
            rec[f'{e}_nan'] = int(np.isnan(F).sum()); rec[f'{e}_std'] = round(float(F.std(0).mean()), 4)
            with h5py.File(f'{FEAT[e]}/{name}.h5.tmp', 'w') as h:
                h['features'] = F.astype(np.float16); h['coords'] = coords
                h.attrs.update({'read_px_level0': rp, 'mpp': TARGET_MPP, 'n_tissue_tiles': n_all, 'max_tiles': MAXT, 'source': 'SurGen CZI', 'oversample': OVER})
            os.replace(f'{FEAT[e]}/{name}.h5.tmp', f'{FEAT[e]}/{name}.h5')
    except Exception as ex:
        rec['error'] = f'{type(ex).__name__}: {str(ex)[:300]}'
    finally:
        if os.path.exists(dst): os.remove(dst)
    rec['total_s'] = round(time.time() - t0, 1); timing.append(rec); print(json.dumps(rec), flush=True); flush()
    np.savez_compressed(f'{OUT}/surgen_pool_bench.npz', **pooled)
    if k == 0:   # after the first slide: try reading a slide over HTTP without downloading it
        try:
            p = subprocess.run([sys.executable, 'stream_probe.py', BASE + open('probe_slide.txt').read().strip()], capture_output=True, text=True, timeout=300)
            meta['stream_probe'] = (p.stdout.strip().splitlines() or ['no output'])[-1][:600] + (' | ERR ' + p.stderr[-300:] if p.returncode else '')
        except subprocess.TimeoutExpired: meta['stream_probe'] = 'timeout after 300 s'
        except Exception as ex: meta['stream_probe'] = f'{type(ex).__name__}: {ex}'
        print('stream probe:', meta['stream_probe'], flush=True); flush()
print('DONE', flush=True)
