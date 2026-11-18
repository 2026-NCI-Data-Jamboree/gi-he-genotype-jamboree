"""Embed SurGen CZI slides (already on the volume) with Phikon + UNI, same tiling rules as embed_stream2.py
(0.5 um/px, 224 px tiles, HSV tissue mask, >=50% tissue, random cap 4,000 tiles, seed 0).
Usage: python embed_surgen.py gpuK.csv   (csv: slide)"""
import os, sys, time, json, csv, shutil, subprocess
import numpy as np, cv2, torch, h5py, timm
from torch.utils.data import Dataset, DataLoader
from transformers import ViTModel
from pylibCZIrw import czi as pyczi

BASE = 'https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/285/S-BIAD1285/Files/'
OUT = 'out'; os.makedirs(OUT, exist_ok=True)
SRC = '/data/surgen_czi'; POOL = '/data/features_surgen/pooled'
FEAT = {'phikon': '/data/features_surgen/phikon_max4k', 'uni': '/data/features_surgen/uni_max4k'}
for p in (POOL, *FEAT.values()): os.makedirs(p, exist_ok=True)
NCPU, BS, MAXT, TARGET_MPP, TILE, OVER = 4, 64, 4000, 0.5, 224, 2
MEAN = np.array([0.485, 0.456, 0.406], np.float32); STD = np.array([0.229, 0.224, 0.225], np.float32)
timing = []; meta = {'gpu': torch.cuda.get_device_name(0), 'list': sys.argv[1]}
def flush(): json.dump({'meta': meta, 'slides': timing}, open(f'{OUT}/timing.json', 'w'), indent=1)

t = time.time()
phikon = ViTModel.from_pretrained('/data/models/phikon', add_pooling_layer=False).to('cuda').half().eval()
uni = timm.create_model('vit_large_patch16_224', img_size=224, patch_size=16, init_values=1e-5, num_classes=0, dynamic_img_size=True)
uni.load_state_dict(torch.load('/data/models/uni/pytorch_model.bin', map_location='cpu'), strict=True)
uni = uni.to('cuda').half().eval()
ENC = {'phikon': lambda x: phikon(pixel_values=x).last_hidden_state[:, 0], 'uni': lambda x: uni(x)}
meta['model_load_s'] = round(time.time() - t, 1); print('models ready', meta, flush=True); flush()

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


TAG = os.path.basename(sys.argv[1]).replace('.csv', ''); pooled = {}; first = True
for r in csv.DictReader(open(sys.argv[1])):
    name = r['slide']; src = f'{SRC}/{name}.czi'; pz = f'{POOL}/{name}.npz'
    rec = {'slide': name}; t0 = time.time()
    try:
        if os.path.exists(pz):
            d = np.load(pz); pooled.update({f'{name}__{e}': d[e] for e in ENC}); pooled[f'{name}__n'] = d['n']; rec['status'] = 'already embedded'
        else:
            if not os.path.exists(src): raise FileNotFoundError('slide not on volume')
            cz = pyczi.CziReader(src); coords, rp, mpp, thumb, wh = tissue_grid(cz); cz.close(); n_all = len(coords)
            if n_all > MAXT: coords = coords[np.sort(np.random.default_rng(0).choice(n_all, MAXT, replace=False))]
            rec.update(grid_s=round(time.time() - t0, 1), n_tissue_tiles=int(n_all), n_tiles=int(len(coords)), read_px=rp, base_mpp=round(mpp, 4), size_gb=round(os.path.getsize(src) / 1e9, 2))
            if len(coords) == 0: raise RuntimeError('no tissue tiles')
            ds = Tiles(src, coords, rp)
            if first:
                f = 900 / max(thumb.shape[:2]); cv2.imwrite(f'{OUT}/check_thumb_{TAG}.jpg', cv2.resize(thumb, None, fx=f, fy=f)[..., ::-1], [cv2.IMWRITE_JPEG_QUALITY, 80])
                ii = np.random.default_rng(1).choice(len(coords), 12, replace=False)
                mont = np.vstack([np.hstack([ds.raw(int(i)) for i in ii[a:a + 4]]) for a in (0, 4, 8)])
                cv2.imwrite(f'{OUT}/check_tiles_{TAG}.jpg', mont[..., ::-1], [cv2.IMWRITE_JPEG_QUALITY, 85]); ds.cz.close(); ds.cz = None; first = False
            t2 = time.time(); feats = {e: [] for e in ENC}
            with torch.no_grad():
                for xb in DataLoader(ds, batch_size=BS, num_workers=NCPU, pin_memory=True, prefetch_factor=4):
                    xb = xb.to('cuda', non_blocking=True).half()
                    for e, fn in ENC.items(): feats[e].append(fn(xb).float().cpu().numpy())
            rec['embed_s'] = round(time.time() - t2, 1); means = {}
            for e in ENC:
                F = np.concatenate(feats[e]); means[e] = F.mean(0).astype(np.float32)
                rec[f'{e}_nan'] = int(np.isnan(F).sum()); rec[f'{e}_std'] = round(float(F.std(0).mean()), 4)
                with h5py.File(f'{FEAT[e]}/{name}.h5.tmp', 'w') as h:
                    h['features'] = F.astype(np.float16); h['coords'] = coords
                    h.attrs.update({'read_px_level0': rp, 'mpp': TARGET_MPP, 'n_tissue_tiles': n_all, 'max_tiles': MAXT, 'source': 'SurGen CZI', 'oversample': OVER})
                os.replace(f'{FEAT[e]}/{name}.h5.tmp', f'{FEAT[e]}/{name}.h5')
            np.savez(pz + '.tmp.npz', n=np.array([n_all, len(coords)]), **means); os.replace(pz + '.tmp.npz', pz)
            pooled.update({f'{name}__{e}': means[e] for e in ENC}); pooled[f'{name}__n'] = np.array([n_all, len(coords)]); rec['status'] = 'embedded'
    except Exception as ex:
        rec['status'] = 'FAILED'; rec['error'] = f'{type(ex).__name__}: {str(ex)[:300]}'
    rec['total_s'] = round(time.time() - t0, 1); timing.append(rec); print(json.dumps(rec), flush=True); flush()
    np.savez_compressed(f'{OUT}/surgen_pool_{TAG}.npz', **pooled)
os.replace(f'{OUT}/timing.json', f'{OUT}/timing_{TAG}.json'); subprocess.run(['sync']); time.sleep(20)
print(f'DONE {sum(x["status"] != "FAILED" for x in timing)}/{len(timing)}', flush=True)
