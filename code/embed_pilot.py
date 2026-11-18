import os, sys, time, json, glob
import numpy as np, cv2, torch, h5py
from torch.utils.data import Dataset, DataLoader
from wsidicom import WsiDicom
from transformers import ViTModel

OUT = 'out'; FEAT = '/data/features/phikon'; os.makedirs(OUT, exist_ok=True); os.makedirs(FEAT, exist_ok=True)
NCPU = int(os.environ.get('NCPU', '4')); BS = 128
TARGET_MPP, TILE = 0.5, 224
MEAN = np.array([0.485, 0.456, 0.406], np.float32); STD = np.array([0.229, 0.224, 0.225], np.float32)
dev = 'cuda'
model = ViTModel.from_pretrained('/data/models/phikon', add_pooling_layer=False).to(dev).half().eval()
print('gpu', torch.cuda.get_device_name(0), flush=True)

def tissue_grid(w):
    base_mpp = w.levels[0].mpp.width * 1000 if w.levels[0].mpp.width < 0.01 else w.levels[0].mpp.width
    lo = w.levels[-1]; lo_mpp = lo.mpp.width
    thumb = np.asarray(w.read_region((0, 0), lo.level, (lo.size.width, lo.size.height)).convert('RGB'))
    hsv = cv2.cvtColor(thumb, cv2.COLOR_RGB2HSV)
    mask = ((hsv[..., 1] > 20) & (hsv[..., 2] > 30) & (hsv[..., 2] < 250)).astype(np.uint8)
    mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN, np.ones((3, 3), np.uint8))
    mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, np.ones((5, 5), np.uint8))
    read_px = int(round(TILE * TARGET_MPP / base_mpp))          # level-0 pixels per tile side
    W0, H0 = w.levels[0].size.width, w.levels[0].size.height
    scale = base_mpp / lo_mpp                                     # level-0 px -> thumbnail px
    coords = []
    for y in range(0, H0 - read_px, read_px):
        for x in range(0, W0 - read_px, read_px):
            x0, y0 = int(x * scale), int(y * scale); s = max(1, int(read_px * scale))
            if mask[y0:y0 + s, x0:x0 + s].mean() >= 0.5: coords.append((x, y))
    return np.array(coords, np.int64), read_px, base_mpp, thumb, mask, scale

class Tiles(Dataset):
    def __init__(self, path, coords, read_px): self.path, self.coords, self.read_px, self.w = path, coords, read_px, None
    def __len__(self): return len(self.coords)
    def __getitem__(self, i):
        if self.w is None: self.w = WsiDicom.open(self.path)
        x, y = self.coords[i]
        im = np.asarray(self.w.read_region((int(x), int(y)), 0, (self.read_px, self.read_px)).convert('RGB'))
        im = cv2.resize(im, (TILE, TILE), interpolation=cv2.INTER_AREA).astype(np.float32) / 255.
        return torch.from_numpy(((im - MEAN) / STD).transpose(2, 0, 1))

timing = []
slides = sorted(glob.glob('/data/slides/TCGA-*'))
for sp in slides:
    name = os.path.basename(sp); rec = {'slide': name}; t0 = time.time()
    try:
        w = WsiDicom.open(sp); rec['open_s'] = round(time.time() - t0, 1)
        t = time.time(); coords, read_px, base_mpp, thumb, mask, scale = tissue_grid(w)
        rec.update(tissue_s=round(time.time() - t, 1), n_tiles=int(len(coords)), read_px=read_px, base_mpp=round(base_mpp, 4))
        w.close()
        if name in [os.path.basename(s) for s in slides[:3]]:
            qc = thumb.copy(); s = max(1, int(read_px * scale))
            for x, y in coords[::1]: cv2.rectangle(qc, (int(x*scale), int(y*scale)), (int(x*scale)+s, int(y*scale)+s), (0, 200, 0), 1)
            h = 900; qc = cv2.resize(qc, (int(qc.shape[1]*h/qc.shape[0]), h))
            cv2.imwrite(f'{OUT}/qc_{name}.jpg', cv2.cvtColor(qc, cv2.COLOR_RGB2BGR), [cv2.IMWRITE_JPEG_QUALITY, 80])
        dl = DataLoader(Tiles(sp, coords, read_px), batch_size=BS, num_workers=NCPU, pin_memory=True, persistent_workers=False, prefetch_factor=4)
        feats = []; t = time.time(); gpu_s = 0.0
        with torch.no_grad():
            for xb in dl:
                g = time.time(); xb = xb.to(dev, non_blocking=True).half()
                feats.append(model(pixel_values=xb).last_hidden_state[:, 0].float().cpu().numpy()); torch.cuda.synchronize(); gpu_s += time.time() - g
        F = np.concatenate(feats).astype(np.float16)
        rec.update(embed_wall_s=round(time.time() - t, 1), gpu_busy_s=round(gpu_s, 1), tiles_per_s=round(len(coords) / (time.time() - t), 1),
                   feat_std=round(float(F.astype(np.float32).std(0).mean()), 4), feat_nan=int(np.isnan(F.astype(np.float32)).sum()))
        with h5py.File(f'{FEAT}/{name}.h5', 'w') as f:
            f['features'] = F; f['coords'] = coords; f.attrs['read_px_level0'] = read_px; f.attrs['mpp'] = TARGET_MPP; f.attrs['model'] = 'owkin/phikon CLS'
    except Exception as e:
        rec['error'] = f'{type(e).__name__}: {str(e)[:300]}'
    rec['total_s'] = round(time.time() - t0, 1); timing.append(rec)
    print(json.dumps(rec), flush=True)
    json.dump(timing, open(f'{OUT}/timing.json', 'w'), indent=1)
print('DONE', flush=True)
