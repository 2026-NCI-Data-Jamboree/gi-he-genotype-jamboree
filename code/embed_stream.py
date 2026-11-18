"""Stream slides from IDC (s3) to local disk, tile at 0.5 um/px, embed with Phikon CLS, keep <=MAX_TILES random tissue tiles.
Usage: python embed_stream.py slides.csv   (csv: slide_barcode,s3_url). Resumable: skips slides whose .h5 exists."""
import os, sys, time, json, glob, csv, shutil, subprocess, threading, queue
import numpy as np, cv2, torch, h5py
from torch.utils.data import Dataset, DataLoader
from wsidicom import WsiDicom
from transformers import ViTModel
OUT = 'out'; FEAT = os.environ.get('FEAT', '/data/features/phikon_max4k'); TMP = '/tmp/slides'
for p in (OUT, FEAT, TMP): os.makedirs(p, exist_ok=True)
NCPU = int(os.environ.get('NCPU', '4')); BS = 128; MAXT = int(os.environ.get('MAX_TILES', '4000')); TARGET_MPP, TILE = 0.5, 224
MEAN = np.array([0.485, 0.456, 0.406], np.float32); STD = np.array([0.229, 0.224, 0.225], np.float32)
rows = [r for r in csv.DictReader(open(sys.argv[1])) if not os.path.exists(f"{FEAT}/{r['slide_barcode']}.h5")]
print(f'{len(rows)} slides to do', flush=True)
model = ViTModel.from_pretrained('/data/models/phikon', add_pooling_layer=False).to('cuda').half().eval()

def fetch(r):
    b = r['slide_barcode']; vol = f'/data/slides/{b}'
    if glob.glob(vol + '/*.dcm'): return vol, 0.0, False
    dst = f'{TMP}/{b}'; os.makedirs(dst, exist_ok=True); t = time.time()
    for attempt in range(3):
        p = subprocess.run(['s5cmd', '--no-sign-request', '--endpoint-url', 'https://s3.amazonaws.com', 'cp', '--show-progress=false',
                            r['s3_url'].rstrip('/') + '/*', dst + '/'], capture_output=True, text=True)
        if p.returncode == 0 and glob.glob(dst + '/*.dcm'): break
        time.sleep(5)
    return dst, time.time() - t, True

q = queue.Queue(maxsize=3)
from concurrent.futures import ThreadPoolExecutor
from collections import deque
def _safe(r):
    try: return (r,) + fetch(r)
    except Exception as e: return (r, None, 0.0, f'{type(e).__name__}: {e}')
def producer():                      # 2 downloads in flight, results queued in order, queue bounded at 3
    ex = ThreadPoolExecutor(int(os.environ.get('DL_THREADS', '2'))); inflight = deque()
    for r in rows:
        inflight.append(ex.submit(_safe, r))
        if len(inflight) >= int(os.environ.get('DL_THREADS', '2')): q.put(inflight.popleft().result())
    while inflight: q.put(inflight.popleft().result())
    q.put(None)
threading.Thread(target=producer, daemon=True).start()

def specimen_desc(path):
    try:
        import pydicom
        ds = pydicom.dcmread(sorted(glob.glob(path + '/*.dcm'))[0], stop_before_pixels=True)
        out = []
        for s in getattr(ds, 'SpecimenDescriptionSequence', []):
            for k in ('SpecimenShortDescription', 'SpecimenDetailedDescription', 'SpecimenIdentifier'):
                if hasattr(s, k): out.append(f'{k}={getattr(s, k)}')
        for k in ('SeriesDescription', 'ContainerIdentifier'):
            if hasattr(ds, k): out.append(f'{k}={getattr(ds, k)}')
        return ' | '.join(map(str, out))[:400]
    except Exception as e: return f'ERR {e}'

def tissue_grid(w):
    base_mpp = w.levels[0].mpp.width; lo = w.levels[-1]
    thumb = np.asarray(w.read_region((0, 0), lo.level, (lo.size.width, lo.size.height)).convert('RGB'))
    hsv = cv2.cvtColor(thumb, cv2.COLOR_RGB2HSV)
    mask = ((hsv[..., 1] > 20) & (hsv[..., 2] > 30) & (hsv[..., 2] < 250)).astype(np.uint8)
    mask = cv2.morphologyEx(cv2.morphologyEx(mask, cv2.MORPH_OPEN, np.ones((3, 3), np.uint8)), cv2.MORPH_CLOSE, np.ones((5, 5), np.uint8))
    read_px = int(round(TILE * TARGET_MPP / base_mpp)); W0, H0 = w.levels[0].size.width, w.levels[0].size.height; sc = base_mpp / lo.mpp.width
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

timing = []; T0 = time.time(); n_done = 0
while True:
    item = q.get()
    if item is None: break
    r, path, dl_s, tmp = item; name = r['slide_barcode']; rec = {'slide': name, 'dl_s': round(dl_s, 1)}; t0 = time.time()
    try:
        if path is None: raise RuntimeError(tmp)
        rec['specimen'] = specimen_desc(path) if name.startswith('CPTAC') else ''
        w = WsiDicom.open(path); coords, rp, mpp = tissue_grid(w); w.close(); n_all = len(coords)
        if n_all > MAXT: coords = coords[np.sort(np.random.default_rng(0).choice(n_all, MAXT, replace=False))]
        rec.update(n_tissue_tiles=int(n_all), n_tiles=int(len(coords)), read_px=rp, base_mpp=round(mpp, 4))
        if len(coords) == 0: raise RuntimeError('no tissue tiles')
        feats = []
        with torch.no_grad():
            for xb in DataLoader(Tiles(path, coords, rp), batch_size=BS, num_workers=NCPU, pin_memory=True, prefetch_factor=4):
                feats.append(model(pixel_values=xb.to('cuda', non_blocking=True).half()).last_hidden_state[:, 0].float().cpu().numpy())
        F = np.concatenate(feats).astype(np.float16)
        rec.update(feat_std=round(float(F.astype(np.float32).std(0).mean()), 4), feat_nan=int(np.isnan(F.astype(np.float32)).sum()))
        with h5py.File(f'{FEAT}/{name}.h5.tmp', 'w') as f:
            f['features'] = F; f['coords'] = coords
            f.attrs.update({'read_px_level0': rp, 'mpp': TARGET_MPP, 'model': 'owkin/phikon CLS', 'n_tissue_tiles': n_all, 'max_tiles': MAXT, 'specimen': rec['specimen']})
        os.replace(f'{FEAT}/{name}.h5.tmp', f'{FEAT}/{name}.h5'); n_done += 1
    except Exception as e:
        rec['error'] = f'{type(e).__name__}: {str(e)[:300]}'
    finally:
        if tmp is True: shutil.rmtree(path, ignore_errors=True)
    rec['proc_s'] = round(time.time() - t0, 1); rec['t_s'] = round(time.time() - T0, 1); timing.append(rec)
    print(json.dumps(rec), flush=True)
    if len(timing) % 10 == 0 or q.empty(): json.dump(timing, open(f'{OUT}/timing.json', 'w'))
json.dump(timing, open(f'{OUT}/timing.json', 'w'))
print(f'DONE {n_done}/{len(rows)} in {time.time() - T0:.0f}s', flush=True)
