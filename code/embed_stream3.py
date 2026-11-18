"""Stream slides from IDC (s3), tile at 0.5 um/px, embed with any subset of
Phikon / UNI / UNI2-h in ONE pass on identical tiles.

ENCODERS env var selects the encoders, e.g. ENCODERS=uni2h (backfill pass over slides
that already have phikon+uni) or ENCODERS=phikon,uni,uni2h (slides with nothing yet).
Resumable: a slide is skipped only if EVERY selected encoder's .h5 already exists.

Tile-alignment guarantee: tiling is deterministic (tissue_grid + rng seed 0), so a later
single-encoder pass reproduces the same coords. The script does not assume this -- when an
earlier feature file exists it reads its coords and compares. On mismatch the slide is
recorded with coords_match=false and NOT written, so misaligned features can never be
concatenated downstream.

Usage: python embed_stream3.py slides.csv   (csv: slide_barcode,s3_url)
"""
import os, sys, time, json, glob, csv, shutil, subprocess, threading, queue
import numpy as np, cv2, torch, h5py, timm
from torch.utils.data import Dataset, DataLoader
from wsidicom import WsiDicom
from transformers import ViTModel

OUT = 'out'; TMP = os.environ.get('SLIDE_TMP', '/tmp/slides')
ALL_FEAT = {'phikon': '/data/features/phikon_max4k',
            'uni':    '/data/features/uni_max4k',
            'uni2h':  '/data/features/uni2h_max4k'}
SEL = [e.strip() for e in os.environ.get('ENCODERS', 'phikon,uni,uni2h').split(',') if e.strip()]
assert set(SEL) <= set(ALL_FEAT), f'unknown encoder in {SEL}'
FEAT = {k: ALL_FEAT[k] for k in SEL}
for p in (OUT, TMP, *ALL_FEAT.values()): os.makedirs(p, exist_ok=True)
NCPU = int(os.environ.get('NCPU', '4')); BS = int(os.environ.get('BS', '64'))
MAXT = int(os.environ.get('MAX_TILES', '4000')); TARGET_MPP, TILE = 0.5, 224
MEAN = np.array([0.485, 0.456, 0.406], np.float32); STD = np.array([0.229, 0.224, 0.225], np.float32)

rows = [r for r in csv.DictReader(open(sys.argv[1]))
        if not all(os.path.exists(f"{d}/{r['slide_barcode']}.h5") for d in FEAT.values())]
print(f'{len(rows)} slides to do | encoders {SEL} | batch {BS} | max_tiles {MAXT}', flush=True)

ENC, MODEL_TAG = {}, {}
if 'phikon' in SEL:
    _ph = ViTModel.from_pretrained('/data/models/phikon', add_pooling_layer=False).to('cuda').half().eval()
    ENC['phikon'] = lambda x: _ph(pixel_values=x).last_hidden_state[:, 0]
    MODEL_TAG['phikon'] = 'owkin/phikon CLS'
if 'uni' in SEL:
    _un = timm.create_model('vit_large_patch16_224', img_size=224, patch_size=16,
                            init_values=1e-5, num_classes=0, dynamic_img_size=True)
    _un.load_state_dict(torch.load('/data/models/uni/pytorch_model.bin', map_location='cpu'), strict=True)
    _un = _un.to('cuda').half().eval()
    ENC['uni'] = lambda x: _un(x)
    MODEL_TAG['uni'] = 'MahmoodLab/UNI vit_large_patch16_224'
if 'uni2h' in SEL:
    # kwargs verbatim from embed_bench3.py, which loaded this checkpoint strict=True on a T4
    _u2_kwargs = dict(model_name='vit_giant_patch14_224', img_size=224, patch_size=14, depth=24,
                      num_heads=24, init_values=1e-5, embed_dim=1536, mlp_ratio=2.66667 * 2,
                      num_classes=0, no_embed_class=True, mlp_layer=timm.layers.SwiGLUPacked,
                      act_layer=torch.nn.SiLU, reg_tokens=8, dynamic_img_size=True)
    _u2 = timm.create_model(pretrained=False, **_u2_kwargs)
    _u2.load_state_dict(torch.load('/data/models/uni2h/pytorch_model.bin', map_location='cpu'), strict=True)
    _u2 = _u2.to('cuda').half().eval()
    ENC['uni2h'] = lambda x: _u2(x)
    MODEL_TAG['uni2h'] = 'MahmoodLab/UNI2-h vit_giant_patch14_224 SwiGLUPacked reg8'
print('models ready:', list(ENC), flush=True)
def fetch(r):
    b = r['slide_barcode']; vol = f'/data/slides/{b}'
    if glob.glob(vol + '/*.dcm'): return vol, 0.0, False
    dst = f'{TMP}/{b}'; os.makedirs(dst, exist_ok=True); t = time.time()
    for _ in range(3):
        p = subprocess.run(['s5cmd', '--no-sign-request', '--endpoint-url', 'https://s3.amazonaws.com',
                            'cp', '--show-progress=false', r['s3_url'].rstrip('/') + '/*', dst + '/'],
                           capture_output=True, text=True)
        if p.returncode == 0 and glob.glob(dst + '/*.dcm'): break
        time.sleep(5)
    return dst, time.time() - t, True

q = queue.Queue(maxsize=int(os.environ.get('PREFETCH', '3')))
from concurrent.futures import ThreadPoolExecutor
from collections import deque
def _safe(r):
    try: return (r,) + fetch(r)
    except Exception as e: return (r, None, 0.0, f'{type(e).__name__}: {e}')
def producer():
    n = int(os.environ.get('DL_THREADS', '2')); ex = ThreadPoolExecutor(n); inflight = deque()
    for r in rows:
        inflight.append(ex.submit(_safe, r))
        if len(inflight) >= n: q.put(inflight.popleft().result())
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

timing = []; T0 = time.time(); n_done = 0
while True:
    item = q.get()
    if item is None: break
    r, path, dl_s, tmp = item; name = r['slide_barcode']
    rec = {'slide': name, 'dl_s': round(dl_s, 1)}; t0 = time.time()
    try:
        if path is None: raise RuntimeError(tmp)
        rec['specimen'] = specimen_desc(path) if name.startswith('CPTAC') else ''
        w = WsiDicom.open(path); coords, rp, mpp = tissue_grid(w); w.close(); n_all = len(coords)
        if n_all > MAXT:
            coords = coords[np.sort(np.random.default_rng(0).choice(n_all, MAXT, replace=False))]
        rec.update(n_tissue_tiles=int(n_all), n_tiles=int(len(coords)), read_px=rp, base_mpp=round(mpp, 4))
        if len(coords) == 0: raise RuntimeError('no tissue tiles')
        ref = None
        for k, d in ALL_FEAT.items():
            fp = f'{d}/{name}.h5'
            if k not in SEL and os.path.exists(fp):
                with h5py.File(fp, 'r') as f: ref = (k, f['coords'][:])
                break
        if ref is not None:
            same = ref[1].shape == coords.shape and bool((ref[1] == coords).all())
            rec['coords_ref'] = ref[0]; rec['coords_match'] = same
            if not same:
                raise RuntimeError(f'coord mismatch vs existing {ref[0]} features '
                                   f'({ref[1].shape} vs {coords.shape}) -- not writing')
        feats = {k: [] for k in ENC}
        with torch.no_grad():
            for xb in DataLoader(Tiles(path, coords, rp), batch_size=BS, num_workers=NCPU,
                                 pin_memory=True, prefetch_factor=4):
                xb = xb.to('cuda', non_blocking=True).half()
                for k, fn in ENC.items(): feats[k].append(fn(xb).float().cpu().numpy())
        for k in ENC:
            F = np.concatenate(feats[k]).astype(np.float16)
            rec[f'{k}_std'] = round(float(F.astype(np.float32).std(0).mean()), 4)
            rec[f'{k}_nan'] = int(np.isnan(F.astype(np.float32)).sum())
            rec[f'{k}_dim'] = int(F.shape[1])
            with h5py.File(f'{FEAT[k]}/{name}.h5.tmp', 'w') as f:
                f['features'] = F; f['coords'] = coords
                f.attrs.update({'read_px_level0': rp, 'mpp': TARGET_MPP, 'model': MODEL_TAG[k],
                                'n_tissue_tiles': n_all, 'max_tiles': MAXT, 'specimen': rec['specimen']})
            os.replace(f'{FEAT[k]}/{name}.h5.tmp', f'{FEAT[k]}/{name}.h5')
        n_done += 1
    except Exception as e:
        rec['error'] = f'{type(e).__name__}: {str(e)[:300]}'
    finally:
        if tmp is True: shutil.rmtree(path, ignore_errors=True)
    rec['proc_s'] = round(time.time() - t0, 1); rec['t_s'] = round(time.time() - T0, 1)
    timing.append(rec); print(json.dumps(rec), flush=True)
    if len(timing) % 5 == 0 or q.empty(): json.dump(timing, open(f'{OUT}/timing.json', 'w'))
json.dump(timing, open(f'{OUT}/timing.json', 'w'))
print(f'DONE {n_done}/{len(rows)} in {time.time() - T0:.0f}s', flush=True)
