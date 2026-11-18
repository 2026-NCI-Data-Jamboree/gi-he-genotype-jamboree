"""Extract high-power fields for the blinded KRAS reader test.
Per slide: download the IDC series, open with wsidicom, build a tissue map on the ~4 um/px level,
sample N_SITES spatially separated cellular fields (random among the upper 60% by haematoxylin density,
excluding pigment/ink), and cut 517 um fields at level 0 resampled to 1024 px (0.505 um/px)."""
import os, sys, json, subprocess, xml.etree.ElementTree as ET
import numpy as np, pandas as pd, requests
from PIL import Image
from wsidicom import WsiDicom
from skimage.color import rgb2hed
from skimage.morphology import remove_small_objects, binary_closing, disk
from scipy.ndimage import uniform_filter

BASE = "https://idc-open-data.storage.googleapis.com/"
NS = {'s': 'http://doc.s3.amazonaws.com/2006-03-01'}
TILE_UM, OUT_PX, N_SITES = 517.0, 1024, 18

def objects(uuid):
    r = requests.get(BASE, params={"prefix": uuid + "/"}, timeout=120); r.raise_for_status()
    root = ET.fromstring(r.text)
    return [(k.findtext('s:Key', namespaces=NS), int(k.findtext('s:Size', namespaces=NS)))
            for k in root.findall('s:Contents', NS)]

def download(uuid, d):
    os.makedirs(d, exist_ok=True)
    for key, size in objects(uuid):
        p = os.path.join(d, os.path.basename(key))
        if os.path.exists(p) and os.path.getsize(p) == size: continue
        subprocess.run(['curl', '-s', '-f', '-o', p, BASE + key], check=True)

def one(slide_id, uuid, outdir, seed):
    d = f'wsik/{slide_id}'
    download(uuid, d)
    s = WsiDicom.open(d)
    lv0 = s.levels[0]; mpp0 = lv0.mpp.width
    low = min(s.levels, key=lambda l: abs(l.mpp.width - 4.0))
    sp = low.mpp.width
    img = np.asarray(s.read_region((0, 0), low.level, (low.size.width, low.size.height)).convert('RGB'))
    rgb = img.astype(np.float32) / 255.
    od = -np.log10(np.clip(rgb, 1/255, 1).mean(2))
    tis = binary_closing(remove_small_objects(od > 0.15, 400), disk(3))
    pigment = (img.max(2) < 60)                         # ink / tattoo / dense pigment
    hema = np.clip(rgb2hed(rgb)[..., 0], 0, None)
    side = int(round(TILE_UM / sp)); half = side // 2
    tf = uniform_filter(tis.astype(np.float32), side)
    pf = uniform_filter(pigment.astype(np.float32), side)
    hd = uniform_filter(np.where(tis, hema, 0).astype(np.float32), side)
    H, W = tis.shape
    cand = [(float(hd[y, x]), y, x) for y in range(half, H - half, max(1, side // 4))
            for x in range(half, W - half, max(1, side // 4)) if tf[y, x] > 0.88 and pf[y, x] < 0.01]
    if len(cand) < N_SITES: raise RuntimeError(f'only {len(cand)} candidate windows')
    thr = np.percentile([c[0] for c in cand], 40)
    pool = [c for c in cand if c[0] >= thr]
    rng = np.random.default_rng(seed); rng.shuffle(pool)
    chosen = []
    for sep in (1.6, 1.3, 1.0, 0.8):
        chosen = []
        for c in pool:
            if all((c[1]-y)**2 + (c[2]-x)**2 > (sep*side)**2 for _, y, x in chosen): chosen.append(c)
            if len(chosen) >= N_SITES: break
        if len(chosen) >= N_SITES: break
    W0, H0 = lv0.size.width, lv0.size.height
    side0 = int(round(TILE_UM / mpp0)); scale = sp / mpp0
    os.makedirs(outdir, exist_ok=True)
    meta = []
    for k, (sc, y, x) in enumerate(chosen):
        cx, cy = int(x * scale), int(y * scale)
        x0 = min(max(0, cx - side0 // 2), W0 - side0); y0 = min(max(0, cy - side0 // 2), H0 - side0)
        reg = s.read_region((x0, y0), 0, (side0, side0)).convert('RGB').resize((OUT_PX, OUT_PX), Image.LANCZOS)
        reg.save(f'{outdir}/{slide_id}_f{k:02d}.jpg', quality=93)
        meta.append(dict(slide=slide_id, field=k, x0=x0, y0=y0, side0=side0, hema=sc))
    f = min(1.0, 1800 / max(img.shape[:2]))
    ov = Image.fromarray(img).resize((int(img.shape[1]*f), int(img.shape[0]*f)), Image.LANCZOS)
    ov.save(f'{outdir}/{slide_id}_overview.jpg', quality=90)
    s.close()
    tiles = np.stack([np.asarray(Image.open(f'{outdir}/{slide_id}_f{k:02d}.jpg')) for k in range(len(chosen))])
    return dict(slide=slide_id, n_fields=len(chosen), sep=sep, mpp0=mpp0, low_mpp=sp,
                tissue_mm2=float(tis.sum() * sp**2 / 1e6), n_candidates=len(cand),
                mean_rgb=tiles.reshape(-1, 3).mean(0).round(1).tolist(), meta=meta)

if __name__ == '__main__':
    from multiprocessing import Pool
    sel = pd.read_csv('kras_pairs_UNBLINDED.csv')
    jobs = []
    for _, r in sel.iterrows():
        jobs.append((f"P{r.pair}_mut", r.mut_uuid, 'ktiles', 100 + r.pair))
        jobs.append((f"P{r.pair}_wt", r.wt_uuid, 'ktiles', 200 + r.pair))
    out = []
    def run(j):
        try: return one(*j)
        except Exception as e: return dict(slide=j[0], error=f'{type(e).__name__}: {e}')
    with Pool(3) as p:
        for res in p.imap_unordered(run, jobs):
            out.append(res)
            print(res['slide'], res.get('error') or f"fields {res['n_fields']} sep {res['sep']} mpp0 {res['mpp0']:.3f} "
                  f"tissue {res['tissue_mm2']:.0f} mm2 rgb {res['mean_rgb']}", flush=True)
            json.dump(out, open('ktiles_summary.json', 'w'))
    print('DONE', len(out), 'errors', sum('error' in o for o in out), flush=True)
