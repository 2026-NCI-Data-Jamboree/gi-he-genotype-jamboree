"""Full-cohort field extraction for the blinded KRAS reader study.
For each slide: download the IDC series to a temp dir, open with wsidicom, build a tissue map on the ~4 um/px level,
score each 517-um window by mean red-channel optical density over tissue (high for haematoxylin-dense epithelium, low for
muscle/stroma), keep windows with tissue coverage > 0.35 and < 2% near-black unsaturated pixels (ink/tattoo/folds), then draw
N_CAND spatially separated fields at random from the top 20% by score. Each is cut at level 0 -> 1024 px JPEG (0.505 um/px);
an overview is saved and the DICOMs are deleted. Resumable: slides with fields/<slide>/meta.json are skipped.
Usage: python extract_kras_full.py <pairs.csv> [n_workers] [limit]"""
import os, sys, json, shutil, subprocess, time, xml.etree.ElementTree as ET, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd, requests
from PIL import Image
from wsidicom import WsiDicom
from skimage.color import rgb2hed
from skimage.morphology import remove_small_objects, binary_closing, disk
from scipy.ndimage import uniform_filter

BASE = "https://idc-open-data.storage.googleapis.com/"
NS = {'s': 'http://doc.s3.amazonaws.com/2006-03-01'}
TILE_UM, OUT_PX, N_CAND = 517.0, 1024, 16
OUT, TMP = 'kfull/fields', 'kfull/tmp'

def objects(uuid):
    for attempt in range(6):
        try:
            r = requests.get(BASE, params={"prefix": uuid + "/"}, timeout=120); r.raise_for_status(); break
        except Exception:
            if attempt == 5: raise
            time.sleep(3 * 2 ** attempt)
    root = ET.fromstring(r.text)
    return [(k.findtext('s:Key', namespaces=NS), int(k.findtext('s:Size', namespaces=NS))) for k in root.findall('s:Contents', NS)]

def one(job):
    slide, uuid, seed = job
    od_ = f'{OUT}/{slide}'
    if os.path.exists(f'{od_}/meta.json'):
        return json.load(open(f'{od_}/meta.json'))
    d = f'{TMP}/{slide}'; os.makedirs(d, exist_ok=True); os.makedirs(od_, exist_ok=True)
    t0 = time.time()
    try:
        for key, size in objects(uuid):
            p = os.path.join(d, os.path.basename(key))
            if not (os.path.exists(p) and os.path.getsize(p) == size):
                subprocess.run(['curl', '-s', '-f', '--retry', '3', '-o', p, BASE + key], check=True)
        t_dl = time.time() - t0
        s = WsiDicom.open(d)
        lv0 = s.levels[0]; mpp0 = lv0.mpp.width
        low = min(s.levels, key=lambda l: abs(l.mpp.width - 4.0)); sp = low.mpp.width
        img = np.asarray(s.read_region((0, 0), low.level, (low.size.width, low.size.height)).convert('RGB'))
        rgb = img.astype(np.float32) / 255.
        odc = -np.log10(np.clip(rgb, 1/255, 1)); od = odc.mean(2)
        tis = binary_closing(remove_small_objects(od > 0.15, 400), disk(3))
        mx = img.max(2).astype(int); mn = img.min(2).astype(int)
        pig = (mx < 45) & ((mx - mn) < 25)                     # near-black unsaturated: ink, tattoo, folds
        side = max(8, int(round(TILE_UM / sp))); half = side // 2
        tf = uniform_filter(tis.astype(np.float32), side); pf = uniform_filter(pig.astype(np.float32), side)
        # score = mean red-channel optical density over tissue: high for haematoxylin-dense epithelium, low for muscle/stroma
        score = uniform_filter(np.where(tis, odc[..., 0], 0).astype(np.float32), side) / np.maximum(tf, 1e-3)
        H, W = tis.shape; step = max(1, side // 4)
        ys, xs = np.mgrid[half:H - half:step, half:W - half:step]; ys = ys.ravel(); xs = xs.ravel()
        ok = (tf[ys, xs] > 0.35) & (pf[ys, xs] < 0.02); ys, xs = ys[ok], xs[ok]; v = score[ys, xs]
        if len(v) == 0: raise RuntimeError('no tissue windows')
        k = min(len(v), max(4 * N_CAND, int(0.20 * len(v))))
        top = np.argsort(-v)[:k]
        rng = np.random.default_rng(seed); order = rng.permutation(top)
        for sep in (1.6, 1.2, 0.9, 0.7, 0.5):
            chosen = []
            for i in order:
                y, x = int(ys[i]), int(xs[i])
                if all((y - yy)**2 + (x - xx)**2 > (sep * side)**2 for _, yy, xx in chosen): chosen.append((float(v[i]), y, x))
                if len(chosen) >= N_CAND: break
            if len(chosen) >= N_CAND: break
        cand = v
        W0, H0 = lv0.size.width, lv0.size.height
        side0 = int(round(TILE_UM / mpp0)); scale = sp / mpp0
        fields = []
        for k, (sc, y, x) in enumerate(chosen):
            cx, cy = int(x * scale), int(y * scale)
            x0 = min(max(0, cx - side0 // 2), W0 - side0); y0 = min(max(0, cy - side0 // 2), H0 - side0)
            reg = s.read_region((x0, y0), 0, (side0, side0)).convert('RGB').resize((OUT_PX, OUT_PX), Image.LANCZOS)
            fn = f'{od_}/f{k:02d}.jpg'; reg.save(fn, quality=90)
            fields.append(dict(field=k, x0=x0, y0=y0, side0=side0, hema=round(sc, 4)))
        f = min(1.0, 1600 / max(img.shape[:2]))
        Image.fromarray(img).resize((max(1, int(img.shape[1]*f)), max(1, int(img.shape[0]*f))), Image.LANCZOS).save(f'{od_}/overview.jpg', quality=88)
        s.close()
        meta = dict(slide=slide, uuid=uuid, ok=True, n_fields=len(fields), sep=sep, mpp0=mpp0, low_mpp=sp,
                    tissue_mm2=round(float(tis.sum() * sp**2 / 1e6), 1), n_cand=len(cand),
                    t_download=round(t_dl, 1), t_total=round(time.time() - t0, 1), fields=fields)
    except Exception as e:
        meta = dict(slide=slide, uuid=uuid, ok=False, error=f'{type(e).__name__}: {e}', t_total=round(time.time() - t0, 1))
    finally:
        shutil.rmtree(d, ignore_errors=True)
    if meta.get('ok'): json.dump(meta, open(f'{od_}/meta.json', 'w'))
    return meta

if __name__ == '__main__':
    pairs = pd.read_csv(sys.argv[1]); nw = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    limit = int(sys.argv[3]) if len(sys.argv) > 3 else None
    jobs = []
    for r in pairs.itertuples():
        jobs.append((f'pair{r.pair:03d}_mut', r.mut_crdc_series_uuid, 1000 + r.pair))
        jobs.append((f'pair{r.pair:03d}_wt', r.wt_crdc_series_uuid, 5000 + r.pair))
    # interleave big and small slides to keep disk/memory bounded
    if limit: jobs = jobs[:limit]
    os.makedirs(OUT, exist_ok=True); os.makedirs(TMP, exist_ok=True)
    from multiprocessing import Pool
    done = 0; t0 = time.time()
    with Pool(nw, maxtasksperchild=4) as p, open('kfull/extract_log.jsonl', 'a') as log:
        for m in p.imap_unordered(one, jobs):
            done += 1
            log.write(json.dumps({k: v for k, v in m.items() if k != 'fields'}) + '\n'); log.flush()
            if done % 10 == 0 or not m.get('ok'):
                print(f"{done}/{len(jobs)} {time.time()-t0:.0f}s | {m['slide']} {'OK' if m.get('ok') else m.get('error')}", flush=True)
    print('DONE', done, 'in', round(time.time() - t0), 's', flush=True)
