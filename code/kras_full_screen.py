"""Per-slide tumour screen + blinded packaging for the full KRAS reader study.
Run inside the python kernel (needs `host`): exec(open('kras_full_screen.py').read())
screen_all() -> kfull/field_screen.csv ; package_all() -> kfull/pkg/pairNNN.zip + kfull/KEY_full.json"""
import os, glob, json, zipfile, secrets
import numpy as np, pandas as pd
from PIL import Image, ImageDraw

CELL = 256
SCREEN_TOOL = {"name": "label_fields", "input_schema": {"type": "object", "required": ["fields"], "properties": {"fields": {"type": "array", "items": {
    "type": "object", "required": ["idx", "label", "carcinoma_fraction"], "properties": {
        "idx": {"type": "integer"},
        "label": {"type": "string", "enum": ["carcinoma", "normal_mucosa", "adenoma_or_dysplasia", "stroma_or_muscle", "lymphoid",
                                             "necrosis_or_mucin_only", "artifact_or_blur", "other"]},
        "carcinoma_fraction": {"type": "number", "minimum": 0, "maximum": 1}}}}}}}
SCREEN_PROMPT = ("Grid of numbered H&E fields from ONE colorectal resection slide. Each cell is a 517 um square field shown at ~2 um/px. "
    "For every numbered field give the dominant tissue label and the fraction of the field that is invasive carcinoma epithelium. "
    "Normal colonic mucosa (regular test-tube crypts, abundant goblet cells, small basal nuclei) is normal_mucosa, not carcinoma. "
    "Carcinoma = crowded stratified enlarged nuclei, complex/fused/cribriform glands, solid nests, or glands infiltrating stroma or muscle. "
    "Use artifact_or_blur when the field is out of focus, folded or crushed so that it cannot be read. Return one entry per numbered field.")

def grid_for(slide):
    fs = sorted(glob.glob(f'kfull/fields/{slide}/f*.jpg')); n = len(fs); cols = 4; rows = -(-n // cols)
    can = Image.new('RGB', (cols * CELL, rows * CELL), 'white'); d = ImageDraw.Draw(can)
    for i, f in enumerate(fs):
        r_, c_ = divmod(i, cols); can.paste(Image.open(f).resize((CELL, CELL), Image.LANCZOS), (c_ * CELL, r_ * CELL))
        d.rectangle([c_ * CELL, r_ * CELL, c_ * CELL + 30, r_ * CELL + 22], fill='white'); d.text((c_ * CELL + 4, r_ * CELL + 4), str(i), fill='black')
    os.makedirs('kfull/screen_grids', exist_ok=True)
    p = f'kfull/screen_grids/{slide}.jpg'; can.save(p, quality=90); return p, fs

def screen_all(slides=None, out='kfull/field_screen.csv'):
    done = pd.read_csv(out) if os.path.exists(out) else pd.DataFrame(columns=['slide'])
    if slides is None:
        slides = sorted(d for d in os.listdir('kfull/fields') if os.path.exists(f'kfull/fields/{d}/meta.json'))
    todo = [s for s in slides if s not in set(done.slide)]
    reqs, meta = [], []
    for s in todo:
        p, fs = grid_for(s); meta.append((s, fs))
        reqs.append({"prompt": SCREEN_PROMPT, "images": [p], "tools": [SCREEN_TOOL], "tool_choice": {"type": "tool", "name": "label_fields"},
                     "max_tokens": 2000, "model": host.reasoning_model()})
    rows = []
    for i in range(0, len(reqs), 64):
        outs = host.llm(reqs[i:i + 64], max_concurrency=16)
        for (s, fs), o in zip(meta[i:i + 64], outs):
            if 'error' in o:
                print('screen error', s, o['error'][:120]); continue
            for f in o['tool_use']['input']['fields']:
                if 0 <= f['idx'] < len(fs): rows.append(dict(slide=s, file=fs[f['idx']], **f))
        pd.concat([done, pd.DataFrame(rows)], ignore_index=True).to_csv(out, index=False)
    return pd.read_csv(out)

def choose_fields(scr, slide, n=8, seed=0):
    g = scr[scr.slide == slide].copy()
    g['carcinoma_fraction'] = pd.to_numeric(g.carcinoma_fraction, errors='coerce').fillna(0)
    good = g[(g.label == 'carcinoma') & (g.carcinoma_fraction >= 0.25)]
    rng = np.random.default_rng(seed)
    pick = list(good.file.sample(frac=1, random_state=int(rng.integers(1e9))))
    if len(pick) < n:   # top up with the most carcinoma-rich remaining fields (not blur/artifact)
        rest = g[~g.file.isin(pick) & (g.label != 'artifact_or_blur')].sort_values('carcinoma_fraction', ascending=False)
        pick += list(rest.file)[: n - len(pick)]
    return pick[:n], int(len(good))

def package_all(pairs_csv='kras_full_pairs_UNBLINDED.csv', scr_csv='kfull/field_screen.csv', n=8):
    P = pd.read_csv(pairs_csv); scr = pd.read_csv(scr_csv); os.makedirs('kfull/pkg', exist_ok=True)
    key, rows = {}, []
    for r in P.itertuples():
        mut, wt = f'pair{r.pair:03d}_mut', f'pair{r.pair:03d}_wt'
        if not all(os.path.exists(f'kfull/fields/{s}/meta.json') for s in (mut, wt)):
            rows.append(dict(pair=r.pair, packaged=False, reason='extraction missing')); continue
        A, B = (mut, wt) if secrets.randbelow(2) == 0 else (wt, mut)
        key[int(r.pair)] = {'A': A, 'B': B}
        zp = f'kfull/pkg/pair{r.pair:03d}.zip'
        with zipfile.ZipFile(zp, 'w', zipfile.ZIP_STORED) as z:
            for slot, s in (('A', A), ('B', B)):
                fs, ngood = choose_fields(scr, s, n=n, seed=int(r.pair) * 10 + (slot == 'B'))
                for j, f in enumerate(fs): z.write(f, f'slide{slot}_field{j + 1}.jpg')
                z.write(f'kfull/fields/{s}/overview.jpg', f'slide{slot}_overview.jpg')
                rows.append(dict(pair=r.pair, slot=slot, slide=s, n_fields=len(fs), n_carcinoma_screened=ngood, packaged=True))
    json.dump(key, open('kfull/KEY_full.json', 'w'))
    pd.DataFrame(rows).to_csv('kfull/package_manifest.csv', index=False)
    return key, pd.DataFrame(rows)
