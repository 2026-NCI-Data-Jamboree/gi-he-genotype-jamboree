"""Per-gene gated-ABMIL on pre-extracted tile features (Phikon 768-d, UNI 1024-d).
Hospital-grouped 5-fold CV on the dev partition (locked cv_fold) -> out-of-fold scores;
then one model per task trained on all dev -> scores on the locked test partition.
Fixed recipe, no early stopping, no tuning (so OOF estimates are not optimistic).
Outputs (./out): oof_preds.csv, test_preds.csv, pooled_means.npz, run_log.csv
"""
import os, json, time, numpy as np, pandas as pd, torch, torch.nn as nn, h5py
torch.backends.cudnn.benchmark = True
DEV = 'cuda' if torch.cuda.is_available() else 'cpu'
ROOT = os.environ.get('FEAT_ROOT', '/data/features'); OUT = 'out'; os.makedirs(OUT, exist_ok=True)
BACK = os.environ.get('BACKUP_DIR', ''); 
if BACK: os.makedirs(BACK, exist_ok=True)
EPOCHS = int(os.environ.get('EPOCHS', '20')); SAMPLE = int(os.environ.get('TRAIN_TILES', '1024'))
ENCS = os.environ.get('ENCS', 'phikon,uni').split(',')
FEATDIR = {'phikon': 'phikon_max4k', 'uni': 'uni_max4k'}
lab = pd.read_csv('mut_labels.csv'); TASKS = json.load(open('tasks.json'))

class ABMIL(nn.Module):
    def __init__(s, d, n_out=1, h=256, a=128, p=0.25):
        super().__init__()
        s.fc = nn.Sequential(nn.Linear(d, h), nn.ReLU(), nn.Dropout(p))
        s.att_a = nn.Sequential(nn.Linear(h, a), nn.Tanh()); s.att_b = nn.Sequential(nn.Linear(h, a), nn.Sigmoid())
        s.att_c = nn.Linear(a, 1); s.cls = nn.Linear(h, n_out)
    def forward(s, x):                       # x: (n_tiles, d)
        h = s.fc(x); w = torch.softmax(s.att_c(s.att_a(h) * s.att_b(h)).squeeze(-1), 0)
        return s.cls((w.unsqueeze(-1) * h).sum(0))

def train_predict(bags, y, train_idx, pred_idx, n_out, seed):
    torch.manual_seed(seed); rng = np.random.default_rng(seed)
    m = ABMIL(bags[0].shape[1], n_out).to(DEV)
    opt = torch.optim.AdamW(m.parameters(), lr=2e-4, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=EPOCHS * len(train_idx))
    ytr = y[train_idx]
    if n_out == 1:
        pw = torch.tensor([(ytr == 0).sum() / max(1, (ytr == 1).sum())], device=DEV)
        lossf = nn.BCEWithLogitsLoss(pos_weight=pw)
    else:
        cnt = np.bincount(ytr.astype(int), minlength=n_out); wt = len(ytr) / (n_out * np.maximum(cnt, 1))
        lossf = nn.CrossEntropyLoss(weight=torch.tensor(wt, dtype=torch.float32, device=DEV))
    m.train()
    for ep in range(EPOCHS):
        for i in rng.permutation(train_idx):
            x = bags[i]
            if len(x) > SAMPLE: x = x[torch.from_numpy(rng.choice(len(x), SAMPLE, replace=False)).to(x.device)]
            with torch.autocast(device_type='cuda', dtype=torch.float16, enabled=DEV == 'cuda'):
                out = m(x.float() if DEV == 'cpu' else x)
            if n_out == 1: loss = lossf(out.float(), torch.tensor([float(y[i])], device=DEV))
            else: loss = lossf(out.float().unsqueeze(0), torch.tensor([int(y[i])], device=DEV))
            opt.zero_grad(set_to_none=True); loss.backward(); opt.step(); sched.step()
    m.eval(); preds = []
    with torch.no_grad():
        for i in pred_idx:
            with torch.autocast(device_type='cuda', dtype=torch.float16, enabled=DEV == 'cuda'):
                o = m(bags[i].float() if DEV == 'cpu' else bags[i]).float()
            preds.append(torch.sigmoid(o).cpu().numpy() if n_out == 1 else torch.softmax(o, 0).cpu().numpy())
    return np.array(preds).reshape(len(pred_idx), -1)

def save(df, name):
    df.to_csv(f'{OUT}/{name}', index=False)
    if BACK: df.to_csv(f'{BACK}/{name}', index=False)

oof, tst, log, pooled = [], [], [], {}
T0 = time.time()
for enc in ENCS:
    slides = lab.slide_barcode.tolist(); bags = []
    t = time.time()
    for s in slides:
        with h5py.File(f'{ROOT}/{FEATDIR[enc]}/{s}.h5', 'r') as f: F = f['features'][:]
        bags.append(torch.from_numpy(F).to(DEV) if DEV == 'cuda' else torch.from_numpy(F.astype(np.float32)))
    pooled[enc] = np.stack([b.float().mean(0).cpu().numpy() for b in bags])
    np.savez_compressed(f'{OUT}/pooled_means.npz', slides=np.array(slides), **pooled)
    print(f'[{enc}] loaded {len(bags)} bags in {time.time()-t:.0f}s', flush=True)
    for group, tasks in TASKS.items():
        gi = np.where(lab.group.values == group)[0]
        jobs = [(tk, 1) for tk in tasks] + ([('crc_KRASxPIK3CA', 4)] if group == 'Colorectal' else [])
        for tk, n_out in jobs:
            y = lab[tk].values; ok = gi[~np.isnan(y[gi])]
            dev = ok[lab.partition.values[ok] == 'dev']; test = ok[lab.partition.values[ok] == 'test']
            t = time.time()
            for fold in range(5):
                tr = dev[lab.cv_fold.values[dev] != fold]; te = dev[lab.cv_fold.values[dev] == fold]
                P = train_predict(bags, y, tr, te, n_out, seed=fold)
                for k, i in enumerate(te):
                    oof.append({'enc': enc, 'task': tk, 'fold': fold, 'slide': slides[i], 'y': y[i],
                                **({'score': P[k, 0]} if n_out == 1 else {f'p{c}': P[k, c] for c in range(n_out)})})
            P = train_predict(bags, y, dev, test, n_out, seed=99)
            for k, i in enumerate(test):
                tst.append({'enc': enc, 'task': tk, 'slide': slides[i], 'y': y[i],
                            **({'score': P[k, 0]} if n_out == 1 else {f'p{c}': P[k, c] for c in range(n_out)})})
            log.append({'enc': enc, 'task': tk, 'n_dev': len(dev), 'n_test': len(test), 'sec': round(time.time()-t, 1)})
            print(json.dumps(log[-1]), f'| elapsed {time.time()-T0:.0f}s', flush=True)
            save(pd.DataFrame(oof), 'oof_preds.csv'); save(pd.DataFrame(tst), 'test_preds.csv'); save(pd.DataFrame(log), 'run_log.csv')
    del bags; torch.cuda.empty_cache() if DEV == 'cuda' else None
print(f'DONE in {time.time()-T0:.0f}s', flush=True)
