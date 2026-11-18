#!/usr/bin/env python
"""Pre-registered KRAS slide model, DEV partition only (preregistration.md, locked 2026-11-15).

Cohort: TCGA-COAD/READ, sequenced, MSS, partition == 'dev', primary_slide == True (1 slide/patient).
Tasks:  label_kras (PRIMARY; negatives include BRAF-mutant KRAS-WT)
        label_kras_sens (SENSITIVITY; negatives also BRAF-WT)
        Both columns come from patient_split_locked and already encode the pre-registered exclusions
        (hypermutated TMB>=10, use_for_training == False, QC fails).
Folds:  the locked TSS-grouped dev folds (cv_fold 0-4). No new folds. Test partition never loaded
        (asserted). For fold k, models train on task-labelled patients of the other 4 folds and score
        EVERY scoreable dev primary slide with cv_fold == k (labelled or not), so dev slides without a
        task label (hypermutated / QC / BRAF-only for the sensitivity task) get a fold-held-out score
        flagged label = NA.
Discovery patients (is_discovery) are kept in training (encoder-only ABMIL is allowed to train on them,
        prereg Section 4); evaluation reports all-dev and dev_excl_discovery.
Repeats: 3 training seeds (11, 22, 33) on the same locked folds (seed-variance only; LR is deterministic).
Encoders: UNI v1 (1024-d, primary available; pre-registered UNI2-h replaced for cost) and Phikon (768-d).
        Only slides with BOTH uni_max4k and phikon_max4k .h5 (tile-aligned dual run) are used.

HYPERPARAMETERS (fixed a priori, not tuned):
  gated ABMIL (Ilse 2018): Linear(d,256)-ReLU-Dropout(0.25); attention V,U: 256->128 (tanh x sigmoid) -> 1;
  classifier Dropout(0.25)-Linear(256,1). AdamW lr 1e-4, weight decay 1e-4, batch = 1 slide,
  random subsample of 2048 tiles per training step (all tiles at inference), max 40 epochs,
  early stopping patience 8 on inner-validation class-weighted BCE; best-epoch weights kept.
  Class weighting: BCE pos_weight = n_neg / n_pos of the inner-training set.
  Inner split: 15% of the training-fold patients, stratified by label x project (seed + fold).
  Mean-pool LR baseline: per-slide mean of tile features -> StandardScaler -> L2 LogisticRegression,
  class_weight='balanced'; C in {1e-3,3e-3,1e-2,3e-2,1e-1,1} chosen by inner-validation AUROC, then refit
  on the full training fold.
OOF predictions are rewritten to out/oof_predictions.csv after every fold.
"""
import argparse, json, os, time, random
from concurrent.futures import ThreadPoolExecutor
import numpy as np, pandas as pd, h5py, torch, torch.nn as nn
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score

ap = argparse.ArgumentParser()
ap.add_argument('--feat_root', default='/data/features')
ap.add_argument('--slides', default='job_slides_dev_primary.csv')
ap.add_argument('--out', default='out')
ap.add_argument('--tasks', default='label_kras,label_kras_sens')
ap.add_argument('--seeds', default='11,22,33')
ap.add_argument('--max_epochs', type=int, default=40)
ap.add_argument('--patience', type=int, default=8)
ap.add_argument('--train_tiles', type=int, default=2048)
ap.add_argument('--lr', type=float, default=1e-4)
ap.add_argument('--wd', type=float, default=1e-4)
ap.add_argument('--device', default='cuda' if torch.cuda.is_available() else 'cpu')
ap.add_argument('--encoders', default='uni,phikon')
args = ap.parse_args()
os.makedirs(args.out, exist_ok=True)
DEV = torch.device(args.device)
ENC_DIR = {'uni': 'uni_max4k', 'phikon': 'phikon_max4k'}
encoders, tasks = args.encoders.split(','), args.tasks.split(',')
seeds = [int(x) for x in args.seeds.split(',')]
t0 = time.time()
def log(*a):
    print(f'[{time.time()-t0:7.1f}s]', *a, flush=True)

# ---------------- data ----------------
tg = pd.read_csv(args.slides)
assert (tg.partition == 'dev').all() and tg.primary_slide.all() and tg.cv_fold.notna().all(), 'non-dev / non-primary slide in input'
assert tg.case_submitter_id.is_unique
have = {e: {f[:-3] for f in os.listdir(os.path.join(args.feat_root, ENC_DIR[e])) if f.endswith('.h5')} for e in ENC_DIR}
both = have['uni'] & have['phikon']
tg['has_features'] = tg.slide_barcode.isin(both)
cov = dict(n_input=len(tg), n_with_both=int(tg.has_features.sum()),
           n_kras_cohort=int(tg.label_kras.notna().sum()), n_kras_cohort_with_both=int((tg.label_kras.notna() & tg.has_features).sum()),
           missing=sorted(tg.loc[~tg.has_features, 'slide_barcode']))
json.dump(cov, open(os.path.join(args.out, 'coverage.json'), 'w'), indent=1)
log('coverage', json.dumps({k: v for k, v in cov.items() if k != 'missing'}))
sl = tg[tg.has_features].reset_index(drop=True)
sl['cv_fold'] = sl.cv_fold.astype(int)

def load(enc, bc):
    with h5py.File(os.path.join(args.feat_root, ENC_DIR[enc], bc + '.h5'), 'r') as f:
        return f['features'][:].astype(np.float16, copy=False), f['coords'][:]

gpu_feats, mean_feats = {}, {}
for enc in encoders:
    with ThreadPoolExecutor(8) as ex:
        res = list(ex.map(lambda b: load(enc, b), sl.slide_barcode))
    if enc == encoders[0]:
        coords0 = [r[1] for r in res]; sl['n_tiles'] = [len(r[0]) for r in res]
    else:
        mis = [b for b, r, c in zip(sl.slide_barcode, res, coords0) if r[1].shape != c.shape or not np.array_equal(r[1], c)]
        log(f'{enc}: coords mismatching {encoders[0]} for {len(mis)} slides', mis[:5])
    mean_feats[enc] = np.stack([r[0].astype(np.float32).mean(0) for r in res])
    gpu_feats[enc] = [torch.from_numpy(r[0]).to(DEV) for r in res]
    log(f'loaded {enc}: {len(res)} slides, dim {res[0][0].shape[1]}, tiles median {int(np.median([len(r[0]) for r in res]))} min {min(len(r[0]) for r in res)}')
    del res
if DEV.type == 'cuda':
    log(f'GPU mem {torch.cuda.memory_allocated()/1e9:.2f} GB on {torch.cuda.get_device_name()}')
for enc in encoders:  # small: lets the LR baseline be re-run locally
    pd.DataFrame(mean_feats[enc], index=sl.slide_barcode).astype(np.float32).to_csv(
        os.path.join(args.out, f'meanpool_{enc}.csv.gz'), float_format='%.5g')
sl.to_csv(os.path.join(args.out, 'slides_used.csv'), index=False)

# ---------------- models ----------------
class ABMIL(nn.Module):
    def __init__(s, d, h=256, a=128, p=0.25):
        super().__init__()
        s.fc = nn.Sequential(nn.Linear(d, h), nn.ReLU(), nn.Dropout(p))
        s.V = nn.Linear(h, a); s.U = nn.Linear(h, a); s.w = nn.Linear(a, 1)
        s.drop = nn.Dropout(p); s.cls = nn.Linear(h, 1)
    def forward(s, x):
        h = s.fc(x)
        A = torch.softmax(s.w(torch.tanh(s.V(h)) * torch.sigmoid(s.U(h))), 0)
        return s.cls(s.drop((A * h).sum(0))).squeeze()

@torch.no_grad()
def predict(m, enc, idx):
    m.eval()
    return np.array([torch.sigmoid(m(gpu_feats[enc][i].float())).item() for i in idx])

def safe_auc(y, p):
    return roc_auc_score(y, p) if len(set(y)) == 2 else float('nan')

def train_abmil(enc, tr, va, y, seed):
    torch.manual_seed(seed); random.seed(seed); np.random.seed(seed)
    m = ABMIL(gpu_feats[enc][0].shape[1]).to(DEV)
    opt = torch.optim.AdamW(m.parameters(), lr=args.lr, weight_decay=args.wd)
    pos = y[tr].sum(); pw = float((len(tr) - pos) / max(pos, 1))
    lossf = nn.BCEWithLogitsLoss(pos_weight=torch.tensor(pw, device=DEV))
    best, best_state, bad, best_ep, yv = 1e9, None, 0, 0, y[va]
    for ep in range(args.max_epochs):
        m.train()
        for i in np.random.permutation(tr):
            x = gpu_feats[enc][i]
            if len(x) > args.train_tiles:
                x = x[torch.randperm(len(x), device=DEV)[:args.train_tiles]]
            loss = lossf(m(x.float()), torch.tensor(float(y[i]), device=DEV))
            opt.zero_grad(); loss.backward(); opt.step()
        pv = predict(m, enc, va)
        vl = float(np.mean(-(yv * np.log(pv + 1e-7) * pw + (1 - yv) * np.log(1 - pv + 1e-7))))
        if vl < best - 1e-4:
            best, bad, best_ep = vl, 0, ep
            best_state = {k: v.detach().clone() for k, v in m.state_dict().items()}
        else:
            bad += 1
            if bad >= args.patience:
                break
    m.load_state_dict(best_state)
    return m, dict(best_epoch=best_ep, epochs_run=ep + 1, val_loss=round(best, 4), val_auc=safe_auc(yv, predict(m, enc, va)))

def fit_lr(enc, tr, va, y):
    X, best = mean_feats[enc], (-1, None)
    for C in [1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 1.0]:
        p = make_pipeline(StandardScaler(), LogisticRegression(C=C, max_iter=5000, class_weight='balanced'))
        p.fit(X[tr], y[tr]); a = safe_auc(y[va], p.predict_proba(X[va])[:, 1])
        if a > best[0]:
            best = (a, C)
    p = make_pipeline(StandardScaler(), LogisticRegression(C=best[1], max_iter=5000, class_weight='balanced'))
    p.fit(X[np.r_[tr, va]], y[np.r_[tr, va]])
    return p, dict(C=best[1], val_auc=best[0])

# ---------------- CV on locked folds ----------------
rows, flog = [], []
def dump():
    pd.DataFrame(rows).to_csv(os.path.join(args.out, 'oof_predictions.csv'), index=False)
    json.dump(flog, open(os.path.join(args.out, 'fold_log.json'), 'w'), indent=1)

for task in tasks:
    y_all = sl[task].astype('float').values  # NaN = not in task cohort
    labelled = ~np.isnan(y_all)
    log(f'{task}: labelled {labelled.sum()} (pos {int(np.nansum(y_all))}), unlabelled scoreable {(~labelled).sum()}')
    for r, seed in enumerate(seeds):
        for k in range(5):
            tr_all = np.where(labelled & (sl.cv_fold.values != k))[0]
            te = np.where(sl.cv_fold.values == k)[0]                      # all scoreable slides of fold k
            strat = (y_all[tr_all].astype(int).astype(str) + '_' + sl.project_id.values[tr_all])
            tr, va = train_test_split(tr_all, test_size=0.15, stratify=strat, random_state=seed + k)
            assert not (set(sl.case_submitter_id[te]) & set(sl.case_submitter_id[tr_all]))
            yl = np.nan_to_num(y_all, nan=-1).astype(int)
            te_lab = te[labelled[te]]
            for enc in encoders:
                for kind in ['abmil', 'meanpool_lr']:
                    if kind == 'meanpool_lr' and r > 0:
                        continue  # deterministic: one repeat is enough
                    mname, ts = f'{kind}_{enc}', time.time()
                    if kind == 'abmil':
                        m, info = train_abmil(enc, tr, va, yl, seed * 100 + k)
                        pte = predict(m, enc, te); del m
                    else:
                        m, info = fit_lr(enc, tr, va, yl)
                        pte = m.predict_proba(mean_feats[enc][te])[:, 1]
                    for i, p in zip(te, pte):
                        rows.append(dict(task=task, model=mname, repeat=r, seed=seed, fold=k, slide_barcode=sl.slide_barcode[i],
                                         case_submitter_id=sl.case_submitter_id[i], project_id=sl.project_id[i], tss=sl.tss[i],
                                         is_discovery=bool(sl.is_discovery[i]), label=y_all[i], prob=float(p)))
                    info.update(task=task, model=mname, repeat=r, fold=k, n_train=len(tr), n_val=len(va), n_test=len(te),
                                n_test_labelled=len(te_lab), test_auc=safe_auc(yl[te_lab], pte[labelled[te]]), sec=round(time.time() - ts, 1))
                    flog.append(info); log(json.dumps(info))
            dump()
log('done')
