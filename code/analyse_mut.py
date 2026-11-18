"""Local, free analysis of the ABMIL outputs + pooled-feature linear/MLP models.
Usage: python analyse_mut.py <job_outdir>"""
import sys, json, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, log_loss
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
import warnings; warnings.filterwarnings('ignore')
D = sys.argv[1]; rng = np.random.default_rng(0)
lab = pd.read_csv('mut_labels.csv'); TASKS = json.load(open('tasks.json'))
oof = pd.read_csv(f'{D}/oof_preds.csv'); tst = pd.read_csv(f'{D}/test_preds.csv')
Z = np.load(f'{D}/pooled_means.npz', allow_pickle=True); zs = list(Z['slides'])
L = lab.set_index('slide_barcode')

def boot_auc(y, s, n=1000):
    y, s = np.asarray(y), np.asarray(s); a = roc_auc_score(y, s); b = []
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if 0 < y[i].sum() < len(i): b.append(roc_auc_score(y[i], s[i]))
    return a, np.percentile(b, 2.5), np.percentile(b, 97.5)

# 1. per-gene performance: ABMIL (OOF dev, locked test) + linear + MLP on mean features (same folds)
rows, linoof = [], []
for group, tasks in TASKS.items():
    for tk in tasks:
        for enc in ('phikon', 'uni'):
            o = oof[(oof.task == tk) & (oof.enc == enc)]
            a, lo, hi = boot_auc(o.y, o.score)
            t = tst[(tst.task == tk) & (tst.enc == enc)]
            ta = roc_auc_score(t.y, t.score) if t.y.sum() >= 15 and (t.y == 0).sum() >= 15 else np.nan
            d = L.loc[o.slide]; X = Z[enc][[zs.index(s) for s in o.slide]]; y = o.y.values; f = d.cv_fold.values
            for name, mk in (('linear', lambda: make_pipeline(StandardScaler(), LogisticRegression(C=0.01, class_weight='balanced', max_iter=2000))),
                             ('mlp', lambda: make_pipeline(StandardScaler(), MLPClassifier((128,), alpha=1.0, max_iter=400, early_stopping=False, random_state=0)))):
                s = np.zeros(len(y))
                for k in range(5):
                    tr, te = f != k, f == k
                    s[te] = mk().fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
                linoof.append(pd.DataFrame({'enc': enc, 'task': tk, 'model': name, 'slide': o.slide.values, 'y': y, 'score': s}))
                rows.append({'group': group, 'task': tk.split('_', 1)[1], 'enc': enc, 'model': name,
                             'n_mut': int(y.sum()), 'n_wt': int((y == 0).sum()), 'auc_dev_cv': round(roc_auc_score(y, s), 3)})
            # site-only baseline: TSS one-hot (folds are hospital-disjoint, so this should sit near 0.5)
            rows.append({'group': group, 'task': tk.split('_', 1)[1], 'enc': enc, 'model': 'abmil',
                         'n_mut': int(o.y.sum()), 'n_wt': int((o.y == 0).sum()), 'auc_dev_cv': round(a, 3),
                         'ci_lo': round(lo, 3), 'ci_hi': round(hi, 3), 'auc_locked_test': None if np.isnan(ta) else round(ta, 3)})
perf = pd.DataFrame(rows); linoof = pd.concat(linoof)
# stomach molecular-subtype baseline (EBV/GS/CIN known from genomics; how much does subtype alone explain?)
sub_rows = []
for tk in TASKS['Stomach']:
    o = oof[(oof.task == tk) & (oof.enc == 'uni')]; st = L.loc[o.slide, 'pancan_molecular_subtype'].fillna('NA')
    X = pd.get_dummies(st).values.astype(float); y = o.y.values; f = L.loc[o.slide, 'cv_fold'].values; s = np.zeros(len(y))
    for k in range(5): s[f == k] = LogisticRegression(class_weight='balanced').fit(X[f != k], y[f != k]).predict_proba(X[f == k])[:, 1]
    sub_rows.append({'task': tk, 'auc_subtype_only': round(roc_auc_score(y, s), 3)})
    print(tk, pd.crosstab(st, y).to_dict())
perf.to_csv('mut_performance.csv', index=False); pd.DataFrame(sub_rows).to_csv('stad_subtype_baseline.csv', index=False)

# 2. co-mutation stratification: gene-A ABMIL score evaluated within gene-B strata, + adjusted logistic
pairs = [('crc_KRAS', 'crc_PIK3CA'), ('crc_PIK3CA', 'crc_KRAS'), ('crc_TP53', 'crc_APC'), ('crc_APC', 'crc_TP53'),
         ('crc_PIK3CA', 'crc_TP53'), ('crc_TP53', 'crc_PIK3CA'), ('crc_SOX9', 'crc_TP53'), ('crc_SOX9', 'crc_PIK3CA'),
         ('crc_KRAS', 'crc_TP53'), ('stad_ARID1A', 'stad_TP53'), ('stad_TP53', 'stad_ARID1A'), ('stad_PIK3CA', 'stad_TP53'),
         ('stad_CDH1', 'stad_TP53'), ('stad_ARID1A', 'stad_PIK3CA')]
srows = []
import statsmodels.api as sm
for A, B in pairs:
    for enc in ('phikon', 'uni'):
        o = oof[(oof.task == A) & (oof.enc == enc)].copy(); o['B'] = L.loc[o.slide, B].values
        r = {'gene_A': A, 'stratified_by': B, 'enc': enc, 'auc_all': round(roc_auc_score(o.y, o.score), 3)}
        for v, nm in ((0, 'B_wt'), (1, 'B_mut')):
            q = o[o.B == v]
            r[f'auc_within_{nm}'] = round(roc_auc_score(q.y, q.score), 3) if q.y.nunique() == 2 and min(q.y.sum(), (q.y == 0).sum()) >= 8 else None
            r[f'n_mut_A_{nm}'] = int(q.y.sum())
        # does the A-score predict A after adjusting for B status? (logit, A ~ score_z + B)
        zsc = (o.score - o.score.mean()) / o.score.std()
        fit = sm.Logit(o.y.values, sm.add_constant(np.column_stack([zsc, o.B.values]))).fit(disp=0)
        r['OR_per_SD_adjusted_for_B'] = round(float(np.exp(fit.params[1])), 2); r['p_adjusted'] = float(fit.pvalues[1])
        srows.append(r)
strat = pd.DataFrame(srows); strat.to_csv('comut_stratified_auc.csv', index=False)

# 3. KRAS x PIK3CA joint genotype: additive (two independent binary models) vs joint 4-class ABMIL
jrows = []
for enc in ('phikon', 'uni'):
    k = oof[(oof.task == 'crc_KRAS') & (oof.enc == enc)].set_index('slide').score
    p = oof[(oof.task == 'crc_PIK3CA') & (oof.enc == enc)].set_index('slide').score
    j = oof[(oof.task == 'crc_KRASxPIK3CA') & (oof.enc == enc)].set_index('slide')
    s = j.index.intersection(k.index).intersection(p.index); y = j.loc[s, 'y'].astype(int).values
    Pj = j.loc[s, ['p0', 'p1', 'p2', 'p3']].values
    # additive: treat the two genes as independent -> product of marginals (each calibrated on OOF via Platt, cross-fitted)
    f = L.loc[s, 'cv_fold'].values; kk, pp = k.loc[s].values, p.loc[s].values
    yk, yp = (y == 1) | (y == 3), (y == 2) | (y == 3)
    ck, cp = np.zeros(len(s)), np.zeros(len(s))
    for fo in range(5):
        tr, te = f != fo, f == fo
        ck[te] = LogisticRegression().fit(kk[tr, None], yk[tr]).predict_proba(kk[te, None])[:, 1]
        cp[te] = LogisticRegression().fit(pp[tr, None], yp[tr]).predict_proba(pp[te, None])[:, 1]
    Pa = np.column_stack([(1 - ck) * (1 - cp), ck * (1 - cp), (1 - ck) * cp, ck * cp])
    # interaction model on the two scores: is there a KRAS x PIK3CA term beyond additive?
    X = np.column_stack([kk, pp]); Xi = np.column_stack([kk, pp, kk * pp]); Pm, Pmi = np.zeros((len(s), 4)), np.zeros((len(s), 4))
    for fo in range(5):
        tr, te = f != fo, f == fo
        Pm[te] = LogisticRegression(max_iter=2000).fit(X[tr], y[tr]).predict_proba(X[te])
        Pmi[te] = LogisticRegression(max_iter=2000).fit(Xi[tr], y[tr]).predict_proba(Xi[te])
    both = (y == 3).astype(int)
    jrows.append({'enc': enc, 'n': len(s), 'logloss_additive_product': round(log_loss(y, Pa, labels=[0, 1, 2, 3]), 4),
                  'logloss_scores_no_interaction': round(log_loss(y, Pm, labels=[0, 1, 2, 3]), 4),
                  'logloss_scores_with_interaction': round(log_loss(y, Pmi, labels=[0, 1, 2, 3]), 4),
                  'logloss_joint_4class_abmil': round(log_loss(y, Pj, labels=[0, 1, 2, 3]), 4),
                  'logloss_prevalence_only': round(log_loss(y, np.tile(np.bincount(y, minlength=4) / len(y), (len(y), 1)), labels=[0, 1, 2, 3]), 4),
                  'auc_double_vs_rest_additive': round(roc_auc_score(both, Pa[:, 3]), 3),
                  'auc_double_vs_rest_joint': round(roc_auc_score(both, Pj[:, 3]), 3),
                  'auc_double_vs_KRASonly_joint': round(roc_auc_score(both[(y == 1) | (y == 3)], (Pj[:, 3] / (Pj[:, 1] + Pj[:, 3]))[(y == 1) | (y == 3)]), 3)})
joint = pd.DataFrame(jrows); joint.to_csv('kras_pik3ca_joint.csv', index=False)
print(perf.to_string()); print(strat.to_string()); print(joint.T.to_string())
