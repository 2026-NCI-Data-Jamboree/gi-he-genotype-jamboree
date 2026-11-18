"""Evaluate pre-registered KRAS dev OOF predictions (preregistration.md Sections 6-7).

Usage: python evaluate_kras_mil.py OOF.csv SLIDES_USED.csv PAIRS116.csv READER116.csv PAIRS129.csv ELIG129.csv OUTDIR

- Section 7 baselines are fitted here on the SAME locked folds and the SAME patients (those with features):
  project-only, TSS-only (target-encoded within training folds; unseen site -> training prior),
  MSI-only (msi_score_mantis), stage-only (AJCC ordinal I-IV + missing indicator),
  technical (objective power 20x/40x/NA + log series size). L2 logistic regression, C = 1.
- Per task x model x subset (all_dev, dev_excl_discovery): pooled OOF AUROC (prob averaged over seeds),
  95% CI from 2,000 patient-level bootstrap draws, Brier score, calibration slope and intercept
  (logistic recalibration: slope = coefficient of logit(p); intercept = calibration-in-the-large with
  logit(p) as offset), per-seed and per-fold AUROC, paired-bootstrap AUROC difference vs the best baseline.
- Pairwise accuracy on dev-only pairs: P(prob mut slide > prob paired control), Clopper-Pearson CI;
  model-vs-reader agreement on the dev-only reader pairs.
"""
import sys, os, json
import numpy as np, pandas as pd, statsmodels.api as sm
from scipy.stats import binomtest
from sklearn.metrics import roc_auc_score, brier_score_loss
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

oof_p, sl_p, p116_p, r116_p, p129_p, e129_p, outd = sys.argv[1:8]
os.makedirs(outd, exist_ok=True)
oof = pd.read_csv(oof_p); sl = pd.read_csv(sl_p)
assert (sl.partition == 'dev').all()
TASKS = list(oof.task.unique())
NB = 2000

# ---------------- Section 7 baselines ----------------
def stage_num(s):
    if not isinstance(s, str): return np.nan
    s = s.upper().replace('STAGE', '').strip()
    for k, v in [('IV', 4), ('III', 3), ('II', 2), ('I', 1)]:
        if s.startswith(k): return v
    return np.nan
sl['stage_num'] = sl.ajcc_pathologic_tumor_stage.map(stage_num)
sl['log_size'] = np.log(sl.series_size_MB)

def design(name, tr, te, y):
    if name == 'project':
        f = lambda d: (d.project_id == 'TCGA-READ').astype(float).values[:, None]
        return f(tr), f(te)
    if name == 'msi':
        med = tr.msi_score_mantis.median()
        f = lambda d: d.msi_score_mantis.fillna(med).values[:, None]
        return f(tr), f(te)
    if name == 'stage':
        med = tr.stage_num.median()
        f = lambda d: np.c_[d.stage_num.fillna(med).values, d.stage_num.isna().astype(float).values]
        return f(tr), f(te)
    if name == 'technical':
        f = lambda d: np.c_[(d.ObjectiveLensPower == 40).astype(float), d.ObjectiveLensPower.isna().astype(float), d.log_size.values]
        return f(tr), f(te)
    raise ValueError(name)

base_rows = []
for task in TASKS:
    # same patients as the models: slides with features that carry a model prediction for this task
    ids = set(oof[oof.task == task].slide_barcode)
    d = sl[sl.slide_barcode.isin(ids)].copy()
    for k in range(5):
        tr = d[(d.cv_fold != k) & d[task].notna()]; te = d[d.cv_fold == k]
        ytr = tr[task].astype(int).values
        # TSS-only target encoding
        prior = ytr.mean(); enc = tr.groupby('tss')[task].mean()
        p_tss = te.tss.map(enc).fillna(prior).values
        preds = {'baseline_tss': p_tss}
        for name in ['project', 'msi', 'stage', 'technical']:
            Xtr, Xte = design(name, tr, te, ytr)
            m = make_pipeline(StandardScaler(), LogisticRegression(C=1.0, class_weight='balanced', max_iter=2000)).fit(Xtr, ytr)
            preds[f'baseline_{name}'] = m.predict_proba(Xte)[:, 1]
        for mname, p in preds.items():
            for (_, r), pi in zip(te.iterrows(), p):
                base_rows.append(dict(task=task, model=mname, repeat=0, seed=0, fold=k, slide_barcode=r.slide_barcode,
                                      case_submitter_id=r.case_submitter_id, project_id=r.project_id, tss=r.tss,
                                      is_discovery=bool(r.is_discovery), label=r[task], prob=float(pi)))
allp = pd.concat([oof, pd.DataFrame(base_rows)], ignore_index=True)
allp.to_csv(os.path.join(outd, 'oof_predictions_with_baselines.csv'), index=False)

# seed-averaged slide (= patient, one primary slide each) predictions
avg = allp.groupby(['task', 'model', 'slide_barcode']).agg(case_submitter_id=('case_submitter_id', 'first'), fold=('fold', 'first'),
        is_discovery=('is_discovery', 'first'), label=('label', 'first'), prob=('prob', 'mean'), n_seeds=('repeat', 'nunique')).reset_index()

def calib(y, p):
    lp = np.log(np.clip(p, 1e-6, 1 - 1e-6) / (1 - np.clip(p, 1e-6, 1 - 1e-6)))
    slope = sm.GLM(y, sm.add_constant(lp), family=sm.families.Binomial()).fit().params[1]
    icpt = sm.GLM(y, np.ones((len(y), 1)), family=sm.families.Binomial(), offset=lp).fit().params[0]
    return slope, icpt

rng = np.random.default_rng(20261116)
metrics, boot_store = [], {}
for task in TASKS:
    for subset in ['all_dev', 'dev_excl_discovery']:
        a = avg[(avg.task == task) & avg.label.notna()]
        if subset == 'dev_excl_discovery': a = a[~a.is_discovery.astype(bool)]
        models = sorted(a.model.unique())
        piv = a.pivot(index='case_submitter_id', columns='model', values='prob')
        y = a.groupby('case_submitter_id').label.first().loc[piv.index].astype(int).values
        n = len(y); idx = [rng.integers(0, n, n) for _ in range(NB)]
        idx = [i for i in idx if 0 < y[i].sum() < len(i)]
        boots = {m: np.array([roc_auc_score(y[i], piv[m].values[i]) for i in idx]) for m in models}
        aucs = {m: roc_auc_score(y, piv[m].values) for m in models}
        bases = [m for m in models if m.startswith('baseline_')]
        best_base = max(bases, key=lambda m: aucs[m])
        for m in models:
            p = piv[m].values
            slope, icpt = calib(y, p)
            lo, hi = np.percentile(boots[m], [2.5, 97.5])
            g = allp[(allp.task == task) & (allp.model == m) & allp.label.notna()]
            if subset == 'dev_excl_discovery': g = g[~g.is_discovery.astype(bool)]
            per_seed = [roc_auc_score(q.label.astype(int), q.prob) for _, q in g.groupby('repeat')]
            per_fold = [roc_auc_score(q.label.astype(int), q.prob) if q.label.nunique() == 2 else np.nan
                        for _, q in a[a.model == m].groupby('fold')]
            diff = boots[m] - boots[best_base]
            metrics.append(dict(task=task, subset=subset, model=m, n_patients=n, n_pos=int(y.sum()), auroc=aucs[m], lo=lo, hi=hi,
                                brier=brier_score_loss(y, p), calib_slope=slope, calib_intercept=icpt,
                                auroc_per_seed=';'.join(f'{v:.4f}' for v in per_seed),
                                auroc_per_fold=';'.join(f'{v:.4f}' for v in per_fold),
                                best_baseline=best_base, best_baseline_auroc=aucs[best_base],
                                delta_vs_best_baseline=aucs[m] - aucs[best_base],
                                delta_lo=np.percentile(diff, 2.5), delta_hi=np.percentile(diff, 97.5),
                                p_delta_le0=float((diff <= 0).mean()),
                                ci_excludes_best_baseline=bool(lo > aucs[best_base]), n_boot=len(idx)))
met = pd.DataFrame(metrics)

# ---------------- pairs ----------------
P = avg.pivot_table(index=['task', 'slide_barcode'], columns='model', values='prob')
lab = sl.set_index('slide_barcode')
def cp(k, n):
    ci = binomtest(int(k), int(n)).proportion_ci(0.95, method='exact'); return ci.low, ci.high
def score_pairs(pairs, idcol):
    rows = []
    for task in TASKS:
        models = [m for m in P.columns if not P.loc[task][m].isna().all()]
        for _, r in pairs.iterrows():
            a_, b_ = r.mut_slide_barcode, r.wt_slide_barcode
            both_lab = all(x in lab.index and pd.notna(lab.at[x, task]) for x in (a_, b_))
            for m in models:
                pa = P.at[(task, a_), m] if (task, a_) in P.index else np.nan
                pb = P.at[(task, b_), m] if (task, b_) in P.index else np.nan
                ok = np.nan if (pd.isna(pa) or pd.isna(pb)) else float(pa > pb) + 0.5 * float(pa == pb)
                rows.append({idcol: r[idcol], 'task': task, 'model': m, 'mut_slide_barcode': a_, 'wt_slide_barcode': b_,
                             'both_in_task_cohort': both_lab, 'p_mut_slide': pa, 'p_wt_slide': pb, 'model_correct': ok})
    return pd.DataFrame(rows)
def summ(df, label, extra_mask=None):
    out = []
    for (task, m), g in df.groupby(['task', 'model']):
        for scope in ['all_scored', 'both_in_task_cohort']:
            h = g.dropna(subset=['model_correct'])
            if scope == 'both_in_task_cohort': h = h[h.both_in_task_cohort]
            kf, n = float(h.model_correct.sum()), len(h)   # ties (equal scores) count 0.5
            if n == 0: continue
            k = int(round(kf)); lo, hi = cp(k, n)
            out.append(dict(pairset=label, scope=scope, task=task, model=m, n_pairs=n, correct=kf, n_ties=int((h.model_correct == 0.5).sum()),
                            acc=kf / n, lo=lo, hi=hi, p_vs_chance=binomtest(k, n, 0.5).pvalue))
    return out

part = lab.partition.to_dict()  # only dev slides are in lab; anything else -> not dev
p116 = pd.read_csv(p116_p); r116 = pd.read_csv(r116_p)
p116 = p116[p116.mut_slide_barcode.isin(lab.index) & p116.wt_slide_barcode.isin(lab.index)]
d116 = score_pairs(p116, 'pair').merge(r116[['pair', 'correct', 'confidence', 'basis']].rename(
    columns={'correct': 'reader_correct', 'confidence': 'reader_confidence', 'basis': 'reader_basis'}), on='pair', how='left')
p129 = pd.read_csv(p129_p); e129 = pd.read_csv(e129_p)[['pair_id', 'eligible']]
p129 = p129[p129.mut_slide_barcode.isin(lab.index) & p129.wt_slide_barcode.isin(lab.index)]
d129 = score_pairs(p129, 'pair_id').merge(e129.rename(columns={'eligible': 'reader_used_63'}), on='pair_id', how='left') \
                                    .merge(p129[['pair_id', 'pik3ca_status', 'project_id', 'stage_matched']], on='pair_id', how='left')
pw = summ(d116, f'reader116_dev{len(p116)}') + summ(d129, f'pik3ca129_dev{len(p129)}') + \
     summ(d129[d129.reader_used_63 == True], 'pik3ca129_dev_reader_used')
# reader on the same dev pairs
rr = r116[r116.pair.isin(p116.pair)]
reader = []
for scope, ids in [('all_scored', set(p116.pair)), ('both_in_task_cohort', set(d116[(d116.task == 'label_kras') & d116.both_in_task_cohort].pair))]:
    q = rr[rr.pair.isin(ids)].correct.astype(float); lo, hi = cp(q.sum(), len(q))
    reader.append(dict(pairset=f'reader116_dev{len(p116)}', scope=scope, task='label_kras', model='reader_single',
                       n_pairs=len(q), correct=int(q.sum()), acc=q.mean(), lo=lo, hi=hi, p_vs_chance=binomtest(int(q.sum()), len(q)).pvalue))
pwdf = pd.DataFrame(pw + reader)
agree = []
for (task, m), g in d116.groupby(['task', 'model']):
    for scope in ['all_scored', 'both_in_task_cohort']:
        h = g.dropna(subset=['model_correct', 'reader_correct'])
        if scope == 'both_in_task_cohort': h = h[h.both_in_task_cohort]
        n_ties = int((h.model_correct == 0.5).sum()); h = h[h.model_correct != 0.5]   # tied model scores carry no A/B choice
        if len(h) == 0:
            agree.append(dict(task=task, model=m, scope=scope, n_pairs=0, n_ties_excluded=n_ties)); continue
        mc, rc = h.model_correct.astype(int).values, h.reader_correct.astype(int).values
        n = len(h); po = (mc == rc).mean(); pm, pr = mc.mean(), rc.mean(); pe = pm * pr + (1 - pm) * (1 - pr)
        b, c = int(((mc == 1) & (rc == 0)).sum()), int(((mc == 0) & (rc == 1)).sum())
        agree.append(dict(task=task, model=m, scope=scope, n_pairs=n, n_ties_excluded=n_ties, same_choice=int((mc == rc).sum()), agreement=po,
                          cohen_kappa=(po - pe) / (1 - pe) if pe < 1 else np.nan, both_correct=int(((mc == 1) & (rc == 1)).sum()),
                          model_only=b, reader_only=c, both_wrong=int(((mc == 0) & (rc == 0)).sum()),
                          mcnemar_exact_p=binomtest(b, b + c, 0.5).pvalue if b + c else 1.0))
agdf = pd.DataFrame(agree)

met.to_csv(os.path.join(outd, 'metrics_auroc_calibration.csv'), index=False)
pwdf.to_csv(os.path.join(outd, 'metrics_pairwise.csv'), index=False)
agdf.to_csv(os.path.join(outd, 'model_reader_agreement.csv'), index=False)
d116.to_csv(os.path.join(outd, 'pairs116_dev_model_vs_reader.csv'), index=False)
d129.to_csv(os.path.join(outd, 'pairs129_dev_model.csv'), index=False)
avg.to_csv(os.path.join(outd, 'oof_seed_averaged.csv'), index=False)
pd.set_option('display.width', 250)
print(met[met.task == 'label_kras'][['subset', 'model', 'n_patients', 'auroc', 'lo', 'hi', 'brier', 'calib_slope', 'calib_intercept', 'delta_vs_best_baseline']].round(3).to_string())
print(pwdf[(pwdf.task == 'label_kras')][['pairset', 'scope', 'model', 'n_pairs', 'correct', 'acc', 'lo', 'hi']].round(3).to_string())
