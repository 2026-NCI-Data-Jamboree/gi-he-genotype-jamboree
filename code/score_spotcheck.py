#!/usr/bin/env python3
"""Score a pathologist's filled spot-check sheet against the blinded AI (sub-agent) grades.

Usage:
  python score_spotcheck.py --sheet spotcheck_grading_sheet_FILLED.xlsx \
      --key spotcheck_key_PRIVATE.csv --ai kras_full_slide_grades.csv \
      [--out-prefix spotcheck_agreement] [--n-boot 2000] [--seed 1]

Inputs
  sheet : the grading workbook (sheet 'grading'; one row per P code; one column per item).
  key   : private key CSV with columns p_code, s_code, slide, arm (arm = mut / wt).
  ai    : AI grades CSV (column slide_id = S code, plus one column per item).

Outputs (printed and written)
  <prefix>_report.txt         human-readable report
  <prefix>_per_item.csv       one row per item: n, kappa_qw (+ bootstrap 95% CI), spearman rho/p,
                              exact agreement, within-1 (ordinal) or within-10-points (% items)
  <prefix>_mut_vs_wt.csv      Mann-Whitney U on the pathologist's own scores, mutant vs wild-type

Dependencies: pandas, numpy, scipy, openpyxl (no scikit-learn needed).
"""
import argparse
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

ORDINAL = {
    "differentiation": ["well", "moderate", "poor"],
    "signet_ring": [0, 1, 2, 3],
    "villous_or_papillary": [0, 1, 2, 3],
    "serrated": [0, 1, 2, 3],
    "tumour_budding": [0, 1, 2, 3],
    "poorly_diff_clusters": [0, 1, 2, 3],
    "dirty_necrosis": [0, 1, 2, 3],
    "intraepithelial_lymphocytes": [0, 1, 2, 3],
    "peritumoral_lymphoid": [0, 1, 2, 3],
    "nuclear_grade": [1, 2, 3],
    "stroma_desmoplasia": [0, 1, 2, 3],
    "adenoma_component": [0, 1, 2, 3],
    "section_quality": ["good", "fair", "poor"],
}
PERCENT = ["gland_formation_pct", "mucinous_pct"]
ITEMS = ["differentiation", "gland_formation_pct", "mucinous_pct", "signet_ring", "villous_or_papillary",
         "serrated", "tumour_budding", "poorly_diff_clusters", "dirty_necrosis", "intraepithelial_lymphocytes",
         "peritumoral_lymphoid", "nuclear_grade", "stroma_desmoplasia", "adenoma_component", "section_quality"]
ARM_TESTS = ["mucinous_pct", "serrated", "dirty_necrosis"]
PCT_TOL = 10  # percentage points counted as "close" for % items


def to_codes(series, levels):
    """Map raw values to integer category codes 0..K-1 (NaN if missing / invalid)."""
    if isinstance(levels[0], str):
        m = {v: i for i, v in enumerate(levels)}
        return series.astype(str).str.strip().str.lower().map(m)
    num = pd.to_numeric(series, errors="coerce")
    m = {v: i for i, v in enumerate(levels)}
    return num.map(lambda x: m.get(int(x)) if pd.notna(x) and float(x).is_integer() else np.nan)


def kappa_qw(a, b, k):
    """Quadratically weighted Cohen kappa for integer codes 0..k-1 on a fixed k-level scale."""
    a = np.asarray(a, int); b = np.asarray(b, int)
    O = np.zeros((k, k))
    np.add.at(O, (a, b), 1)
    n = O.sum()
    E = np.outer(O.sum(1), O.sum(0)) / n
    i, j = np.indices((k, k))
    W = (i - j) ** 2 / (k - 1) ** 2
    den = (W * E).sum()
    if den == 0:
        return np.nan  # both raters constant on the same category
    return 1.0 - (W * O).sum() / den


def boot_ci(fn, a, b, n_boot, rng):
    n = len(a); vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        v = fn(a[idx], b[idx])
        if np.isfinite(v):
            vals.append(v)
    if len(vals) < n_boot * 0.5:
        return (np.nan, np.nan)
    return tuple(np.percentile(vals, [2.5, 97.5]))


def safe_spearman(a, b):
    if len(a) < 3 or np.std(a) == 0 or np.std(b) == 0:
        return np.nan, np.nan
    r = spearmanr(a, b)
    return float(r.statistic if hasattr(r, "statistic") else r[0]), float(r.pvalue if hasattr(r, "pvalue") else r[1])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--key", required=True)
    ap.add_argument("--ai", required=True)
    ap.add_argument("--out-prefix", default="spotcheck_agreement")
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args(argv)
    rng = np.random.default_rng(args.seed)

    path = pd.read_excel(args.sheet, sheet_name="grading", dtype=object)
    path.columns = [str(c).strip() for c in path.columns]
    path["p_code"] = path["p_code"].astype(str).str.strip().str.upper()
    key = pd.read_csv(args.key)
    ai = pd.read_csv(args.ai)
    missing_cols = [c for c in ITEMS if c not in path.columns or c not in ai.columns]
    if missing_cols:
        sys.exit(f"ERROR: columns missing from sheet or AI file: {missing_cols}")

    d = key.merge(path, on="p_code", how="left").merge(
        ai[["slide_id"] + ITEMS], left_on="s_code", right_on="slide_id", how="left", suffixes=("_path", "_ai"))
    lines = []
    out = lines.append
    out("Pathologist vs AI-grader agreement - colorectal spot-check")
    out(f"slides in key: {len(key)} | rows in sheet: {len(path)} | key P codes absent from sheet: "
        f"{sorted(set(key.p_code) - set(path.p_code))} | AI rows matched: {d['slide_id'].notna().sum()}")
    out(f"arms: {key.arm.value_counts().to_dict()}")
    out("")

    rows = []
    for it in ITEMS:
        pa, aa = d[f"{it}_path"], d[f"{it}_ai"]
        if it in ORDINAL:
            lv = ORDINAL[it]; k = len(lv)
            pc, ac = to_codes(pa, lv), to_codes(aa, lv)
            ok = pc.notna() & ac.notna()
            x, y = pc[ok].astype(int).to_numpy(), ac[ok].astype(int).to_numpy()
            n_invalid = int((pa.notna() & pc.isna()).sum())
            kq = kappa_qw(x, y, k) if len(x) else np.nan
            lo, hi = boot_ci(lambda a, b: kappa_qw(a, b, k), x, y, args.n_boot, rng) if len(x) >= 5 else (np.nan, np.nan)
            rho, p = safe_spearman(x, y)
            exact = float(np.mean(x == y)) if len(x) else np.nan
            close = float(np.mean(np.abs(x - y) <= 1)) if len(x) else np.nan
            rows.append(dict(item=it, type="ordinal", n=len(x), n_invalid_entries=n_invalid, kappa_qw=kq,
                             kappa_ci_lo=lo, kappa_ci_hi=hi, spearman_rho=rho, spearman_p=p,
                             exact_agree=exact, close_agree=close, close_rule="|diff|<=1 category",
                             path_mean=float(np.mean(x)) if len(x) else np.nan,
                             ai_mean=float(np.mean(y)) if len(y) else np.nan,
                             note=("AI constant on these slides - kappa uninformative; use exact agreement"
                                   if len(y) and np.unique(y).size == 1 else
                                   "pathologist constant - kappa uninformative" if len(x) and np.unique(x).size == 1 else "")))
        else:
            pn, an = pd.to_numeric(pa, errors="coerce"), pd.to_numeric(aa, errors="coerce")
            ok = pn.notna() & an.notna() & pn.between(0, 100)
            x, y = pn[ok].to_numpy(float), an[ok].to_numpy(float)
            rho, p = safe_spearman(x, y)
            rows.append(dict(item=it, type="percent", n=len(x),
                             n_invalid_entries=int((pa.notna() & ~(pn.between(0, 100))).sum()),
                             kappa_qw=np.nan, kappa_ci_lo=np.nan, kappa_ci_hi=np.nan,
                             spearman_rho=rho, spearman_p=p,
                             exact_agree=float(np.mean(x == y)) if len(x) else np.nan,
                             close_agree=float(np.mean(np.abs(x - y) <= PCT_TOL)) if len(x) else np.nan,
                             close_rule=f"|diff|<={PCT_TOL} points",
                             path_mean=float(np.mean(x)) if len(x) else np.nan,
                             ai_mean=float(np.mean(y)) if len(y) else np.nan,
                             mean_abs_diff=float(np.mean(np.abs(x - y))) if len(x) else np.nan))
    res = pd.DataFrame(rows)
    out("Per-item agreement (pathologist vs AI). kappa_qw = quadratically weighted Cohen kappa with bootstrap 95% CI;")
    out(f"close = within 1 category (ordinal) or within {PCT_TOL} percentage points (% items).")
    fmt = res.copy()
    for c in ["kappa_qw", "kappa_ci_lo", "kappa_ci_hi", "spearman_rho", "spearman_p", "exact_agree", "close_agree",
              "path_mean", "ai_mean"]:
        fmt[c] = fmt[c].map(lambda v: "" if pd.isna(v) else f"{v:.3f}")
    out(fmt[["item", "type", "n", "n_invalid_entries", "kappa_qw", "kappa_ci_lo", "kappa_ci_hi", "spearman_rho",
             "spearman_p", "exact_agree", "close_agree", "path_mean", "ai_mean"]].to_string(index=False))
    for _, r in res[res.get("note", pd.Series(dtype=str)).fillna("").astype(str) != ""].iterrows():
        out(f"  note - {r['item']}: {r['note']}")
    out("")

    lv = ORDINAL["differentiation"]
    pc = pd.Categorical(d["differentiation_path"].astype(str).str.strip().str.lower(), categories=lv)
    ac = pd.Categorical(d["differentiation_ai"].astype(str).str.strip().str.lower(), categories=lv)
    ct = pd.crosstab(pd.Series(pc, name="pathologist"), pd.Series(ac, name="AI"), dropna=False)
    ct = ct.reindex(index=lv, columns=lv, fill_value=0)
    out("Differentiation confusion table (rows = pathologist, columns = AI):")
    out(ct.to_string())
    out("")

    arm_rows = []
    for it in ARM_TESTS:
        for src in ("path", "ai"):
            v = pd.to_numeric(d[f"{it}_{src}"], errors="coerce")
            m, w = v[(d.arm == "mut") & v.notna()], v[(d.arm == "wt") & v.notna()]
            if len(m) and len(w) and (pd.concat([m, w]).nunique() > 1):
                u = mannwhitneyu(m, w, alternative="two-sided")
                U, pv = float(u.statistic), float(u.pvalue)
            else:
                U, pv = np.nan, np.nan
            arm_rows.append(dict(item=it, rater="pathologist" if src == "path" else "AI (same 30 slides, reference)",
                                 n_mut=len(m), n_wt=len(w), mean_mut=m.mean(), mean_wt=w.mean(),
                                 median_mut=m.median(), median_wt=w.median(),
                                 frac_pos_mut=float((m > 0).mean()) if len(m) else np.nan,
                                 frac_pos_wt=float((w > 0).mean()) if len(w) else np.nan,
                                 U=U, p_mannwhitney=pv, rank_biserial=(2 * U / (len(m) * len(w)) - 1) if np.isfinite(U) else np.nan))
    arm = pd.DataFrame(arm_rows)
    out("Mutant vs wild-type, two-sided Mann-Whitney U (pathologist's own scores; AI rows on the same slides for reference).")
    out("rank_biserial > 0 means higher in mutant.")
    out("CAVEAT: the 30 slides were chosen with per-arm quotas on mucinous_pct>0, serrated>=1 and dirty_necrosis>=2")
    out("(by AI grade), so arm differences are partly balanced away by design. Read these tests as a check that the")
    out("pathologist's scores move in the same direction as the AI's on the same slides, not as a cohort-level estimate.")
    out(arm.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    res.to_csv(f"{args.out_prefix}_per_item.csv", index=False)
    arm.to_csv(f"{args.out_prefix}_mut_vs_wt.csv", index=False)
    ct.to_csv(f"{args.out_prefix}_differentiation_confusion.csv")
    text = "\n".join(lines)
    with open(f"{args.out_prefix}_report.txt", "w") as fh:
        fh.write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
