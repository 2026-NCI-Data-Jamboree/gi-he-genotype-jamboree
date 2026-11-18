# Blinded KRAS test: 6 matched pairs, 18 readers

## Design
Six KRAS-mutant vs KRAS/BRAF-wild-type pairs from TCGA COAD/READ, drawn at random. Both slides in each pair are MSS
(non-hypermutated, MANTIS < 0.4), 40x, single-DX-slide patients, from the **same tissue source site** and the **same AJCC
stage group** — so hospital, staining batch, MSI and stage cannot separate them. Slot A/B assignment was randomised per pair
and sealed (../../reader_study/keys/KEY_kras.json).

Fields: 517 um at 0.505 um/px, cut from the 0.25 um/px base level with wsidicom, sampled at random from cellular tissue,
then screened field-by-field (blind to arm) for epithelial content; stroma/muscle-only fields were dropped. Three readers
per pair, five fields per slide each (non-overlapping where the slide allowed; **pair 6 readers all saw the same 5 fields**
because the slides are small). Readers were told KRAS is hard, given the reported weak priors (mucinous, serrated, villous,
budding), and asked for a calibrated probability.

## Result
| pair | site | stage | KRAS-mutant | wild-type | mean p(true mutant) | readers correct |
|---|---|---|---|---|---|---|
| 1 | TCGA-AY | III | TCGA-AY-6386 (G12D) | TCGA-AY-4070 | 0.64 | 3/3 |
| 2 | TCGA-AZ | IV | TCGA-AZ-6600 (G13D) | TCGA-AZ-4682 | 0.60 | 3/3 |
| 3 | TCGA-AD | III | TCGA-AD-6965 (G12A) | TCGA-AD-6899 | 0.51 | 2/3 |
| 4 | TCGA-D5 | II | TCGA-D5-5537 (A146T) | TCGA-D5-6924 | 0.57 | 3/3 |
| 5 | TCGA-CA | II | TCGA-CA-5796 (G12V) | TCGA-CA-5256 | 0.66 | 3/3 |
| 6 | TCGA-EI | IV | TCGA-EI-7002 (G12V) | TCGA-EI-7004 | 0.54 | 3/3 |

- **Pair level: 6/6 correct** (mean reader probability > 0.5 on the true mutant). Two-sided binomial p = 0.031.
- Reader level: 17/18 correct; readers within a pair are not independent, so the pair count is the honest n.
- Calibration: mean probability on the true mutant **0.59**; Brier 0.174 vs 0.25 for an uninformative reader. Readers
  were correct but not confident — every call was 52–67%.

## What the readers used
| feature (0–3, mean of 18) | KRAS-mut | WT | delta | pairs mut>wt – wt>mut |
|---|---|---|---|---|
| mucinous differentiation | 1.44 | 0.22 | +1.22 | 5–1 |
| serrated architecture | 0.39 | 0.11 | +0.28 | 3–0 |
| villous or papillary | 1.61 | 0.78 | +0.83 | 4–1 |
| tumour budding or poorly diff clusters | 0.67 | 1.50 | -0.83 | 1–5 |
| intraepithelial lymphocytes | 0.39 | 0.44 | -0.06 | 1–2 |
| dirty necrosis | 0.89 | 1.22 | -0.33 | 2–4 |
| nuclear pleomorphism | 1.56 | 2.00 | -0.44 | 1–3 |
| eosinophilic cytoplasm | 1.28 | 1.50 | -0.22 | 2–2 |

Gland formation: KRAS-mutant 6 well / 12 moderate / 0 poor; WT 0 well / 10 moderate / 8 poor.

Every rationale applied the same rule: **mucinous and villous/low-grade, adenoma-like = KRAS; solid, pleomorphic, budding,
poorly differentiated = wild-type.** In these six pairs that rule happened to point the right way every time.

## Why this is not yet evidence that KRAS is readable from H&E
1. **n = 6 pairs.** 6/6 has p = 0.031, but the 95% CI on pair accuracy is 54–100%. A real accuracy of 65–70% — the
   published ceiling — produces 6/6 about 8–12% of the time (0.65^6 = 0.075, 0.70^6 = 0.118).
2. **Two pairs were decided by tissue quality, not tumour biology.** In pair 4 the WT slide had severe crush/chatter
   artefact; in pair 6 the WT slide was out of focus with freezing artefact and readers said so explicitly. Readers defaulted to
   the readable slide. Excluding those two leaves **4/4** (p = 0.125).
3. **The readers were primed.** I supplied the literature priors; the swarm applied them consistently. What was tested is
   whether *grade and mucinous/villous phenotype* track KRAS in these pairs, which is a known weak population association —
   not whether the model sees something KRAS-specific.
4. **Grade is a plausible confounder.** The WT arm here skews poorly differentiated (8/18 vs 0/18 'poor'). Whether that
   holds across the 1,152-patient cohort is checkable and is the obvious next step.
5. Most readers viewed some fields downsampled to 700–900 px because of image-size limits; nuclear-detail scores are
   provisional. Architecture and mucin calls were less affected.

## Next step that would make this a result
20–30 pairs, matched additionally on **differentiation grade** (or with grade recorded and regressed out), with section-quality
QC before pairing, and readers given **no** priors. If calls stay above chance with grade matched, that is a KRAS-linked
morphologic signal; if they collapse to 50%, the swarm was reading grade.
