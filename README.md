# GI H&E → molecular phenotype (NCI Data Jamboree 2026)

Predicting molecular phenotypes of gastrointestinal adenocarcinoma from routine H&E
whole-slide images with frozen pathology foundation-model encoders, and testing how far
those predictions travel outside the cohort they were fitted on.

Built for the NCI Data Jamboree, November 2026. This repository holds the code, the
derived result tables, the pre-registration, the blinded reader-study materials and the
figures. It holds no whole-slide images and no model weights.

## What was done

Slides from three NCI commons (TCGA, CPTAC, IDC) were tiled at 0.5 µm/px into 224 px
tiles, embedded with two **frozen** encoders — UNI v1 and Phikon — and aggregated to a
slide vector by mean pooling or attention-MIL. Linear heads were fitted on a
hospital-grouped development partition with a **locked** test partition and a
pre-registered analysis plan (`docs/preregistration.md`), then applied unchanged to two
external cohorts: CPTAC-COAD and SurGen SR386.

Cohort: 1,396 slides, 928 patients (519 colorectal, 409 gastric), 7 endpoints.

## Headline results

Patient- or slide-level AUROC with 95% confidence intervals. Full tables in `results/`.

| Endpoint | Evaluation | Encoder | n | AUROC (95% CI) |
|---|---|---|---|---|
| MSI-H, all GI | locked test | UNI | 191 | 0.826 (0.732–0.905) |
| MSI-H, all GI | locked test | Phikon | 196 | 0.760 (0.667–0.838) |
| MSI-H, colorectal | CPTAC external | UNI | 140 | 0.794 (0.702–0.878) |
| MMR loss (IHC), colorectal | SurGen external | Phikon | 126 | 0.760 (0.638–0.862) |
| CDH1, gastric | dev CV (attention-MIL) | Phikon | 241 | 0.844 (0.775–0.903) |
| TP53, colorectal | dev CV (attention-MIL) | UNI | 324 | 0.772 (0.718–0.824) |
| TP53, colorectal | CPTAC external | Phikon | 80 | 0.642 (0.513–0.770) |
| KRAS, colorectal MSS | dev CV (attention-MIL) | Phikon | 324 | 0.659 (0.596–0.718) |
| KRAS, colorectal MSS | SurGen external | UNI | 90 | 0.579 (0.455–0.713) |

The pattern that matters is the second column, not the first. MSI survives the move to
another hospital system; TP53 and KRAS lose most of their signal. A tile-count-only
confounder control is reported alongside every external evaluation and sits at chance
(0.39–0.51), so the external numbers are not an artefact of how much tissue each slide had.

## The reader study — read this before citing it

`reader_study/` contains a blinded forced-choice experiment on hospital-matched slide
pairs. **The readers were language-model sub-agents shown image fields, not human
pathologists.** They are not independent observers — they share one underlying model — so
the confidence intervals understate the true uncertainty and the result is not a
substitute for a pathologist reader study. Reported as 74/116 pairs called correctly
(63.8%, 95% CI 54.4–72.5%) for colorectal KRAS, with blinded grading showing higher
mucinous and serrated content in the mutant arm. Details and caveats:
`results/mutation/kras_full_cohort_report.md`.

Blinding keys and unblinded pair tables are included in `reader_study/keys/` at the
author's choice, for full reproducibility. Anyone wanting to re-run these packets as a
genuinely blind study should ignore that directory.

## Layout

```
code/            embedding, training, evaluation and analysis scripts
docs/            pre-registration, reports, Jamboree decks and dossier
results/
  msi/           MSI/MMR performance, attribution maps, tile maps
  mutation/      per-gene mutation results, KRAS/CDH1/TP53 analyses
  external/      CPTAC-COAD and SurGen SR386 external validation
  cohort/        slide manifests, locked patient split, label coverage
  literature/    prior-art tables and the literature synthesis
reader_study/    manifests, scored calls, image fields, blinding keys
figures/         figure PNGs used in the deck and reports
data/            DOWNLOAD.md — the large assets kept outside git
notes/           run logs and session state
```

## Reproducing

The pipeline needs the encoders, which this repository does not and cannot ship:

- **Phikon** (`owkin/phikon`) is ungated on Hugging Face.
- **UNI v1** (`MahmoodLab/UNI`) is access-gated; you must accept the terms on the model
  card with an institutional account. `docs/uni_access_request.md` documents that process.

Slides come from the GDC (TCGA, CPTAC), the Imaging Data Commons, and the SurGen
deposition. Embedding was run on cloud GPUs; `code/embed_stream*.py` and
`code/embed_surgen.py` are the entry points, `code/train_kras_mil.py` and
`code/mut_abmil.py` fit the heads, `code/analyse_mut.py` and `code/evaluate_kras_mil.py`
produce the tables in `results/`.

## Data use and attribution

Derived data only; no patient-identifiable information. TCGA and CPTAC data are used under
the NCI open-access tier — cite the GDC and the relevant TCGA/CPTAC marker papers.
SurGen (SR386/SR1482) is distributed as open data under CC BY 4.0; cite
Myles et al., *GigaScience* 2025 (doi:10.1093/gigascience/giaf086). The SurGen label CSVs
are referenced rather than re-hosted — download them from the source deposition.

## Licence

Code in `code/` is MIT (see `LICENSE`). Derived result tables and figures are CC BY 4.0,
subject to the upstream attribution requirements above.

## Status

Work in progress, presented at the NCI Data Jamboree November 2026. Nothing here is
validated for clinical use.
