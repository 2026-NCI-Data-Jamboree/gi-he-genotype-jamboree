# Pre-registration: Reading the Genome Off the Slide (NCI Jamboree, Nov 2026)

**Status:** LOCKED before any model training, feature extraction review or test-set inspection. Date locked: 2026-11-15.
**Input:** TCGA GI FFPE diagnostic H&E slide manifest (1,225 slides, 1,188 patients; COAD, READ, STAD, ESCA; IDC DICOM).
**Split file (authoritative):** [patient_split_locked.csv](../results/cohort/patient_split_locked.csv) / [patient_split_locked.parquet](../results/cohort/patient_split_locked.parquet)
**Balance report:** [split_balance_report.csv](../results/cohort/split_balance_report.csv)
**Gene count table:** [tier1_gene_counts_by_stratum.csv](../results/mutation/tier1_gene_counts_by_stratum.csv)

Any deviation from this document after locking must be logged in a "Deviations" section with its date and reason.

---

## 1. Unit of analysis, eligibility, slide choice

- **Unit = patient.** Every slide of a patient sits in the same partition and fold. Labels are patient-level.
- **Label-usable patients:** `use_for_training == True` (excludes 13 patients without an MSI call because GDC WXS and WGS disagree, and 13 with a MANTIS/GDC MSI discordance) **and** slide QC pass. That leaves 1,159 of 1,188 patients.
- **QC exclusions** (`qc_flag`): TCGA-D5-6924 (crush artefact), TCGA-EI-7004 (out of focus, ~13 mm² tissue), TCGA-EI-7002 (~27 mm² tissue, below the 30 mm² floor). Each of these patients has one slide only, so they carry NA for every label. Every other slide has `qc_flag = 'unchecked'`; slides failing tissue segmentation at feature extraction will be logged as deviations.
- **Sequenced** = `mutation_count` not null. 105 patients are unsequenced (COAD 66, READ 31, STAD 7, ESCA 1). In 64 of them the source `mut_*` flags read 0, which is false; they are NA for every gene task (including KRAS).
- **Hypermutation:** gene tasks exclude `hypermutated_tmb_ge10 == True`. Among sequenced MSS patients this removes COAD 10, READ 4, STAD 15 and ESCA 4 (POLE or other non-MSI hypermutators).
- **Primary slide** (`primary_slide == True`, one per patient): DX1 from the primary tumour (sample type 01), else the smallest series. TCGA-NH-A6GC has a DX1 from both the primary (01Z) and a metastasis (06Z); the 01Z slide is primary. Training and evaluation use the primary slide. Non-primary slides stay in the file, in the same partition, and are used only in a pre-specified slide-aggregation sensitivity analysis.
- **Mutation flag definition:** `mut_<GENE> = True` means any non-silent call in the TCGA PanCancer Atlas MAF (via cBioPortal). Classes observed in the protein-change strings: missense, nonsense, frameshift, splice and in-frame indel; no synonymous calls. Flags are **not** filtered to oncogenic or hotspot variants, so TP53, APC, SMAD4 and CDH1 positives include some VUS and passenger missense calls. An OncoKB-filtered relabel is a possible later sensitivity analysis and is not primary.

## 2. Tasks and exact label definitions

Counts are label-usable patients, positives/negatives, from [split_balance_report.csv](../results/cohort/split_balance_report.csv).

| Task column | Definition | Status | All | Locked test | Dev | Dev excl. discovery |
|---|---|---|---|---|---|---|
| `label_msi` | MSI-H (1) vs MSS (0); COAD+READ+STAD; ESCA excluded (only 2 MSI-H) | primary | 162/842 | 31/165 | 131/677 | 101/602 |
| `label_kras` | COAD+READ, MSS, sequenced, non-hypermutated; 1 = `mut_KRAS`; 0 = KRAS not mutated (includes 14 BRAF-mutant) | primary | 176/243 | 35/49 | 141/194 | 111/164 |
| `label_kras_sens` | as `label_kras` but negatives must also be BRAF-WT | sensitivity | 176/229 | 35/44 | 141/185 | 111/155 |
| `label_tp53_coadread` | COAD+READ, MSS, seq, non-hyper; 1 = `mut_TP53` | primary | 270/149 | 55/29 | 215/120 | 172/103 |
| `label_tp53_stad` | STAD, MSS, seq, non-hyper; 1 = `mut_TP53` | primary | 161/156 | 35/33 | 126/123 | 115/120 |
| `label_pik3ca_coadread` | COAD+READ, MSS, seq, non-hyper; 1 = `mut_PIK3CA` | primary | 98/321 | 20/64 | 78/257 | 65/210 |
| `label_apc_coadread` | COAD+READ, MSS, seq, non-hyper; 1 = `mut_APC` | primary | 321/98 | 64/20 | 257/78 | 202/73 |
| `label_pik3ca_stad` | STAD, MSS, seq, non-hyper; 1 = `mut_PIK3CA` | exploratory | 31/286 | 8/60 | 23/226 | 21/214 |
| `label_smad4_coadread` | COAD+READ, MSS, seq, non-hyper; 1 = `mut_SMAD4` | exploratory | 47/372 | 11/73 | 36/299 | 28/247 |
| `label_cdh1_stad` | STAD, MSS, seq, non-hyper; 1 = `mut_CDH1`; **pre-declared positive control** | exploratory (positive control) | 30/287 | 7/61 | 23/226 | 21/214 |

MSI by cohort (MSI-H/MSS): COAD 76/352, READ 8/152, STAD 78/338.

**KRAS hotspot class** (`kras_hotspot_class`, recorded and not modelled as a primary endpoint): G12V 44, G12D 43, G13D 28, other G12/G13 22, G12C 13, A146 10, Q61 6, other 7, multiple 3; WT 229, WT_BRAFmut 14.

## 3. Inclusion rule for gene tasks and outcome

Rule (fixed before counting): smaller arm ≥ 50 patients → **primary**; 30–49 → **exploratory**; < 30 → **dropped**. Eligibility: label-usable, sequenced, MSS, non-hypermutated.

| Gene × cohort | n mut | n WT | % MSI-H among mutants (unrestricted, sequenced) | Cohort MSI-H baseline | Outcome |
|---|---|---|---|---|---|
| TP53 COADREAD | 270 | 149 | 7.1 | 14.1 | primary |
| TP53 STAD | 161 | 156 | 13.4 | 18.8 | primary |
| PIK3CA COADREAD | 98 | 321 | 20.4 | 14.1 | primary |
| APC COADREAD | 321 | 98 | 8.5 | 14.1 | primary |
| PIK3CA STAD | 31 | 286 | 47.0 | 18.8 | exploratory |
| SMAD4 COADREAD | 47 | 372 | 20.0 | 14.1 | exploratory |
| CDH1 STAD | 30 | 287 | 23.1 | 18.8 | exploratory (positive control) |
| SMAD4 STAD | 22 | 295 | 20.0 | 18.8 | **dropped** |

Screened but not planned, all dropped: APC STAD (20 mutant), CDH1 COADREAD (2), and every gene in ESCA adenocarcinoma (n = 55 eligible) and ESCA squamous (n = 82). In ESCA, TP53 has only 8 and 5 WT patients respectively. ESCA histology from the PanCancer subtype: squamous (ESCC) 85, adenocarcinoma (CIN/GS/MSI) 58, missing 13. ESCA patients are in the split but carry no task label. They are not used for any primary task.

The % MSI-H among mutants shows that PIK3CA STAD (47%), and to a lesser degree PIK3CA, SMAD4 and CDH1, are MSI-enriched in the unrestricted cohort. That is why every gene task is restricted to MSS and benchmarked against an MSI-only baseline.

## 4. Split (one patient-level split shared by all tasks)

- **Grouping:** tissue source site (TSS = barcode field 2; 77 TSSs; no TSS spans projects). No TSS appears in both test and dev, and no TSS spans two dev folds.
- **Composite stratification label** (`strat_composite`): project × MSI (MSI-H/MSS/NA) × KRAS (mut/WT/NA) for COAD and READ; project × MSI for STAD; project × histology for ESCA.
- **Locked TEST (seed 20261115).**
  1. 20,000 random TSS subsets drawn with `numpy.random.default_rng(20261115)` and filled to 18–22% of patients.
  2. Feasibility constraints: every task's test positive fraction within ±5 pp of its overall fraction; every primary task with ≥15 positives and ≥15 negatives in test; every project's test share within 20 ± 8%. 18 of 20,000 draws were feasible.
  3. Starting from the best-scoring feasible draw, 3,000 local swap/add/remove moves (rng seed 20261116) were accepted only if they stayed feasible and lowered the score. Score = squared deviations of task positive fractions + per-stratum test-share deviations + per-task test-share deviation from 20%.
  4. Result: **29 TSSs, 228 patients (19.2%; 262 slides)**: COAD 79, READ 31, STAD 87, ESCA 31. Test TSSs: 3L, 4N, 4T, AD, AH, AY, AZ, B7, BM, CA, CG, CK, D7, DC, DY, G5, HJ, IG, L7, M9, QG, QL, RE, SS, T9, V5, WS, X8, ZR.
  5. Largest deviation of a test positive fraction from overall: 2.0 pp (PIK3CA STAD).
- **DEV 5-fold CV (960 patients, 48 TSSs).** sklearn `StratifiedGroupKFold` (v1.9.1) gave identical poor balance for all 300 seeds tried (worst fold deviation 17 pp), so it was replaced by a seeded assignment search.
  1. Greedy size-balanced initialisation, then 6,000 random move/swap steps × 20 restarts (rng seed 20261117).
  2. Score = squared deviation of fold positive fractions from the dev fraction (primary tasks weighted 3×) + per-task eligible share per fold + fold size + composite-stratum share.
  3. Fold sizes 173/176/213/184/214.
  4. Fold TSSs:
     - Fold 0: AF, DM, DT, EQ, G4, HF, KH, L5, VQ, XP, Z6
     - Fold 1: 2H, A6, AM, CD, D5, EI, MX, RD, S8, VR
     - Fold 2: AA, EF, F5, FP, IP, Q9, R5, SW, ZA, ZQ
     - Fold 3: 3M, 5M, AG, AU, CI, F1, HU, IN, KB, LN, R6
     - Fold 4: BR, CL, CM, F4, JY, NH
  5. Maximum fold deviation from the dev positive fraction: MSI 3.8 pp, KRAS 5.6 pp, TP53 STAD 8.7 pp, PIK3CA COADREAD 3.3 pp, APC 5.8 pp, **TP53 COADREAD 14.6 pp**. Every primary task has ≥10 patients in its smaller arm in every fold.
- **Known structural imbalances (accepted, and the reason for the TSS-only baseline):**
  - **TSS AA** (173 patients, all COAD, fold 2) has a TP53 mutation rate of 44.2% among eligible patients (n = 95), against 70.4% in all other CRC sites (Fisher p = 5×10⁻⁶). APC is also lower (68.4% vs 79.0%, p = 0.039). Median mutation counts are similar, and the cause is not established. Pre-specified sensitivity analysis: TP53 COADREAD and APC COADREAD re-evaluated excluding TSS AA from the pooled out-of-fold predictions.
  - **TSS BR** (136 STAD patients, fold 4) carries 15 of the 23 dev CDH1-mutant STAD patients. Fold 0 has 0 CDH1 positives, and fold 1 has 1 PIK3CA-STAD positive. A re-optimisation with a minimum-positives penalty could not remove this without degrading primary-task balance. Exploratory STAD tasks are therefore evaluated only by **pooled out-of-fold AUROC**; per-fold AUROC is not reported for them.
- **DISCOVERY subset** (`is_discovery`, `discovery_group`; dev only; 105 patients; 89.5% scanned at 40×; QC-pass only; drawn per fold with rng seed 20261118 and 40× preferred):
  - 30 MSI-H: 15 COAD + 15 STAD.
  - 30 MSS: 15 COAD + 15 STAD.
  - 30 KRAS-mutant and 30 KRAS/BRAF-WT (MSS COAD+READ).
  - The 15 MSS COAD patients are a subset of the KRAS discovery patients (8 mutant, 7 WT), which limits training loss.
  - **Rule:** discovery patients serve only morphology hypothesis generation and rubric development (LLM rubric drafting, tile review). They are **excluded from training and evaluation of the fused (encoder + rubric) model** in every task. Encoder-only ABMIL models may train on them, and every dev metric is also reported with them excluded (`dev_excl_discovery` rows).

## 5. Features and models

- **Encoders (frozen, no fine-tuning):** UNI2-h (primary) and Phikon (comparison). Same tissue mask, tiling (fixed µm/px target, documented at extraction), and tile set for both encoders.
- **Model:** gated attention-based MIL (gated ABMIL, Ilse et al. 2018), one model per task, trained on the 5 site-grouped dev folds from this file. Hyperparameters are fixed a priori or tuned only by inner splits within training folds. Class imbalance is handled by class weighting.
- **Fusion:** encoder ABMIL score + frozen blinded LLM morphology rubric scores, combined by a logistic regression fitted out-of-fold on dev, excluding discovery patients.

## 6. Metrics

- **Primary metric:** AUROC on pooled out-of-fold dev predictions, with a 95% CI from **2,000 patient-level bootstrap draws**; the same for the locked test set.
- **Calibration:** reliability diagram, calibration slope and intercept, Brier score.
- **Secondary:** AUPRC, sensitivity at 90% specificity (MSI), per-cohort AUROC for MSI (COADREAD vs STAD).
- **Multiplicity:** Benjamini–Hochberg FDR across all gene tasks: KRAS, TP53 ×2, PIK3CA ×2, APC, SMAD4 COADREAD, CDH1 STAD; the BH family includes the exploratory tasks. The null is AUROC = 0.5, tested by permutation or bootstrap p-value on the locked test set. MSI is reported separately, not in the gene family.
- **Positive control:** CDH1 STAD (diffuse / signet-ring morphology). If CDH1 fails to beat its baselines, the pipeline's sensitivity is questioned before negative gene results are interpreted.

## 7. Mandatory confounder baselines

Reported alongside every task, on the same folds and the same patients:
1. **Project-only** (cancer type).
2. **TSS-only** (site identity, target-encoded within training folds; new sites in a held-out fold fall back to the training prior).
3. **MSI-only** (MANTIS score) for every gene task. This is relevant even inside the MSS-restricted cohorts, as a borderline-MSI check.
4. **Stage-only** (AJCC pathologic stage).
5. **LLM-scored differentiation grade** (blinded rubric score) for KRAS.
6. Technical: scanner objective power (20× vs 40×) and series size, which should sit near AUROC 0.5.

A gene-task model is reported as adding information only if its AUROC CI excludes the best single baseline's point estimate, or it improves on that baseline in a nested comparison (DeLong or paired bootstrap).

## 8. Locked test and external validation (single look)

- **The locked TEST partition (228 patients) and CPTAC-COAD are evaluated exactly once, at the end**, after every model, fusion weight, threshold and rubric is frozen and recorded with a hash. No test-set metric is computed beforehand, including for debugging.
- **CPTAC-COAD external cohort:** see the CPTAC manifest (built separately). Its counts and label definitions are recorded there. The MSI and KRAS label definitions must match Section 2 as closely as the available data allow, and any mismatch is logged.

## 9. Reproducibility

Seeds: test search 20261115, test refinement 20261116, dev folds 20261117, discovery per-fold draws 20261118. The MSS-COAD discovery subsample from the KRAS discovery set and the one MSI-H COAD top-up use pandas `sample(random_state=20261115)`. The top-up was needed because fold 3 has only 2 dev MSI-H COAD patients. Environment: pandas, numpy and scikit-learn 1.9.1. The split file is the authoritative record; regenerating it must give identical partitions.


---

## Deviations (logged after locking)

**D1 — 2026-11-15/16 — Blinded reader and grading studies looked at locked-TEST slides (Section 8).**
Exploratory LLM-reader studies were run before this pre-registration was checked, and they drew on all colorectal MSS patients, regardless of partition:
- KRAS 6-pair pilot (TSS-matched): 4 of 6 pairs come from test TSSs (AY, AZ, AD, CA).
- MSI two-slide demo: the MSS slide TCGA-4T-AA8H is from test TSS 4T.
- Single-reader KRAS study (`kras_full_reader_results.csv`): 22 of 116 pairs contain a test slide.
- Three-reader PIK3CA-matched KRAS study (`kras3_reader_results.csv`): 16 of 63 pairs contain a test slide.
- Blinded per-slide grading (`kras_full_slide_grades.csv`): 44 of 232 slides are test slides.
- Pathologist spot-check packet: 5 of 30 slides are test slides.
No encoder model was trained on, or produced predictions for, any test slide. These reader/grader outputs are not the pre-registered model, but under a strict reading of Section 8 they are test-set predictions.
Handling: every reported reader or grader number is restated on dev-only pairs/slides (single-reader 61/94 = 64.9%; three-reader majority vote 23/47 = 48.9%). Test-slide rows are quarantined from all analyses and from rubric development. The LLM morphology rubric (Section 5, fusion) must not be tuned on any of these outputs. Resolution (investigator decision, 2026-11-16): the locked test set **remains valid for encoder-only (ABMIL) models**, because no encoder model was trained on, tuned on, or produced predictions for any test slide. For the **fused encoder + LLM-rubric model**: the rubric is developed **only on the discovery subset**; no existing reader or grader output for a test slide is reused; and the test slides are re-graded fresh with the frozen rubric at the single final look (Section 8). The pathologist spot-check packet is rebuilt with dev-partition slides only (the 5 test slides are swapped for dev replacements under the same stratification).

**D2 — 2026-11-15 — Primary encoder.** UNI2-h (Section 5) was replaced by **UNI (v1, vit_large_patch16_224)**, and Phikon is kept as the comparison. Reason: UNI2-h embedding of the cohort was projected at about $43 on T4, over the $40 compute cap. UNI v1 and Phikon features come from one tiling pass with identical tile coordinates (embed_stream2.py: 224 px at 0.5 µm/px, at most 4,000 tiles per slide, seed 0).

**D3 — 2026-11-16 — Count check.** Among the 443 colorectal MSS sequenced slides, the split assigns 92 slides / 87 patients to test (85 slides carry the KRAS task label). All are excluded from training, OOF prediction and pair scoring.


**D4 — 2026-11-16 — Prior leakage into the three-reader KRAS study.** The platform automatically injects relevant project memory into each new sub-agent. All 189 three-reader KRAS readers (`kras3_reader_results.csv`) therefore received a memory entry summarising the earlier full-cohort study: mucinous component higher in KRAS-mutant tumours, serrated architecture higher, dirty necrosis lower, and the PIK3CA co-mutation imbalance. That study is a *primed* replication, not the no-hint design it was meant to be. The 116 single-reader study readers and the 8 graders received no KRAS morphology associations (checked transcript by transcript; one grader saw the MSI demo summary). No memory entry contained pair identities or labels. Rule going forward: before reporting any blinded LLM-reader study, audit the injected recall in every reader transcript, and keep result summaries out of the memory that blinded readers can be injected with.

**D5 — 2026-11-16 — CPTAC-COAD opened for two encoder-only tasks before the final look (Section 8).**
CPTAC-COAD was evaluated for MSI (mean-pooled features + logistic regression; results in msi_performance.csv) and then for colorectal TP53, each once. No fused or rubric model has been evaluated on CPTAC-COAD, and no CPTAC label was used for training or tuning.
TP53 look: specification frozen before any CPTAC TP53 label was read (tp53_external_spec.json, hash 95a60aef1afefcbc). A first draft specification (C=1.0, label_tp53_coadread) failed to reproduce the published development result and was replaced, before the look, by the recipe already fixed in analyse_mut.py (StandardScaler + LogisticRegression C=0.01, mut_labels.csv crc_TP53, 324 dev patients; dev CV reproduced at Phikon 0.762 / UNI 0.765).
Result, 80 patients (47 mutant / 33 wild-type), patient-level AUROC: UNI 0.592 [0.465-0.717], Phikon 0.642 [0.513-0.770]. Restricted to the 71 patients with GDC-concordant TP53 calls: UNI 0.641 [0.496-0.780], Phikon 0.703 [0.555-0.829]. The model evaluated is the linear head, not attention-MIL. These numbers stand as the external TP53 result and are not to be re-run with a different recipe.

**D6 — 2026-11-16/17 — Second external cohort for MSI: SurGen SR386 (not in the original plan).**
The MSI model evaluated on CPTAC-COAD under D5 was frozen (msi_model_frozen.npz) and applied once to SurGen (BioStudies S-BIAD1285), cohort SR386: all 32 mismatch-repair-deficient cases plus 96 proficient cases drawn at random with a recorded seed (surgen_msi_cohort.csv). Specification frozen before any SurGen feature existed (surgen_msi_spec.json, hash c678b34ef0fb3b94); the frozen model reproduced the CPTAC AUROCs exactly before the look.
Two proficient slides (SR386_40X_HE_T264_01, SR386_40X_HE_T353_01) were lost to a storage fault (4 KB stub files on the volume) and excluded as the specification prescribes, leaving 126 slides (32 / 94).
Result, AUROC: UNI 0.665 [0.562-0.765], Phikon 0.760 [0.638-0.862]; tissue tile count alone 0.512 [0.384-0.641]; UNI minus Phikon -0.095 [-0.223-0.034].
Differences from the TCGA setting: the label is immunohistochemistry (not sequencing-based MSI), slides are Zeiss CZI at 0.1112 um/px read with a different library at 2x oversampling, and the cohort is a case-control sample. These numbers stand as the SurGen MSI result and are not to be re-run with a different recipe.

**D7 — 2026-11-17 — SurGen SR386 opened for colorectal KRAS (microsatellite-stable only).**
The KRAS arm of the SurGen slides already embedded under D6 was scored once. No new slides were downloaded or embedded; the 126 D6 slides were restricted to those that are mismatch-repair-proficient by IHC and carry a determinate KRAS call. Recipe is the one fixed in analyse_mut.py (StandardScaler + LogisticRegression C=0.01, class_weight balanced, mut_labels.csv crc_KRAS, 324 TCGA dev slides, 136 mutant / 188 wild-type); its development CV was reproduced exactly before the look (Phikon 0.661, UNI 0.632, matching mut_performance.csv).

**Label-rule correction, same day, before any further analysis.** The first frozen specification (hash 418f482d36950b53) defined wild-type as kras_ex_2 and kras_ex_3 both 'WT' and did not consult kras_codon_117. Three proficient cases (327, 428, 551) carry kras_codon_117 == 'FAIL' with no mutation on the other two assays: their KRAS status is indeterminate, not wild-type, and they were wrongly retained in the wild-type arm (n was 28 / 65). The rule was corrected to require all three assays determinate — mutant if any begins with 'M', wild-type only if all three read 'WT', otherwise excluded — and the look was recomputed under spec hash 28556c6fc25e5c0d on 90 slides (28 mutant / 62 wild-type). Cases 239 and 459 also carry a codon-117 FAIL but are mutant on exon 2, so their status is established and they remain. This is a correction of a label-definition error, not a change of recipe or model; the superseded numbers were Phikon 0.578 [0.434-0.705] and UNI 0.583 [0.437-0.712], and they are recorded here so the correction is auditable.

**Result (authoritative), slide-level AUROC on 90 slides:** Phikon 0.570 [0.431-0.692], UNI 0.579 [0.455-0.713]. Tissue tile count alone 0.389 [0.258-0.530]. Both encoder intervals include 0.5.
Interpretation: colorectal KRAS-from-H&E, already weak within TCGA (dev CV 0.63-0.66), is indistinguishable from chance in an independent cohort once microsatellite-unstable cases are removed. The arm is small (28 mutant), so the interval is wide by construction and this result bounds rather than estimates the effect.
The same differences from the TCGA setting as D6 apply (IHC-based mismatch-repair status, Zeiss CZI at 0.1112 um/px, case-control sampling). These numbers stand as the SurGen KRAS result and are not to be re-run with a different recipe.
