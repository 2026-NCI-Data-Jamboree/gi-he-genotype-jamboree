# Visual deck — what is on each slide, and where the number comes from

File: NCI_Jamboree_ReportOut_GI_HE_molecular.pptx (14 slides, built on the organizers' REVISED template of 17 Nov; template slides 1, 2 and 5, their frame for the rest)

| # | Slide | What the audience sees | Source of every number |
|---|---|---|---|
| 1 | Project Description | Title + category from your Sept 8 lightning deck; author line (lead first, D. Pasnoor second, D. Castillo last as senior author), affiliations, roles — no table | Jamboree_LightningTalk_MolinaRodriguez_Sept8.pptx |
| 2 | Project Rationale and Overall Approach | 3 one-line statements: problem / gap / approach | preregistration.md |
| 3 | How a slide becomes a number | Six-box chain: slide -> tiles -> frozen encoder -> pooled vector -> logistic -> prediction | fig_pipeline.png panel B |
| 4 | The cohort, in five numbers | 1,396 slides · 928 patients · 3 commons · 2 external cohorts · 7 endpoints | feature_coverage.csv (1,396 rows); patient_split_locked.csv (519 CRC + 409 gastric) |
| 5 | Data and Tools Used | Template 5-column table, 5 rows | feature_coverage.csv, three_commons_table.csv, external_cohort_options.csv |
| 6 | Does MSI detection hold up? | 4 evaluations x 2 encoders with 95% intervals, against the Kather 2019 band | msi_performance.csv; surgen_msi_external.csv; msi_sota_prior_art.csv (Kather 0.77 TCGA-CRC / 0.84 DACHS) |
| 7 | A check that does not use the image | MSH2 and MSH6 protein, MSI-H vs MMR-intact, with p-values | msi_protein_validation.csv; three_commons_table.csv (94 patients) |
| 8 | What this would mean in clinic | 140 squares: 35 MSI-H, 18 ruled out, 87 still tested | msi_external_operating_points.csv (Phikon, 100% sensitivity) |
| 9 | Every endpoint we tested | 11 endpoint x evaluation rows, external cohorts shaded | msi_performance.csv, mut_performance.csv (attention-MIL), tp53_cptac_external.csv, surgen_kras_external.csv, surgen_msi_external.csv |
| 10 | High burden is mostly just MSI | Stacked bars + the drop when MSI-H is removed | tmb_cause_decomposition.csv; tmb_from_msi_score.csv |
| 11 | What travels to another hospital | MSI / TP53 / KRAS, inside TCGA -> outside TCGA | msi_performance.csv, mut_performance.csv (linear), tp53_cptac_external.csv, surgen_kras_external.csv |
| 12 | Feedback to the commons | 3 cards: data gaps / pain points / what we give back | commons_gap_catalogue.csv; pdc000109_mmr_coverage.csv; PDC API notes |
| 13 | What we suggest, where it goes, what is next | 3 cards: suggestions to the commons / deposition plan / next steps | commons_gap_catalogue.csv; external_cohort_options.csv (SurGen 301 remaining cases, CI half-width 0.13 -> 0.06); pdc000109_mmr_coverage.csv |
| 14 | In one sentence | Three closing lines | — |

## Why slide 13 exists
The template's own speaker notes on its results slide ask for three things its visible bullets do not: suggestion of improvements, documentation and tools deposition, and future directions. Slide 13 answers all three. Deposition is worded as planned, because no public repository exists yet.

## Speaker notes are in the file
Slides 9 and 11 carry the extra numbers (CDH1 0.84, Ding 2022 comparison, MSIntuit comparison, POLE/POLD1 counts) in the notes pane, not on the slide.

## Not on the slides (ready if asked)
- Gastric CDH1 beyond Lauren histotype: 0.844 / 0.826 on 241 slides, 0.896 / 0.867 in non-diffuse only, histotype alone 0.659 (cdh1_vs_lauren_histotype.csv)
- Blinded reader arm, multimodal AI readers: MSI 6/6 pairs, CDH1 14/19 (msi_blinded_reader_results.csv, cdh1_swarm_majority.csv)
- Nuclear-morphometry fusion did not transfer: 0.701->0.722 and 0.794->0.807, both spanning zero (colorectal_nuclear_fusion_results.csv)
- TP53 flat across variant classes (tp53_by_variant_class.csv); per-tile attribution maps (tilemap_*.png)
- Cloud cost: about $30 of GPU time for the whole project (Modal ledger) — cut from the slides, mention only if asked
