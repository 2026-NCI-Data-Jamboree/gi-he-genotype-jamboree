# CPTAC-COAD external validation manifest — notes (IDC v25, built Nov 15 2026)

## Slides
- IDC collection `cptac_coad`, Modality SM: **372 series / 178 patients / 68.6 GB**. All are H&E, Aperio (converted by PixelMed), 40x objective, base pixel spacing 0.25 µm/px.
- Download: `s3://idc-open-data/<crdc_series_uuid>/` (column `s3_url`).
- **Embedding medium and fixative are not recorded** for any series (sm_index, sm_instance_index and DICOM SpecimenPreparationSequence contain only the staining step). So FFPE vs frozen **cannot be verified from metadata**: frozen 0 / FFPE 0 / unknown 372.
  - The CPTAC pan-cancer histology papers describe CPTAC H&E images as FFPE "mirror-image" sections. That description is documented for the CPTAC-3 cancer types (CCRCC, HNSCC, LSCC, LUAD, PDA, UCEC), **not specifically for this CPTAC-2 prospective colon cohort**.
  - The *molecular* specimens are frozen (GDC preservation_method = Frozen; cBioPortal SPECIMEN_PRESERVATION = Frozen Tissue).
  - Morphology in the few slides inspected looked paraffin-quality, but this is not a formal frozen/FFPE call.
  - Treat CPTAC as a **potential domain shift** relative to TCGA FFPE DX slides (different institutions, Aperio software versions v12.0.15 on 347 slides and vFS90 on 25, and a different tissue-procurement protocol), with FFPE status presumed rather than confirmed.

## Tumour vs normal — the DICOM tag is unreliable
- DICOM SpecimenDetailedDescription / SeriesDescription: 82 `tumor_tissue`, 114 `normal_tissue`, 176 untagged (identical to the sm_index modifier; untagged slides carry no specimen metadata at all).
- The DICOM ContainerIdentifier equals the GDC CPTAC-2 **sample submitter_id** (truncated UUID). This join covered 108/372 slides, and all 108 match a GDC *Primary Tumor* sample of the same case (0 case mismatches). GDC holds only tumour and blood samples for these cases.
- **Conflict:** 54 slides tagged `normal_tissue` in DICOM have a container ID equal to the profiled *Primary Tumor* sample. Thumbnail and 1 µm/px review of 11 sampled slides (see QC figure) showed tumour-like glandular architecture on most `normal_tissue`-tagged slides, sometimes alongside normal crypts. This is a visual impression, not a pathology review.
- **Eligibility rule used (strict):** `eligible_tumour_slide_strict` = GDC tumour-sample match OR DICOM `tumor_tissue` tag. This gives **148 slides / 108 patients**, all in cBioPortal.
  - Breakdown: 42 slides are both GDC-matched and DICOM-tumour, 54 GDC-matched but DICOM-normal, 12 GDC-matched and untagged, 40 DICOM-tumour only.
  - The other 224 slides (164 untagged, 60 normal-tag-only) have unknown tumour content. They are flagged `eligible_slide_permissive` only and must pass a tumour-tile QC before use.
  - **Warning for the embedding/analysis plan:** "filter out normal slides" by DICOM tissue_type would wrongly discard 54 slides from the profiled tumour specimen. Use `eligible_tumour_slide_strict` instead.
- The 68 IDC patients without cBioPortal data (139 slides) have no GDC-matched or tumour-tagged slide and no labels.

## ID mapping
- IDC PatientID (e.g. `01CO005`) = cBioPortal `coad_cptac_2019` patientId = sampleId = GDC case submitter_id. All **110/110** cBioPortal patients are in IDC; 68 IDC patients are unmatched; 0 cBioPortal patients are missing from IDC.
- `coad_cptac_gdc` sampleIds are GDC sample submitter_ids: 109 samples, 105 patients, 107 of them equal an IDC ContainerIdentifier (`in_cbio_gdc_study`).

## Labels
- MSI: `coad_cptac_2019` MSI_STATUS. 105 called (24 MSI-H / 81 MSS), 5 not called; one uncalled patient has TMB 315/Mb (ultramutated, excluded). MUTATION_PHENOTYPE is identical to MSI_STATUS.
- Hypermutation: TMB_NONSYNONYMOUS ≥10/Mb. MSS range 2.2–9.4, so no MSS tumour is excluded; all MSI-H tumours are ≥23.6.
- Mutations: `coad_cptac_2019_mutations`, sample list `coad_cptac_2019_sequenced` (106 patients). Gene = mutant if any protein-altering variant (missense, nonsense, frameshift, in-frame, splice-site); Splice_Region is excluded (1 CDH1 call).
  - KRAS: 35 mutant overall, 33 at hotspots (codons 12/13/61/117/146); `KRAS_hotspot` is given for sensitivity analyses.
  - KRAS-WT arm = KRAS-WT and BRAF-WT (BRAF 16 mutant, 13 V600E).
- Cross-check against `coad_cptac_gdc` (GDC harmonised calls, 104 shared sequenced patients). Discordant patients: KRAS 1, BRAF 2, TP53 10, PIK3CA 5, APC 3, SMAD4 0, CDH1 1. Per-gene `<GENE>_gdc_concordant` flags are in the manifest.

## Externally testable tasks (strict tumour slide, ≥15 patients per arm)
- MSI-H vs MSS: 24 vs 80 — **testable**
- KRAS-mut vs KRAS/BRAF-WT (MSS, TMB<10): 34 vs 44 — **testable** (34 vs 42 if restricted to GDC-concordant calls)
- TP53 (MSS): 47 vs 33 — **testable** (45 vs 26 concordant-only)
- APC (MSS): 67 vs 13 — not testable (WT arm 13)
- PIK3CA (MSS): 12 vs 68 — not testable
- SMAD4 (MSS): 8 vs 72 — not testable
- CDH1 (MSS): 0 vs 80 — not testable
- The permissive (all-slides) tier changes these by at most 1 patient per arm and does not change any testability call.
- External sets are small. Hanley–McNeil 95% CI half-width at AUROC 0.75 / 0.85: MSI (24 vs 80) ±0.12 / ±0.10, KRAS (34 vs 44) ±0.11 / ±0.09, TP53 (47 vs 33) ±0.11 / ±0.08.
