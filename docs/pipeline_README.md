# Datasets 1 + 2: IDC FFPE diagnostic slides joined to MSI labels
Project: Reading the Genome Off the Slide (NCI Data Jamboree, Nov 2026). Built 2026-10-24.
Sources:
- IDC data release v25, via idc-index 0.13.4.
- GDC API, file-level `msi_status`.
- cBioPortal TCGA PanCancer Atlas 2018 clinical data. `MSI_SCORE_MANTIS` comes from Bonneville et al. 2017, doi:10.1200/PO.17.00073.

## Coverage
| project_id   |   patients |   slides |   labelled |   msi_h |   trainable |   size_GB |   pct_msi_h |
|:-------------|-----------:|---------:|-----------:|--------:|------------:|----------:|------------:|
| TCGA-COAD    |        451 |      459 |        442 |      88 |         429 |     290.1 |        19.9 |
| TCGA-ESCA    |        156 |      158 |        155 |       2 |         155 |     192.4 |         1.3 |
| TCGA-READ    |        165 |      166 |        162 |       8 |         162 |      80.4 |         4.9 |
| TCGA-STAD    |        416 |      442 |        416 |      78 |         416 |     366.5 |        18.8 |
| TOTAL        |       1188 |     1225 |       1175 |     176 |        1162 |     929.4 |        15   |

## Slide selection (dataset 1)
- Started from collections tcga_stad, tcga_coad, tcga_read and tcga_esca with Modality = SM (3,565 series).
- Kept slides with embedding medium = paraffin wax (formalin-fixed). Frozen sections (tissue freezing medium) were excluded.
- That filter matches the TCGA `-DX` barcode suffix exactly: 1,225 slides from 1,188 patients. 36 patients have more than one DX slide.
- All slides are H&E. Sample type 01 (primary tumour) covers 1224 slides and 06 (metastatic) covers 1.
- Objective power: 40x for 860 slides, 20x for 298, missing for 67.
- `gcs_url` points at the series folder in the public bucket `gs://idc-open-data/<crdc_series_uuid>/`. Read it in place from GCP; there is no need to copy it.
- On 2026-10-24, listing the bucket and a byte-range read (valid DICOM preamble) both worked.
- Size: 929.4 GB in decimal units (`series_size_MB` is decimal MB). The deck's "908 GB" came from dividing MB by 1024.

## MSI label (dataset 2)
- `msi_label` uses MANTIS first: MSI-H means a score of 0.4 or more, the Bonneville threshold. 1,108 patients have a MANTIS score.
- 67 patients without MANTIS fall back to GDC `msi_status`, where the WXS and WGS calls agree.
- 13 patients have no label (GDC WXS and WGS calls disagree, and there is no MANTIS score).
- MANTIS and GDC agree on 1,071 of 1,084 patients with single calls from both (98.8%). The 13 discordant patients are flagged `mantis_gdc_discordant` and set to `use_for_training = False`.
- Per-patient extras: MSIsensor score, PanCan molecular subtype, age, sex and stage. Join to slides on `case_submitter_id`, the 12-character TCGA patient barcode.
- MSI-H makes up 15.0% of labelled patients, so report AUROC/AUPRC rather than accuracy.

## Files
- idc_ffpe_dx_slide_manifest.csv / .parquet: one row per slide (1,225). Includes GCS/S3 URLs, DICOM UIDs, the IDC viewer link and the MSI label.
- patient_msi_labels.csv: one row per patient (1,188).

## Reading a slide in place on GCP (example only; not run in the build sandbox)
```python
# pip install gcsfs wsidicom
import gcsfs, pandas as pd
fs = gcsfs.GCSFileSystem(token="anon")
m = pd.read_parquet("idc_ffpe_dx_slide_manifest.parquet")
files = fs.ls(m.gcs_url.iloc[0].removeprefix("gs://"))   # DICOM instances (pyramid levels) of one slide
```
