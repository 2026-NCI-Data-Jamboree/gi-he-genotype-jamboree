# Blinded MSI test: can a reader swarm separate one MSI-H slide from one MSS slide?

## Design
Two TCGA-COAD diagnostic slides, matched on project, scanner objective (40x) and stage (IIA):

| | case | MANTIS | MSIsensor | truth |
|---|---|---|---|---|
| Slide A | TCGA-A6-5665 | 1.18 | 46.8 | **MSI-H** |
| Slide B | TCGA-4T-AA8H | 0.31 | 0.17 | **MSS** |

Slide/slot assignment was drawn at random and sealed until after all calls were in.

36 fields per slide were sampled from cellular tissue (ranked by haematoxylin density on the 4 um/px level, spatially
separated), cut from the 0.25 um/px level and resampled to **0.505 um/px, 1024 x 1024 px = a 517 um field** (~20x).
Six readers each received a private, non-overlapping set of 6 fields per slide plus both ~8 um/px whole-slide overviews;
12 fields per reader, 72 fields total. Each scored six dMMR-associated features 0-3 per slide, then called A or B with
a confidence and a declared basis. No abstention allowed; lookups forbidden.

## Result
**6 of 6 readers called the MSI-H slide correctly** (one-sided binomial p = 0.016), mean confidence
75% (range 65-88).
All six declared **morphology**, not stain colour, as the basis.

Mean feature scores (0-3), averaged over the six readers:

| feature | MSI-H slide | MSS slide | delta |
|---|---|---|---|
| intraepithelial lymphocyte density | 2.50 | 0.83 | +1.67 |
| peritumoral / Crohn-like lymphoid reaction | 2.33 | 1.50 | +0.83 |
| mucinous or signet-ring | 0.67 | 1.50 | -0.83 |
| medullary / solid poorly-differentiated growth | 2.67 | 0.17 | +2.50 |
| intratumoural heterogeneity | 2.17 | 1.00 | +1.17 |
| gland formation | poor (6/6) | well (6/6) | - |

The separation ran almost entirely on **solid/medullary growth with poor gland formation** and **dense TILs** - the two
features with the largest and most consistent gaps. Mucin ran the wrong way and was provisional in four readers, who
scored it off the 8 um/px overview because no high-power field of the MSS slide sampled a mucin lake.

## Example (Reader 3's field set)
MSI-H slide (../../reader_study/fields/reader3_slideA.jpg) and MSS slide (../../reader_study/fields/reader3_slideB.jpg); overviews ../../reader_study/fields/overview_slideA.jpg / ../../reader_study/fields/overview_slideB.jpg.

## What this does and does not show
- n = 1 slide pair, chosen at the **extremes** of the MANTIS/MSIsensor distribution and matched on stage. This is the
  easiest possible instance of the task, not an accuracy estimate. A two-slide forced choice also lets readers use
  *relative* comparison, which a per-slide classifier will not have.
- The six readers are not independent of each other in the way human readers would be: same model, same prompt,
  same feature checklist. Agreement here is weaker evidence than 6 independent pathologists agreeing.
- Tile selection was cellularity-based, not tumour-restricted: two of the six MSS fields in one reader's set contained
  no tumour (normal mucosa; deep wall with tattoo pigment). A tumour-region filter is needed before this scales.
- The two slides differ measurably in colour (mean RGB {"A": [140.6, 80.2, 119.3], "B": [137.9, 86.2, 142.0]}); readers
  were told to ignore it and all declared morphology as the basis, but stain/scanner batch remains a confound in TCGA
  and is not controlled here.
- Nothing here is a validated diagnostic. MSI status for care is established by IHC for MMR proteins or PCR/NGS testing.

## Next step if this is worth pursuing
Run the same forced-choice protocol over ~20-30 randomly drawn MSI-H/MSS pairs spanning the full MANTIS range, with
tumour-restricted tile selection and stain normalisation, and report accuracy against the MANTIS label.
