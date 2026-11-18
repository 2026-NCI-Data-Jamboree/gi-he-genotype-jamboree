# KRAS signal: PIK3CA entanglement and tumour-sidedness checks (dev partition only)
No new training or reader runs. Test-partition predictions were not read.

## 1. Slide model (ABMIL-Phikon, label_kras, n=334)
| Cell | n | mean OOF prob |
|---|---|---|
| KRAS only | 91 | 0.550 |
| KRAS + PIK3CA | 50 | 0.546 |
| PIK3CA only | 28 | 0.457 |
| Neither | 165 | 0.451 |

KRAS AUROC within PIK3CA-WT 0.658 [0.582-0.720]; within PIK3CA-mut 0.651 [0.521-0.776].
Site-clustered OLS on logit score: KRAS +0.468 (p=0.001), PIK3CA +0.031 (p=0.74), interaction -0.062 (p=0.71). **The model tracks KRAS, not PIK3CA.**

## 2. Readers (116-pair study, 94 dev pairs; blinded, no priors)
| Pair type (mut arm vs WT arm) | n | reader correct | Phikon correct |
|---|---|---|---|
| KRAS only vs neither | 50 | 27 (54%, p=0.67) | 35 (70%, p=0.007) |
| KRAS only vs PIK3CA only | 11 | 9 (82%) | 9 |
| KRAS+PIK3CA vs neither | 27 | 20 (74%, p=0.019) | 16 |
| KRAS+PIK3CA vs PIK3CA only | 6 | 5 | 3 |

Readers succeed mainly when the two slides differ in PIK3CA. On clean KRAS-only vs neither pairs they are at chance and the model is not (McNemar p=0.12, so the difference between them is not itself significant). Readers did not just pick the PIK3CA-mutant slide (22/38 discordant dev pairs). The primed 3-reader PIK3CA-matched study agrees: majority vote 18/34 in both-WT pairs, 5/13 in both-mut pairs. Phikon on the same PIK3CA-matched dev pairs scores 47/73 (WT, p=0.019) and 11/24 (mut).

## 3. Independent PIK3CA model (other chat's multi-gene ABMIL, fixed 20-epoch recipe, dev OOF n=321)
| Model | AUROC all | within KRAS-WT | within KRAS-mut |
|---|---|---|---|
| PIK3CA-Phikon | 0.507 [0.42-0.59] | 0.501 | 0.519 |
| PIK3CA-UNI | 0.497 [0.42-0.57] | 0.520 | 0.471 |

The same run's KRAS-Phikon head: AUROC 0.657 [0.595-0.712]; within PIK3CA-WT 0.678; within PIK3CA-mut 0.641. **PIK3CA is not learnable from these H&E features. The KRAS result reproduces with a second, untuned recipe.**

## 4. Tumour sidedness baseline (GDC site_of_resection_or_biopsy; right = caecum to transverse)
Site is known for 275/334 dev KRAS-cohort slides. KRAS is mutant in 55/94 right-sided vs 55/181 left-sided tumours (OR 3.23, Fisher p<1e-4).
- **Sidedness alone: AUROC 0.632 [0.574-0.688]**. This is stronger than the pre-registered best baseline (technical, 0.582).
- Phikon on the same 275: 0.659 [0.589-0.724]. Its score partly encodes side (AUROC 0.604 for right vs left).
- Phikon adjusted for side: OR 1.59 per SD (TSS-clustered p=0.005; LRT p=0.0006). Within left 0.626, within right 0.674.
- UNI adjusted for side: OR 1.36/SD (clustered p=0.18; LRT p=0.02).

Reader pairs by side (94 dev): **same-side pairs reader 27/39 (69%, p=0.024), Phikon 26/39**. KRAS-right/WT-left 12/17; KRAS-left/WT-right 3/6; unknown side 19/32.

## Interpretation
1. The slide model's KRAS signal is not a PIK3CA artefact. It holds within side and adds to sidedness, but about half its apparent margin over chance is shared with location. Sidedness should be reported as the relevant clinical baseline (post-hoc; not in the pre-registration, so it is logged as a deviation, not a change to the primary).
2. The readers' signal is concentrated in PIK3CA-discordant pairs. It is not explained by side (holds in same-side pairs). On clean KRAS-only vs neither pairs it is at chance. Their cues probably track a broader mucinous/serrated-type phenotype, not KRAS itself.
3. Do not launch a PIK3CA reader study. An H&E PIK3CA model with two encoders is null, so the prior for readers is low.
