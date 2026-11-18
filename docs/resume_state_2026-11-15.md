# Resume state — TCGA/CPTAC GI slide embedding (paused 2026-11-15 ~08:50 UTC)

## Spend state
- All Modal containers terminated. Nothing is computing; compute billing is $0/h.
- Posted spend: $5.62 through 07:59 UTC (Modal's report lags ~1 h; the 08:00 hour was not posted at pause time).
- Only remaining charge is volume storage: 8.75 GiB on volume `idc-gi-slides`
  (features only; the `/slides` download cache is empty, 0 files).

## Work completed and stored (volume `idc-gi-slides`)
- `/data/features/uni_max4k/`     550 files, 1024-d, 3.66 GiB
- `/data/features/phikon_max4k/`  980 files, 768-d, 4.86 GiB
- Slides with BOTH encoders (the usable, tile-aligned set): **550**
- The 430 phikon-only files came from the retired script and are NOT tile-aligned with UNI;
  they will be overwritten when the dual run reaches those slides. Do not concatenate
  768+1024 features unless both files came from embed_stream2.py.

## Work remaining
- Dual-encoder cohort target: 1,396 slides (1,175 TCGA + 221 CPTAC) -> **846 remaining**
- KRAS pair top-up list: 113 slides, 19 done -> 94 remaining (94 of these overlap the cohort)
- Union of outstanding slides: **865**
- Estimated cost to finish: ~$13 (T4 + 4 cores + 16 GiB = $1.54/h per container;
  measured ~35 s/slide/container -> ~$0.015/slide). About 2 h wall clock on 4 containers.

## How to resume
1. Say "resume the embedding run" in the MSI Slide Swarm chat.
2. The run is resumable by design: embed_stream2.py skips a slide only when BOTH .h5 exist,
   so nothing already done is repeated. Only the slide in flight at pause time is redone.
3. Rebuild the remaining list from the volume (dual target minus slides with both files),
   quote the cost, then relaunch after approval.

## Do not
- Do not launch `embed_stream.py` (Phikon-only). It is retired: it emits no UNI features and,
  because the dual script requires both files, its output saves the dual run nothing.
  It was relaunched twice by a now-closed session and cost roughly $6 in duplicated work.
- Closing a chat does NOT stop its containers. Audit with modal.Sandbox.list() and read each
  sandbox's `claude-science-session` tag.

## Correction to an earlier statement
The dual run does NOT include all 372 CPTAC-COAD manifest slides. It includes exactly the
**221 CPTAC slides that carry an MSI label** (53 MSI-H, 168 MSS) - verified by set fingerprint
against the manifest after barcode normalisation (the manifest pads the last field with '08').
The 151 excluded CPTAC slides have no MSI call and are unusable for external validation.
All 140 slides flagged `eligible_MSI` are inside the included 221.
