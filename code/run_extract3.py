"""Runs extract_kras_full.one() (unchanged) with N_CAND=48 over every slide of the 3-reader pair set.
Slide directory = TCGA slide barcode (no arm in names). Seeds follow the original scheme (1000+pair mut, 5000+pair wt).
Usage: python run_extract3.py [limit]"""
import os, sys, json, time
import pandas as pd
sys.path.append(os.getcwd())
import extract_kras_full as E
E.N_CAND = 48
if __name__ == '__main__':
    P = pd.read_csv('kras3_pairs_UNBLINDED.csv'); limit = int(sys.argv[1]) if len(sys.argv) > 1 else None
    jobs = []
    for r in P.itertuples():
        jobs.append((r.mut_slide_barcode, r.mut_crdc_series_uuid, 1000 + r.pair))
        jobs.append((r.wt_slide_barcode, r.wt_crdc_series_uuid, 5000 + r.pair))
    if limit: jobs = jobs[:limit]
    os.makedirs(E.OUT, exist_ok=True); os.makedirs(E.TMP, exist_ok=True)
    from multiprocessing import Pool
    done = 0; t0 = time.time()
    with Pool(3, maxtasksperchild=4) as p, open('kfull/extract_log.jsonl', 'a') as log:
        for m in p.imap_unordered(E.one, jobs):
            done += 1
            log.write(json.dumps({k: v for k, v in m.items() if k != 'fields'}) + '\n'); log.flush()
            if done % 10 == 0 or not m.get('ok'):
                print(f"{done}/{len(jobs)} {time.time()-t0:.0f}s | {m['slide']} {'OK' if m.get('ok') else m.get('error')}", flush=True)
    print('DONE', done, 'in', round(time.time() - t0), 's', flush=True)
