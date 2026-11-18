"""Download SurGen CZI slides onto the mounted volume. Resumable: a slide with a final .czi of the expected size is skipped.
Usage: python download_surgen.py dlK.csv   (csv: slide,rel_path,size_gb)"""
import os, sys, csv, json, time, subprocess, urllib.request
BASE = 'https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/285/S-BIAD1285/Files/'
DST = '/data/surgen_czi'; os.makedirs(DST, exist_ok=True); os.makedirs('out', exist_ok=True)
log = []; T0 = time.time()
def remote_size(url):
    for _ in range(5):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, method='HEAD'), timeout=60) as r: return int(r.headers['Content-Length'])
        except Exception: time.sleep(5)
    return None
for r in csv.DictReader(open(sys.argv[1])):
    name = r['slide']; url = BASE + r['rel_path']; fin = f'{DST}/{name}.czi'; part = fin + '.part'
    rec = {'slide': name}; t = time.time(); want = remote_size(url); rec['remote_bytes'] = want
    try:
        if os.path.exists(fin) and want and os.path.getsize(fin) == want: rec['status'] = 'already present'
        else:
            for attempt in range(3):
                p = subprocess.run(['aria2c', '-x', '16', '-s', '16', '-k', '8M', '--max-tries=20', '--retry-wait=5', '--file-allocation=none', '--allow-overwrite=true',
                                    '--console-log-level=error', '--summary-interval=0', '-d', DST, '-o', os.path.basename(part), url], capture_output=True, text=True)
                if p.returncode == 0 and os.path.exists(part) and (want is None or os.path.getsize(part) == want): break
                rec['last_error'] = (p.stderr or p.stdout)[-200:]; time.sleep(10)
            else: raise RuntimeError('download failed after 3 attempts: ' + rec.get('last_error', ''))
            os.replace(part, fin); rec['status'] = 'downloaded'; rec['attempts'] = attempt + 1
        rec['bytes'] = os.path.getsize(fin); rec['s'] = round(time.time() - t, 1); rec['MBps'] = round(rec['bytes'] / 1e6 / max(rec['s'], 0.1), 1)
    except Exception as ex:
        rec['status'] = 'FAILED'; rec['error'] = str(ex)[:300]
        if os.path.exists(part): os.remove(part)
    log.append(rec); print(json.dumps(rec), flush=True); json.dump(log, open(f'out/{os.path.basename(sys.argv[1])}.log.json', 'w'))
subprocess.run(['sync']); time.sleep(20)
print(f'DONE {sum(x["status"] != "FAILED" for x in log)}/{len(log)} in {time.time() - T0:.0f}s', flush=True)
