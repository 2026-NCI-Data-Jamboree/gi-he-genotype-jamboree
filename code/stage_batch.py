import csv, json, os, subprocess, time, glob, sys
from concurrent.futures import ThreadPoolExecutor
rows = list(csv.DictReader(open(sys.argv[1])))
def dl(r):
    dst = f"/data/slides/{r['slide_barcode']}"; os.makedirs(dst, exist_ok=True); t = time.time()
    if len(glob.glob(dst + '/*.dcm')) == 0:
        p = subprocess.run(['s5cmd', '--no-sign-request', '--endpoint-url', 'https://s3.amazonaws.com', 'cp', '--show-progress=false',
                            r['s3_url'].rstrip('/') + '/*', dst + '/'], capture_output=True, text=True); rc, err = p.returncode, p.stderr[-200:]
    else: rc, err = 0, 'already present'
    fs = glob.glob(dst + '/*.dcm')
    return {'slide': r['slide_barcode'], 'n_files': len(fs), 'GB': round(sum(map(os.path.getsize, fs))/1e9, 3), 'dl_s': round(time.time()-t, 1), 'rc': rc, 'err': err}
t0 = time.time()
with ThreadPoolExecutor(4) as ex:
    for rec in ex.map(dl, rows): print(json.dumps(rec), flush=True)
print('TOTAL_S', round(time.time()-t0, 1), flush=True)
