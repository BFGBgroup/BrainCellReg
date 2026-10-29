"""在超算登录节点下载一个 Zenodo 记录的全部文件（可续传，下完核 md5）。
用法: python zenodo_download.py <record_id> <out_dir>
"""
import json, sys, os, hashlib, subprocess, time, urllib.request

rid, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)


def get(url):
    with urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'curl/8'}), timeout=120) as r:
        return json.loads(r.read().decode('utf-8', 'replace'))


rec = get(f'https://zenodo.org/api/records/{rid}')
files = rec.get('files', [])
json.dump(rec.get('metadata', {}), open(os.path.join(out, '_metadata.json'), 'w'), ensure_ascii=False, indent=1)
print('record', rid, rec.get('metadata', {}).get('title'), 'files', len(files), 'GB %.2f' % (sum(f['size'] for f in files) / 1e9), flush=True)


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''):
            h.update(b)
    return h.hexdigest()


bad = []
for i, f in enumerate(sorted(files, key=lambda x: x['key']), 1):
    key, size = f['key'], f['size']
    want = (f.get('checksum') or '').replace('md5:', '')
    dst = os.path.join(out, key)
    url = f.get('links', {}).get('self') or f'https://zenodo.org/api/records/{rid}/files/{key}/content'
    ok = False
    for attempt in range(6):
        if os.path.exists(dst) and os.path.getsize(dst) == size and (not want or md5(dst) == want):
            ok = True; break
        if os.path.exists(dst) and os.path.getsize(dst) > size:
            os.remove(dst)
        subprocess.run(['wget', '-q', '-c', '--tries=5', '--timeout=60', '-O', dst, url])
        time.sleep(2)
    if not ok:
        bad.append(key)
    if i % 10 == 0 or not ok:
        print(time.strftime('%H:%M:%S'), f'{i}/{len(files)}', key, 'OK' if ok else 'FAILED', flush=True)
print('DONE files', len(files), 'failed', len(bad), bad[:10], flush=True)
open(os.path.join(out, '_download.done' if not bad else '_download.failed'), 'w').write('\n'.join(bad))
