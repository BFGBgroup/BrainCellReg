"""在超算登录节点从 EBI GWAS Catalog 的 FTP（http）下载指定 GCST 的汇总数据（优先 harmonised）。
用法: python gwas_fetch.py <out_dir> GCSTxxxx[:标签] ...
"""
import os, re, sys, subprocess, urllib.request

out = sys.argv[1]
BASE = 'http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics'


def listing(url):
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'curl/8'}), timeout=90) as r:
            html = r.read().decode('utf-8', 'replace')
    except Exception as e:
        return None
    return [m for m in re.findall(r'href="([^"?/][^"]*)"', html) if not m.startswith('/')]


for arg in sys.argv[2:]:
    gid, _, tag = arg.partition(':')
    num = gid[4:]
    n = int(num); lo = ((n - 1) // 1000) * 1000 + 1; hi = lo + 999
    bucket = f'GCST{str(lo).zfill(len(num))}-GCST{str(hi).zfill(len(num))}'
    root = f'{BASE}/{bucket}/{gid}'
    d = os.path.join(out, f'{gid}_{tag}' if tag else gid)
    os.makedirs(d, exist_ok=True)
    top = listing(root + '/')
    if top is None:
        print(gid, 'LIST_FAILED', root, flush=True); continue
    harm = listing(root + '/harmonised/') if any(x.rstrip('/') == 'harmonised' for x in top) else None
    todo = []
    if harm:
        todo = [(root + '/harmonised/' + f, f) for f in harm if f.endswith(('.h.tsv.gz', '.yaml', '.log', '.tbi')) or 'README' in f]
    main = [f for f in top if re.search(r'\.(tsv|txt|csv)(\.gz)?$|\.gz$|\.yaml$|\.zip$', f)]
    if not any(f.endswith('.h.tsv.gz') for _, f in todo):
        todo += [(root + '/' + f, f) for f in main]
    else:
        todo += [(root + '/' + f, f) for f in main if f.endswith('.yaml')]
    print(gid, tag, '| top:', ', '.join(top)[:200], '| harmonised:', (', '.join(harm)[:200] if harm else 'none'), flush=True)
    for url, f in todo:
        dst = os.path.join(d, f)
        subprocess.run(['wget', '-q', '-c', '--tries=8', '--timeout=90', '-O', dst, url])
        print('   ', f, os.path.getsize(dst) if os.path.exists(dst) else 'MISSING', flush=True)
print('FETCH_DONE', flush=True)
