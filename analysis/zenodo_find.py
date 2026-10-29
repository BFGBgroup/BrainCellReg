"""用 Zenodo 官方接口找单核 eQTL 的公开汇总数据（只查询、只打印）。"""
import json, sys, urllib.request, urllib.parse

QUERIES = sys.argv[1:] or [
    'Bryois eQTL brain cell types',
    'cell-type-specific cis-eQTLs eight human brain cell types',
    'SingleBrain single-nucleus eQTL meta-analysis',
    'single-nucleus eQTL brain summary statistics',
]
seen = set()
for q in QUERIES:
    url = 'https://zenodo.org/api/records?' + urllib.parse.urlencode({'q': q, 'size': 8, 'sort': 'bestmatch'})
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'curl/8'}), timeout=60) as r:
            d = json.loads(r.read().decode('utf-8', 'replace'))
    except Exception as e:
        print('ERR', q, e); continue
    print('\n### query:', q, '| hits', d.get('hits', {}).get('total'))
    for h in d.get('hits', {}).get('hits', []):
        rid = h.get('id')
        if rid in seen:
            continue
        seen.add(rid)
        md = h.get('metadata', {})
        files = h.get('files', []) or []
        tot = sum(f.get('size', 0) for f in files)
        creators = ', '.join(c.get('name', '') for c in md.get('creators', [])[:3])
        print(f"- id {rid} | {md.get('publication_date')} | {md.get('title', '')[:110]} | {creators} | doi {h.get('doi')} | files {len(files)} | {tot / 1e9:.2f} GB")
        for f in files[:14]:
            print(f"    {f.get('size', 0) / 1e6:10.1f} MB  {f.get('key')}")
