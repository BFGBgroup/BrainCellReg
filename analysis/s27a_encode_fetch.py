"""Fetch ENCODE peaks for the prefrontal cortex (Fig 4e-f of the v2 figure plan, as in the reference manuscript's Figure 4c/e):
histone ChIP-seq for 7 marks and DNase-seq, GRCh38, released. Picks one peak file per experiment
(histone: replicated / pseudoreplicated peaks; DNase: IDR thresholded / conservative peaks, else any bed narrowPeak),
writes the manifest and downloads to scbrain/external/encode_pfc/. Runs on the login node (download only).
Usage: python s27a_encode_fetch.py [--list-only]"""
import os, sys, json, urllib.request, urllib.parse, hashlib

OUT = '/public/home/yucong/scbrain/external/encode_pfc'
os.makedirs(OUT, exist_ok=True)
MARKS = ['H3K27me3', 'H3K9me3', 'H3K27ac', 'H3K4me1', 'H3K4me3', 'H3K9ac', 'H3K36me3']
TISSUES = ['dorsolateral prefrontal cortex', 'prefrontal cortex', 'Brodmann (1909) area 9', 'Brodmann (1909) area 46']
BASE = 'https://www.encodeproject.org'


def get(url):
    req = urllib.request.Request(url, headers={'Accept': 'application/json'})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)


def experiments(assay, target=None):
    q = [('type', 'Experiment'), ('status', 'released'), ('assay_title', assay), ('replicates.library.biosample.donor.organism.scientific_name', 'Homo sapiens'),
         ('format', 'json'), ('limit', 'all')]
    for t in TISSUES:
        q.append(('biosample_ontology.term_name', t))
    if target:
        q.append(('target.label', target))
    return get(BASE + '/search/?' + urllib.parse.urlencode(q))['@graph']


def pick_file(exp):
    e = get(BASE + exp['@id'] + '?format=json')
    files = [f for f in e.get('files', []) if f.get('status') == 'released' and f.get('assembly') == 'GRCh38'
             and f.get('file_format') == 'bed']
    order = ['replicated peaks', 'pseudoreplicated peaks', 'conservative IDR thresholded peaks', 'IDR thresholded peaks', 'peaks']
    for o in order:
        c = [f for f in files if f.get('output_type') == o]
        if c:
            return e, c[0]
    return e, None


rows = []
for assay, targets in (('Histone ChIP-seq', MARKS), ('DNase-seq', [None])):
    for t in targets:
        for x in experiments(assay, t):
            e, f = pick_file(x)
            bs = e.get('biosample_summary', '')
            donor = ';'.join(sorted({r['library']['biosample']['donor'].get('accession', '') for r in e.get('replicates', [])
                                     if isinstance(r.get('library', {}).get('biosample', {}).get('donor'), dict)}))
            rows.append(dict(assay=assay, feature=t or 'DNase', experiment=e['accession'], biosample=bs, donor=donor,
                             file=f['accession'] if f else '', output_type=f.get('output_type', '') if f else '',
                             md5=f.get('md5sum', '') if f else '', href=(BASE + f['href']) if f else ''))
            print(rows[-1]['feature'], rows[-1]['experiment'], rows[-1]['file'], rows[-1]['output_type'], bs[:70], flush=True)

with open(f'{OUT}/manifest.tsv', 'w') as fh:
    keys = list(rows[0].keys())
    fh.write('\t'.join(keys) + '\n')
    for r in rows:
        fh.write('\t'.join(str(r[k]) for k in keys) + '\n')
print('experiments', len(rows), '| with file', sum(1 for r in rows if r['file']))
if '--list-only' in sys.argv:
    sys.exit(0)

bad = 0
for r in rows:
    if not r['file']:
        continue
    dst = f"{OUT}/{r['feature']}__{r['experiment']}__{r['file']}.bed.gz"
    if not os.path.exists(dst):
        urllib.request.urlretrieve(r['href'], dst)
    md5 = hashlib.md5(open(dst, 'rb').read()).hexdigest()
    if md5 != r['md5']:
        bad += 1
        print('MD5 MISMATCH', dst)
print('downloaded; md5 mismatches', bad)
print('FETCH_DONE')
