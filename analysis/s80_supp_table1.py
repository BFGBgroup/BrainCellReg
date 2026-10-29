# Supplementary Table 1: donor characteristics per brain region, recomputed from GEO GSE60862 sample annotation
# for exactly the samples that entered deconvolution / eQTL mapping (estimated_frac files), to check and extend
# Heyijie's Han-Supplementary Table S1 (N, age at death mean (SD), sex ratio).
# python s80_supp_table1.py   (admin2; reads only the header block of the series matrix)
import gzip, re, statistics as st, sys, io, os, collections
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

GEO = '/public/home/yucong/scbrain/supp/geo/GSE60862_series_matrix.txt.gz'
FR = '/public/home/yucong/scbrain/local_handoff/scbrain2/BayesPrism_InstaPrism_results/estimated_frac'
OUT = '/public/home/yucong/scbrain/supp/table1'
REG = ['FCTX', 'TCTX', 'OCTX', 'HIPP', 'PUTM', 'THAL', 'SNIG', 'MEDU', 'WHMT', 'CRBL']
os.makedirs(OUT, exist_ok=True)

rows = collections.defaultdict(list)
with gzip.open(GEO, 'rt', encoding='utf-8', errors='replace') as f:
    for line in f:
        if line.startswith('!series_matrix_table_begin'):
            break
        if line.startswith('!Sample_characteristics_ch1') or line.startswith('!Sample_geo_accession'):
            parts = [p.strip().strip('"') for p in line.rstrip('\n').split('\t')]
            if parts[0] == '!Sample_geo_accession':
                rows['gsm'] = parts[1:]
                continue
            key = parts[1].split(':', 1)[0].strip().lower() if ':' in parts[1] else parts[1]
            rows[key] = [p.split(':', 1)[1].strip() if ':' in p else p for p in parts[1:]]
n = len(rows['gsm'])
print('GEO samples', n, '| fields', [k for k in rows if k != 'gsm'])
print('regions in GEO:', collections.Counter(rows['brain region']))
print('banks in GEO:', collections.Counter(rows['brain bank']))

REGMAP = {'cerebellar cortex': 'CRBL', 'frontal cortex': 'FCTX', 'hippocampus': 'HIPP', 'medulla': 'MEDU',
          'occipital cortex': 'OCTX', 'putamen': 'PUTM', 'substantia nigra': 'SNIG', 'temporal cortex': 'TCTX',
          'thalamus': 'THAL', 'white matter': 'WHMT'}
def reg_of(s):
    s = s.lower()
    hits = [v for k, v in REGMAP.items() if k in s]
    return hits[0] if len(hits) == 1 else None
def num(x):
    try:
        return float(x)
    except ValueError:
        return None
def norm_id(x):
    return re.sub(r'^X', '', x).replace('-', '_').replace('/', '_')   # GEO writes Edinburgh donors as 001/06

geo = {}
for i in range(n):
    r = reg_of(rows['brain region'][i])
    iid = norm_id(rows['individual id'][i])
    geo[(r, iid)] = dict(sex=rows['sex'][i], age=num(rows['age at death (in years)'][i]),
                         bank=rows['brain bank'][i], pmi=num(rows['post-mortem interval (in hours)'][i]),
                         rin=num(rows['rin'][i]), anc=rows['ancestry'][i])
unmapped = [rows['brain region'][i] for i in range(n) if reg_of(rows['brain region'][i]) is None]
print('unmapped region strings:', collections.Counter(unmapped))

def msd(v):
    v = [x for x in v if x is not None]
    return (st.mean(v), st.stdev(v), min(v), max(v), len(v)) if len(v) > 1 else None
def mediqr(v):
    v = sorted(x for x in v if x is not None)
    if not v:
        return None
    q = st.quantiles(v, n=4)
    return (st.median(v), q[0], q[2], len(v))

out = []
indiv = {}
for r in REG:
    with open(os.path.join(FR, f'{r}.estimated_frac.txt'), encoding='utf-8') as f:
        f.readline()
        ids = [norm_id(l.split()[0].strip('"')) for l in f if l.strip()]
    recs = [geo.get((r, i)) for i in ids]
    miss = [i for i, x in zip(ids, recs) if x is None]
    recs = [x for x in recs if x]
    for i, x in zip(ids, [geo.get((r, i)) for i in ids]):
        if x:
            indiv.setdefault(i, x)
    a = msd([x['age'] for x in recs]); p = mediqr([x['pmi'] for x in recs]); q = mediqr([x['rin'] for x in recs])
    fem = sum(x['sex'].lower() == 'female' for x in recs); mal = sum(x['sex'].lower() == 'male' for x in recs)
    banks = collections.Counter(x['bank'] for x in recs); anc = collections.Counter(x['anc'] for x in recs)
    out.append([r, len(ids), len(recs), len(miss), f'{a[0]:.1f} ({a[1]:.1f})', f'{a[2]:.0f}-{a[3]:.0f}', a[4], fem, mal,
                f'{fem / mal:.4f}', f'{p[0]:.1f} ({p[1]:.1f}-{p[2]:.1f})' if p else 'NA', p[3] if p else 0,
                f'{q[0]:.1f} ({q[1]:.1f}-{q[2]:.1f})' if q else 'NA', q[3] if q else 0,
                '; '.join(f'{k}={v}' for k, v in banks.items()), '; '.join(f'{k}={v}' for k, v in anc.items()),
                ','.join(miss[:5])])
recs = list(indiv.values())
a = msd([x['age'] for x in recs]); p = mediqr([x['pmi'] for x in recs]); q = mediqr([x['rin'] for x in recs])
fem = sum(x['sex'].lower() == 'female' for x in recs); mal = sum(x['sex'].lower() == 'male' for x in recs)
banks = collections.Counter(x['bank'] for x in recs)
out.append(['Donors', len(indiv), len(recs), 0, f'{a[0]:.1f} ({a[1]:.1f})', f'{a[2]:.0f}-{a[3]:.0f}', a[4], fem, mal, f'{fem / mal:.4f}',
            f'{p[0]:.1f} ({p[1]:.1f}-{p[2]:.1f})' if p else 'NA', p[3] if p else 0, 'per-sample', 0,
            '; '.join(f'{k}={v}' for k, v in banks.items()), '', ''])
hdr = ['region', 'n_samples', 'n_matched_geo', 'n_unmatched', 'age_mean_sd', 'age_range', 'age_n', 'female', 'male',
       'f_to_m', 'pmi_median_iqr', 'pmi_n', 'rin_median_iqr', 'rin_n', 'banks', 'ancestry', 'unmatched_ids']
with open(os.path.join(OUT, 'table1_recomputed.tsv'), 'w', encoding='utf-8') as f:
    f.write('\t'.join(hdr) + '\n')
    for o in out:
        f.write('\t'.join(str(x) for x in o) + '\n')
for o in out:
    print('\t'.join(str(x) for x in o))
