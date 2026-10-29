"""External replication of the SMR genes (spec: progress/分析规格_独立数据验证_v1.md, part 2).
For the 253 disease x gene pairs (SMR q<0.05 in the main analysis): same GWAS, same reference panel, same SMR settings
(--diff-freq 1, other options default), eQTL side replaced by BrainMeta v2 and the 13 GTEx v8 brain tissues.
  python s25_smr_external.py run      -> indep/smr_ext/<dataset>/<GWAS>.smr
  python s25_smr_external.py stats    -> indep/smr_ext_summary.tsv, smr_ext_pairs.tsv
"""
import os, sys, glob, subprocess
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd

B = '/public/home/yucong/scbrain'
SMR = '/public/home/yucong/scbrain_in/SMR_inputs/smr-1.3.1/smr'
BF = B + '/local_handoff/scbrain1/g1000_eur/g1000_eur'
GW = B + '/gwas10'
EXT = B + '/external/indep/unz'
OUT = B + '/indep/smr_ext'
os.makedirs(OUT, exist_ok=True)
NPROC = int(os.environ.get('NPROC', '16'))


def datasets():
    d = {}
    for f in glob.glob(f'{EXT}/**/*.besd', recursive=True):
        name = os.path.basename(f)[:-5]          # BrainMeta may come split by chromosome: one dataset per file is fine
        if 'BrainMeta' in f and 'BrainMeta' not in name:
            name = 'BrainMeta_' + name
        d[name] = f[:-5]
    return d


def run():
    ge = pd.read_csv(B + '/integrate10_fix/gene_evidence.tsv', sep='\t')
    ds = datasets()
    print('datasets', len(ds), sorted(ds))
    jobs = []
    for name, pre in ds.items():
        epi = pd.read_csv(pre + '.epi', sep=r'\s+', header=None, usecols=[1], names=['probe'])
        epi['gene'] = epi.probe.str.split('.').str[0]
        od = f'{OUT}/{name}'; os.makedirs(od, exist_ok=True)
        for g, x in ge.groupby('GWAS'):
            pr = epi[epi.gene.isin(x.gene)].probe
            if pr.empty:
                continue
            pf = f'{od}/{g}.probes'; pr.to_csv(pf, index=False, header=False)
            jobs.append([SMR, '--bfile', BF, '--gwas-summary', f'{GW}/{g}-smr.tsv', '--beqtl-summary', pre, '--extract-probe', pf,
                         '--out', f'{od}/{g}', '--diff-freq', '1', '--thread-num', '1'])
    print('SMR runs', len(jobs), flush=True)

    def go(c):
        r = subprocess.run(c, capture_output=True, text=True)
        ok = os.path.exists(c[c.index('--out') + 1] + '.smr')
        return ok, c[c.index('--out') + 1], r.stdout[-300:]
    with ThreadPoolExecutor(NPROC) as ex:
        res = list(ex.map(go, jobs))
    bad = [o for ok, o, _ in res if not ok]
    print('done', len(res) - len(bad), '| no .smr output (no probe passing the eQTL threshold, or error):', len(bad))
    for o in bad[:20]:
        print('  ', o)


def stats():
    ge = pd.read_csv(B + '/integrate10_fix/gene_evidence.tsv', sep='\t')
    tr = pd.read_csv(B + '/integrate10_fix/smr_coloc_triplets.tsv', sep='\t')
    all3 = tr[tr.heidi_ok_001 & tr.h4_08][['disease', 'Gene']].drop_duplicates().rename(columns={'Gene': 'gene'})
    all3['all3'] = True
    main = tr.loc[tr.groupby(['disease', 'Gene']).p_SMR.idxmin(), ['disease', 'Gene', 'b_SMR', 'p_SMR']].rename(columns={'Gene': 'gene', 'b_SMR': 'b_main', 'p_SMR': 'p_main'})
    rows = []
    for f in glob.glob(f'{OUT}/*/*.smr'):
        d = pd.read_csv(f, sep='\t')
        d['dataset'] = os.path.basename(os.path.dirname(f)); d['GWAS'] = os.path.basename(f)[:-4]
        rows.append(d)
    E = pd.concat(rows, ignore_index=True)
    E['gene'] = E.probeID.str.split('.').str[0]
    E = E.merge(ge[['disease', 'GWAS', 'gene', 'symbol']], on=['GWAS', 'gene'], how='inner')
    E = E.merge(main, on=['disease', 'gene'], how='left').merge(all3, on=['disease', 'gene'], how='left')
    E['all3'] = E.all3.eq(True)
    E['same_sign'] = np.sign(E.b_SMR) == np.sign(E.b_main)
    E['rep'] = (E.p_SMR < 0.05 / 32) & E.same_sign
    E.to_csv(OUT + '/../smr_ext_pairs.tsv', sep='\t', index=False)
    P = E.groupby(['disease', 'gene', 'symbol', 'all3']).agg(datasets_tested=('dataset', 'nunique'), datasets_rep=('rep', 'sum'),
                                                             best_p=('p_SMR', 'min'), brainmeta_p=('p_SMR', lambda v: v[E.loc[v.index, 'dataset'].str.contains('BrainMeta')].min()),
                                                             any_heidi_ok=('p_HEIDI', lambda v: bool((v >= 0.01).any()))).reset_index()
    P['replicated'] = P.datasets_rep > 0
    P = ge[['disease', 'gene', 'symbol']].merge(P, on=['disease', 'gene', 'symbol'], how='left')
    P['all3'] = P.all3.eq(True) | P.set_index(['disease', 'gene']).index.isin(all3.set_index(['disease', 'gene']).index)
    P['replicated'] = P.replicated.eq(True)
    P.to_csv(OUT + '/../smr_ext_summary.tsv', sep='\t', index=False)
    for lab, x in (('32 pairs (SMR + HEIDI + coloc)', P[P.all3]), ('221 pairs (SMR only)', P[~P.all3])):
        t = x.datasets_tested.fillna(0) > 0
        print(f'{lab}: {len(x)} | testable in >=1 external dataset {int(t.sum())} | replicated {int(x.replicated.sum())} '
              f'({100 * x.replicated.sum() / max(1, t.sum()):.0f}% of testable)')
    print(P[P.all3].sort_values(['disease', 'replicated'], ascending=[True, False])[
        ['disease', 'symbol', 'datasets_tested', 'datasets_rep', 'best_p', 'brainmeta_p', 'any_heidi_ok']].to_string(index=False))


if __name__ == '__main__':
    run() if sys.argv[1] == 'run' else stats()
