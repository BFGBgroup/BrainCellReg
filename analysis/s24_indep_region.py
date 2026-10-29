"""Independent region-level replication in GTEx v8 brain tissues (spec: progress/分析规格_独立数据验证_v1.md, part 1).
  python s24_indep_region.py query <tissue>   -> indep/query/<tissue>.txt   (smr --query of our lead + control pairs)
  python s24_indep_region.py stats            -> indep/region_pi1.tsv, region_matrix.tsv, criteria.txt
"""
import os, sys, glob, subprocess
import numpy as np
import pandas as pd

B = '/public/home/yucong/scbrain'
SMR = '/public/home/yucong/scbrain_in/SMR_inputs/smr-1.3.1/smr'
EXT = B + '/external/indep'
OUT = B + '/indep'
os.makedirs(OUT + '/query', exist_ok=True)
MAIN = {'FCTX': 'Brain_Frontal_Cortex_BA9', 'PUTM': 'Brain_Putamen_basal_ganglia', 'HIPP': 'Brain_Hippocampus',
        'SNIG': 'Brain_Substantia_nigra', 'CRBL': 'Brain_Cerebellum'}
ALT = {'FCTX': 'Brain_Cortex', 'CRBL': 'Brain_Cerebellar_Hemisphere'}
REG5 = list(MAIN)
SEED = 20261015


def besd_prefix(t):
    f = glob.glob(f'{EXT}/unz/{t}/**/*.besd', recursive=True)
    if not f:
        raise SystemExit(f'no besd for {t}')
    return f[0][:-5]


def pairs():
    lead = pd.read_csv(B + '/sn_replication/pairs_lead.tsv.gz', sep='\t')
    ctrl = pd.read_csv(B + '/sn_replication/pairs_ctrl.tsv.gz', sep='\t')
    lead['kind'] = 'lead'; ctrl['kind'] = 'control'
    return pd.concat([lead, ctrl], ignore_index=True)


def query(t):
    pre = besd_prefix(t)
    epi = pd.read_csv(pre + '.epi', sep=r'\s+', header=None, usecols=[1], names=['probe'])
    epi['gene'] = epi.probe.str.split('.').str[0]
    P = pairs()
    pr = epi[epi.gene.isin(P.gene)].probe.drop_duplicates()
    sn = P.snp.drop_duplicates()
    pr.to_csv(f'{OUT}/query/{t}.probes', index=False, header=False)
    sn.to_csv(f'{OUT}/query/{t}.snps', index=False, header=False)
    cmd = [SMR, '--beqtl-summary', pre, '--extract-probe', f'{OUT}/query/{t}.probes', '--extract-snp', f'{OUT}/query/{t}.snps',
           '--query', '1', '--out', f'{OUT}/query/{t}']
    print(' '.join(cmd), flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True)
    print(r.stdout[-1500:], r.stderr[-1500:])
    q = pd.read_csv(f'{OUT}/query/{t}.txt', sep=r'\s+')
    q['gene'] = q.Probe.str.split('.').str[0]
    m = P.merge(q[['SNP', 'gene', 'A1', 'A2', 'b', 'p']].rename(columns={'SNP': 'snp', 'A1': 'g_a1', 'A2': 'g_a2', 'b': 'g_b', 'p': 'g_p'}),
                on=['snp', 'gene'], how='inner')
    m.to_csv(f'{OUT}/query/{t}.pairs.tsv.gz', sep='\t', index=False)
    print(t, '| probes', len(pr), '| snps', len(sn), '| query rows', len(q), '| our pairs found', m.groupby('kind').size().to_dict())


def pi1(p):
    p = np.asarray(p, float); p = p[~np.isnan(p)]
    if len(p) < 30:
        return np.nan
    return float(min(1.0, max(0.0, 1 - np.mean(p > 0.5) / 0.5)))


def stats():
    at = pd.read_csv(B + '/allele_fix/allele_table_v2.tsv.gz', sep='\t', usecols=['snp', 'eff', 'oth', 'status'])
    at = at[at.status.isin(['ok', 'ok_ld'])].drop_duplicates('snp')
    tissues = sorted(set(MAIN.values()) | set(ALT.values()))
    Q = {}
    for t in tissues:
        f = f'{OUT}/query/{t}.pairs.tsv.gz'
        if os.path.exists(f):
            d = pd.read_csv(f, sep='\t')
            d = d.merge(at, on='snp', how='left')
            # direction: our slope refers to the counted allele (eff); GTEx b refers to g_a1
            ok = d.eff.notna() & (((d.g_a1 == d.eff) & (d.g_a2 == d.oth)) | ((d.g_a1 == d.oth) & (d.g_a2 == d.eff)))
            d['g_b_al'] = np.where(d.g_a1 == d.eff, d.g_b, -d.g_b)
            d.loc[~ok, 'g_b_al'] = np.nan
            Q[t] = d
    rows = []
    for kind in ('lead', 'control'):
        for r in REG5:
            for t in tissues:
                if t not in Q:
                    continue
                d = Q[t][(Q[t].kind == kind) & (Q[t].region == r)].drop_duplicates(['gene', 'snp'])
                s = d[(d.g_p < 0.05) & d.g_b_al.notna()]
                rows.append(dict(kind=kind, region=r, tissue=t, n=len(d), pi1=pi1(d.g_p),
                                 n_p05=len(s), concordance=float((np.sign(s.g_b_al) == np.sign(s.slope)).mean()) if len(s) else np.nan))
    R = pd.DataFrame(rows)
    R.to_csv(OUT + '/region_pi1.tsv', sep='\t', index=False)

    def S_stat(lead, main):
        M = pd.DataFrame({r: {t: pi1(lead[(lead.region == r) & (lead.tissue == t)].g_p) for t in main.values()} for r in REG5}).T
        D = M - M.mean(axis=0)
        diag = np.array([D.loc[r, main[r]] for r in REG5])
        rowmax = sum(D.loc[r].idxmax() == main[r] for r in REG5)
        return M, D, float(diag.mean()), int(rowmax)

    lines = []
    for label, main in (('main', MAIN), ('alt_cortex_cerebellar_hemisphere', {**MAIN, **ALT})):
        if not all(t in Q for t in main.values()):
            lines.append(f'[{label}] missing tissues, skipped'); continue
        L = pd.concat([Q[t][Q[t].kind == 'lead'].assign(tissue=t) for t in main.values()], ignore_index=True)
        L = L[L.region.isin(REG5)].drop_duplicates(['region', 'gene', 'snp', 'tissue'])
        M, D, S, rowmax = S_stat(L, main)
        M.to_csv(f'{OUT}/region_matrix_{label}.tsv', sep='\t'); D.to_csv(f'{OUT}/region_matrix_centered_{label}.tsv', sep='\t')
        # permutation: shuffle the region label of the (gene, snp) pairs, keeping the number of pairs per region
        keys = L[['region', 'gene', 'snp']].drop_duplicates(['region', 'gene', 'snp'])
        rng = np.random.default_rng(SEED)
        null = []
        for _ in range(1000):
            k = keys.copy(); k['region_perm'] = rng.permutation(k.region.values)
            Lp = L.merge(k, on=['region', 'gene', 'snp']).drop(columns='region').rename(columns={'region_perm': 'region'})
            null.append(S_stat(Lp, main)[2])
        p = (1 + sum(x >= S for x in null)) / 1001
        verdict = 'SUPPORTED' if (S > 0 and p < 0.05 and rowmax >= 3) else 'NOT SUPPORTED'
        lines.append(f'[{label}] S (mean centred diagonal pi1) = {S:.4f} | permutation P = {p:.4f} | rows where the matched tissue is the row maximum = {rowmax}/5 -> {verdict}')
        lines.append('pi1 matrix (rows = our regions, columns = GTEx tissues):\n' + M.round(3).to_string())
        lines.append('column-centred:\n' + D.round(3).to_string())
    open(OUT + '/criteria.txt', 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    print(R[R.tissue.isin(MAIN.values())].pivot_table(index=['kind', 'region'], columns='tissue', values='concordance').round(3).to_string())


if __name__ == '__main__':
    if sys.argv[1] == 'query':
        query(sys.argv[2])
    else:
        stats()
