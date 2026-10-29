# -*- coding: utf-8 -*-
"""Replication of the lead eSNP of every eGene in the single-nucleus eQTL of Bryois 2022 (spec:
progress/分析规格_单核eQTL复现检验_v1.md). Three steps, run in this order by s14_snrep.slurm:
  prep   lead pair of every eGene (q<0.05) in every context + allele bookkeeping
  match  one matched non-eQTL pair per lead pair (same context, 10-kb TSS-distance bin, 0.05 MAF bin, P>0.05)
  stats  look-ups in Bryois -> pi1 (lambda 0.5), sign concordance, Spearman; matrices and the two criteria
"""
import glob
import os
import sys

import numpy as np
import pandas as pd

IN = '/public/home/yucong/scbrain_in/paper2/single-cell-eqtl/data'
BRY = '/public/home/yucong/scbrain/external/bryois2022'
OUT = '/public/home/yucong/scbrain/sn_replication'
REG = ['CRBL', 'FCTX', 'HIPP', 'MEDU', 'OCTX', 'PUTM', 'SNIG', 'TCTX', 'THAL', 'WHMT']
CELL = ['Astro', 'Endo', 'ExN', 'Gran', 'InN', 'Micro', 'NPC', 'Olig', 'OPC', 'Perc', 'Purk']
MAIN = ['FCTX', 'TCTX', 'WHMT']
MAP = {'Astro': 'Astrocytes', 'Endo': 'Endothelial.cells', 'ExN': 'Excitatory.neurons', 'InN': 'Inhibitory.neurons',
       'Micro': 'Microglia', 'Olig': 'Oligodendrocytes', 'OPC': 'OPCs...COPs', 'Perc': 'Pericytes'}
SETS = list(MAP.values()) + ['pb']
SEED = 20261009
LAMBDA = 0.5


def prep():
    os.makedirs(OUT, exist_ok=True)
    bim = pd.read_csv(IN + '/genotype/out.bim', sep='\t', header=None, usecols=[1, 4, 5], names=['snp', 'a1', 'a2'])
    bim = bim.drop_duplicates('snp', keep='last')
    frq = pd.read_csv(IN + '/genotype/MAF.txt', sep=r'\s+', usecols=[1, 3, 5], names=['snp', 'alt', 'alt_frq'], header=0)
    frq = frq.drop_duplicates('snp', keep='last')
    rows, summ = [], []
    for r in REG:
        for c in CELL:
            f = '%s/run/new2/%s.%s.cis_qtl.txt' % (IN, r, c)
            if not os.path.exists(f):
                summ.append((r, c, 0, 0))
                continue
            d = pd.read_csv(f, sep='\t')
            e = d[d.qval < 0.05]
            summ.append((r, c, len(d), len(e)))
            if len(e):
                e = e[['phenotype_id', 'variant_id', 'start_distance', 'af', 'pval_nominal', 'slope']].copy()
                e.columns = ['gene', 'snp', 'dist', 'af', 'p', 'slope']
                e.insert(0, 'cell', c)
                e.insert(0, 'region', r)
                rows.append(e)
    lead = pd.concat(rows, ignore_index=True)
    lead = lead.merge(bim, on='snp', how='left').merge(frq, on='snp', how='left')
    # which allele does the slope refer to: compare the af column with the frequency of bim A1
    k = lead.dropna(subset=['a1', 'alt', 'alt_frq'])
    same = k[k.alt == k.a1]
    print('lead pairs', len(lead), '| contexts with eGenes', lead.groupby(['region', 'cell']).ngroups,
          '| in bim', int(lead.a1.notna().sum()))
    print('afreq ALT equals bim A1: %d of %d' % (len(same), len(k)))
    print('corr(af, freq of bim A1) = %.4f ; corr(af, 1 - freq) = %.4f ; mean |af - freq| = %.4f' % (
        np.corrcoef(same.af, same.alt_frq)[0, 1], np.corrcoef(same.af, 1 - same.alt_frq)[0, 1],
        float(np.abs(same.af - same.alt_frq).mean())))
    lead['maf'] = np.minimum(lead.af, 1 - lead.af)
    lead.drop(columns=['alt', 'alt_frq']).to_csv(OUT + '/pairs_lead.tsv.gz', sep='\t', index=False)
    pd.DataFrame(summ, columns=['region', 'cell', 'genes_tested', 'egenes']).to_csv(OUT + '/context_counts.tsv', sep='\t', index=False)
    with open(OUT + '/contexts.txt', 'w') as fh:
        for (r, c), _ in lead.groupby(['region', 'cell']):
            fh.write('%s.%s\n' % (r, c))


def match():
    lead = pd.read_csv(OUT + '/pairs_lead.tsv.gz', sep='\t')
    bim = pd.read_csv(IN + '/genotype/out.bim', sep='\t', header=None, usecols=[1, 4, 5], names=['snp', 'a1', 'a2'])
    bim = bim.drop_duplicates('snp', keep='last')
    rng = np.random.default_rng(SEED)
    out, log = [], []
    for (r, c), L in lead.groupby(['region', 'cell'], sort=True):
        f = '%s/cand/%s.%s.tsv' % (OUT, r, c)
        if not os.path.exists(f):
            log.append((r, c, len(L), 0, 0, len(L)))
            continue
        C = pd.read_csv(f, sep='\t')
        C['maf'] = np.minimum(C.af, 1 - C.af)
        C['db'] = np.floor(C.dist / 10000.0).astype(int)
        C['mb'] = np.minimum(np.floor(C.maf / 0.05).astype(int), 9)
        own = set(zip(L.gene, L.snp))
        C = C[[(g, s) not in own for g, s in zip(C.gene, C.snp)]]
        groups = {k: list(v) for k, v in C.groupby(['db', 'mb']).indices.items()}
        used, n_exact, n_near, n_none = set(), 0, 0, 0
        L = L.sort_values(['gene', 'snp'])
        for g, s, dist, maf in zip(L.gene, L.snp, L.dist, L.maf):
            db, mb = int(np.floor(dist / 10000.0)), min(int(np.floor(maf / 0.05)), 9)
            pick = None
            for step in range(0, 6):
                for d2 in ([db] if step == 0 else [db - step, db + step]):
                    pool = [i for i in groups.get((d2, mb), []) if i not in used]
                    if pool:
                        pick = pool[int(rng.integers(len(pool)))]
                        break
                if pick is not None:
                    n_exact += step == 0
                    n_near += step > 0
                    break
            if pick is None:
                n_none += 1
                continue
            used.add(pick)
            x = C.iloc[pick]
            out.append((r, c, x.gene, x.snp, x.dist, x.af, x.p, x.slope, x.maf, g, s))
        log.append((r, c, len(L), n_exact, n_near, n_none))
    ctrl = pd.DataFrame(out, columns=['region', 'cell', 'gene', 'snp', 'dist', 'af', 'p', 'slope', 'maf', 'lead_gene', 'lead_snp'])
    ctrl = ctrl.merge(bim, on='snp', how='left')
    ctrl.to_csv(OUT + '/pairs_ctrl.tsv.gz', sep='\t', index=False)
    lg = pd.DataFrame(log, columns=['region', 'cell', 'lead', 'same_bin', 'neighbour_bin', 'unmatched'])
    lg.to_csv(OUT + '/match_log.tsv', sep='\t', index=False)
    print('controls', len(ctrl), '| same bin', int(lg.same_bin.sum()), '| neighbour bin', int(lg.neighbour_bin.sum()),
          '| unmatched', int(lg.unmatched.sum()))
    keys = pd.concat([lead[['gene', 'snp']], ctrl[['gene', 'snp']]]).drop_duplicates()
    keys.to_csv(OUT + '/keys.txt', sep=' ', index=False, header=False)
    print('keys', len(keys))


def pi1(p):
    p = np.asarray(p, dtype=float)
    if len(p) < 20:
        return np.nan
    return float(np.clip(1 - (p > LAMBDA).mean() / (1 - LAMBDA), 0, 1))


def metrics(d):
    """d: aligned rows with p_b, beta_b, slope (already on the same allele)."""
    n = len(d)
    if n == 0:
        return dict(n=0, pi1=np.nan, n_p05=0, concord=np.nan, rho=np.nan)
    sig = d[d.p_b < 0.05]
    conc = float((np.sign(sig.slope) == np.sign(sig.beta_b)).mean()) if len(sig) >= 10 else np.nan
    rho = float(d[['slope', 'beta_b']].corr(method='spearman').iloc[0, 1]) if n >= 20 else np.nan
    return dict(n=n, pi1=pi1(d.p_b), n_p05=len(sig), concord=conc, rho=rho)


def stats():
    global RES
    RES = os.environ.get('SNREP_RES', OUT)
    os.makedirs(RES, exist_ok=True)
    from scipy.stats import wilcoxon
    pos = pd.read_csv(BRY + '/snp_pos.txt', sep='\t', usecols=['SNP', 'effect_allele', 'other_allele'])
    pos.columns = ['snp', 'ea', 'oa']
    pos = pos.drop_duplicates('snp', keep='last')
    bry = []
    for f in sorted(glob.glob(OUT + '/bry/*.tsv')):
        if os.path.getsize(f) == 0:
            continue
        s = os.path.basename(f)[:-4].rsplit('.', 1)[0]
        d = pd.read_csv(f, sep='\t', header=None, names=['gene', 'snp', 'dist_b', 'p_b', 'beta_b'])
        d['set'] = s
        bry.append(d)
    bry = pd.concat(bry, ignore_index=True).drop_duplicates(['set', 'gene', 'snp'])
    print('bryois rows kept', len(bry), '| sets', sorted(bry.set.unique()))
    qc = []

    # The allele letters in out.bim do not tell which allele the genotypes count (s14b_allele_diag.py: for about
    # three quarters of the variants the frequency of the allele written in column 5 is 1 - its 1000G EUR frequency).
    # The counted allele is therefore identified per variant by frequency: of the two letters, the one whose
    # 1000G EUR frequency is closer to the af reported by tensorQTL. Variants where this is not clear are dropped.
    g = pd.read_csv(OUT + '/g1000.frq', sep=r'\s+', usecols=['SNP', 'A1', 'A2', 'MAF'])
    g.columns = ['snp', 'g1', 'g2', 'gf']
    g = g.drop_duplicates('snp', keep='last')

    AT = os.environ.get('ALLELE_TABLE')     # sensitivity run: counted allele from allele_table_v2 (frequency + LD), not re-inferred here
    if AT:
        at = pd.read_csv(AT, sep='\t', usecols=['snp', 'eff', 'oth', 'status'])
        at = at[at.status.isin(['ok', 'ok_ld'])].drop_duplicates('snp')
        print('allele table', AT, '| usable variants', len(at))

    def align(P, kind):
        n0 = len(P)
        if AT:
            P = P.merge(at, on='snp', how='left').merge(pos, on='snp', how='left')
            known = P.eff.notna()
            set_ok = known & (((P.a1 == P.eff) & (P.a2 == P.oth)) | ((P.a1 == P.oth) & (P.a2 == P.eff)))
            nopos = P.ea.isna()
            same = (P.eff == P.ea) & (P.oth == P.oa)
            flip = (P.eff == P.oa) & (P.oth == P.ea)
            use = set_ok & (same | flip)
            qc.append((kind, n0, int((~known).sum()), int((known & ~set_ok).sum()), 0, int((set_ok & nopos).sum()),
                       int((set_ok & ~nopos & ~same & ~flip).sum()), int(use.sum()), float((P.eff[use] != P.a1[use]).mean())))
            P = P[use].copy()
            P['slope'] = np.where(P.eff == P.ea, P.slope, -P.slope)
            return P.drop(columns=['ea', 'oa', 'status'])
        P = P.merge(g, on='snp', how='left').merge(pos, on='snp', how='left')
        no_ref = P.gf.isna()
        f1 = np.where(P.a1 == P.g1, P.gf, np.where(P.a1 == P.g2, 1 - P.gf, np.nan))
        set_ok = ((P.a1 == P.g1) & (P.a2 == P.g2)) | ((P.a1 == P.g2) & (P.a2 == P.g1))
        d1, d2 = np.abs(P.af - f1), np.abs(P.af - (1 - f1))
        clear = set_ok & (np.abs(d1 - d2) > 0.2)
        P['eff'] = np.where(d1 <= d2, P.a1, P.a2)
        P['oth'] = np.where(d1 <= d2, P.a2, P.a1)
        nopos = P.ea.isna()
        same = (P.eff == P.ea) & (P.oth == P.oa)
        flip = (P.eff == P.oa) & (P.oth == P.ea)
        use = clear & (same | flip)
        qc.append((kind, n0, int(no_ref.sum()), int((~no_ref & ~set_ok).sum()), int((set_ok & ~clear).sum()),
                   int((clear & nopos).sum()), int((clear & ~nopos & ~same & ~flip).sum()), int(use.sum()),
                   float((P.eff[use] != P.a1[use]).mean())))
        P = P[use].copy()
        P['slope'] = np.where(P.eff == P.ea, P.slope, -P.slope)
        return P.drop(columns=['ea', 'oa', 'g1', 'g2', 'gf'])

    lead = align(pd.read_csv(OUT + '/pairs_lead.tsv.gz', sep='\t'), 'lead')
    ctrl = align(pd.read_csv(OUT + '/pairs_ctrl.tsv.gz', sep='\t'), 'control')
    pd.DataFrame(qc, columns=['kind', 'pairs', 'snp_not_in_1000G', 'letters_differ_from_1000G', 'frequency_not_decisive',
                              'snp_not_in_bryois_table', 'letters_differ_from_bryois', 'usable',
                              'share_counted_allele_is_bim_col6']).to_csv(RES + '/align_qc.tsv', sep='\t', index=False)

    def look(P, s):
        return P.merge(bry[bry.set == s][['gene', 'snp', 'p_b', 'beta_b']], on=['gene', 'snp'], how='inner')

    # 1) every context x every Bryois set
    rows = []
    for kind, P in (('lead', lead), ('control', ctrl)):
        for (r, c), d in P.groupby(['region', 'cell']):
            for s in SETS:
                m = metrics(look(d, s))
                rows.append(dict(kind=kind, region=r, cell=c, set=s, pairs=len(d), **m))
    pd.DataFrame(rows).to_csv(RES + '/rep_by_context.tsv', sep='\t', index=False)

    # 2) region groups, pairs pooled over regions and de-duplicated (smallest P kept)
    def pooled(P, regs, cells=None):
        d = P[P.region.isin(regs)]
        if cells is not None:
            d = d[d.cell.isin(cells)]
        return d.sort_values('p').drop_duplicates(['gene', 'snp'])

    groups = {'main_FCTX_TCTX_WHMT': MAIN, 'other_7_regions': [r for r in REG if r not in MAIN], 'all_10_regions': REG}
    rows = []
    for gname, regs in groups.items():
        for kind, P in (('lead', lead), ('control', ctrl)):
            for c in CELL:
                d = pooled(P, regs, [c])
                for s in SETS:
                    rows.append(dict(group=gname, kind=kind, cell=c, set=s, pairs=len(d), **metrics(look(d, s))))
            d = pooled(P, regs)
            rows.append(dict(group=gname, kind=kind, cell='ALL', set='pb', pairs=len(d), **metrics(look(d, 'pb'))))
    M = pd.DataFrame(rows)
    M.to_csv(RES + '/rep_matrix.tsv', sep='\t', index=False)

    # 3) the two criteria fixed in the spec
    lines = []
    for gname in groups:
        a = M[(M.group == gname) & (M.cell == 'ALL')].set_index('kind')
        L, C = a.loc['lead'], a.loc['control']
        if L.pi1 >= 0.5 and L.concord >= 0.80:
            va = 'HOLDS'
        elif L.pi1 < 0.2 or L.concord < 0.60:
            va = 'FAILS'
        else:
            va = 'IN BETWEEN'
        lines.append('[%s] A (vs pb): lead n=%d pi1=%.3f concordance=%.3f rho=%.3f | control n=%d pi1=%.3f concordance=%.3f -> %s' % (
            gname, L.n, L.pi1, L.concord, L.rho, C.n, C.pi1, C.concord, va))
        X = M[(M.group == gname) & (M.kind == 'lead') & M.cell.isin(list(MAP)) & M.set.isin(list(MAP.values()))]
        W = X.pivot(index='cell', columns='set', values='pi1').reindex(index=list(MAP), columns=list(MAP.values()))
        W.to_csv(RES + '/pi1_8x8_%s.tsv' % gname, sep='\t')
        diag, off, top = [], [], 0
        for c, s in MAP.items():
            row = W.loc[c]
            if pd.isna(row[s]):
                continue
            diag.append(row[s])
            off.append(row.drop(s).mean())
            top += int(row[s] >= row.max() - 1e-12)
        if len(diag) >= 5:
            pw = float(wilcoxon(np.array(diag) - np.array(off), alternative='greater').pvalue)
        else:
            pw = np.nan
        vb = 'SUPPORTED' if (top >= 5 and pw < 0.05) else 'NOT SUPPORTED'
        lines.append('[%s] B (8x8): rows with pi1 available %d; diagonal is the row maximum in %d; mean diagonal %.3f vs off-diagonal %.3f; '
                     'one-sided paired Wilcoxon P=%.4g -> %s' % (gname, len(diag), top, float(np.mean(diag)) if diag else np.nan,
                                                               float(np.mean(off)) if off else np.nan, pw, vb))
    open(RES + '/criteria.txt', 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    print(pd.read_csv(RES + '/align_qc.tsv', sep='\t').to_string(index=False))


if __name__ == '__main__':
    {'prep': prep, 'match': match, 'stats': stats}[sys.argv[1]]()
