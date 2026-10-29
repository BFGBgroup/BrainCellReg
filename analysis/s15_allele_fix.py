# -*- coding: utf-8 -*-
"""Which allele do the genotypes count? out.bim letters are unreliable (see s14b_allele_diag.py), so the counted
allele is identified per variant from its frequency against 1000G EUR.
  table                 genome-wide table of counted / other allele  -> allele_fix/allele_table.tsv.gz, exclude.txt
  esi <ctx>             rewrite the SMR .esi of one context with these alleles (besd / epi are symlinked)
  compare <old> <new>   compare two SMR result files
"""
import os
import sys

import numpy as np
import pandas as pd

IN = '/public/home/yucong/scbrain_in'
GEN = IN + '/paper2/single-cell-eqtl/data/genotype'
EQ = IN + '/SMR_inputs/SMR'
FIX = '/public/home/yucong/scbrain/allele_fix'
G1000 = '/public/home/yucong/scbrain/sn_replication/g1000.frq'


def table():
    os.makedirs(FIX + '/eqtl', exist_ok=True)
    bim = pd.read_csv(GEN + '/out.bim', sep='\t', header=None, usecols=[1, 4, 5], names=['snp', 'a1', 'a2'])
    dup = int(bim.snp.duplicated().sum())
    bim = bim.drop_duplicates('snp', keep='last')
    frq = pd.read_csv(GEN + '/MAF.txt', sep=r'\s+', usecols=[1, 5], names=['snp', 'f1'], header=0).drop_duplicates('snp', keep='last')
    g = pd.read_csv(G1000, sep=r'\s+', usecols=['SNP', 'A1', 'A2', 'MAF'])
    g.columns = ['snp', 'g1', 'g2', 'gf']
    g = g.drop_duplicates('snp', keep='last')
    d = bim.merge(frq, on='snp', how='left').merge(g, on='snp', how='left')
    fa1 = np.where(d.a1 == d.g1, d.gf, np.where(d.a1 == d.g2, 1 - d.gf, np.nan))
    set_ok = ((d.a1 == d.g1) & (d.a2 == d.g2)) | ((d.a1 == d.g2) & (d.a2 == d.g1))
    d1, d2 = np.abs(d.f1 - fa1), np.abs(d.f1 - (1 - fa1))
    clear = set_ok & (np.abs(d1 - d2) > 0.2)
    d['status'] = np.select([d.gf.isna(), ~set_ok, ~clear], ['not_in_1000G', 'letters_differ_from_1000G', 'frequency_not_decisive'], 'ok')
    col5 = (d1 <= d2) | (d.status != 'ok')
    d['eff'] = np.where(col5, d.a1, d.a2)
    d['oth'] = np.where(col5, d.a2, d.a1)
    d[['snp', 'eff', 'oth', 'f1', 'status', 'a1', 'a2']].to_csv(FIX + '/allele_table.tsv.gz', sep='\t', index=False)
    d.loc[d.status != 'ok', 'snp'].to_csv(FIX + '/exclude.txt', index=False, header=False)
    print('variants', len(d), '| duplicated ids dropped', dup)
    print(d.status.value_counts().to_string())
    ok = d[d.status == 'ok']
    print('among usable variants: counted allele is the letter in bim column 6 for %.4f ; column 5 for %.4f' % (
        float((ok.eff == ok.a2).mean()), float((ok.eff == ok.a1).mean())))


TABLES = {'v1': ('allele_table.tsv.gz', 'eqtl', 'v2'),          # frequency only
          'ld': ('allele_table_v2.tsv.gz', 'eqtl_ld', 'ld')}    # frequency + LD (s17_ld_resolve.py)


def esi(ctx, tag='ld'):
    # SMR 1.3.1 hangs on --exclude-snp with a list of 2.3 million ids, so variants whose counted allele cannot be
    # identified are taken out here instead: their id gets the prefix "unres_" and no longer matches the GWAS or the LD reference.
    tab, sub, mark = TABLES[tag]
    os.makedirs('%s/%s' % (FIX, sub), exist_ok=True)
    t = pd.read_csv('%s/%s' % (FIX, tab), sep='\t', usecols=['snp', 'eff', 'oth', 'f1', 'status'])
    e = pd.read_csv('%s/%s.esi' % (EQ, ctx), sep=r'\s+', header=None, names=['chr', 'snp', 'gd', 'bp', 'o1', 'o2', 'of'], dtype=str)
    n = len(e)
    m = e.merge(t, on='snp', how='left')
    assert len(m) == n and (m.snp.values == e.snp.values).all()
    miss = m.eff.isna()
    bad = miss | ~m.status.isin(['ok', 'ok_ld'])
    same_old = int(((m.o1 == m.eff) & (m.o2 == m.oth) & ~bad).sum())
    swap_old = int(((m.o1 == m.oth) & (m.o2 == m.eff) & ~bad).sum())
    m.loc[miss, 'eff'] = m.o1[miss]
    m.loc[miss, 'oth'] = m.o2[miss]
    m['f1'] = m.f1.astype(float).round(6).fillna(0.5)
    m['snp'] = np.where(bad, 'unres_' + m.snp, m.snp)
    m[['chr', 'snp', 'gd', 'bp', 'eff', 'oth', 'f1']].to_csv('%s/%s/%s.esi' % (FIX, sub, ctx), sep='\t', header=False, index=False)
    for ext in ('besd', 'epi'):
        dst = '%s/%s/%s.%s' % (FIX, sub, ctx, ext)
        if not os.path.lexists(dst):
            os.symlink('%s/%s.%s' % (EQ, ctx, ext), dst)
    open('%s/%s/%s.esi.%s' % (FIX, sub, ctx, mark), 'w').close()
    nb = int(bad.sum())
    print('%s esi rows %d | counted allele not identified, id prefixed %d (not in allele table %d) | usable %d: old letters equal to counted/other %d, '
          'reversed %d, a different pair %d' % (ctx, n, nb, int(miss.sum()), n - nb, same_old, swap_old, n - nb - same_old - swap_old))


def compare(old, new):
    a = pd.read_csv(old, sep='\t')
    b = pd.read_csv(new, sep='\t')
    print('probes old %d new %d' % (len(a), len(b)))
    m = a.merge(b, on='probeID', suffixes=('_o', '_n'))
    print('probes in both', len(m), '| same top SNP: %.4f' % float((m.topSNP_o == m.topSNP_n).mean()))
    tab = FIX + '/allele_table_v2.tsv.gz'
    st = pd.read_csv(tab if os.path.exists(tab) else FIX + '/allele_table.tsv.gz', sep='\t', usecols=['snp', 'status']).set_index('snp').status
    print('old top SNPs by status in the allele table:', a.topSNP.map(st).fillna('absent').value_counts().to_dict())
    lost = a[~a.probeID.isin(b.probeID)]
    print('probes only in old %d (of which p_SMR<1e-4: %d) | only in new %d' % (len(lost), int((lost.p_SMR < 1e-4).sum()), int((~b.probeID.isin(a.probeID)).sum())))
    s = m[m.topSNP_o == m.topSNP_n]
    hh = s.dropna(subset=['p_HEIDI_o', 'p_HEIDI_n'])
    print('same top SNP, HEIDI in both %d: verdict at 0.01 agrees %.4f | old pass %d, new pass %d' % (
        len(hh), float(((hh.p_HEIDI_o >= 0.01) == (hh.p_HEIDI_n >= 0.01)).mean()), int((hh.p_HEIDI_o >= 0.01).sum()), int((hh.p_HEIDI_n >= 0.01).sum())))
    print('same top SNP: corr of -log10 p_SMR = %.5f | b_SMR same sign: %.4f | median |b| ratio new/old = %.4f' % (
        np.corrcoef(-np.log10(s.p_SMR_o), -np.log10(s.p_SMR_n))[0, 1], float((np.sign(s.b_SMR_o) == np.sign(s.b_SMR_n)).mean()),
        float((s.b_SMR_n.abs() / s.b_SMR_o.abs()).median())))
    for tag, x in (('old', a), ('new', b)):
        h = x.dropna(subset=['p_HEIDI'])
        sig = x[x.p_SMR < 1e-4]
        hs = sig.dropna(subset=['p_HEIDI'])
        print('%s: HEIDI available for %d probes, median nsnp_HEIDI %.0f, share p_HEIDI>=0.01: %.4f | probes with p_SMR<1e-4: %d, of which '
              'HEIDI available %d and p_HEIDI>=0.01: %d' % (tag, len(h), h.nsnp_HEIDI.median(), float((h.p_HEIDI >= 0.01).mean()), len(sig),
                                                             len(hs), int((hs.p_HEIDI >= 0.01).sum())))
    print('p_SMR<1e-4 in old %d, in new %d, in both %d' % (int((a.p_SMR < 1e-4).sum()), int((b.p_SMR < 1e-4).sum()),
                                                         int(((m.p_SMR_o < 1e-4) & (m.p_SMR_n < 1e-4)).sum())))


if __name__ == '__main__':
    {'table': table, 'esi': lambda: esi(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else 'ld'), 'compare': lambda: compare(sys.argv[2], sys.argv[3])}[sys.argv[1]]()
