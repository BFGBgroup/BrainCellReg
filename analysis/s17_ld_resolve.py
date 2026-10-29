# -*- coding: utf-8 -*-
"""Counted allele of the variants that allele frequency cannot settle (frequency near 0.5), resolved from LD.

s15 identifies the allele counted by the genotypes from its frequency against 1000G EUR; about 0.94 million variants
are left open because both letters fit. For such a variant j and a nearby variant k whose counted allele is known,
the sign of the genotype correlation in our 134 donors and in 1000G EUR tells whether j is coded the same way in the
two data sets: s_j = sign(r_ours(j,k)) * sign(r_1000G(j,k)) * s_k. Neighbours vote; a call needs >= MIN_INF informative
neighbours (|r| >= RMIN in both data sets) and >= AGREE of them on the same side.

  resolve   LD calls for the open variants + leave-one-out check on every frequency-resolved variant
            -> allele_fix/allele_table_v2.tsv.gz ; exits 1 if the leave-one-out accuracy is below ACC_MIN
  bryois    sign concordance with Bryois 2022 (pooled nuclei, P < 0.05) of lead pairs, by how the allele was resolved;
            exits 1 if the LD-resolved pairs fall below CONC_MIN
"""
import os
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

IN = '/public/home/yucong/scbrain_in'
GEN = IN + '/paper2/single-cell-eqtl/data/genotype'
REF = '/public/home/yucong/scbrain/local_handoff/scbrain1/g1000_eur/g1000_eur'
FIX = '/public/home/yucong/scbrain/allele_fix'
SNR = '/public/home/yucong/scbrain/sn_replication'
BRY = '/public/home/yucong/scbrain/external/bryois2022'
N1, N2 = 134, 503
K, W, RMIN, MIN_INF, AGREE = 20, 500000, 0.3, 3, 0.9
ACC_MIN, CONC_MIN = 0.99, 0.85
CH = 3000
LUT = np.array([2, np.nan, 1, 0], dtype=np.float32)   # 00 hom A1 (bim column 5), 01 missing, 10 het, 11 hom A2
D = None


def bed_rows(path, n, idx):
    """Counts of the bim column 5 allele for the variants in idx (0-based rows of the bim)."""
    with open(path, 'rb') as fh:
        assert fh.read(3) == b'\x6c\x1b\x01', 'not a SNP-major PLINK bed: ' + path
    nb = (n + 3) // 4
    mm = np.memmap(path, dtype=np.uint8, mode='r', offset=3).reshape(-1, nb)
    raw = np.asarray(mm[idx])
    g = np.empty((len(idx), nb * 4), dtype=np.uint8)
    for k in range(4):
        g[:, k::4] = (raw >> (2 * k)) & 3
    return LUT[g[:, :n]]


def zs(G):
    mu = np.nanmean(G, axis=1, keepdims=True)
    G = np.where(np.isnan(G), mu, G) - mu
    sd = G.std(axis=1, keepdims=True)
    ok = sd[:, 0] > 0
    G[ok] /= sd[ok]
    return G.astype(np.float32), ok, (mu[:, 0] / 2).astype(np.float32)


def run_chr(ch):
    d = D[D.chr == ch].sort_values('bp').reset_index(drop=True)
    Z1, ok1, fr1 = zs(bed_rows(GEN + '/out.bed', N1, d.oi.values))
    Z2, ok2, fr2 = zs(bed_rows(REF + '.bed', N2, d.ri.values))
    bp = d.bp.values
    s = d.s.values.astype(np.int8)
    pool = np.where((s != 0) & ok1 & ok2)[0]
    T = np.where(ok1 & ok2)[0]
    off = np.arange(-K, K + 1)
    n_inf = np.zeros(len(d), dtype=np.int16)
    n_pos = np.zeros(len(d), dtype=np.int16)
    for a in range(0, len(T), CH):
        t = T[a:a + CH]
        ci = np.searchsorted(pool, t)[:, None] + off[None, :]
        valid = (ci >= 0) & (ci < len(pool))
        nb = pool[np.clip(ci, 0, len(pool) - 1)]
        valid &= (nb != t[:, None]) & (np.abs(bp[nb] - bp[t][:, None]) <= W)
        r1 = np.matmul(Z1[nb], Z1[t][:, :, None])[:, :, 0] / N1
        r2 = np.matmul(Z2[nb], Z2[t][:, :, None])[:, :, 0] / N2
        inf = valid & (np.abs(r1) >= RMIN) & (np.abs(r2) >= RMIN)
        v = np.sign(r1 * r2) * s[nb]
        n_inf[t] = inf.sum(1)
        n_pos[t] = ((v > 0) & inf).sum(1)
    d['ld_n'] = n_inf
    d['ld_pos'] = n_pos
    d['fr_ours'] = fr1
    d['fr_ref'] = fr2
    return d[['snp', 'chr', 's', 'ld_n', 'ld_pos', 'fr_ours', 'fr_ref']]


def resolve():
    global D
    t = pd.read_csv(FIX + '/allele_table.tsv.gz', sep='\t')
    ob = pd.read_csv(GEN + '/out.bim', sep='\t', header=None, usecols=[1], names=['snp'])
    ob['oi'] = np.arange(len(ob))
    ob = ob.drop_duplicates('snp', keep='last')
    rb = pd.read_csv(REF + '.bim', sep='\t', header=None, usecols=[0, 1, 3, 4, 5], names=['chr', 'snp', 'bp', 'g1', 'g2'], dtype={'chr': str})
    rb['ri'] = np.arange(len(rb))
    rb = rb.drop_duplicates('snp', keep='last')
    d = t[t.status.isin(['ok', 'frequency_not_decisive'])].merge(ob, on='snp').merge(rb, on='snp')
    # s = +1 when the counted allele is the 1000G bim column 5 allele, -1 when it is the column 6 allele, 0 when open
    d['s'] = np.where(d.status != 'ok', 0, np.where(d.eff == d.g1, 1, np.where(d.eff == d.g2, -1, 0)))
    print('variants entering: frequency-resolved %d, open %d | chromosomes %s' % (
        int((d.status == 'ok').sum()), int((d.status != 'ok').sum()), ' '.join(sorted(d.chr.unique(), key=lambda x: (len(x), x)))))
    assert int(((d.status == 'ok') & (d.s == 0)).sum()) == 0
    D = d[['snp', 'chr', 'bp', 'oi', 'ri', 's']]
    chs = sorted(D.chr.unique(), key=lambda x: -int((D.chr == x).sum()))
    with Pool(6) as p:
        res = pd.concat(p.map(run_chr, chs, chunksize=1), ignore_index=True)
    d = d.merge(res[['snp', 'ld_n', 'ld_pos', 'fr_ours', 'fr_ref']], on='snp')

    # the bed parser must reproduce the frequencies on file
    c = d[d.f1.notna() & d.fr_ours.notna()]
    print('parser check, ours: corr(freq from bed, f1 on file) = %.5f, median |diff| = %.5f' % (
        np.corrcoef(c.fr_ours, c.f1)[0, 1], float(np.abs(c.fr_ours - c.f1).median())))
    assert np.corrcoef(c.fr_ours, c.f1)[0, 1] > 0.999
    g = pd.read_csv(SNR + '/g1000.frq', sep=r'\s+', usecols=['SNP', 'A1', 'MAF']).drop_duplicates('SNP', keep='last')
    g.columns = ['snp', 'fa1', 'maf']
    c = d.merge(g, on='snp')
    c = c[c.fr_ref.notna() & c.maf.notna()]
    exp = np.where(c.fa1 == c.g1, c.maf, 1 - c.maf)
    print('parser check, 1000G: corr(freq from bed, plink --freq) = %.5f, median |diff| = %.5f' % (
        np.corrcoef(c.fr_ref, exp)[0, 1], float(np.abs(c.fr_ref - exp).median())))
    assert np.corrcoef(c.fr_ref, exp)[0, 1] > 0.999

    frac = np.where(d.ld_n > 0, d.ld_pos / np.maximum(d.ld_n, 1), np.nan)
    called = (d.ld_n >= MIN_INF) & ((frac >= AGREE) | (frac <= 1 - AGREE))
    d['ld_s'] = np.where(called, np.where(frac >= 0.5, 1, -1), 0)
    d['ld_agree'] = np.round(np.maximum(frac, 1 - frac), 3)

    # leave-one-out on the frequency-resolved variants
    k = d[d.status == 'ok']
    kc = k[k.ld_s != 0]
    acc = float((kc.ld_s == kc.s).mean())
    print('leave-one-out on frequency-resolved variants: %d, LD call made for %d (%.4f), LD call equals frequency call %.5f' % (
        len(k), len(kc), len(kc) / len(k), acc))
    fa = np.where(k.a1 == k.g1, k.fr_ref, 1 - k.fr_ref)
    margin = np.abs(np.abs(k.f1 - fa) - np.abs(k.f1 - (1 - fa)))
    for lo, hi in ((0.2, 0.3), (0.3, 0.5), (0.5, 0.8), (0.8, 2.1)):
        b = k[(margin > lo) & (margin <= hi)]
        bc = b[b.ld_s != 0]
        print('  frequency margin %.1f-%.1f: n %d, called %.4f, agree %.5f' % (lo, hi, len(b), len(bc) / max(len(b), 1),
                                                                         float((bc.ld_s == bc.s).mean()) if len(bc) else float('nan')))
    conflict = (d.status == 'ok') & (d.ld_s != 0) & (d.ld_s != d.s) & (d.ld_n >= 5)
    o = d[d.status != 'ok']
    print('open variants: %d, LD call made for %d (%.4f); no informative neighbour %d, fewer than %d %d, neighbours split %d' % (
        len(o), int((o.ld_s != 0).sum()), float((o.ld_s != 0).mean()), int((o.ld_n == 0).sum()),
        MIN_INF, int(((o.ld_n > 0) & (o.ld_n < MIN_INF)).sum()), int(((o.ld_n >= MIN_INF) & (o.ld_s == 0)).sum())))
    print('frequency-resolved variants contradicted by LD (>= 5 informative neighbours, >= %.0f%% against): %d -> taken out' % (
        AGREE * 100, int(conflict.sum())))

    new = d.status.copy()
    new[(d.status != 'ok') & (d.ld_s != 0)] = 'ok_ld'
    new[conflict] = 'frequency_ld_conflict'
    use_ld = new == 'ok_ld'
    d.loc[use_ld, 'eff'] = np.where(d.ld_s[use_ld] > 0, d.g1[use_ld], d.g2[use_ld])
    d.loc[use_ld, 'oth'] = np.where(d.ld_s[use_ld] > 0, d.g2[use_ld], d.g1[use_ld])
    d['status'] = new
    out = t.merge(d[['snp', 'eff', 'oth', 'status', 'ld_n', 'ld_agree']], on='snp', how='left', suffixes=('', '_n'))
    hit = out.status_n.notna()
    for col in ('eff', 'oth', 'status'):
        out.loc[hit, col] = out.loc[hit, col + '_n']
    out = out.drop(columns=['eff_n', 'oth_n', 'status_n'])
    assert len(out) == len(t)
    out.to_csv(FIX + '/allele_table_v2.tsv.gz', sep='\t', index=False)
    print(out.status.value_counts().to_string())
    ok = out[out.status.isin(['ok', 'ok_ld'])]
    print('usable variants %d of %d | among LD-resolved, counted allele is bim column 6 for %.4f' % (
        len(ok), len(out), float((out[out.status == 'ok_ld'].eff == out[out.status == 'ok_ld'].a2).mean())))
    if acc < ACC_MIN:
        print('LD_CHECK_FAILED accuracy %.5f < %.2f' % (acc, ACC_MIN))
        sys.exit(1)
    print('LD_RESOLVE_DONE')


def bryois():
    t = pd.read_csv(FIX + '/allele_table_v2.tsv.gz', sep='\t', usecols=['snp', 'eff', 'oth', 'status'])
    pos = pd.read_csv(BRY + '/snp_pos.txt', sep='\t', usecols=['SNP', 'effect_allele', 'other_allele'])
    pos.columns = ['snp', 'ea', 'oa']
    pos = pos.drop_duplicates('snp', keep='last')
    fs = [f for f in os.listdir(SNR + '/bry') if f.startswith('pb.') and os.path.getsize(SNR + '/bry/' + f) > 0]
    b = pd.concat([pd.read_csv(SNR + '/bry/' + f, sep='\t', header=None, names=['gene', 'snp', 'dist_b', 'p_b', 'beta_b']) for f in fs])
    b = b.drop_duplicates(['gene', 'snp'])
    rc = 0
    for kind in ('lead', 'ctrl'):
        P = pd.read_csv('%s/pairs_%s.tsv.gz' % (SNR, kind), sep='\t', usecols=['region', 'cell', 'gene', 'snp', 'slope'])
        P = P.merge(t, on='snp').merge(pos, on='snp').merge(b[['gene', 'snp', 'p_b', 'beta_b']], on=['gene', 'snp'])
        same = (P.eff == P.ea) & (P.oth == P.oa)
        flip = (P.eff == P.oa) & (P.oth == P.ea)
        P = P[same | flip].copy()
        P['slope'] = np.where(P.eff == P.ea, P.slope, -P.slope)
        for st in ('ok', 'ok_ld', 'frequency_not_decisive'):
            q = P[(P.status == st) & (P.p_b < 0.05)]
            conc = float((np.sign(q.slope) == np.sign(q.beta_b)).mean()) if len(q) else float('nan')
            print('%s pairs, allele status %s: in Bryois pooled %d, P<0.05 %d, same direction %.4f' % (
                kind, st, int((P.status == st).sum()), len(q), conc))
            if kind == 'lead' and st == 'ok_ld' and len(q) >= 100 and conc < CONC_MIN:
                rc = 1
    if rc:
        print('BRYOIS_CHECK_FAILED: LD-resolved lead pairs below %.2f' % CONC_MIN)
        sys.exit(1)
    print('BRYOIS_CHECK_DONE')


if __name__ == '__main__':
    {'resolve': resolve, 'bryois': bryois}[sys.argv[1]]()
