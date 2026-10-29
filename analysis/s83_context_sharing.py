# -*- coding: utf-8 -*-
"""Sharing of eQTLs between the 110 region-cell type contexts, and expression of eGenes vs non-eGenes in each context.
  pair_sharing.tsv        every pair of contexts: Jaccard of eGenes, Jaccard of eSNPs (official set), rb, pair class
  egene_expr_by_ctx.tsv   expression of eGenes vs non-eGenes in every context (medians, Wilcoxon P)
"""
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

SETS = '/public/home/yucong/scbrain/esnp_eg/sets'
LI = '/public/home/yucong/scbrain/local_inputs'
HO = '/public/home/yucong/scbrain/local_handoff/scbrain2/eqtl'
OUT = '/public/home/yucong/scbrain/results/context_sharing'
FIG5 = '/public/home/yucong/scbrain/figures/fig5'

m = np.load(SETS + '/membership.npz', allow_pickle=True)
ctx = [str(x) for x in m['ctx']]
# official eSNP set of a context = independent (clumped) eSNPs of its eGenes (q < 0.05): e3 & eg, as in s4_stats.py
E = m['e3'] & np.load(SETS + '/membership_eg.npz')['eg']
E = E[E.any(1)].astype(np.float32)
n_s = E.sum(0)
I_s = E.T @ E
print('eSNPs in any context', E.shape[0], '| contexts with eSNPs', int((n_s > 0).sum()))
assert E.shape[0] == 37412, 'whole-brain eSNP count differs from the official 37,412'

eg = pd.read_csv(FIG5 + '/egene_by_ctx.tsv', sep='\t', header=None, names=['ctx', 'gene', 'biotype'])
G = pd.crosstab(eg.gene, eg.ctx).reindex(columns=ctx, fill_value=0).values.astype(np.float32)
n_g = G.sum(0)
I_g = G.T @ G
print('eGenes in any context', G.shape[0], '| contexts with eGenes', int((n_g > 0).sum()))

rb = pd.read_csv(LI + '/rALLmatr.txt', sep=r'\s+', header=0, names=['pair', 'rb'], quotechar='"')
se = pd.read_csv(LI + '/r_seALLmatr.txt', sep=r'\s+', header=0, names=['pair', 'rb_se'], quotechar='"')
rb = rb.merge(se, on='pair')
rb[['a', 'b']] = rb.pair.str.split('vs', expand=True)
rbd = {}
for a, b, r, s in zip(rb.a, rb.b, rb.rb, rb.rb_se):
    rbd[(a, b)] = (r, s)
    rbd[(b, a)] = (r, s)
print('rb pairs', len(rb), '| missing', int(rb.rb.isna().sum()), '| range %.3f to %.3f' % (rb.rb.min(), rb.rb.max()))

rows = []
for i, a in enumerate(ctx):
    for j, b in enumerate(ctx):
        if j <= i:
            continue
        ra, ca = a.split('.')
        rb_, cb = b.split('.')
        cls = 'Same region' if ra == rb_ else ('Same cell type' if ca == cb else 'Different region and cell type')
        r, s = rbd.get((a, b), (np.nan, np.nan))
        ug, us = n_g[i] + n_g[j] - I_g[i, j], n_s[i] + n_s[j] - I_s[i, j]
        rows.append((a, b, cls, int(n_g[i]), int(n_g[j]), int(I_g[i, j]), I_g[i, j] / ug if ug > 0 and min(n_g[i], n_g[j]) > 0 else np.nan,
                     int(n_s[i]), int(n_s[j]), int(I_s[i, j]), I_s[i, j] / us if us > 0 and min(n_s[i], n_s[j]) > 0 else np.nan, r, s))
P = pd.DataFrame(rows, columns=['ctx1', 'ctx2', 'pair_class', 'egene1', 'egene2', 'egene_shared', 'egene_jaccard',
                                'esnp1', 'esnp2', 'esnp_shared', 'esnp_jaccard', 'rb', 'rb_se'])
P.to_csv(OUT + '/pair_sharing.tsv', sep='\t', index=False)
print(P.groupby('pair_class')[['egene_jaccard', 'esnp_jaccard', 'rb']].agg(['count', 'median']).round(3).to_string())

x = pd.read_csv(HO + '/plot_input_mea.txt', sep=r'\s+', quotechar='"')
x.columns = ['gene', 'region', 'cell', 'is_egene', 'expr']
print('expression rows', len(x), '| quantiles of expr', x.expr.quantile([0, 0.01, 0.25, 0.5, 0.75, 0.99, 1]).round(3).to_dict())
out = []
for (r, c), d in x.groupby(['region', 'cell']):
    e, n = d.expr[d.is_egene == 'y'], d.expr[d.is_egene == 'n']
    p = mannwhitneyu(e, n, alternative='two-sided').pvalue if len(e) >= 5 and len(n) >= 5 else np.nan
    out.append((r + '.' + c, r, c, len(e), len(n), e.median() if len(e) else np.nan, n.median(), p))
X = pd.DataFrame(out, columns=['ctx', 'region', 'cell', 'n_egene', 'n_non', 'median_egene', 'median_non', 'p_wilcoxon'])
X.to_csv(OUT + '/egene_expr_by_ctx.tsv', sep='\t', index=False)
ok = X.dropna(subset=['p_wilcoxon'])
print('contexts tested', len(ok), '| eGene median higher in', int((ok.median_egene > ok.median_non).sum()), '| P<0.05 in', int((ok.p_wilcoxon < 0.05).sum()))
# pooled densities for the violin: a thinned copy of the long table (every non-eGene row is not needed)
thin = pd.concat([x[x.is_egene == 'y'], x[x.is_egene == 'n'].sample(n=200000, random_state=1)])
thin[['region', 'cell', 'is_egene', 'expr']].to_csv(OUT + '/egene_expr_thin.tsv.gz', sep='\t', index=False)
print('FIG3_PREP_DONE')
