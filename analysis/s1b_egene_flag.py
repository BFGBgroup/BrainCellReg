"""补一张成员矩阵：每个组合里，该变异是否作为 q<0.05 的 eGene 的 eSNP 出现。
背景：交接来的 <组合>-esnp.txt 只按「nominal P < 该基因阈值」取对，没有再限 q<0.05，
所以含有一部分来自非 eGene 的对。这里不改原集合，只加标记，后面两套都算。
"""
import os, gzip, time
import numpy as np
import pandas as pd

P = '/public/home/yucong/scbrain_in/paper2'
D = f'{P}/single-cell-eqtl/data'
OUT = os.environ.get('SCB_OUT', '/public/home/yucong/scbrain/esnp_func') + '/sets'
ESNP = os.environ.get('ESNP_DIR', f'{D}/esnp')
t0 = time.time()

m = np.load(f'{OUT}/membership.npz', allow_pickle=True)
U = [str(x) for x in m['rsid']]
CTX = [str(x) for x in m['ctx']]
e3 = m['e3']
idx = {r: i for i, r in enumerate(U)}
eg = np.zeros(e3.shape, dtype=bool)

rows = []
with gzip.open(f'{OUT}/esnp_pairs_q.tsv.gz', 'wt') as w:
    w.write('ctx\tgene\trsid\tdist\tpval\tqval\n')
    for j, c in enumerate(CTX):
        path = f'{ESNP}/{c}-esnp.txt'
        n_eg = n_non = 0; g_eg = set(); g_non = set(); s_eg = set(); s_non = set()
        if os.path.exists(path):
            with open(path) as f:
                for line in f:
                    a = line.split()
                    if len(a) < 27:
                        continue
                    q = float(a[25])
                    w.write(f'{c}\t{a[0]}\t{a[1]}\t{a[2]}\t{a[6]}\t{a[25]}\n')
                    if q < 0.05:
                        n_eg += 1; g_eg.add(a[0]); s_eg.add(a[1])
                        i = idx.get(a[1])
                        if i is not None:
                            eg[i, j] = True
                    else:
                        n_non += 1; g_non.add(a[0]); s_non.add(a[1])
        col = e3[:, j]
        rows.append(dict(ctx=c, pairs_egene=n_eg, pairs_non_egene=n_non, genes_egene=len(g_eg), genes_non_egene=len(g_non),
                         snps_egene=len(s_eg), snps_only_non_egene=len(s_non - s_eg),
                         esnp3=int(col.sum()), esnp3_from_egene=int((col & eg[:, j]).sum()),
                         esnp3_only_non_egene=int((col & ~eg[:, j]).sum())))
df = pd.DataFrame(rows)
df.to_csv(f'{OUT}/egene_restriction_summary.tsv', sep='\t', index=False)
np.savez_compressed(f'{OUT}/membership_eg.npz', eg=eg)
tot = df.sum(numeric_only=True)
print('pairs: eGene %d, non-eGene %d (%.1f%%)' % (tot.pairs_egene, tot.pairs_non_egene, 100 * tot.pairs_non_egene / (tot.pairs_egene + tot.pairs_non_egene)))
print('independent eSNP (sum over contexts): %d, from eGene %d, only non-eGene %d (%.1f%%)' % (
    tot.esnp3, tot.esnp3_from_egene, tot.esnp3_only_non_egene, 100 * tot.esnp3_only_non_egene / tot.esnp3))
wb_all = e3.any(1); wb_eg = (e3 & eg).any(1)
print('whole-brain union: all %d, from eGene %d' % (wb_all.sum(), wb_eg.sum()))
print('done %.0fs' % (time.time() - t0))
