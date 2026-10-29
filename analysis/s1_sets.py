"""第 5 步的集合整理：把贺奕捷现成的 eSNP / non-eSNP 名单读成统一的表和成员矩阵。
只读 scbrain_in，产物写 esnp_func/sets/。
"""
import os, sys, gzip, time
import numpy as np
import pandas as pd

P = '/public/home/yucong/scbrain_in/paper2'
D = f'{P}/single-cell-eqtl/data'
OUT = os.environ.get('SCB_OUT', '/public/home/yucong/scbrain/esnp_func') + '/sets'
ESNP = os.environ.get('ESNP_DIR', f'{D}/esnp')   # <组合>-esnp.txt 和 <组合>-esnp3.txt 所在目录
os.makedirs(OUT, exist_ok=True)
print('ESNP_DIR', ESNP, 'OUT', OUT, flush=True)

REGIONS = ['CRBL', 'FCTX', 'HIPP', 'MEDU', 'OCTX', 'PUTM', 'SNIG', 'TCTX', 'THAL', 'WHMT']
CELLS = ['Astro', 'Endo', 'ExN', 'Gran', 'InN', 'Micro', 'NPC', 'Olig', 'OPC', 'Perc', 'Purk']
CTX = [f'{r}.{c}' for r in REGIONS for c in CELLS]
t0 = time.time()


def log(*a):
    print(f'[{time.time() - t0:7.1f}s]', *a, flush=True)


def read_ids(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return []
    with open(path) as f:
        ids = [l.split()[0] for l in f if l.strip()]
    return [x for x in ids if x != 'SNP']


# 基因类型
gt = pd.read_csv(f'{D}/expression/Homo_sapiens.GRCh38.82.chr.gtf.txt', sep='\t', dtype=str)
biotype = dict(zip(gt['a5'], gt['a10']))
log('genes with biotype', len(biotype), pd.Series(list(biotype.values())).value_counts().head(8).to_dict())

# 第一遍：全集 U
U = set()
n_e3, n_n4 = {}, {}
for c in CTX:
    e3 = read_ids(f'{ESNP}/{c}-esnp3.txt')
    n4 = read_ids(f'{D}/non-esnp/{c}-non-esnp4.txt')
    n_e3[c], n_n4[c] = len(e3), len(n4)
    U.update(e3); U.update(n4)
U = sorted(U)
idx = {r: i for i, r in enumerate(U)}
log('universe', len(U))

nU, nC = len(U), len(CTX)
mem_e3 = np.zeros((nU, nC), dtype=bool)
mem_n4 = np.zeros((nU, nC), dtype=bool)
mem_pre = np.zeros((nU, nC), dtype=bool)   # 该组合里是（clump 前）eSNP
mem_lnc = np.zeros((nU, nC), dtype=bool)   # 该组合里关联到 lncRNA eGene
mem_pc = np.zeros((nU, nC), dtype=bool)    # 该组合里关联到蛋白编码 eGene

for j, c in enumerate(CTX):
    for r in read_ids(f'{ESNP}/{c}-esnp3.txt'):
        mem_e3[idx[r], j] = True
    for r in read_ids(f'{D}/non-esnp/{c}-non-esnp4.txt'):
        mem_n4[idx[r], j] = True
log('membership filled')

# eSNP—基因对
rows_summary = []
pre_any = set()
with gzip.open(f'{OUT}/esnp_pairs.tsv.gz', 'wt') as w:
    w.write('ctx\tgene\trsid\tdist\tpval\tslope\tse\taf\tthreshold\tbiotype\n')
    for j, c in enumerate(CTX):
        path = f'{ESNP}/{c}-esnp.txt'
        npair = 0; genes = set(); snps = set(); lncg = set()
        if os.path.exists(path):
            with open(path) as f:
                for line in f:
                    a = line.split()
                    if len(a) < 27:
                        continue
                    g, r = a[0], a[1]
                    bt = biotype.get(g, 'NA')
                    w.write(f'{c}\t{g}\t{r}\t{a[2]}\t{a[6]}\t{a[7]}\t{a[8]}\t{a[3]}\t{a[26]}\t{bt}\n')
                    npair += 1; genes.add(g); snps.add(r)
                    is_lnc = bt != 'protein_coding'
                    if is_lnc:
                        lncg.add(g)
                    i = idx.get(r)
                    if i is not None:
                        mem_pre[i, j] = True
                        if is_lnc:
                            mem_lnc[i, j] = True
                        else:
                            mem_pc[i, j] = True
        pre_any |= snps
        rows_summary.append(dict(ctx=c, region=c.split('.')[0], cell=c.split('.')[1], n_pairs=npair,
                                 n_egene_with_esnp=len(genes), n_lnc_egene=len(lncg), n_esnp_preclump=len(snps),
                                 n_esnp3=n_e3[c], n_nonesnp4=n_n4[c],
                                 n_esnp3_lnc=int((mem_e3[:, j] & mem_lnc[:, j]).sum()),
                                 n_esnp3_pc=int((mem_e3[:, j] & mem_pc[:, j]).sum())))
log('pairs written')
pd.DataFrame(rows_summary).to_csv(f'{OUT}/context_summary.tsv', sep='\t', index=False)

# 位置和频率
bim = pd.read_csv(f'{D}/genotype/out.bim', sep='\t', header=None, names=['chr', 'rsid', 'cm', 'pos', 'a1', 'a2'],
                  dtype={'chr': str, 'rsid': str, 'pos': np.int64, 'a1': str, 'a2': str})
maf = pd.read_csv(f'{D}/genotype/MAF.txt', sep=r'\s+', usecols=['ID', 'MAF'])
bim = bim.merge(maf.rename(columns={'ID': 'rsid'}), on='rsid', how='left')
bim['chr'] = bim['chr'].str.replace('chr', '', regex=False)
bim.drop_duplicates('rsid', inplace=True)
var = pd.DataFrame({'rsid': U}).merge(bim[['rsid', 'chr', 'pos', 'a1', 'a2', 'MAF']], on='rsid', how='left')
var['n_ctx_e3'] = mem_e3.sum(1)
var['n_ctx_n4'] = mem_n4.sum(1)
var['pre_any'] = var['rsid'].isin(pre_any)
e_any = mem_e3.any(1)
var['wb_esnp'] = e_any
var['wb_non'] = mem_n4.any(1) & ~e_any & ~var['pre_any'].values
var['wb_esnp_lnc'] = (mem_e3 & mem_lnc).any(1)
var['wb_esnp_pc'] = (mem_e3 & mem_pc).any(1)
var.to_csv(f'{OUT}/variants.tsv.gz', sep='\t', index=False)
log('variants', len(var), 'no position:', int(var['pos'].isna().sum()),
    'wb_esnp', int(var.wb_esnp.sum()), 'wb_non', int(var.wb_non.sum()),
    'wb_esnp_lnc', int(var.wb_esnp_lnc.sum()), 'wb_esnp_pc', int(var.wb_esnp_pc.sum()))

np.savez_compressed(f'{OUT}/membership.npz', rsid=np.array(U), ctx=np.array(CTX),
                    e3=mem_e3, n4=mem_n4, pre=mem_pre, lnc=mem_lnc, pc=mem_pc)

ok = var.dropna(subset=['pos']).copy()
ok['pos'] = ok['pos'].astype(np.int64)
ok['start'] = ok['pos'] - 1
ok['c'] = 'chr' + ok['chr']
ok[['c', 'start', 'pos', 'rsid']].to_csv(f'{OUT}/universe_hg38.bed', sep='\t', header=False, index=False)
(ok['chr'] + '-' + ok['pos'].astype(str)).to_csv(f'{OUT}/universe_chrpos.txt', header=False, index=False)
ok['rsid'].to_csv(f'{OUT}/universe_rsid.txt', header=False, index=False)

# 全部被检验变异的位置和频率（MAF 分箱的分母用）
bim[['rsid', 'chr', 'pos', 'MAF']].to_csv(f'{OUT}/all_variants_maf.tsv.gz', sep='\t', index=False)
log('done')
