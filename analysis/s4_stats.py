"""eSNP 功能特征：第 5–10、12 步 + MAF 分箱的统计表。
读 esnp_func/sets、annot、tss，写 esnp_func/results/*.tsv。每步给 all / lnc / pc 三套 eSNP。
"""
import os, sys, glob, time
from multiprocessing import Pool
import numpy as np
import pandas as pd
from scipy import stats

S = '/public/home/yucong/scbrain'
O = os.environ.get('SCB_OUT', f'{S}/esnp_func')
R = f'{O}/results'
P = '/public/home/yucong/scbrain_in/paper2'
os.makedirs(R, exist_ok=True)
NPROC = int(os.environ.get('NPROC', '24'))
REGIONS = ['CRBL', 'FCTX', 'HIPP', 'MEDU', 'OCTX', 'PUTM', 'SNIG', 'TCTX', 'THAL', 'WHMT']
CELLS = ['Astro', 'Endo', 'ExN', 'Gran', 'InN', 'Micro', 'NPC', 'Olig', 'OPC', 'Perc', 'Purk']
t0 = time.time()


def log(*a):
    print(f'[{time.time() - t0:7.1f}s]', *a, flush=True)


def bh(p):
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    n = ok.sum()
    if n == 0:
        return q
    o = np.argsort(p[ok])
    r = p[ok][o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(r, 1)
    q[ok] = out
    return q


def fisher(a, b, c, d):
    """[[a, b], [c, d]] = [[E 在类内, E 不在], [N 在类内, N 不在]]"""
    a, b, c, d = int(a), int(b), int(c), int(d)
    if (a + b) == 0 or (c + d) == 0 or (a + c) == 0:
        return np.nan, np.nan, np.nan, np.nan
    _, p = stats.fisher_exact([[a, b], [c, d]])
    aa, bb, cc, dd = (a + .5, b + .5, c + .5, d + .5) if min(a, b, c, d) == 0 else (a, b, c, d)
    lor = np.log(aa * dd / (bb * cc)); se = np.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    return np.exp(lor), np.exp(lor - 1.96 * se), np.exp(lor + 1.96 * se), p


# ---------------- 集合 ----------------
var = pd.read_csv(f'{O}/sets/variants.tsv.gz', sep='\t', dtype={'chr': str})
m = np.load(f'{O}/sets/membership.npz', allow_pickle=True)
ctx = [str(x) for x in m['ctx']]
e3, n4, pre, lnc, pc = m['e3'], m['n4'], m['pre'], m['lnc'], m['pc']
eg = np.load(f'{O}/sets/membership_eg.npz')['eg']   # 该组合里是 q<0.05 的 eGene 的 eSNP
nU, nC = e3.shape
assert nU == len(var) and eg.shape == e3.shape
log('universe', nU, 'contexts', nC)

GROUPS = {'WB': list(range(nC))}
for c in CELLS:
    GROUPS['cell:' + c] = [j for j, x in enumerate(ctx) if x.split('.')[1] == c]
for r in REGIONS:
    GROUPS['region:' + r] = [j for j, x in enumerate(ctx) if x.split('.')[0] == r]
CTXG = {'ctx:' + x: [j] for j, x in enumerate(ctx)}
MASK = {}


def masks(g):
    if g in MASK:
        return MASK[g]
    J = (GROUPS.get(g) or CTXG[g])
    e = e3[:, J]
    E = e.any(1)
    El = (e & lnc[:, J]).any(1)
    Ep = (e & pc[:, J]).any(1)
    Eg = (e & eg[:, J]).any(1)
    N = n4[:, J].any(1) & ~pre[:, J].any(1) & ~E
    # all = 交接来的集合；eg = 只留 q<0.05 的 eGene 的 eSNP；lnc / pc = 按关联基因的类型
    MASK[g] = {'all': E, 'eg': Eg, 'lnc': El, 'pc': Ep, 'N': N}
    return MASK[g]


rows = []
for g in list(GROUPS) + list(CTXG):
    mk = masks(g)
    rows.append(dict(group=g, n_esnp_all=int(mk['all'].sum()), n_esnp_eg=int(mk['eg'].sum()), n_esnp_lnc=int(mk['lnc'].sum()),
                     n_esnp_pc=int(mk['pc'].sum()), n_non=int(mk['N'].sum())))
sets_tab = pd.DataFrame(rows)
sets_tab.to_csv(f'{R}/step5_sets_by_group.tsv', sep='\t', index=False)
log('step5 written', sets_tab.iloc[0].to_dict())

var['idx'] = np.arange(nU)
rs_index = pd.Index(var['rsid'])

# ---------------- 注释载入 ----------------
APC = ['apc_conservation', 'apc_conservation_v2', 'apc_epigenetics', 'apc_epigenetics_active',
       'apc_epigenetics_repressed', 'apc_epigenetics_transcription', 'apc_local_nucleotide_diversity',
       'apc_local_nucleotide_diversity_v2', 'apc_local_nucleotide_diversity_v3', 'apc_mappability',
       'apc_micro_rna', 'apc_mutation_density', 'apc_protein_function', 'apc_proximity_to_coding',
       'apc_proximity_to_coding_v2', 'apc_proximity_to_tsstes', 'apc_transcription_factor']
FCOLS = ['chr', 'pos', 'ref', 'alt'] + APC + ['rsid_f', 'category', 'linsight', 'cadd_phred']
FEAT = {}

fl = sorted(glob.glob(f'{O}/annot/favor/chr*.tsv'))
fav = pd.concat([pd.read_csv(f, sep='\t', header=None, names=FCOLS, na_values=['\\N'], dtype={'chr': str, 'ref': str, 'alt': str})
                 for f in fl if os.path.getsize(f) > 0], ignore_index=True)
log('favor rows', len(fav), 'chromosomes', sorted(fav['chr'].unique(), key=lambda x: int(x)))
fav = fav.merge(var[['idx', 'chr', 'pos', 'a1', 'a2']], on=['chr', 'pos'], how='inner')
fav['match'] = ((fav.ref == fav.a1) & (fav.alt == fav.a2)) | ((fav.ref == fav.a2) & (fav.alt == fav.a1))
fav = fav.sort_values(['idx', 'match'], ascending=[True, False]).drop_duplicates('idx')
log('favor matched variants', len(fav), 'allele-matched', int(fav['match'].sum()))
for c in APC + ['linsight', 'cadd_phred']:
    a = np.full(nU, np.nan, dtype=np.float32); a[fav['idx'].values] = fav[c].values.astype(np.float32); FEAT[c] = a
cat = np.full(nU, None, dtype=object); cat[fav['idx'].values] = fav['category'].values
has_favor = np.zeros(nU, bool); has_favor[fav['idx'].values] = True
del fav

el = sorted(glob.glob(f'{O}/annot/enformer/chr*.tsv'))
PCS = [f'PC{i}' for i in range(20)]
enf = pd.concat([pd.read_csv(f, sep='\t', header=None, names=['chrom', 'pos19', 'rsid', 'ref', 'alt'] + PCS, dtype={'ref': str, 'alt': str})
                 for f in el if os.path.getsize(f) > 0], ignore_index=True)
enf = enf.merge(var[['idx', 'rsid', 'a1', 'a2']], on='rsid', how='inner')
enf['match'] = ((enf.ref == enf.a1) & (enf.alt == enf.a2)) | ((enf.ref == enf.a2) & (enf.alt == enf.a1))
enf = enf.sort_values(['idx', 'match'], ascending=[True, False]).drop_duplicates('idx')
log('enformer matched variants', len(enf), 'allele-matched', int(enf['match'].sum()))
for c in PCS:
    a = np.full(nU, np.nan, dtype=np.float32); a[enf['idx'].values] = enf[c].values.astype(np.float32); FEAT['enformer_' + c] = a
l2 = np.full(nU, np.nan, dtype=np.float32); l2[enf['idx'].values] = np.sqrt((enf[PCS].values ** 2).sum(1)).astype(np.float32)
FEAT['enformer_L2'] = l2
del enf

cc = pd.read_csv(f'{O}/annot/ccre_hits.tsv', sep='\t', header=None, names=['rsid', 'cls']).drop_duplicates('rsid')
ci = rs_index.get_indexer(cc['rsid'])
CCRE = {}
for name, pat in [('any_cCRE', None), ('PLS', 'PLS'), ('pELS', 'pELS'), ('dELS', 'dELS'),
                  ('DNase-H3K4me3', 'DNase-H3K4me3'), ('CTCF-only', 'CTCF-only'), ('CTCF-bound', 'CTCF-bound')]:
    a = np.zeros(nU, bool)
    sel = np.ones(len(cc), bool) if pat is None else cc['cls'].str.contains(pat, regex=False).values
    a[ci[sel & (ci >= 0)]] = True
    CCRE[name] = a
log('ccre variants', int(CCRE['any_cCRE'].sum()))

# ---------------- 通用检验 ----------------
ESETS = ['all', 'eg', 'lnc', 'pc']


def fisher_block(flag_dict, groups, valid=None, esets=ESETS):
    out = []
    for g in groups:
        mk = masks(g)
        N = mk['N'] if valid is None else (mk['N'] & valid)
        nN = int(N.sum())
        for es in esets:
            E = mk[es] if valid is None else (mk[es] & valid)
            nE = int(E.sum())
            for name, f in flag_dict.items():
                a = int((E & f).sum()); c = int((N & f).sum())
                orr, lo, hi, p = fisher(a, nE - a, c, nN - c)
                out.append(dict(group=g, eset=es, feature=name, n_esnp=nE, n_non=nN, esnp_in=a, non_in=c,
                                pct_esnp=100 * a / nE if nE else np.nan, pct_non=100 * c / nN if nN else np.nan,
                                OR=orr, OR_lo=lo, OR_hi=hi, p=p))
    df = pd.DataFrame(out)
    df['q'] = df.groupby(['group', 'eset'])['p'].transform(bh)
    return df


def mwu_task(t):
    g, es, feat = t
    mk = MASK[g]
    v = FEAT[feat]
    x = v[mk[es]]; y = v[mk['N']]
    x = x[~np.isnan(x)]; y = y[~np.isnan(y)]
    if len(x) < 5 or len(y) < 5:
        return dict(group=g, eset=es, feature=feat, n_esnp=len(x), n_non=len(y))
    u, p = stats.mannwhitneyu(x, y, alternative='two-sided')
    return dict(group=g, eset=es, feature=feat, n_esnp=len(x), n_non=len(y), median_esnp=float(np.median(x)),
                median_non=float(np.median(y)), mean_esnp=float(x.mean()), mean_non=float(y.mean()),
                auc=float(u) / (len(x) * len(y)), p=p)


def mwu_block(feats, groups, esets=ESETS):
    for g in groups:
        masks(g)
    tasks = [(g, es, f) for g in groups for es in esets for f in feats]
    with Pool(NPROC) as pool:
        res = pool.map(mwu_task, tasks, chunksize=4)
    df = pd.DataFrame(res)
    df['q'] = df.groupby(['group', 'eset'])['p'].transform(bh)
    return df


G_MAIN = list(GROUPS)
G_CTX = list(CTXG)

# ---------------- 第 7 步：区域 11 类 ----------------
CLASSES = ['ncRNA_intronic', 'ncRNA_exonic', 'ncRNA_splicing', 'upstream', 'downstream', 'UTR5', 'UTR3',
           'exonic', 'intronic', 'splicing', 'intergenic']
cat_s = pd.Series(cat)
single = cat_s.isin(CLASSES).values
pd.DataFrame({'category': cat_s[has_favor].value_counts()}).to_csv(f'{R}/step7_category_counts_universe.tsv', sep='\t')
flags7 = {c: (cat_s == c).values for c in CLASSES}
s7 = fisher_block(flags7, G_MAIN + G_CTX, valid=single)
s7.to_csv(f'{R}/step7_region_class.tsv', sep='\t', index=False)
log('step7 done; multi-class or unannotated among FAVOR-annotated:', int((has_favor & ~single).sum()))

# ---------------- 第 10 步：cCRE ----------------
s10c = fisher_block(CCRE, G_MAIN + G_CTX)
s10c.to_csv(f'{R}/step10_ccre.tsv', sep='\t', index=False)
log('step10 ccre done')

# ---------------- 第 8、9 步：Enformer、FAVOR aPC ----------------
enf_feats = ['enformer_' + c for c in PCS] + ['enformer_L2']
s8 = mwu_block(enf_feats, G_MAIN)
s8.to_csv(f'{R}/step8_enformer.tsv', sep='\t', index=False)
s8c = mwu_block(['enformer_L2'], G_CTX, esets=['all'])
s8c.to_csv(f'{R}/step8_enformer_L2_by_ctx.tsv', sep='\t', index=False)
log('step8 done')
s9 = mwu_block(APC + ['linsight', 'cadd_phred'], G_MAIN)
s9.to_csv(f'{R}/step9_favor_apc.tsv', sep='\t', index=False)
log('step9 done')

# 画分布图用的分位数（全脑）
qs = [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]
qrows = []
mk = masks('WB')
for f in enf_feats + APC + ['linsight', 'cadd_phred']:
    for es in ESETS + ['N']:
        v = FEAT[f][mk[es]]; v = v[~np.isnan(v)]
        if len(v):
            qrows.append(dict(feature=f, set=es, n=len(v), **{f'q{int(q * 100):02d}': float(np.quantile(v, q)) for q in qs}))
pd.DataFrame(qrows).to_csv(f'{R}/step8_9_quantiles_WB.tsv', sep='\t', index=False)

# ---------------- 第 10 步：ReMap ----------------
tf_codes, tf_idx = [], []
tf_names = {}
for ch in pd.read_csv(f'{O}/annot/remap_hits.tsv.gz', sep='\t', header=None, names=['rsid', 'tf'], chunksize=5_000_000):
    i = rs_index.get_indexer(ch['rsid'])
    ok = i >= 0
    codes = ch['tf'].map(lambda x: tf_names.setdefault(x, len(tf_names))).values
    tf_idx.append(i[ok].astype(np.int32)); tf_codes.append(codes[ok].astype(np.int32))
tf_idx = np.concatenate(tf_idx); tf_codes = np.concatenate(tf_codes)
TF = [None] * len(tf_names)
for k, v in tf_names.items():
    TF[v] = k
nTF = len(TF)
log('remap pairs', len(tf_idx), 'TFs', nTF)
any_tf = np.zeros(nU, bool); any_tf[tf_idx] = True
out = []
for g in G_MAIN:
    mk = masks(g)
    cN = np.bincount(tf_codes[mk['N'][tf_idx]], minlength=nTF); nN = int(mk['N'].sum())
    for es in (ESETS if g == 'WB' else ['all']):
        E = mk[es]; nE = int(E.sum())
        cE = np.bincount(tf_codes[E[tf_idx]], minlength=nTF)
        a = int((E & any_tf).sum()); c = int((mk['N'] & any_tf).sum())
        orr, lo, hi, p = fisher(a, nE - a, c, nN - c)
        out.append(dict(group=g, eset=es, feature='any_TF_peak', n_esnp=nE, n_non=nN, esnp_in=a, non_in=c, OR=orr, OR_lo=lo, OR_hi=hi, p=p))
        for k in range(nTF):
            orr, lo, hi, p = fisher(cE[k], nE - cE[k], cN[k], nN - cN[k])
            out.append(dict(group=g, eset=es, feature=TF[k], n_esnp=nE, n_non=nN, esnp_in=int(cE[k]), non_in=int(cN[k]),
                            OR=orr, OR_lo=lo, OR_hi=hi, p=p))
s10t = pd.DataFrame(out)
s10t['pct_esnp'] = 100 * s10t.esnp_in / s10t.n_esnp; s10t['pct_non'] = 100 * s10t.non_in / s10t.n_non
s10t['q'] = s10t.groupby(['group', 'eset'])['p'].transform(bh)
s10t.to_csv(f'{R}/step10_remap_tf.tsv', sep='\t', index=False)
log('step10 remap done')

# 脑相关生物型（含锚文用的 brain-prefrontal-cortex SIN3A / OLIG2）
bpath = f'{O}/annot/remap_brain_hits.tsv'
if os.path.exists(bpath) and os.path.getsize(bpath) > 0:
    bh_ = pd.read_csv(bpath, sep='\t', header=None, names=['rsid', 'tf', 'biotype'])
    bh_['biotype'] = bh_['biotype'].str.split(',')
    bh_ = bh_.explode('biotype')
    bh_ = bh_[bh_['biotype'].str.lower().str.contains('brain|cortex|neuro|astro|glia|hippocamp|cerebell|microgl|oligodend')]
    bh_['key'] = bh_['tf'] + ':' + bh_['biotype']
    bh_ = bh_.drop_duplicates(['rsid', 'key'])
    flagsb = {}
    for key, sub in bh_.groupby('key'):
        ii = rs_index.get_indexer(sub['rsid'])
        ii = ii[ii >= 0]
        if len(ii) < 50:
            continue
        a = np.zeros(nU, bool); a[ii] = True
        flagsb[key] = a
    s10b = fisher_block(flagsb, G_MAIN)
    s10b.to_csv(f'{R}/step10_remap_brain_biotypes.tsv', sep='\t', index=False)
    log('step10 brain biotypes done', len(flagsb))

# ---------------- 第 6 步：TSS 距离 ----------------
tested = pd.read_csv(f'{O}/tss/tested_pairs_bins.tsv', sep='\t')
tested['cell'] = tested['ctx'].str.split('.').str[1]; tested['region'] = tested['ctx'].str.split('.').str[0]
pairs = pd.read_csv(f'{O}/sets/esnp_pairs.tsv.gz', sep='\t', usecols=['ctx', 'dist', 'pval', 'biotype'])
pairs['bin'] = np.floor(pairs['dist'] / 10000.0).astype(int)
pairs['cell'] = pairs['ctx'].str.split('.').str[1]; pairs['region'] = pairs['ctx'].str.split('.').str[0]
pairs['lnc'] = pairs['biotype'] != 'protein_coding'
log('tested pairs total', int(tested['n'].sum()), 'eSNP pairs', len(pairs))


def tss_scope(tt, pp, scope):
    t = tt.groupby('bin')['n'].sum().rename('n_tested')
    e = pp.groupby('bin').agg(n_esnp_pairs=('dist', 'size'), n_esnp_pairs_lnc=('lnc', 'sum'))
    d = pd.concat([t, e], axis=1).fillna(0).reset_index()
    d['n_esnp_pairs_pc'] = d['n_esnp_pairs'] - d['n_esnp_pairs_lnc']
    d['pct_esnp'] = 100 * d['n_esnp_pairs'] / d['n_tested'].replace(0, np.nan)
    d.insert(0, 'scope', scope)
    return d


tss = [tss_scope(tested, pairs, 'WB')]
for c in CELLS:
    tss.append(tss_scope(tested[tested.cell == c], pairs[pairs.cell == c], 'cell:' + c))
for r in REGIONS:
    tss.append(tss_scope(tested[tested.region == r], pairs[pairs.region == r], 'region:' + r))
pd.concat(tss).to_csv(f'{R}/step6_tss_bins.tsv', sep='\t', index=False)
tested.groupby(['bin', 'pbin'])['n'].sum().reset_index().to_csv(f'{R}/step6_tss_density_WB.tsv', sep='\t', index=False)

# 按信号强度：P 小于阈值的对里，落在 TSS 近旁的比例，对全部被检验的对
srows = []
WIN = {'10kb': (-1, 0), '50kb': (-5, 4), '100kb': (-10, 9), '250kb': (-25, 24)}
THR = {'1e-2': 4, '1e-3': 6, '1e-4': 8, '1e-5': 10, '1e-6': 12, '1e-8': 16, '1e-10': 20}
for scope, tt, pp in [('WB', tested, pairs)] + [('cell:' + c, tested[tested.cell == c], pairs[pairs.cell == c]) for c in CELLS]:
    tot = int(tt['n'].sum())
    for wn, (lo_, hi_) in WIN.items():
        inw = tt['bin'].between(lo_, hi_)
        all_in = int(tt.loc[inw, 'n'].sum())
        for tn, pb in THR.items():
            sig = tt['pbin'] >= pb
            a = int(tt.loc[sig & inw, 'n'].sum()); b = int(tt.loc[sig & ~inw, 'n'].sum())
            c_ = all_in - a; d_ = tot - all_in - b
            orr = (a * d_) / (b * c_) if b and c_ else np.nan
            srows.append(dict(scope=scope, signal='P<' + tn, window=wn, sig_in=a, sig_out=b, rest_in=c_, rest_out=d_,
                              pct_sig_in=100 * a / (a + b) if a + b else np.nan, pct_all_in=100 * all_in / tot, OR=orr))
        ein = int(pp['bin'].between(lo_, hi_).sum()); eout = len(pp) - ein
        c_ = all_in - ein; d_ = tot - all_in - eout
        srows.append(dict(scope=scope, signal='eSNP', window=wn, sig_in=ein, sig_out=eout, rest_in=c_, rest_out=d_,
                          pct_sig_in=100 * ein / len(pp) if len(pp) else np.nan, pct_all_in=100 * all_in / tot,
                          OR=(ein * d_) / (eout * c_) if eout and c_ else np.nan))
pd.DataFrame(srows).to_csv(f'{R}/step6_tss_by_strength.tsv', sep='\t', index=False)
log('step6 done')

# ---------------- MAF 分箱 ----------------
allv = pd.read_csv(f'{O}/sets/all_variants_maf.tsv.gz', sep='\t', dtype={'chr': str})
gb = pd.read_csv(f'{P}/single-cell-eqtl/data/expression/110/CRBL.Astro.bed', sep='\t', usecols=[0, 1, 2, 3])
gb.columns = ['chr', 'start', 'end', 'gene']
gb['chr'] = gb['chr'].astype(str).str.replace('chr', '', regex=False)
cis = np.zeros(len(allv), bool)
for ch, sub in allv.groupby('chr'):
    tssp = np.sort(gb.loc[gb.chr == ch, 'end'].values)
    if len(tssp) == 0:
        continue
    pos = sub['pos'].values
    k = np.searchsorted(tssp, pos)
    dl = np.abs(pos - tssp[np.clip(k - 1, 0, len(tssp) - 1)]); dr = np.abs(tssp[np.clip(k, 0, len(tssp) - 1)] - pos)
    cis[sub.index.values] = np.minimum(dl, dr) <= 1_000_000
edges = np.round(np.arange(0.05, 0.5001, 0.05), 2)
labels = [f'{edges[i]:.2f}-{edges[i + 1]:.2f}' for i in range(len(edges) - 1)]


def mafbin(x):
    return np.clip(np.digitize(x, edges[1:-1]), 0, len(labels) - 1)


cis_b = np.bincount(mafbin(allv.loc[cis, 'MAF'].values), minlength=len(labels))
vb = mafbin(var['MAF'].fillna(-1).values); vok = var['MAF'].notna().values
mrows = []
for g in G_MAIN + G_CTX:
    mk = masks(g)
    for es in ESETS:
        E = mk[es] & vok; N = mk['N'] & vok
        eb = np.bincount(vb[E], minlength=len(labels)); nb = np.bincount(vb[N], minlength=len(labels))
        for i, lab in enumerate(labels):
            orr, lo, hi, p = fisher(eb[i], eb.sum() - eb[i], nb[i], nb.sum() - nb[i])
            mrows.append(dict(group=g, eset=es, maf_bin=lab, n_cis=int(cis_b[i]), n_esnp=int(eb[i]), n_non=int(nb[i]),
                              prop_esnp_in_cis=eb[i] / cis_b[i], prop_non_in_cis=nb[i] / cis_b[i],
                              OR=orr, OR_lo=lo, OR_hi=hi, p=p))
pd.DataFrame(mrows).to_csv(f'{R}/maf_bins.tsv', sep='\t', index=False)
log('maf done; cis variants', int(cis.sum()), 'of', len(allv))

# ---------------- 第 12 步：丰度 × 富集 ----------------
frac = pd.read_csv(f'{S}/local_inputs/cell_frac_mean.tsv', sep='\t')
fr_ctx = dict(zip('ctx:' + frac['ctx'], frac['mean_frac']))
fr_cell = ('cell:' + frac['cell']).to_frame('group').assign(mean_frac=frac['mean_frac']).groupby('group')['mean_frac'].mean().to_dict()
long = pd.concat([s7.assign(step='7_region'), s10c.assign(step='10_cCRE')], ignore_index=True)
long = long[long.eset == 'all']
enfc = s8c.rename(columns={'auc': 'OR'}).assign(step='8_enformer')
arows = []
for level, fr in [('ctx', fr_ctx), ('cell', fr_cell)]:
    sub = long[long.group.isin(fr.keys())]
    for (step, feat), d in sub.groupby(['step', 'feature']):
        d = d.dropna(subset=['p', 'OR'])
        d = d[(d.esnp_in + d.non_in) > 0]
        if len(d) < 6:
            continue
        x = d.group.map(fr).values
        y1 = -np.log10(np.clip(d.p.values, 1e-300, 1)); y2 = np.log(np.clip(d.OR.values, 1e-6, 1e6)); ne = d.n_esnp.values
        r1, p1 = stats.spearmanr(x, y1); r2, p2 = stats.spearmanr(x, y2); r3, p3 = stats.spearmanr(ne, y1); r4, p4 = stats.spearmanr(x, ne)
        arows.append(dict(level=level, step=step, feature=feat, n_points=len(d), rho_frac_vs_neglog10p=r1, p1=p1,
                          rho_frac_vs_logOR=r2, p2=p2, rho_nesnp_vs_neglog10p=r3, p3=p3, rho_frac_vs_nesnp=r4, p4=p4))
if len(enfc):
    d = enfc[enfc.group.isin(fr_ctx.keys())].dropna(subset=['p'])
    if len(d) >= 6:
        x = d.group.map(fr_ctx).values; y1 = -np.log10(np.clip(d.p.values, 1e-300, 1))
        r1, p1 = stats.spearmanr(x, y1); r2, p2 = stats.spearmanr(x, d.OR.values); r3, p3 = stats.spearmanr(d.n_esnp.values, y1)
        arows.append(dict(level='ctx', step='8_enformer', feature='enformer_L2(AUC)', n_points=len(d), rho_frac_vs_neglog10p=r1, p1=p1,
                          rho_frac_vs_logOR=r2, p2=p2, rho_nesnp_vs_neglog10p=r3, p3=p3))
pd.DataFrame(arows).to_csv(f'{R}/step12_abundance_vs_enrichment.tsv', sep='\t', index=False)
log('step12 done')

cov = dict(universe=nU, favor_annotated=int(has_favor.sum()), favor_single_class=int(single.sum()),
           enformer_annotated=int((~np.isnan(FEAT['enformer_L2'])).sum()), ccre_any=int(CCRE['any_cCRE'].sum()),
           remap_any=int(any_tf.sum()), wb_esnp=int(masks('WB')['all'].sum()), wb_non=int(masks('WB')['N'].sum()),
           wb_esnp_favor=int((masks('WB')['all'] & has_favor).sum()), wb_esnp_enformer=int((masks('WB')['all'] & ~np.isnan(FEAT['enformer_L2'])).sum()),
           wb_esnp_chr6=int((masks('WB')['all'] & (var['chr'] == '6').values).sum()), wb_esnp_nopos=int((masks('WB')['all'] & var['pos'].isna().values).sum()))
pd.Series(cov).to_csv(f'{R}/coverage.tsv', sep='\t', header=False)
log('ALL DONE', cov)
