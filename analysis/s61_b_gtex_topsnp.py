"""分析 B 第二步：18 对「只在皮层下过 SMR+HEIDI+共定位」的病–基因，在 GTEx v8（eQTL Catalogue r7 *.all.tsv.gz，远程 tabix）
里查同一个 top SNP 对同一个基因的效应：对应脑组织 + cortex + frontal_cortex。方向：本稿 b_eQTL（对 A1）对齐到 GTEx alt。
admin2 上跑（只发几十次 HTTP range 请求）。in gtex_rep/B_subcortical_pairs.tsv（本地上传）；out gtex_rep/B_topsnp_gtex.tsv"""
import sys, urllib.request, zlib
import duckdb
import pandas as pd
import pysam

B = '/public/home/yucong/scbrain'
FTP = 'https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/{q}/{q}.all.tsv.gz'
Q = {'substantia_nigra': 'QTD000206', 'putamen': 'QTD000196', 'cerebellum': 'QTD000166', 'cerebellar_hemisphere': 'QTD000161',
     'spinal_cord': 'QTD000201', 'cortex': 'QTD000171', 'frontal_cortex': 'QTD000176'}
MAP = {'SNIG': ['substantia_nigra'], 'PUTM': ['putamen'], 'CRBL': ['cerebellum', 'cerebellar_hemisphere'], 'MEDU': ['spinal_cord'],
       'THAL': [], 'WHMT': []}
COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

pr = pd.read_csv(B + '/gtex_rep/B_subcortical_pairs.tsv', sep='\t')
smr = duckdb.sql(f"SELECT disease, probeID, ctx, topSNP, topSNP_chr, topSNP_bp, A1, A2, b_eQTL, p_eQTL FROM '{B}/db/tables/smr.parquet'").df()
v = duckdb.sql(f"SELECT rsid, chr, pos FROM '{B}/db/tables_allele_v3/variants.parquet'").df()
tbx = {}
hdr = None
rows = []
for r in pr.itertuples():
    s = smr[(smr.disease == r.disease) & (smr.probeID == r.gene_id) & (smr.ctx == r.best_ctx)].iloc[0]
    pos = v[v.rsid == s.topSNP]
    if pos.empty:
        print('no position', s.topSNP); continue
    c, p = str(pos.chr.iat[0]).replace('chr', ''), int(pos.pos.iat[0])
    tissues = sorted(set(t for reg in r.subcortical_ctx.split(',') for t in MAP[reg.split('.')[0]])) + ['cortex', 'frontal_cortex']
    for t in tissues:
        if t not in tbx:
            tbx[t] = pysam.TabixFile(FTP.format(q=Q[t]), index=FTP.format(q=Q[t]) + '.tbi')
            if hdr is None:                    # header line has no '#': read it from the first BGZF block
                req = urllib.request.Request(FTP.format(q=Q[t]), headers={'Range': 'bytes=0-65535'})
                raw = urllib.request.urlopen(req, timeout=120).read()
                hdr = zlib.decompressobj(31).decompress(raw).decode().split('\n')[0].split('\t')
                print('columns', hdr)
        hit = None
        for ln in tbx[t].fetch(c, p - 1, p):
            f = dict(zip(hdr, ln.split('\t')))
            if f['molecular_trait_id'].split('.')[0] == r.gene_id and int(f['position']) == p:
                hit = f; break
        row = dict(disease=r.disease, gene=r.gene, best_ctx=r.best_ctx, topSNP=s.topSNP, A1=s.A1, b_eQTL_ours=s.b_eQTL, tissue=t,
                   matched=t not in ('cortex', 'frontal_cortex'))
        if hit is None:
            row.update(found=False)
        else:
            ref, alt, beta, pv = hit['ref'], hit['alt'], float(hit['beta']), float(hit['pvalue'])
            a1 = s.A1
            if a1 == alt or COMP.get(a1) == alt:
                b_al = beta
            elif a1 == ref or COMP.get(a1) == ref:
                b_al = -beta
            else:
                b_al = float('nan')
            row.update(found=True, ref=ref, alt=alt, beta_gtex_alt=beta, beta_gtex_on_A1=b_al, p_gtex=pv,
                       same_direction=(b_al * s.b_eQTL > 0) if b_al == b_al else None)
        rows.append(row)
out = pd.DataFrame(rows)
out.to_csv(B + '/gtex_rep/B_topsnp_gtex.tsv', sep='\t', index=False)
pd.set_option('display.width', 250)
print(out.to_string(index=False))
