"""scGWAS 结果汇总：每个（组合、GWAS）读 4.R 写出的 sig_module.txt，数显著模块、涉及的 DER-22 细胞类型和基因。
同时把贺奕捷原来的结果按同样办法汇总一遍，放在旁边对照。用法: python s12_scgwas_summary.py
"""
import os, glob
import pandas as pd

NEW = '/public/home/yucong/scbrain/scgwas10/out'
HIS = '/public/home/yucong/scbrain_in/paper2/eqtl-gwas/scGWAS/run/out'
OUT = '/public/home/yucong/scbrain/scgwas10'
NAME = {'GCST90432156': 'MG', 'GCST003566': 'MS', 'GCST90726344': 'MDD', 'GCST009325': 'PD', 'GCST90018844': 'FebrileSeizures',
        'GCST90128471': 'SCZ', 'GCST011179': 'CUD', 'GCST013196': 'AD', 'GCST90014122': 'LacunarStroke', 'GCST90477562': 'Narcolepsy',
        'GCST90043769': 'MG(old)', 'GCST90473316': 'MS(old)', 'GCST90473258': 'SCZ(old)', 'GCST90225526': 'Depression(old)',
        'GCST90473307': 'PD(old)', 'GCST90435857': 'MoodDisorders(old)', 'GCST90104541': 'CardioembolicStroke(old)'}


def collect(root, tag):
    rows, mods = [], []
    for f in sorted(glob.glob(f'{root}/*/*.sig_module.txt')):
        b = os.path.basename(f)[:-len('.sig_module.txt')]
        ctx, g = b.rsplit('-', 1)
        n = 0; cells = set(); genes = set()
        if os.path.getsize(f) > 1:
            try:
                d = pd.read_csv(f, sep='\t')
            except Exception:
                d = pd.DataFrame()
            if len(d) and 'module_genes' in d.columns:
                n = len(d); cells = set(d.cell_type); genes = set(x for m in d.module_genes for x in str(m).split(':'))
                d.insert(0, 'GWAS', g); d.insert(0, 'ctx', ctx); mods.append(d)
        rows.append(dict(set=tag, ctx=ctx, GWAS=g, disease=NAME.get(g, g), n_sig_modules=n, n_cell_types=len(cells),
                         cell_types=';'.join(sorted(cells)), n_genes=len(genes)))
    return pd.DataFrame(rows), (pd.concat(mods, ignore_index=True) if mods else pd.DataFrame())


new, newm = collect(NEW, 'new'); his, hism = collect(HIS, 'handover')
new.to_csv(f'{OUT}/scgwas_summary_new.tsv', sep='\t', index=False)
his.to_csv(f'{OUT}/scgwas_summary_handover.tsv', sep='\t', index=False)
newm.to_csv(f'{OUT}/scgwas_sig_modules_new.tsv', sep='\t', index=False)
for tag, d in (('new', new), ('handover', his)):
    print(f'== {tag}: sig_module files {len(d)}; with >=1 significant module {(d.n_sig_modules > 0).sum()}')
    if len(d):
        s = d.groupby(['GWAS', 'disease']).agg(files=('ctx', 'size'), ctx_with_modules=('n_sig_modules', lambda x: int((x > 0).sum())),
                                                modules=('n_sig_modules', 'sum')).reset_index()
        print(s.to_string(index=False))
n_cfg = len(glob.glob('/public/home/yucong/scbrain/scgwas10/configure/*-configure.txt'))
n_cz = len(glob.glob(f'{NEW}/*/*.combined-z.txt'))
print('configs', n_cfg, 'combined-z', n_cz, 'sig_module files', len(new))
print('SCGWAS_SUMMARY_DONE')
