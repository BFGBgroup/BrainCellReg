"""Genes per disease in the significant scGWAS modules, and the size of the scGWAS network (GO background for Fig 6, PI #85)."""
import pandas as pd
B = '/public/home/yucong/scbrain'
NAME = {'GCST90432156': 'MG', 'GCST003566': 'MS', 'GCST90726344': 'DEP', 'GCST009325': 'PD', 'GCST90018844': 'FS',
        'GCST90128471': 'SCZ', 'GCST011179': 'CD', 'GCST013196': 'AD', 'GCST90014122': 'LS', 'GCST90477562': 'NAR'}
sg = pd.read_csv(B + '/scgwas10/scgwas_sig_modules_new.tsv', sep='\t')
sg['dz'] = sg.GWAS.map(NAME)
for dz, d in sg.groupby('dz'):
    g = set(x for m in d.module_genes for x in str(m).split(':'))
    print(dz, 'modules', len(d), 'genes', len(g))
print('all', len(set(x for m in sg.module_genes for x in str(m).split(':'))))
net = pd.read_csv('/public/home/yucong/scbrain_in/scGWAS_inputs/run/network.tsv', sep='\t', header=None, nrows=5)
print(net)
net = pd.read_csv('/public/home/yucong/scbrain_in/scGWAS_inputs/run/network.tsv', sep='\t', header=None, usecols=[0, 1])
print('network nodes', len(set(net[0]) | set(net[1])), 'edges', len(net))
