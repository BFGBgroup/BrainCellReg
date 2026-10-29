# -*- coding: utf-8 -*-
"""Pieces around the SMR rerun with corrected alleles.
  merge     coloc10_fix/coloc_all.tsv = existing coloc results of the triplets that are still q_SMR < 0.05
            + the results newly computed for triplets that had none (coloc does not use allele labels)
  compare   integrate10 (handover .esi) vs integrate10_fix: disease x gene pairs with SMR + coloc, with and without HEIDI
"""
import os
import sys

import pandas as pd

B = '/public/home/yucong/scbrain'
KEY = ['ctx', 'gwas', 'gene']


def merge():
    sig = pd.read_csv(B + '/smr10_fix/smr-sig-gene.txt', sep='\t', header=None, names=KEY)
    old = pd.read_csv(B + '/coloc10/coloc_all.tsv', sep='\t')
    keep = old.merge(sig, on=KEY)
    f = B + '/coloc10_fix/coloc_new_only.tsv'
    new = pd.read_csv(f, sep='\t') if os.path.exists(f) and os.path.getsize(f) > 0 else old.iloc[:0]
    out = pd.concat([keep, new], ignore_index=True).drop_duplicates(KEY)
    out.to_csv(B + '/coloc10_fix/coloc_all.tsv', sep='\t', index=False)
    miss = sig.merge(out[KEY], on=KEY, how='left', indicator=True)
    print('q_SMR<0.05 triplets %d | coloc reused %d, newly run %d, without a result %d | PP.H4>=0.8: %d' % (
        len(sig), len(keep), len(new), int((miss._merge == 'left_only').sum()), int((out['PP.H4.abf'] >= 0.8).sum())))
    print('COLOC_MERGE_DONE')


def compare():
    o = pd.read_csv(B + '/integrate10/gene_evidence.tsv', sep='\t')
    n = pd.read_csv(B + '/integrate10_fix/gene_evidence.tsv', sep='\t')
    k = ['disease', 'gene']
    for tag, col in (('SMR (q<0.05)', None), ('SMR + PP.H4>=0.8', 'n_ctx_h4_ge_0_8'), ('SMR + HEIDI>=0.01 + PP.H4>=0.8', 'n_ctx_smr_heidi_h4')):
        a = o if col is None else o[o[col] > 0]
        b = n if col is None else n[n[col] > 0]
        both = a.merge(b[k], on=k)
        print('%s: disease x gene pairs old %d, new %d, in both %d' % (tag, len(a), len(b), len(both)))
        if col is not None:
            lost = a.merge(b[k], on=k, how='left', indicator=True)
            gain = b.merge(a[k], on=k, how='left', indicator=True)
            print('   only old:', '; '.join('%s-%s' % (r.disease, r.symbol) for r in lost[lost._merge == 'left_only'].itertuples()) or 'none')
            print('   only new:', '; '.join('%s-%s' % (r.disease, r.symbol) for r in gain[gain._merge == 'left_only'].itertuples()) or 'none')
    b = n[n.n_ctx_h4_ge_0_8 > 0].sort_values(['disease', 'n_ctx_smr_heidi_h4'], ascending=[True, False])
    cols = ['disease', 'symbol', 'MHC', 'n_ctx_smr', 'n_ctx_heidi_ge_0_01', 'n_ctx_h4_ge_0_8', 'n_ctx_smr_heidi_h4', 'max_h4', 'scgwas_module']
    pd.set_option('display.width', 250); pd.set_option('display.max_rows', 200)
    print(b[cols].to_string(index=False))
    print('FIX_COMPARE_DONE')


if __name__ == '__main__':
    {'merge': merge, 'compare': compare}[sys.argv[1]]()
