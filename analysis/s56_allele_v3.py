"""Counted allele for the variants allele_table_v2 left unresolved (PI 10-19 #135「不应该有这种不确定的」).
Same frequency rule as s15_allele_fix.py: of the two letters, the counted allele is the one whose reference frequency is closer to
the frequency in our genotypes (f1 = frequency of out.bim column 5), and only when the two distances differ by more than 0.2.
  step 1  letters_differ_from_1000G: the letters on the other strand (A<->T, C<->G) match 1000G -> status ok_strand
          (A/T and C/G variants are skipped: their strand cannot be told from the letters)
  step 2  still unresolved (mostly not_in_1000G): dbSNP ALFA European frequencies (SAMN10492695) -> ok_alfa / ok_alfa_strand
  check   the ALFA rule is applied blind to variants already resolved by 1000G frequency (status ok) and must agree >= 99%;
          otherwise nothing is written.
Inputs  allele_fix/allele_table_v2.tsv.gz, sn_replication/g1000.frq, external/alfa/freq.vcf.gz (filtered to our rs ids first)
Output  allele_fix/allele_table_v3.tsv.gz (same columns as v2; only previously unresolved rows change), check lines printed."""
import gzip, os, subprocess, sys
import numpy as np
import pandas as pd

B = '/public/home/yucong/scbrain'
FIX = B + '/allele_fix'
ALFA = B + '/external/alfa/freq.vcf.gz'
EUR = 'SAMN10492695'
COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
RES = ('ok', 'ok_ld')

t = pd.read_csv(FIX + '/allele_table_v2.tsv.gz', sep='\t', dtype={'a1': str, 'a2': str})
g = pd.read_csv(B + '/sn_replication/g1000.frq', sep=r'\s+', usecols=['SNP', 'A1', 'A2', 'MAF'])
g.columns = ['snp', 'g1', 'g2', 'gf']; g = g.drop_duplicates('snp', keep='last')
print('allele_table_v2', len(t)); print(t.status.value_counts().to_string())


def rule(f1, fa1):
    """-> (decisive, a1_is_counted)"""
    d1, d2 = np.abs(f1 - fa1), np.abs(f1 - (1 - fa1))
    return np.abs(d1 - d2) > 0.2, d1 <= d2


# ---- step 1: other strand vs 1000G
s1 = t.status == 'letters_differ_from_1000G'
x = t[s1].merge(g, on='snp', how='left')
ca1, ca2 = x.a1.map(COMP), x.a2.map(COMP)
pal = x.a1.map(COMP) == x.a2
match = (~pal) & (((ca1 == x.g1) & (ca2 == x.g2)) | ((ca1 == x.g2) & (ca2 == x.g1)))
fa1 = np.where(ca1 == x.g1, x.gf, 1 - x.gf)
dec, a1c = rule(x.f1.values, fa1)
ok1 = match.values & dec
idx = t.index[s1][ok1]
t.loc[idx, 'eff'] = np.where(a1c[ok1], ca1[ok1], ca2[ok1])
t.loc[idx, 'oth'] = np.where(a1c[ok1], ca2[ok1], ca1[ok1])
t.loc[idx, 'status'] = 'ok_strand'
print(f'step 1: letters_differ {s1.sum():,}; other strand matches 1000G {int(match.sum()):,} (palindromic skipped {int(pal.sum()):,}); decided {int(ok1.sum()):,}')

# ---- ALFA: keep only our rs ids (zcat | awk), then parse
ids = B + '/allele_fix/v3_ids.txt'
t.snp.to_csv(ids, index=False, header=False)
sub = B + '/allele_fix/v3_alfa_subset.tsv'
if not os.path.exists(sub):
    hdr = None
    with gzip.open(ALFA, 'rt') as fh:
        for ln in fh:
            if ln.startswith('#CHROM'):
                hdr = ln.rstrip('\n').split('\t'); break
    col = hdr.index(EUR) + 1
    cmd = (f"zcat {ALFA} | awk -F'\\t' 'NR==FNR{{a[$1];next}} !/^#/ && ($3 in a) {{print $3\"\\t\"$4\"\\t\"$5\"\\t\"${col}}}' {ids} - > {sub}.tmp && mv {sub}.tmp {sub}")
    subprocess.run(cmd, shell=True, check=True)
al = pd.read_csv(sub, sep='\t', header=None, names=['snp', 'ref', 'alt', 'eur'], dtype=str).drop_duplicates('snp')
print('ALFA rows for our ids', len(al))


def alfa_freqs(r):
    """allele -> EUR frequency, or None"""
    try:
        an, ac = r.eur.split(':')[:2]           # FORMAT AN:AC:HWEP:GR:GV:GA
        an = int(an)
        acs = [int(v) for v in ac.split(',')]
    except Exception:
        return None
    if an < 200:
        return None
    alts = r.alt.split(',')
    fr = {a: c / an for a, c in zip(alts, acs)}
    fr[r.ref] = 1 - sum(acs) / an
    return fr


al['fr'] = [alfa_freqs(r) for r in al.itertuples()]
al = al[al.fr.notna()]


def alfa_call(df):
    """frequency rule with ALFA for the rows of df (a1, a2, f1) -> (decided mask, eff, oth, status)"""
    m = df.merge(al[['snp', 'fr']], on='snp', how='left')
    dec = np.zeros(len(m), bool); eff = np.array([None] * len(m), object); oth = eff.copy(); st = eff.copy()
    for i, r in enumerate(m.itertuples()):
        fr = r.fr
        if not isinstance(fr, dict):
            continue
        for a1, a2, tag in ((r.a1, r.a2, 'ok_alfa'), (COMP.get(r.a1), COMP.get(r.a2), 'ok_alfa_strand')):
            if tag == 'ok_alfa_strand' and COMP.get(r.a1) == r.a2:
                break                           # palindromic: no strand flip test
            if a1 in fr and a2 in fr and fr[a1] + fr[a2] > 0.9:
                fa1 = fr[a1] / (fr[a1] + fr[a2])
                d, c = rule(np.array([r.f1]), np.array([fa1]))
                if d[0]:
                    dec[i] = True; eff[i], oth[i] = (a1, a2) if c[0] else (a2, a1); st[i] = tag
                break
    return dec, eff, oth, st


# ---- check: ALFA rule blind on variants resolved by 1000G frequency
okv = t[t.status == 'ok'].sample(n=min(200000, int((t.status == 'ok').sum())), random_state=1)
dec, eff, oth, st = alfa_call(okv)
agree = float((eff[dec] == okv.eff.values[dec]).mean()) if dec.any() else float('nan')
print(f'check: ALFA rule on {len(okv):,} variants already resolved by 1000G frequency: decided {int(dec.sum()):,}, agreement {agree:.5f}')
if not agree >= 0.99:
    sys.exit('ALFA rule disagrees with 1000G calls: table not written')

# ---- step 2: everything still unresolved
s2 = ~t.status.isin(RES + ('ok_strand',))
dec, eff, oth, st = alfa_call(t[s2])
idx = t.index[s2][dec]
t.loc[idx, 'eff'] = eff[dec]; t.loc[idx, 'oth'] = oth[dec]
before = t.loc[idx, 'status'].value_counts()
t.loc[idx, 'status'] = st[dec]
print(f'step 2: unresolved {int(s2.sum()):,}; ALFA decided {int(dec.sum()):,}; by previous status:'); print(before.to_string())
t.to_csv(FIX + '/allele_table_v3.tsv.gz', sep='\t', index=False)
print('allele_table_v3'); print(t.status.value_counts().to_string())
res = t.status.isin(RES + ('ok_strand', 'ok_alfa', 'ok_alfa_strand')).sum()
print(f'usable {res:,} of {len(t):,} ({res / len(t):.1%}); still unresolved {len(t) - res:,}')
