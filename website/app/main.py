"""scBRAIN results database: API (FastAPI + DuckDB over the Parquet tables) and the single-page front end.
Run:  SCBRAIN_DATA=<dir with tables/ and demo/> uvicorn app.main:app --port 8765
Every figure is drawn in the browser from these JSON answers (PI 10-12: the plots are drawn on demand).
"""
import os, re, glob, math
from functools import lru_cache
import duckdb
import pandas as pd
from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

DATA = os.environ.get('SCBRAIN_DATA', 'D:/超算/scBRAIN/db')
TB = DATA + '/tables'
DEMO = DATA + '/demo'
# full cis results of every tested variant–gene pair, repacked by server_scripts/s52_repack_nominal.py:
# <NOMINAL_DIR>/chr<N>/<ctx>.parquet with pos, phenotype_id, variant_id, af, slope, slope_se, mlog10p, sorted by pos
NOMINAL_DIR = os.environ.get('SCBRAIN_NOMINAL_DIR')
CIS = 1_000_000  # tensorQTL cis window around the TSS
HERE = os.path.dirname(os.path.abspath(__file__))
STATIC = os.path.join(os.path.dirname(HERE), 'static')

REG = ['FCTX', 'TCTX', 'OCTX', 'HIPP', 'PUTM', 'THAL', 'SNIG', 'MEDU', 'WHMT', 'CRBL']
CELL = ['ExN', 'InN', 'Astro', 'Olig', 'OPC', 'Micro', 'Endo', 'Perc', 'NPC', 'Gran', 'Purk']
REG_NAME = {'FCTX': 'Frontal cortex', 'TCTX': 'Temporal cortex', 'OCTX': 'Occipital cortex', 'HIPP': 'Hippocampus', 'PUTM': 'Putamen',
            'THAL': 'Thalamus', 'SNIG': 'Substantia nigra', 'MEDU': 'Medulla', 'WHMT': 'Intralobular white matter', 'CRBL': 'Cerebellum'}
CELL_NAME = {'ExN': 'Excitatory neurons', 'InN': 'Inhibitory neurons', 'Astro': 'Astrocytes', 'Olig': 'Oligodendrocytes', 'OPC': 'OPCs',
             'Micro': 'Microglia', 'Endo': 'Endothelial cells', 'Perc': 'Pericytes', 'NPC': 'Neural progenitor cells', 'Gran': 'Granule cells',
             'Purk': 'Purkinje cells'}
# the 10 diseases of the integration (SMR / colocalisation / scGWAS) and their GWAS Catalog trait in the Fisher enrichment
DISEASES = [
    ('PD', "Parkinson's disease", 'GCST009325', 'Parkinsons disease'),
    ('AD', "Alzheimer's disease", 'GCST013196', 'Alzheimers disease'),
    ('MS', 'Multiple sclerosis', 'GCST003566', 'Multiple sclerosis'),
    ('MG', 'Myasthenia gravis', 'GCST90432156', 'Myasthenia gravis'),
    ('SCZ', 'Schizophrenia', 'GCST90128471', 'Schizophrenia'),
    ('MDD', 'Depression', 'GCST90726344', 'Depression'),
    ('Narcolepsy', 'Narcolepsy', 'GCST90477562', 'Narcolepsy'),
    ('LacunarStroke', 'Lacunar stroke', 'GCST90014122', 'Lacunar stroke'),
    ('CUD', 'Cannabis dependence', 'GCST011179', 'Cannabis dependence'),
    ('FS', 'Febrile seizures', 'GCST90018844', 'Febrile seizures'),
]
DZ = {d[0]: d for d in DISEASES}
GWAS2DZ = {d[2]: d[0] for d in DISEASES}

# same tables stored sorted by their lookup key (server_scripts/s62_sort_expression.py; identical rows), plus copies sorted by rsid;
# a file here replaces the one of the same name in tables/ (unset the variable to go back)
SORTED_DIR = os.environ.get('SCBRAIN_SORTED_DIR')
FILES = {os.path.basename(t)[:-8]: t for t in glob.glob(TB + '/*.parquet')}
if SORTED_DIR and os.path.isdir(SORTED_DIR):
    FILES.update({os.path.basename(t)[:-8]: t for t in glob.glob(SORTED_DIR + '/*.parquet')})
con = duckdb.connect(read_only=False)
con.execute("SET threads TO 4")
for name, t in FILES.items():
    con.execute(f"CREATE VIEW {name} AS SELECT * FROM read_parquet('{t}')")
HAVE = set(FILES)


def q(sql, *args):
    cur = con.cursor()
    return cur.execute(sql, list(args)).df()


def qd(sql, **frames):
    """SQL over DataFrames given as keyword arguments (each registered under its name), next to the table views."""
    cur = con.cursor()
    for k, v in frames.items():
        cur.register(k, v)
    return cur.execute(sql).df()


class RowGroups:
    """A parquet file sorted by one key column: key -> the row groups that can hold it, from the file's own min/max statistics
    (metadata only, read once). DuckDB 0.10 does not skip row groups on string statistics, so a lookup in the 807 MB expression
    table read the whole file (3.3 s); reading only the matching row group takes milliseconds."""

    def __init__(self, path, col):
        import pyarrow.parquet as pq
        self.pq, self.path, self.col = pq, path, col
        self.md = pq.ParquetFile(path).metadata
        # pyarrow ignores the string statistics DuckDB writes (no column order), DuckDB's own metadata reader returns them
        st = duckdb.execute(f"""SELECT row_group_id, coalesce(stats_min_value, stats_min) lo, coalesce(stats_max_value, stats_max) hi
                                FROM parquet_metadata('{path}') WHERE path_in_schema = '{col}' ORDER BY row_group_id""").df()
        if len(st) != self.md.num_row_groups or st.lo.isna().any() or st.hi.isna().any():
            raise ValueError('no statistics')
        self.lo, self.hi = st.lo.tolist(), st.hi.tolist()
        if any(self.lo[i] > self.lo[i + 1] or self.hi[i] > self.lo[i + 1] for i in range(len(self.lo) - 1)):
            raise ValueError('not sorted')   # overlapping ranges: not sorted by this column

    def rows(self, val):
        import bisect
        i = bisect.bisect_left(self.hi, val)
        rgs = []
        while i < len(self.lo) and self.lo[i] <= val:
            rgs.append(i); i += 1
        f = self.pq.ParquetFile(self.path, metadata=self.md)
        t = f.read_row_groups(rgs) if rgs else f.schema_arrow.empty_table()
        d = t.to_pandas()
        return d[d[self.col] == val].reset_index(drop=True)


RG = {}
for name, col in [('expression', 'gene'), ('egene_stats', 'gene'), ('esnp_pairs', 'gene'), ('geno_box', 'gene'), ('expr_summary', 'gene'),
                  ('esnp_pairs_by_rsid', 'rsid'), ('geno_box_by_rsid', 'rsid'), ('egene_stats_by_lead', 'lead_rsid')]:
    if name in FILES:
        try:
            RG[name] = RowGroups(FILES[name], col)
        except (ValueError, ImportError, KeyError):
            pass                      # unsorted (e.g. a local copy): plain SQL below


def keyed(table, col, val):
    """Every row of a table whose column col equals val (val already checked); by row groups when the file is sorted by col."""
    ix = RG.get(table)
    if ix is not None and ix.col == col:
        return ix.rows(val)
    return q(f"SELECT * FROM {table} WHERE {col} = '{val}'")


# Lookups by gene / variant / context are written into the SQL as checked literals: DuckDB 0.10 skips parquet row groups by their
# min/max statistics only for constants, and the big tables are stored sorted by these keys (server_scripts/s62_sort_expression.py)
CTXS = {f'{r}.{c}' for r in REG for c in CELL}
_PAT = {'gene': r'ENSG\d+(\.\d+)?', 'rsid': r'rs\d+'}


def lit(x, kind):
    """x as a SQL string literal after checking it is a gene id, an rs number or one of the 110 contexts (404 otherwise)."""
    ok = (x in CTXS) if kind == 'ctx' else bool(re.fullmatch(_PAT[kind], str(x)))
    if not ok:
        raise HTTPException(404, f'{x} is not a valid {kind}.')
    return f"'{x}'"


# copies sorted by rsid for the variant pages, when present (otherwise the gene-sorted tables are scanned)
ESNP_R = 'esnp_pairs_by_rsid' if 'esnp_pairs_by_rsid' in HAVE else 'esnp_pairs'
GENO_R = 'geno_box_by_rsid' if 'geno_box_by_rsid' in HAVE else 'geno_box'
EGENE_L = 'egene_stats_by_lead' if 'egene_stats_by_lead' in HAVE else 'egene_stats'


def rec(df):
    """DataFrame -> JSON-safe records (NaN -> null)."""
    out = df.astype(object).where(pd.notna(df), None).to_dict('records')
    for r in out:
        for k, v in r.items():
            if isinstance(v, float) and (math.isinf(v) or math.isnan(v)):
                r[k] = None
    return out


app = FastAPI(title='BrainCellReg', docs_url='/api/docs', openapi_url='/api/openapi.json')

# --- lookup index -------------------------------------------------------------------------------------------------
GENES = q("SELECT gene, symbol, biotype, chr, start, \"end\", strand, is_egene_anywhere FROM genes")
SYM2ID = {}
for r in GENES.itertuples():
    SYM2ID.setdefault(str(r.symbol).upper(), r.gene)
    SYM2ID[r.gene.upper()] = r.gene
RSIDS = set(q("SELECT DISTINCT rsid FROM esnp_pairs").rsid)
FULL = bool(NOMINAL_DIR) and 'variants' in HAVE and os.path.isdir(NOMINAL_DIR or '')
SYM = dict(zip(GENES.gene, GENES.symbol))


def variant_info(rs):
    """chr / pos / alleles of any tested variant (variants table, all 6.4 million), or None."""
    if 'variants' not in HAVE:
        return None
    # chrX variants are genotyped but were not in the cis-eQTL scan (autosomes only), so they count as not tested
    v = q(f"SELECT * FROM variants WHERE rsid = {lit(rs, 'rsid')} AND chr <> 'chrX'")     # in_1000g only in the allele v3 table
    return None if v.empty else rec(v)[0]


def nominal(chrom, lo, hi, rsid=None, gene=None, ctx=None):
    """Rows of the full cis results on one chromosome between two positions, all 110 contexts (context from the file name),
    or only the file of one context when ctx is given (the locus plot: one file instead of 110, 5-7 s -> 0.05 s).
    The filter is written into the SQL as literals (ints / checked ids): DuckDB 0.10 skips row groups by their pos statistics
    only for constants, not for ? parameters (one variant took 16-23 s with parameters)."""
    c = 'chr' + str(chrom).replace('chr', '')
    empty = pd.DataFrame(columns=['ctx', 'pos', 'gene', 'rsid', 'af', 'slope', 'se', 'p'])
    if ctx is not None and ctx not in CTXS:
        return empty
    f = f'{NOMINAL_DIR}/{c}/{ctx or "*"}.parquet'
    if not glob.glob(f):
        return empty
    where = f'pos BETWEEN {int(lo)} AND {int(hi)}'
    if rsid is not None:
        if not re.fullmatch(r'rs\d+', rsid):
            return empty
        where += f" AND variant_id = '{rsid}'"
    if gene is not None:
        if not re.fullmatch(r'ENSG\d+(\.\d+)?', gene):
            return empty
        where += f" AND phenotype_id = '{gene}'"
    return q(f"""SELECT regexp_extract(filename, '([^/\\\\]+)\\.parquet$', 1) ctx, pos, phenotype_id::VARCHAR gene, variant_id rsid,
                     af, slope, slope_se se, power(10, -mlog10p::DOUBLE) p
              FROM read_parquet('{f}', filename = true) WHERE {where}""")


@app.get('/api/meta')
@lru_cache(maxsize=1)            # the tables do not change while the site runs; every page load asks for this
def meta():
    cx = q("SELECT * FROM contexts")
    tot = {
        'contexts': len(cx), 'contexts_with_egenes': int((cx.egenes > 0).sum()),
        'samples': int(q("SELECT count(*) n FROM samples").n[0]), 'individuals': int(q("SELECT count(DISTINCT individual) n FROM samples").n[0]),
        'genes_tested': len(GENES), 'egenes': int(GENES.is_egene_anywhere.sum()),
        'independent_esnps': int(q("SELECT count(DISTINCT rsid) n FROM esnp_pairs WHERE independent").n[0]),
        'esnp_pairs': int(q("SELECT count(*) n FROM esnp_pairs").n[0]),
    }
    return {'regions': [{'id': r, 'name': REG_NAME[r]} for r in REG], 'cells': [{'id': c, 'name': CELL_NAME[c]} for c in CELL],
            'contexts': rec(cx), 'totals': tot,
            'diseases': [{'id': d[0], 'name': d[1], 'gwas': d[2]} for d in DISEASES],
            'geno_box': 'geno_box' in HAVE, 'expression': 'expression' in HAVE, 'full': FULL}


@app.get('/api/search')
def search(text: str = Query('', alias='q'), limit: int = 12):
    s = text.strip()
    if not s:
        return []
    out = []
    S = s.upper()
    if re.fullmatch(r'RS\d+', S):
        rs = s.lower()
        if rs in RSIDS:
            out.append({'type': 'variant', 'id': rs, 'label': rs, 'sub': 'eSNP'})
        elif FULL and variant_info(rs):
            out.append({'type': 'variant', 'id': rs, 'label': rs, 'sub': 'tested, not an eSNP'})
    for d in DISEASES:
        if S in d[0].upper() or S in d[1].upper():
            out.append({'type': 'disease', 'id': d[0], 'label': d[1]})
    m = GENES[GENES.symbol.str.upper().str.startswith(S) | (GENES.gene == S)]
    m = m.assign(exact=m.symbol.str.upper() == S).sort_values(['exact', 'is_egene_anywhere', 'symbol'], ascending=[False, False, True])
    out += [{'type': 'gene', 'id': r.gene, 'label': r.symbol, 'sub': r.gene + (' · eGene' if r.is_egene_anywhere else '')}
            for r in m.head(limit).itertuples()]
    return out[:limit]


def gene_id(g):
    gid = SYM2ID.get(g.upper())
    if not gid:
        raise HTTPException(404, f'{g} is not among the 15,661 genes tested.')
    return gid


@app.get('/api/gene/{g}')
def gene(g: str):
    gid = gene_id(g)
    info = rec(GENES[GENES.gene == gid])[0]
    G = lit(gid, 'gene')
    eq = qd("""SELECT ctx, region, cell, lead_rsid, effect_allele, other_allele, allele_status, af, pval_nominal, slope, slope_se, qval, is_egene, num_var
               FROM t ORDER BY ctx""", t=keyed('egene_stats', 'gene', gid))
    ex = qd("SELECT ctx, n, mean, sd, q1, median, q3 FROM t ORDER BY ctx", t=keyed('expr_summary', 'gene', gid)) if 'expr_summary' in HAVE else pd.DataFrame()
    ev = q(f"SELECT * FROM gene_evidence WHERE gene = {G}")
    ev['disease_name'] = ev.disease.map(lambda d: DZ.get(d, (d, d))[1])
    smr = q(f"""SELECT disease, GWAS, ctx, topSNP, p_SMR, q_SMR, b_SMR, se_SMR, p_HEIDI, PP_H4_abf, smr_heidi_coloc
                FROM smr WHERE probeID = {G} AND q_SMR < 0.05 ORDER BY disease, p_SMR, ctx""")
    smr['disease'] = smr.disease.fillna(smr.GWAS.map(GWAS2DZ))
    n_pairs = len(keyed('esnp_pairs', 'gene', gid))
    leads = eq[eq.is_egene].sort_values('qval')[['ctx', 'lead_rsid']]
    return {'info': info, 'eqtl': rec(eq), 'expr': rec(ex), 'evidence': rec(ev), 'smr': rec(smr), 'n_esnp_pairs': n_pairs,
            'leads': rec(leads), 'locus': locus_available(info['symbol'])}


@app.get('/api/gene/{g}/groups')
def gene_groups(g: str, ctx: str):
    """Expression by sex and by age group in one context (box statistics, computed on the fly)."""
    gid = gene_id(g)
    if 'expression' not in HAVE:
        raise HTTPException(404, 'Per-sample expression is not in this build.')
    C = lit(ctx, 'ctx')
    region = ctx.split('.')[0]
    lit(gid, 'gene')
    d = qd(f"""SELECT e.expr, s.sex, s.age_group FROM (SELECT sample, expr FROM t WHERE ctx = {C}) e
               JOIN samples s ON s.individual = e.sample AND s.region = '{region}'""", t=keyed('expression', 'gene', gid))
    d['expr'] = d.expr.astype(float)

    def box(by):
        out = []
        for k, v in d.groupby(by):
            x = v.expr
            q1, med, q3 = x.quantile([.25, .5, .75])
            iqr = q3 - q1
            out.append({'group': k, 'n': len(x), 'q1': float(q1), 'median': float(med), 'q3': float(q3),
                        'lo': float(x[x >= q1 - 1.5 * iqr].min()), 'hi': float(x[x <= q3 + 1.5 * iqr].max()),
                        'out': x[(x < q1 - 1.5 * iqr) | (x > q3 + 1.5 * iqr)].round(4).tolist()})
        return out
    return {'ctx': ctx, 'sex': box('sex'), 'age': box('age_group')}


@app.get('/api/gene/{g}/esnps')
def gene_esnps(g: str, ctx: str):
    """Every eSNP of the gene in one context, strongest first (for the genotype panel's variant list)."""
    gid = gene_id(g)
    lit(gid, 'gene')
    d = qd(f"""SELECT rsid, pval, slope, se, independent, effect_allele, other_allele FROM t
               WHERE ctx = {lit(ctx, 'ctx')} ORDER BY pval, rsid""", t=keyed('esnp_pairs', 'gene', gid))
    return {'gene': gid, 'ctx': ctx, 'rows': rec(d)}


@app.get('/api/genotype')
def genotype(gene: str, rsid: str):
    """Expression by genotype of one variant for one gene, every context where the pair is a lead or independent eSNP."""
    if 'geno_box' not in HAVE:
        raise HTTPException(404, 'Genotype groups are not in this build.')
    gid = gene_id(gene)
    G, R = lit(gid, 'gene'), lit(rsid.lower(), 'rsid')
    d = qd(f"SELECT * FROM t WHERE rsid = {R} ORDER BY ctx, n_effect", t=keyed('geno_box', 'gene', gid))
    st = qd(f"SELECT ctx, pval, slope, se FROM t WHERE rsid = {R} ORDER BY ctx", t=keyed('esnp_pairs', 'gene', gid))
    return {'boxes': rec(d), 'stats': rec(st)}


@app.get('/api/variant/{rs}')
def variant(rs: str):
    rs = rs.lower()
    R = lit(rs, 'rsid')
    pr = qd("""SELECT p.ctx, p.region, p.cell, p.gene, g.symbol, p.biotype, p.chr, p.pos, p.dist, p.MAF, p.effect_allele, p.other_allele, p.allele_status,
                      p.pval, p.slope, p.se, p.independent
               FROM t p LEFT JOIN genes g USING (gene) ORDER BY p.pval, p.ctx, p.gene""", t=keyed(ESNP_R, 'rsid', rs))
    info = variant_info(rs)
    if pr.empty and not (FULL and info):
        raise HTTPException(404, f'{rs} is not an eSNP of any eGene in the 110 contexts.' if not FULL else f'{rs} was not tested (not among the genotyped variants).')
    lead = qd("SELECT ctx, gene FROM t WHERE is_egene ORDER BY ctx, gene", t=keyed(EGENE_L, 'lead_rsid', rs))
    gb = qd("SELECT DISTINCT gene FROM t ORDER BY gene", t=keyed(GENO_R, 'rsid', rs)).gene.tolist() if 'geno_box' in HAVE else []
    smr = q(f"""SELECT disease, GWAS, ctx, Gene gene, p_SMR, q_SMR, p_HEIDI, PP_H4_abf, smr_heidi_coloc FROM smr WHERE topSNP = {R} AND q_SMR < 0.05
                ORDER BY disease, ctx, gene""")
    smr['disease'] = smr.disease.fillna(smr.GWAS.map(GWAS2DZ))
    return {'rsid': rs, 'info': info, 'pairs': rec(pr), 'lead_of': rec(lead), 'geno_genes': gb, 'smr': rec(smr), 'full': FULL}


@app.get('/api/variant/{rs}/tested')
def variant_tested(rs: str):
    """Every gene tested with this variant in every context, significant or not (full cis results)."""
    rs = rs.lower()
    if not FULL:
        raise HTTPException(404, 'The full results of every tested pair are not in this build.')
    info = variant_info(rs)
    if not info:
        raise HTTPException(404, f'{rs} was not tested (not among the genotyped variants).')
    d = nominal(info['chr'], info['pos'], info['pos'], rsid=rs)
    d['symbol'] = d.gene.map(SYM)
    lit(rs, 'rsid')
    sig = keyed(ESNP_R, 'rsid', rs)[['ctx', 'gene', 'independent']]
    d = d.merge(sig.assign(esnp=True), on=['ctx', 'gene'], how='left')
    d['esnp'] = d.esnp.fillna(False).astype(bool)
    d['independent'] = d.independent.fillna(False).astype(bool)
    return {'rsid': rs, 'info': info, 'rows': rec(d.drop(columns=['pos', 'rsid']).sort_values('p'))}


@app.get('/api/pair')
def pair(gene: str, rsid: str):
    """One variant–gene pair in all 110 contexts (full cis results)."""
    gid = gene_id(gene)
    rs = rsid.lower()
    if not FULL:
        raise HTTPException(404, 'The full results of every tested pair are not in this build.')
    info = variant_info(rs)
    if not info:
        raise HTTPException(404, f'{rs} was not tested (not among the genotyped variants).')
    d = nominal(info['chr'], info['pos'], info['pos'], rsid=rs, gene=gid)
    if d.empty:
        raise HTTPException(404, f'{rs} and {SYM.get(gid, gid)} were not tested together (the variant is outside the gene\'s ±1 Mb cis window).')
    sig = qd(f"SELECT ctx, independent FROM t WHERE rsid = {lit(rs, 'rsid')}", t=keyed('esnp_pairs', 'gene', lit(gid, 'gene')[1:-1]))
    d = d.merge(sig.assign(esnp=True), on='ctx', how='left')
    d['esnp'] = d.esnp.fillna(False).astype(bool)
    d['independent'] = d.independent.fillna(False).astype(bool)
    return {'gene': gid, 'symbol': SYM.get(gid), 'rsid': rs, 'info': info, 'rows': rec(d.drop(columns=['pos', 'rsid', 'gene']))}


@app.get('/api/disease/{d}')
@lru_cache(maxsize=32)           # 10 diseases, 11 cell types, 110 contexts: the answers never change while the site runs
def disease(d: str):
    if d not in DZ:
        raise HTTPException(404, f'{d} is not one of the 10 diseases.')
    code, name, gwas, trait = DZ[d]
    fi = q("""SELECT ctx, E_N, E_NN, NE_N, NE_NN, "OR", OR_lower95, OR_upper95, p_value, q_value, significant FROM fisher
              WHERE disease = ? AND class = 'neurological'""", trait)
    ev = q("SELECT * FROM gene_evidence WHERE GWAS = ? ORDER BY n_ctx_smr_heidi_h4 DESC, min_q_smr", gwas)
    sm = q("""SELECT ctx, probeID gene, Gene, p_SMR, q_SMR, b_SMR, p_HEIDI, PP_H4_abf, smr_heidi_coloc FROM smr WHERE GWAS = ? AND q_SMR < 0.05
              ORDER BY p_SMR, ctx, gene""", gwas)
    sg = q("SELECT ctx, seed, module_genes, module_score, p, p_gwas, p_scrn FROM scgwas_modules WHERE GWAS = ? ORDER BY p, ctx", gwas)
    sym = dict(zip(GENES.gene, GENES.symbol))
    sg['seed_symbol'] = sg.seed.map(sym)
    sg['module_symbols'] = sg.module_genes.map(lambda s: [sym.get(x, x) for x in str(s).split(':')])
    sm['symbol'] = sm.gene.map(sym)
    return {'id': code, 'name': name, 'gwas': gwas, 'trait': trait, 'fisher': rec(fi), 'genes': rec(ev), 'smr': rec(sm), 'scgwas': rec(sg)}


@app.get('/api/context/{ctx}')
@lru_cache(maxsize=128)
def context(ctx: str):
    cx = q("SELECT * FROM contexts WHERE ctx = ?", ctx)
    if cx.empty:
        raise HTTPException(404, f'{ctx} is not one of the 110 contexts.')
    region = ctx.split('.')[0]
    fr = q(f"SELECT {', '.join(CELL)} FROM samples WHERE region = ?", region)
    box = []
    for c in CELL:
        x = fr[c]
        q1, med, q3 = x.quantile([.25, .5, .75])
        box.append({'cell': c, 'q1': float(q1), 'median': float(med), 'q3': float(q3), 'lo': float(x.min()), 'hi': float(x.max())})
    eg = q("""SELECT gene, symbol, biotype, lead_rsid, slope, pval_nominal, qval FROM egene_stats WHERE ctx = ? AND is_egene ORDER BY qval, pval_nominal""", ctx)
    # sharing of eGenes with every other context (Jaccard)
    sh = q("""WITH e AS (SELECT ctx, gene FROM egene_stats WHERE is_egene),
                   a AS (SELECT gene FROM e WHERE ctx = ?)
              SELECT e.ctx, count(*) FILTER (WHERE e.gene IN (SELECT gene FROM a)) inter, count(*) n FROM e GROUP BY e.ctx""", ctx)
    na = len(eg)
    sh['jaccard'] = sh.inter / (sh.n + na - sh.inter).where(lambda s: s > 0)
    return {'context': rec(cx)[0], 'region_name': REG_NAME[region], 'cell_name': CELL_NAME[ctx.split('.')[1]],
            'fractions': box, 'egenes': rec(eg), 'sharing': rec(sh[sh.ctx != ctx].sort_values(['jaccard', 'ctx'], ascending=[False, True]))}


@app.get('/api/celltype/{cell}')
@lru_cache(maxsize=32)
def celltype(cell: str):
    """One cell type's estimated fraction in each of the 10 regions (box statistics over samples), for the home-page grid."""
    if cell not in CELL:
        raise HTTPException(404, f'{cell} is not one of the 11 cell types.')
    fr = q(f"SELECT region, {cell} AS f FROM samples")
    cx = q("SELECT ctx, egenes, independent_esnps FROM contexts WHERE cell = ?", cell).set_index('ctx')
    out = []
    for r in REG:
        x = fr.f[fr.region == r]
        q1, med, q3 = x.quantile([.25, .5, .75])
        k = f'{r}.{cell}'
        ie = cx.independent_esnps.get(k, 0)
        out.append({'region': r, 'region_name': REG_NAME[r], 'ctx': k, 'n': int(len(x)), 'mean': float(x.mean()),
                    'q1': float(q1), 'median': float(med), 'q3': float(q3), 'lo': float(x.min()), 'hi': float(x.max()),
                    'egenes': int(cx.egenes.get(k, 0)), 'independent_esnps': 0 if pd.isna(ie) else int(ie)})
    return {'cell': cell, 'cell_name': CELL_NAME[cell], 'regions': out}


# --- locus plot ---------------------------------------------------------------------------------------------------
# GWAS P of the 10 diseases for every variant (server_scripts/s55_web_gwas_ld.py) and 1000G EUR per chromosome for LD on demand
GWAS_DIR = os.environ.get('SCBRAIN_GWAS_DIR')
LD_DIR = os.environ.get('SCBRAIN_LD_DIR')
PLINK = os.environ.get('SCBRAIN_PLINK', 'plink')


@lru_cache(maxsize=512)
def _ld_cached(chrom, index, lo, hi):
    import subprocess, tempfile
    bf = f'{LD_DIR}/chr{chrom}'
    if not (LD_DIR and os.path.exists(bf + '.bed')):
        return ()
    with tempfile.TemporaryDirectory() as td:
        kb = int((hi - lo) / 1000) + 50
        r = subprocess.run([PLINK, '--bfile', bf, '--ld-snp', index, '--ld-window-kb', str(kb), '--ld-window', '999999',
                            '--ld-window-r2', '0', '--r2', '--threads', '1', '--memory', '1000', '--out', f'{td}/x'],
                           capture_output=True, text=True, timeout=60)
        f = f'{td}/x.ld'
        if r.returncode != 0 or not os.path.exists(f):
            return ()
        d = pd.read_csv(f, sep=r'\s+', usecols=['SNP_B', 'R2'])
    return tuple(zip(d.SNP_B, d.R2.round(4)))


def ld_r2(chrom, index, lo, hi):
    """r2 (1000G EUR) between the index variant and every variant of the window; empty when the index is not in 1000G"""
    rows = _ld_cached(str(chrom).replace('chr', ''), index, lo, hi)
    return pd.DataFrame(list(rows), columns=['rsid', 'r2'])


def locus_available(symbol):
    return FULL or os.path.exists(f'{DEMO}/{symbol}_nominal.parquet')


@app.get('/api/locus/{g}')
def locus(g: str, ctx: str, disease: str = 'PD', index: str = ''):
    gid = gene_id(g)
    info = GENES[GENES.gene == gid].iloc[0]
    sym = info.symbol
    if not FULL and os.path.exists(f'{DEMO}/{sym}_nominal.parquet'):   # the demo files are only the fallback without full results
        e = q(f"SELECT rsid, pos38 pos, p, slope FROM read_parquet('{DEMO}/{sym}_nominal.parquet') WHERE ctx = ?", ctx)
        gw = q(f"SELECT rsid, pos38 pos, p FROM read_parquet('{DEMO}/{sym}_gwas.parquet') WHERE disease = ?", disease)
        ldf = f'{DEMO}/{sym}_ld.parquet'
        have_gwas = True
    elif FULL:
        tss = int(info.start if info.strand == '+' else info['end'])
        e = nominal(info.chr, tss - CIS - 1000, tss + CIS + 1000, gene=gid, ctx=ctx)
        e = e[e.ctx == ctx][['rsid', 'pos', 'p', 'slope']]
        gwf = f'{GWAS_DIR}/{disease}.parquet'
        have_gwas = bool(GWAS_DIR) and os.path.exists(gwf) and len(e) > 0
        if have_gwas:
            chrom = 'chr' + str(info.chr).replace('chr', '')
            gw = q(f"""SELECT rsid, pos, power(10, -mlog10p::DOUBLE) p FROM read_parquet('{gwf}')
                       WHERE chr = '{chrom}' AND pos BETWEEN {int(e.pos.min())} AND {int(e.pos.max())}""")
        else:
            gw = pd.DataFrame(columns=['rsid', 'pos', 'p'])
        ldf = None
    else:
        raise HTTPException(404, f'Locus data for {sym} is not in this build.')
    e = e.dropna(subset=['pos'])
    gw = gw.dropna(subset=['pos'])
    if not index:
        index = e.loc[e.p.idxmin(), 'rsid'] if len(e) and e.p.notna().any() else ''
    if ldf:
        ld = q(f"SELECT rsid, r2 FROM read_parquet('{ldf}') WHERE index_rsid = ?", index)
    elif index and len(e):
        ld = ld_r2(str(info.chr), index, int(e.pos.min()), int(e.pos.max()))
    else:
        ld = pd.DataFrame(columns=['rsid', 'r2'])
    r2 = dict(zip(ld.rsid, ld.r2)); r2[index] = 1.0
    e['r2'] = e.rsid.map(r2); gw['r2'] = gw.rsid.map(r2)
    near = GENES[(GENES.chr == info.chr) & (GENES.end > e.pos.min()) & (GENES.start < e.pos.max())][['symbol', 'start', 'end', 'strand', 'gene']]
    return {'gene': sym, 'ctx': ctx, 'disease': disease, 'index': index, 'chr': str(info.chr),
            'eqtl': rec(e), 'gwas': rec(gw), 'genes': rec(near), 'ld_indexed': len(ld) > 0, 'have_gwas': have_gwas}


@app.get('/api/download/{table}')
def download(table: str):
    p = f'{TB}/{table}.parquet'
    if table not in HAVE or table in ('expression',):
        raise HTTPException(404, 'No such table.')
    return FileResponse(p, filename=f'BrainCellReg_{table}.parquet')


@app.get('/')
def index():
    return FileResponse(os.path.join(STATIC, 'index.html'))


app.mount('/static', StaticFiles(directory=STATIC), name='static')
# Let's Encrypt http-01 challenge files (acme.sh --webroot, see deploy/ssl.sh); nginx passes every path of the site here
ACME = os.environ.get('SCBRAIN_ACME_DIR')
if ACME and os.path.isdir(ACME):
    app.mount('/.well-known/acme-challenge', StaticFiles(directory=ACME), name='acme')
