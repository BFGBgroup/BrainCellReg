# scBRAIN database site

Results database of the scBRAIN study (110 brain region × cell-type eQTL contexts × 10 diseases).
Back end: FastAPI + DuckDB reading Parquet tables directly. Front end: Vue 3 + ECharts (vendored in `static/vendor/`, no build step);
every figure is drawn in the browser and can be saved as SVG or PNG.

## Layout

```
app/main.py          API (/api/...) and the page
static/index.html    page shell
static/app.js        all pages and plotting functions
static/style.css
static/vendor/       vue 3.4.21, echarts 5.5.0
```

Data directory (`SCBRAIN_DATA`):

```
tables/  contexts, samples, genes, egene_stats, esnp_pairs, fisher, smr, gene_evidence, scgwas_modules   (step 1, s22_db_tables.py)
         expression, expr_summary, geno_box                                                                 (step 2, s26_db_plot_tables.py)
demo/    <SYMBOL>_nominal / _gwas / _ld .parquet  locus data for the sample pages (GPNMB, ADAM15, CD38, CLU)
```

## Run

```
pip install fastapi uvicorn duckdb==0.10.3 pandas
SCBRAIN_DATA=/path/to/db uvicorn app.main:app --host 0.0.0.0 --port 8765
```

API reference: `/api/docs`.

## Open decisions (database plan v1)

- Locus plots for every gene need the 176 GB of nominal cis results on the server (`SCBRAIN_NOMINAL_DIR`, not wired yet);
  without them, locus plots are limited to genes exported to `demo/`.
- `geno_box` holds box statistics by genotype group (summary level); no per-person genotype is stored on the site.
- Fonts load from Google Fonts; for a server in mainland China, self-host them under `static/fonts/`.
