# Functional enrichment of the disease-associated genes (PI 10-14 #34: "还差计算出疾病相关基因后的功能解释，比如富集分析").
# Rules fixed before looking at results:
#   gene set   = genes with q_SMR < 0.05 in >= 1 context for the disease (integrate10_fix/gene_evidence.tsv; corrected alleles;
#                p_SMR is identical to the original run, the gene list differs in 5 of 253 disease x gene pairs)
#   universe   = all genes tested by SMR for that disease (smr10_fix/smr_all_results.tsv.gz)
#   databases  = GO biological process (org.Hs.eg.db) and Reactome (ReactomePA + reactome.db); both offline. KEGG not used:
#                it needs online access (KEGGREST) and the msigdbr data package is not installed on the cluster
#   test       = one-sided hypergeometric, BH within disease x database, significant at q < 0.05; gene-set size 10-500
#   diseases   = only those with >= 10 genes; also all diseases pooled; main analysis keeps MHC genes, sensitivity drops them
# Rscript s20_enrich.R  -> $OUT/enrich_all.tsv, enrich_summary.tsv
suppressPackageStartupMessages({ library(data.table); library(clusterProfiler); library(org.Hs.eg.db); library(ReactomePA) })
B <- "/public/home/yucong/scbrain"
OUT <- file.path(B, "enrich10_fix"); dir.create(OUT, showWarnings = FALSE)
ge <- fread(file.path(B, "integrate10_fix/gene_evidence.tsv"))
sm <- fread(file.path(B, "smr10_fix/smr_all_results.tsv.gz"), select = c("GWAS", "probeID", "ProbeChr", "Probe_bp"))
dz <- unique(ge[, .(disease, GWAS)])
uni <- unique(sm[, .(GWAS, gene = probeID, mhc = ProbeChr == 6 & Probe_bp > 25e6 & Probe_bp < 35e6)])
cat("genes per disease:\n"); print(ge[, .N, by = disease][order(-N)])
cat("universe per GWAS:\n"); print(uni[, .N, by = GWAS])

map <- suppressMessages(bitr(unique(uni$gene), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
map <- as.data.table(map)[!duplicated(ENSEMBL)]
cat("Ensembl -> Entrez mapped:", nrow(map), "of", uniqueN(uni$gene), "\n")

run <- function(genes, universe, label, mhc_rule) {
  g <- unique(map[ENSEMBL %in% genes, ENTREZID]); u <- unique(map[ENSEMBL %in% universe, ENTREZID])
  if (length(g) < 10) return(data.table(set = label, mhc = mhc_rule, db = NA, n_genes = length(g), note = "fewer than 10 mapped genes; not tested"))
  res <- list()
  go <- enrichGO(g, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", universe = u, pvalueCutoff = 1, qvalueCutoff = 1,
                 pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, readable = TRUE)
  if (!is.null(go)) res$GO_BP <- as.data.table(go@result)
  rp <- ReactomePA::enrichPathway(g, organism = "human", universe = u, pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH",
                                  minGSSize = 10, maxGSSize = 500, readable = TRUE)
  if (!is.null(rp)) res$Reactome <- as.data.table(rp@result)
  out <- rbindlist(lapply(names(res), function(d) res[[d]][, db := d]), fill = TRUE)
  out[, `:=`(set = label, mhc = mhc_rule, n_genes = length(g), n_universe = length(u))]
  out
}
sym <- as.data.table(suppressMessages(bitr(map$ENTREZID, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)))
map2sym <- function(x) sym$SYMBOL[match(x, sym$ENTREZID)]

all <- list()
for (i in seq_len(nrow(dz))) {
  d <- dz$disease[i]; gw <- dz$GWAS[i]
  u <- uni[GWAS == gw]
  gg <- ge[disease == d, gene]
  all[[paste(d, "all")]] <- run(gg, u$gene, d, "with MHC")
  all[[paste(d, "nomhc")]] <- run(setdiff(gg, ge[disease == d & MHC == TRUE, gene]), u[mhc == FALSE, gene], d, "without MHC")
}
upool <- unique(uni$gene)
all[["pooled all"]] <- run(unique(ge$gene), upool, "All diseases", "with MHC")
all[["pooled nomhc"]] <- run(unique(ge[MHC == FALSE, gene]), unique(uni[mhc == FALSE, gene]), "All diseases", "without MHC")
E <- rbindlist(all, fill = TRUE)
fwrite(E, file.path(OUT, "enrich_all.tsv"), sep = "\t")
S <- E[!is.na(db), .(n_genes = n_genes[1], n_terms_tested = sum(!is.na(pvalue)), n_q05 = sum(p.adjust < 0.05, na.rm = TRUE),
                    top_term = Description[which.min(pvalue)], top_q = min(p.adjust, na.rm = TRUE)), by = .(set, mhc, db)]
fwrite(S, file.path(OUT, "enrich_summary.tsv"), sep = "\t")
print(S[order(set, mhc, db)])
cat("\nnot tested:\n"); print(unique(E[is.na(db), .(set, mhc, n_genes, note)]))
cat("\ntop significant terms (q < 0.05), up to 5 per set x db:\n")
print(E[p.adjust < 0.05][order(set, mhc, db, pvalue)][, head(.SD, 5), by = .(set, mhc, db)][, .(set, mhc, db, Description, GeneRatio, BgRatio, p.adjust)])
cat("S20_DONE\n")
