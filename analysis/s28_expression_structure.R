# Expression structure of the deconvolved profiles, recomputed with Heyijie's code (handoff code.txt) on the deconvolved
# profiles (braineac_combined_seurat.RData, 13,541 "sample x cell type" profiles):
#   c  hierarchical clustering of the profiles (Euclidean distance; his dist_matrix.RData when it matches, else recomputed)
#   d  UMAP: 2,000 HVGs (vst), ScaleData, PCA, UMAP on PCs 1-20 (code.txt l.1040-1073)
#   e  CellChat, Secreted Signaling, grouped by Tissue and by CellType; over-expressed signalling genes from FindAllMarkers
#      (adj P < 0.05, log2FC > 1), triMean, min.cells = 0 (code.txt l.1874-2100). Gene symbols from the local Ensembl 82 GTF
#      (his code used biomaRt online; compute nodes have no internet).
# Outputs (tables only; figures are drawn by fig2cde_plot.R): /public/home/yucong/scbrain/fig2cde/
suppressPackageStartupMessages({ library(Seurat); library(CellChat); library(data.table); library(dplyr) })
set.seed(42)
H <- "/public/home/yucong/scbrain/local_handoff/scbrain2/BayesPrism_InstaPrism_results"
OUT <- "/public/home/yucong/scbrain/fig2cde"; dir.create(OUT, showWarnings = FALSE)
load(file.path(H, "braineac_combined_seurat.RData"))   # -> combined_seurat
so <- combined_seurat; rm(combined_seurat)
md <- so@meta.data
cat("profiles", ncol(so), "genes", nrow(so), "| meta columns:", paste(colnames(md), collapse = ","), "\n")
stopifnot(all(c("Tissue", "CellType") %in% colnames(md)))
fwrite(data.table(id = colnames(so), Tissue = md$Tissue, CellType = md$CellType), file.path(OUT, "profiles.tsv"), sep = "\t")

## c ---------------------------------------------------------------------------------------------------------
X <- as.matrix(GetAssayData(so, layer = "data"))   # SeuratObject 5: layer, not slot
d <- NULL
e <- new.env(); ok <- tryCatch({ load(file.path(H, "dist_matrix.RData"), envir = e); TRUE }, error = function(x) FALSE)
if (ok) {
  obj <- mget(ls(e), envir = e)
  cat("dist_matrix.RData objects:", paste(names(obj), sapply(obj, function(z) paste(class(z), collapse = "/")), collapse = "; "), "\n")
  cand <- Filter(function(z) inherits(z, "dist") && attr(z, "Size") == ncol(so), obj)
  if (length(cand)) d <- cand[[1]]
}
if (is.null(d)) { cat("recomputing Euclidean distance on the expression profiles\n"); d <- dist(t(X)) }
hc <- hclust(d)
ord <- hc$order
fwrite(data.table(id = colnames(so)[ord], pos = seq_along(ord)), file.path(OUT, "hc_order.tsv"), sep = "\t")
# block-averaged distance matrix in clustering order, for a raster heatmap (13,541^2 cells cannot be drawn one by one)
D <- as.matrix(d)[ord, ord]; rm(d)
nb <- 600; br <- cut(seq_len(nrow(D)), nb, labels = FALSE)
B <- matrix(NA_real_, nb, nb)
for (i in seq_len(nb)) { ri <- which(br == i); B[i, ] <- tapply(colMeans(D[ri, , drop = FALSE]), br, mean) }
rm(D)
fwrite(as.data.table(B), file.path(OUT, "hc_dist_blocks.tsv"), sep = "\t", col.names = FALSE)
cat("c done: profiles", length(ord), "blocks", nb, "\n")

## d ---------------------------------------------------------------------------------------------------------
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)
so <- ScaleData(so, verbose = FALSE)
so <- RunPCA(so, features = VariableFeatures(so), npcs = 100, verbose = FALSE)
so <- RunUMAP(so, reduction = "pca", dims = 1:20, seed.use = 42, verbose = FALSE)
um <- Embeddings(so, "umap")
fwrite(data.table(id = rownames(um), UMAP_1 = um[, 1], UMAP_2 = um[, 2], Tissue = so$Tissue, CellType = so$CellType), file.path(OUT, "umap.tsv"), sep = "\t")
fwrite(data.table(pc = 1:100, sdev = so[["pca"]]@stdev), file.path(OUT, "pca_sdev.tsv"), sep = "\t")
cat("d done\n")

## e ---------------------------------------------------------------------------------------------------------
gtf <- fread("/public/home/yucong/scbrain_in/paper2/single-cell-eqtl/data/expression/Homo_sapiens.GRCh38.82.chr.gtf.txt", header = TRUE)
setnames(gtf, 1:7, c("chr", "start", "end", "strand", "gene", "symbol", "biotype"))
map <- setNames(gtf$symbol, gtf$gene)
cnt <- GetAssayData(so, layer = "counts"); dat <- GetAssayData(so, layer = "data")
sym <- map[rownames(dat)]; sym[is.na(sym) | sym == ""] <- rownames(dat)[is.na(sym) | sym == ""]
keep <- !duplicated(sym)
cat("genes", length(sym), "| mapped to symbols", sum(sym != rownames(dat)), "| kept after removing duplicate symbols", sum(keep), "\n")
cnt <- cnt[keep, ]; dat <- dat[keep, ]; rownames(cnt) <- sym[keep]; rownames(dat) <- sym[keep]
so2 <- CreateSeuratObject(counts = cnt, meta.data = so@meta.data[, c("Tissue", "CellType")])
so2 <- SetAssayData(so2, layer = "data", new.data = dat)
rm(cnt, dat, X)
run_cc <- function(group) {
  cc <- createCellChat(object = so2, group.by = group, assay = "RNA")
  cc@DB <- subsetDB(CellChatDB.human, search = "Secreted Signaling", key = "annotation")
  cc <- subsetData(cc)
  Idents(so2) <- group
  sg <- rownames(cc@data.signaling)
  mk <- FindAllMarkers(so2, features = sg, only.pos = TRUE, verbose = FALSE)
  mk <- mk[, c("cluster", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
  names(mk) <- c("clusters", "features", "pvalues", "logFC", "pct.1", "pct.2", "pvalues.adj")
  sel <- mk %>% filter(pvalues.adj < 0.05 & logFC > 1)
  cc@var.features$features <- sel$features
  cc@var.features$features.info <- sel
  cc <- subsetData(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean")
  cc <- filterCommunication(cc, min.cells = 0)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  n_genes <- length(unique(sel$features)); n_lr <- nrow(cc@LR$LRsig)
  cat("CellChat by", group, ": over-expressed signalling genes", n_genes, "| over-expressed L-R interactions", n_lr, "\n")
  for (k in c("count", "weight")) {
    m <- cc@net[[k]]
    fwrite(data.table(source = rep(rownames(m), ncol(m)), target = rep(colnames(m), each = nrow(m)), value = as.vector(m)),
           file.path(OUT, sprintf("cellchat_%s_%s.tsv", group, k)), sep = "\t")
  }
  fwrite(data.table(group = names(table(cc@idents)), n = as.integer(table(cc@idents))), file.path(OUT, sprintf("cellchat_%s_groupsize.tsv", group)), sep = "\t")
  fwrite(data.table(by = group, overexpressed_genes = n_genes, overexpressed_LR = n_lr), file.path(OUT, sprintf("cellchat_%s_summary.tsv", group)), sep = "\t")
  saveRDS(cc, file.path(OUT, sprintf("cellchat_%s.rds", group)))
}
run_cc("Tissue")
run_cc("CellType")
cat("manuscript: by cell type 662 genes / 975 interactions; by brain region 166 genes / 57 interactions\n")
cat("FIG2CDE_DONE\n")
