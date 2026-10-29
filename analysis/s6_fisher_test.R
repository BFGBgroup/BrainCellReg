# 疾病风险位点富集：Fisher 检验 + 筛选 + 清洗 + 疾病汇总，口径照贺奕捷的
# fisher2.R / filter_by_q.R / qingxi.R / keshihua.R。三套计数并列（见 s6_fisher_counts.py）。
suppressMessages(library(parallel))
OUT <- Sys.getenv("FISHER_OUT", "/public/home/yucong/scbrain/fisher_redo")
P <- "/public/home/yucong/scbrain_in/paper2/eqtl-gwas/fisher"
nc <- as.integer(Sys.getenv("NPROC", "20"))
d <- read.table(file.path(OUT, "counts.tsv"), header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
cat("rows", nrow(d), "\n")

ft <- function(i) {
  v <- c(d[i, "E.N"], d[i, "E.NN"], d[i, "NE.N"], d[i, "NE.NN"])
  if (any(v < 0) || (v[1] + v[2]) == 0) return(c(NA, NA, NA, NA))
  t <- tryCatch(fisher.test(matrix(v, ncol = 2)), error = function(e) NULL)
  if (is.null(t)) return(c(NA, NA, NA, NA))
  c(t$p.value, t$estimate[[1]], t$conf.int[1], t$conf.int[2])
}
res <- do.call(rbind, mclapply(seq_len(nrow(d)), ft, mc.cores = nc))
colnames(res) <- c("p.value", "OR", "OR.lower95", "OR.upper95")
d <- cbind(d, res)

summ <- list(); kept_all <- list()
for (v in unique(d$variant)) for (cl in unique(d$class)) {
  x <- d[d$variant == v & d$class == cl & !is.na(d$p.value), ]
  x$q.value <- p.adjust(x$p.value, method = "fdr")
  f <- subset(x, q.value <= 0.05 & OR > 1 & OR.lower95 > 1 & OR.upper95 > 1 & E.N > 0 & E.NN > 0 & NE.N > 0 & NE.NN > 0)
  rm_flag <- !is.finite(f$OR) | ((f$E.NN < 3 | f$NE.N < 3) & f$OR > 1000)
  k <- f[!rm_flag, ]
  write.table(x, file.path(OUT, paste0("fisher_all_", v, "_", cl, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(k, file.path(OUT, paste0("fisher_kept_", v, "_", cl, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  summ[[length(summ) + 1]] <- data.frame(variant = v, class = cl, tested = nrow(x), q_filtered = nrow(f), removed = sum(rm_flag), kept = nrow(k),
                                           diseases_kept = length(unique(k$disease)), diseases_total = length(unique(x$disease)),
                                           median_OR_kept = median(k$OR), stringsAsFactors = FALSE)
  if (nrow(k)) {
    k$region_name <- sub("\\..*$", "", k$ctx); k$cell_type <- sub(".*\\.", "", k$ctx)
    ds <- do.call(rbind, lapply(split(k, k$disease), function(g) data.frame(
      variant = v, class = cl, disease = g$disease[1], n_combos = nrow(g), n_regions = length(unique(g$region_name)),
      n_cell_types = length(unique(g$cell_type)), min_q = min(g$q.value), mean_q = mean(g$q.value), max_or = max(g$OR),
      median_or = median(g$OR), median_EN = median(g$E.N), score = sum(-log10(g$q.value) * log2(g$OR)) / nrow(g), stringsAsFactors = FALSE)))
    kept_all[[length(kept_all) + 1]] <- ds[order(-ds$n_combos, ds$min_q), ]
  }
}
summ <- do.call(rbind, summ)
write.table(summ, file.path(OUT, "summary_by_variant.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, kept_all), file.path(OUT, "disease_summary_all_variants.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summ)

# 和他原来的结果对一遍（神经系统疾病，orig_asis）
his <- tryCatch(read.table(file.path(P, "nervous/nervous-disease-fisher2/fisheroutzong_all.txt"), header = FALSE, fill = TRUE, stringsAsFactors = FALSE), error = function(e) NULL)
if (!is.null(his) && any(d$variant == "orig_asis")) {
  his <- his[his$V1 != "region", ]
  his$disease <- sub("fisherout\\.txt$", "", his$V10)
  for (j in 2:9) his[[j]] <- suppressWarnings(as.numeric(his[[j]]))
  mine <- d[d$variant == "orig_asis" & d$class == "nervous", ]
  m <- merge(his, mine, by.x = c("disease", "V1"), by.y = c("disease", "ctx"))
  cat("对回原结果：他", nrow(his), "行，我", nrow(mine), "行，配上", nrow(m), "行\n")
  cat("四格计数完全一致的行：", sum(m$V2 == m$E.N & m$V3 == m$E.NN & m$V4 == m$NE.N & m$V5 == m$NE.NN, na.rm = TRUE), "\n")
  ok <- !is.na(m$V6) & !is.na(m$p.value)
  cat("P 值相对差 <1e-6 的行：", sum(abs(m$V6[ok] - m$p.value[ok]) <= 1e-6 * pmax(abs(m$V6[ok]), 1e-300)), "/", sum(ok), "\n")
}
cat("FISHER_DONE\n")
