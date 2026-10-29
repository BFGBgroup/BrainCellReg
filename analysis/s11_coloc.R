# 共定位，口径照贺奕捷 1-mk-gene-new.sh + 2-coloc.sh + 2-coloc_analysis_new.R：
#   eQTL 一侧：该组合里这个基因的全部 nominal 结果（基因、SNP、到 TSS 的距离、P、beta）配队列 MAF，type="quant"，N=134；
#   GWAS 一侧：<编号>-coloc.tsv 的 P、b、se²，type="cc"；按 SNP 取交集，少于 10 个不跑；coloc.abf 用默认先验。
# 和他的差别只在跑法：GWAS 表每份只读一次（他是每个基因读一次），结果应当一致，用没换 GWAS 的 3 种病核对。
# 用法: Rscript s11_coloc.R <smr-sig-gene.txt> <eqtl目录> <输出目录> [核数]
suppressMessages({library(coloc); library(data.table); library(parallel)})
a <- commandArgs(trailingOnly = TRUE)
sigf <- a[1]; eqdir <- a[2]; out <- a[3]; nc <- if (length(a) >= 4) as.integer(a[4]) else 8L
GW <- "/public/home/yucong/scbrain/gwas10"
MAFF <- "/public/home/yucong/scbrain_in/paper2/single-cell-eqtl/data/genotype/MAF.txt"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
cat("coloc version", as.character(packageVersion("coloc")), "\n")
sig <- fread(sigf, header = FALSE, col.names = c("ctx", "gwas", "gene"))
sig <- unique(sig)
maftab <- fread(MAFF, select = c(2, 8), col.names = c("snps", "maf"))
maftab <- maftab[!duplicated(snps, fromLast = TRUE)]   # 他的 awk 是后出现的覆盖先出现的

one_gwas <- function(g) {
  s <- sig[gwas == g]
  gw <- fread(file.path(GW, paste0(g, "-coloc.tsv")), select = c("SNP", "b", "se", "p"))
  gw[, varbeta := se^2]
  setkey(gw, SNP)
  res <- list()
  for (cx in unique(s$ctx)) {
    ef <- file.path(eqdir, paste0(cx, ".tsv"))
    if (!file.exists(ef)) { cat("NO_EQTL", cx, "\n"); next }
    eq <- fread(ef)
    # nominal 结果里有少数 SNP 没有 P 值（tensorQTL 对该基因没算出来）；他转出的文本里没有这些行，这里同样去掉
    eq <- eq[!is.na(pvalue)]
    eq[, maf := maftab$maf[match(snps, maftab$snps)]]
    for (ge in s[ctx == cx, gene]) {
      e1 <- eq[genes == ge]
      input <- merge(e1, gw, by.x = "snps", by.y = "SNP", all = FALSE)
      r <- data.table(ctx = cx, gwas = g, gene = ge, n_eqtl = nrow(e1), nsnps = nrow(input), PP.H0.abf = NA_real_, PP.H1.abf = NA_real_,
                      PP.H2.abf = NA_real_, PP.H3.abf = NA_real_, PP.H4.abf = NA_real_, note = "")
      if (nrow(input) < 10) { r$note <- "Too few SNPs for coloc"; res[[length(res) + 1]] <- r; next }
      cr <- tryCatch({
        invisible(capture.output(x <- suppressWarnings(suppressMessages(
          coloc.abf(dataset1 = list(snp = input$snps, pvalues = input$pvalue, beta = input$beta, MAF = input$maf, type = "quant", N = 134),
                    dataset2 = list(snp = input$snps, pvalues = input$p, beta = input$b, varbeta = input$varbeta, type = "cc"))))))
        x
      }, error = function(e) conditionMessage(e))
      if (is.character(cr)) { r$note <- paste("ERROR:", cr) } else {
        sm <- cr$summary
        r$nsnps <- sm[["nsnps"]]; r$PP.H0.abf <- sm[["PP.H0.abf"]]; r$PP.H1.abf <- sm[["PP.H1.abf"]]
        r$PP.H2.abf <- sm[["PP.H2.abf"]]; r$PP.H3.abf <- sm[["PP.H3.abf"]]; r$PP.H4.abf <- sm[["PP.H4.abf"]]
      }
      res[[length(res) + 1]] <- r
    }
  }
  rbindlist(res)
}
rr <- mclapply(unique(sig$gwas), function(g) tryCatch(one_gwas(g), error = function(e) { cat("GWAS_FAIL", g, conditionMessage(e), "\n"); NULL }), mc.cores = nc)
res <- rbindlist(rr)
fwrite(res, file.path(out, "coloc_all.tsv"), sep = "\t")
for (k in unique(paste(res$ctx, res$gwas))) {
  x <- res[paste(ctx, gwas) == k & note == ""]
  if (nrow(x)) fwrite(x[, .(gene, nsnps, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf, PP.H4.abf)],
                      file.path(out, paste0(x$ctx[1], "-", x$gwas[1], "-coloc-result.txt")), sep = "\t")
}
cat("triplets", nrow(sig), "done", nrow(res), "ok", sum(res$note == ""), "too_few", sum(res$note == "Too few SNPs for coloc"),
    "error", sum(grepl("^ERROR", res$note)), "\n")
if (any(grepl("^ERROR", res$note))) print(head(res[grepl("^ERROR", note), .(ctx, gwas, gene, note)], 5))
cat("PP.H4 >= 0.8:", sum(res$PP.H4.abf >= 0.8, na.rm = TRUE), " >= 0.5:", sum(res$PP.H4.abf >= 0.5, na.rm = TRUE), "\n")

# 有他的原结果就对一遍（只在核对模式下有意义）
HIS <- "/public/home/yucong/scbrain_in/paper2/eqtl-gwas/coloc"
cmp <- rbindlist(lapply(unique(paste(res$ctx, res$gwas)), function(k) {
  p <- strsplit(k, " ")[[1]]; f <- file.path(HIS, p[1], paste0(p[1], "-", p[2], "-coloc-result.txt"))
  if (!file.exists(f)) return(NULL)
  h <- tryCatch(fread(f), error = function(e) NULL); if (is.null(h) || !nrow(h)) return(NULL)
  h <- h[!duplicated(gene)]; h[, `:=`(ctx = p[1], gwas = p[2])]; h
}), fill = TRUE)
if (nrow(cmp)) {
  m <- merge(res, cmp, by = c("ctx", "gwas", "gene"), suffixes = c("", ".his"))
  m <- m[!is.na(PP.H4.abf) & !is.na(as.numeric(PP.H4.abf.his))]
  d <- abs(m$PP.H4.abf - as.numeric(m$PP.H4.abf.his))
  cat("对回原结果：配上", nrow(m), "行；nsnps 相同", sum(m$nsnps == m$nsnps.his), "；PP.H4 差 <1e-6 的", sum(d < 1e-6), "；最大差", max(d), "\n")
  fwrite(m, file.path(out, "compare_with_his.tsv"), sep = "\t")
}
cat("COLOC_DONE\n")
