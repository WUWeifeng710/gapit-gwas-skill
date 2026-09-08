# =====================================================================
# summarize.R — GAPIT 结果汇总 (用法: mamba run -n gapit Rscript summarize.R <outdir> [cutOff])
# 输出: 每性状 Bonferroni 显著数 + Top5 + 跨性状共享位点, 并落盘 hits_summary.csv
# =====================================================================
args    <- commandArgs(trailingOnly = TRUE)
outdir  <- if (length(args) >= 1) args[1] else getwd()
cutOff  <- if (length(args) >= 2) as.numeric(args[2]) else 0.01
setwd(outdir)
f <- list.files(".", pattern = "GWAS_Results.*\\.csv$")
if (!length(f)) stop("未找到 GWAS_Results*.csv — 分析是否完成? (tail ../gapit_input/run.log)")
cat("found result files:\n "); cat(paste(basename(f), collapse = "\n  "), "\n\n")
allhits <- list()
for (i in seq_along(f)) {
  tr <- sub(".*GWAS_Results\\..*?\\.([^.]+(\\.[^.]+)*)\\((NYC|Kansas)\\)\\.csv", "\\1", f[i])
  d  <- read.csv(f[i], check.names = FALSE)
  pv <- grep("^P\\.value", colnames(d), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(pv)) { cat("skip:", f[i], "\n"); next }
  thr <- cutOff / nrow(d)
  o <- d[order(d[[pv]]), ]
  n_sig <- sum(o[[pv]] < thr, na.rm = TRUE)
  cat(sprintf("=== [%d/%d] trait=%s | markers=%d | Bonf(%.2f/M)=%.2e | sig=%d ===\n",
              i, length(f), tr, nrow(d), cutOff, thr, n_sig))
  keep <- intersect(c("SNP","Chromosome","Chr","Position","Pos", pv, "MAF","Effect","H&B.P.Value"), colnames(o))
  print(head(o[, keep], 5)); cat("\n")
  if (n_sig > 0) { h <- o[which(o[[pv]] < thr), keep, drop = FALSE]; h$trait <- tr
                   allhits[[tr]] <- h }
}
if (length(allhits)) {
  hits <- do.call(rbind, allhits)
  write.csv(hits, "hits_summary.csv", row.names = FALSE)
  cat("saved hits_summary.csv\n\n")
  if ("SNP" %in% colnames(hits)) {
    tb <- table(hits$SNP)
    sh <- tb[tb > 1]
    if (length(sh)) { cat("== 跨性状共享显著位点 ==\n")
      pcol <- intersect(c("P.value", "p.value", "P.Value"), colnames(hits))[1]
      if (is.na(pcol)) pcol <- grep("^P\\.value", colnames(hits), value = TRUE, ignore.case = TRUE)[1]
      print(hits[hits$SNP %in% names(sh), c("SNP", "trait", pcol)]) }
  }
}
ff <- "GAPIT.Association.Filter_GWAS_results.csv"
if (file.exists(ff)) { cat("\n== GAPIT 内置显著位点汇总 (Filter_GWAS_results.csv) ==\n"); print(read.csv(ff, check.names = FALSE)) }
