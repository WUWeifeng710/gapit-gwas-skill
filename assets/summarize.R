# =====================================================================
# summarize.R — GAPIT 结果汇总 (用法: mamba run -n gapit Rscript summarize.R <outdir> [cutOff])
# 输出: 每性状 Bonferroni 显著数 + Top5 + 跨性状共享位点, 并落盘 hits_summary.csv
# 文件名形如 GAPIT.association.<MODEL>.<memo>.<trait>(NYC|Kansas).csv
# =====================================================================
args    <- commandArgs(trailingOnly = TRUE)
outdir  <- if (length(args) >= 1) args[1] else getwd()
cutOff  <- if (length(args) >= 2) as.numeric(args[2]) else 0.01
setwd(outdir)
f <- list.files(".", pattern = "GWAS_Results.*\\.csv$")
if (!length(f)) stop("未找到 GWAS_Results*.csv — 分析是否完成? (tail ../gapit_input/run.log)")
cat("found result files:\n  "); cat(paste(basename(f), collapse = "\n  "), "\n\n")
allhits <- list()
for (i in seq_along(f)) {
  fn  <- basename(f[i])
  tag <- if (grepl("(Kansas)", fn, fixed = TRUE)) "Kansas" else "NYC"
  tr  <- sub("^.*\\.", "", sub("\\((NYC|Kansas)\\)$", "", sub("\\.csv$", "", fn)))
  d   <- read.csv(f[i], check.names = FALSE)
  pv  <- grep("^P\\.value", colnames(d), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(pv)) { cat("skip:", fn, "\n"); next }
  thr <- cutOff / nrow(d)
  o   <- d[order(d[[pv]]), ]
  n_sig <- sum(o[[pv]] < thr, na.rm = TRUE)
  cat(sprintf("=== [%d/%d] trait=%s type=%s | markers=%d | Bonf(%.2f/M)=%.2e | sig=%d ===\n",
              i, length(f), tr, tag, nrow(d), cutOff, thr, n_sig))
  keep <- intersect(c("SNP","Chromosome","Chr","Position","Pos","P.value","MAF","Effect","effect","H&B.P.Value"), colnames(o))
  print(head(o[, keep], 5)); cat("\n")
  if (n_sig > 0) {
    h <- o[which(o[[pv]] < thr), keep, drop = FALSE]
    names(h)[names(h) == pv] <- "P.value"          # 归一列名, 保证 rbind 与后续引用一致
    names(h)[names(h) == "Chromosome"] <- "Chr"; names(h)[names(h) == "Position"] <- "Pos"
    names(h)[names(h) == "effect"] <- "Effect"     # GAPIT 4.1: Kansas 文件为小写 effect
    h$trait <- tr; h$file_type <- tag
    allhits[[paste(tr, tag, sep = "|")]] <- h      # NYC/Kansas 各自入表, 不互相覆盖
  }
}
if (length(allhits)) {
  # 不同结果文件列集可能不同 (如个别缺 H&B.P.Value) -> 按并集对齐再 rbind
  allcols <- unique(unlist(lapply(allhits, colnames), use.names = FALSE))
  hits <- do.call(rbind, lapply(allhits, function(h) {
    for (cc in setdiff(allcols, colnames(h))) h[[cc]] <- NA
    h[, allcols, drop = FALSE]
  }))
  write.csv(hits, "hits_summary.csv", row.names = FALSE)
  cat("saved hits_summary.csv\n\n")
  if ("SNP" %in% colnames(hits)) {
    u  <- unique(hits[hits$file_type == "NYC", c("SNP", "trait")])   # 共享=跨【不同性状】
    tb <- table(u$SNP); sh <- names(tb[tb > 1])
    if (length(sh)) { cat("== 跨性状共享显著位点 (NYC) ==\n")
      print(hits[hits$SNP %in% sh & hits$file_type == "NYC",
                 intersect(c("SNP", "Chr", "Pos", "trait", "P.value", "MAF", "Effect"), colnames(hits))]) }
    else cat("(无跨性状共享显著位点)\n")
  }
}
ff <- "GAPIT.Association.Filter_GWAS_results.csv"
if (file.exists(ff)) { cat("\n== GAPIT 内置显著位点汇总 (Filter_GWAS_results.csv) ==\n"); print(read.csv(ff, check.names = FALSE)) }
