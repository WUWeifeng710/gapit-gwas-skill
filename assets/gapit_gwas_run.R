# =====================================================================
# gapit_gwas_run.R — GAPIT 4.x GWAS 全流程模板 (mamba env: gapit)
# 依据 2026-09 实测通过的生产脚本泛化而来。
#
# 用法:
#   1) 复制本文件到工作目录, 只改 CONFIG 区
#   2) 预计运行 >15 分钟时必须后台分离运行:
#      cd ${wd} && mkdir -p gapit_input
#      nohup mamba run -n gapit Rscript gapit_gwas_run.R > gapit_input/run.log 2>&1 &
#   3) 监控: tail gapit_input/run.log ; pgrep -f "exec/R --no-echo" ; ps -o etime,%cpu,rss
#      进度标志: "markers passing QC" -> "PC created" -> 模型 banner ->
#                每性状 "has been analyzed" -> "GAPIT has done all analysis!!!" -> "DONE in"
#
# genotype_format 支持两种 (均为实测路径):
#   "rdata_basecall_hmp" : Rdata 内为变异行式 HapMap 布局 data.frame
#                          (样本名在数据第1行, 基因型为单字母碱基 A/C/G/T + '-' 缺失, 列为 factor)
#   "gd_gm_txt"          : 数值型 GD 文本(行=样本, 首列=taxa, 列=SNP, 值0/1/2, 可含NA)
#                          + GM 文本(表头 SNP Chromosome Position, 行序与 GD 标记列一致)
# 标准 HapMap 文本(SNP Chr Pos Allele A/G): 本模板不处理 —
#   小数据集直接写最小脚本 `GAPIT(Y=, G=读入的hmp data.frame, ...)` 交给 GAPIT 原生转换,
#   或先离线转成 gd_gm_txt 格式。
# =====================================================================
options(warn = 1)
t0 <- Sys.time()
suppressPackageStartupMessages(library(GAPIT))

## ============================== CONFIG ==============================
wd    <- "/path/to/workdir"                       # 工作目录 (绝对路径)
model <- "BLINK"          # GLM / MLM / CMLM / SUPER / MLMM / FarmCPU / BLINK ...
pca.total <- 3            # BLINK/FarmCPU: PCA 协变量; MLM 演示常用 1
kinship.algorithm <- "VanRaden"                   # MLM/CMLM 等用; BLINK 计算但不用
MAF.min  <- 0.05          # 位点最小 MAF (GAPIT 调用前 QC)
MISS.max <- 0.20          # 位点最大缺失率 (在填补之前计算!)
CH       <- 250000L       # rdata 分支位点分块大小
cutOff   <- 0.01          # 显著性 = cutOff / n_markers (GAPIT 约定)
MIN_COMMON <- 20L         # taxa 对齐最少共同样本数 (小数据集按需调低; 另要求重叠比例>=50%)
pheno    <- file.path(wd, "phenotype.txt")        # 首列 Taxa, 行=样本, 余列=数值性状
outdir   <- file.path(wd, paste0(tolower(model), "_results"))
memo     <- paste0(toupper(model), "_proj")

geno <- list(
  format = "rdata_basecall_hmp",
  rdata  = file.path(wd, "genotype.Rdata"), obj = "myG",
  meta_n = 11L,                                   # 前11列元信息; 样本列自 12 起
  rs = "V1", alleles = "V2", chr = "V3", pos = "V4",   # 元数据列名(按实际改)
  gd_txt = NULL, gm_txt = NULL                    # format="gd_gm_txt" 时给出
)
## ====================================================================

gdir <- file.path(wd, "gapit_input"); dir.create(gdir, showWarnings = FALSE)
# 缓存键含 QC 参数; 任一输入源比缓存新则自动失效
rds <- file.path(gdir, sprintf("gapit_input_maf%g_miss%g.rds", MAF.min, MISS.max))
rds_valid <- file.exists(rds)
if (rds_valid) {
  srcs <- Filter(function(p) !is.null(p) && file.exists(p),
                 list(pheno, geno$rdata, geno$gd_txt, geno$gm_txt))
  # 只比较实际存在的源文件, 避免 file.mtime(缺失文件)=NA 传播导致 if(!rds_valid) 崩溃
  rds_valid <- all(file.mtime(rds) > file.mtime(unlist(srcs)))
}
dir.create(outdir, showWarnings = FALSE)
log <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
BASES <- c("A","C","G","T")

# ---------- 1. 基因型 -> GD(行=样本,首列Taxa,列=标记) + GM ----------
if (!rds_valid) {  if (isTRUE(file.exists(rds))) log("cache stale (inputs newer) — re-convert")
  if (geno$format == "rdata_basecall_hmp") {
    log("Loading genotype Rdata ...")
    e <- new.env(); load(geno$rdata, envir = e); x <- get(geno$obj, envir = e)
    row1 <- as.character(unlist(x[1, ]))
    sn <- row1[(geno$meta_n + 1L):length(row1)]                  # 第1行样本名
    d <- x[-1, c(geno$rs, geno$alleles, geno$chr, geno$pos), drop = FALSE]   # 去标签行
    G <- x[-1, (geno$meta_n + 1L):ncol(x), drop = FALSE]
    rm(x, e, row1); invisible(gc())
    stopifnot(nrow(d) == nrow(G), length(sn) == ncol(G))
    for (cc in seq_len(ncol(G))) G[[cc]] <- as.character(G[[cc]]) # factor->character(水平序陷阱)

    # 双等位判定: 按 '/' 拆分后恰含两个不同 ACGT ('A/C'、'G/-/A' 均保留; indel/三连碱基剔除)
    tok  <- lapply(strsplit(as.character(d[[geno$alleles]]), "/"), function(s) unique(s[s %in% BASES]))
    okbi <- vapply(tok, length, 1L) == 2L
    log(sprintf("sites total=%d biallelic-ACGT=%d (%.1f%%)", nrow(d), sum(okbi), 100*mean(okbi)))
    d <- d[okbi, , drop = FALSE]; G <- G[okbi, , drop = FALSE]; tok <- tok[okbi]
    snmp  <- paste0(as.character(d[[geno$rs]]), ":", as.character(d[[geno$pos]]))
    chrom <- as.character(d[[geno$chr]]); pos <- suppressWarnings(as.numeric(as.character(d[[geno$pos]])))
    rm(d); invisible(gc())
    ref <- vapply(tok, `[`, "", 1L); alt <- vapply(tok, `[`, "", 2L); rm(tok); invisible(gc())
    rc <- match(ref, BASES); ac <- match(alt, BASES); rm(ref, alt); invisible(gc())

    nr <- nrow(G); D_keep <- list(); meta_keep <- list()
    for (i in seq(1, nr, by = CH)) {
      j <- min(i + CH - 1L, nr)
      gm <- as.matrix(G[i:j, , drop = FALSE])
      M <- match(gm, BASES, 0L); dim(M) <- dim(gm)         # match() 丢维度, 必须恢复!
      rm(gm)
      E  <- matrix(ac[i:j], nrow = j-i+1L, ncol = ncol(M)) # 列优先广播 (byrow=TRUE 是错的)
      Rm <- matrix(rc[i:j], nrow = j-i+1L, ncol = ncol(M))
      dos <- 2L * (M == E); dos[M == Rm] <- 0L
      dos[!(M == Rm | M == E)] <- NA_integer_              # '-' / 第三碱基 -> 缺失
      freq <- rowMeans(dos, na.rm = TRUE) / 2              # 位点在行 -> row 方向统计!
      maf  <- pmin(freq, 1 - freq)
      miss <- rowSums(is.na(dos)) / ncol(dos)
      k    <- !is.na(maf) & maf >= MAF.min & miss <= MISS.max
      D_keep[[length(D_keep) + 1L]]   <- dos[k, , drop = FALSE]
      meta_keep[[length(meta_keep) + 1L]] <-
        data.frame(SNP = snmp[i:j], Chromosome = chrom[i:j], Position = pos[i:j],
                   MAF = maf, MISS = miss, stringsAsFactors = FALSE)[k, , drop = FALSE]
      log(sprintf("sites %s-%s kept %d / %d (cum %d)", i, j, sum(k), j-i+1L,
                  sum(vapply(D_keep, nrow, 1L))))
      rm(M, E, Rm, dos, k); invisible(gc())
    }
    D  <- do.call(rbind, D_keep); MK <- do.call(rbind, meta_keep)
    rm(D_keep, meta_keep, snmp, chrom, pos, rc, ac, G); invisible(gc())
    MK$SNP <- make.unique(MK$SNP)
    log(sprintf("markers passing QC: %d (of %d biallelic)", nrow(MK), nr))

    naim <- which(is.na(D), arr.ind = TRUE)                # 位点(行)均值填补 -> 取整 0/1/2
    if (nrow(naim) > 0) { rowm <- rowMeans(D, na.rm = TRUE); D[naim] <- rowm[naim[, 1]]; D <- round(D) }
    log(sprintf("imputed %d genotype cells (%.2f%%)", nrow(naim), 100*nrow(naim)/(nrow(D)*ncol(D))))

    rownames(D) <- MK$SNP
    GD <- data.frame(Taxa = sn, t(D), check.names = FALSE, stringsAsFactors = FALSE)
    GM <- MK[, c("SNP", "Chromosome", "Position")]
    Y  <- read.table(pheno, header = TRUE, check.names = FALSE)
    saveRDS(list(GD = GD, GM = GM, Y = Y, MK = MK), rds, compress = FALSE)
    log(paste("input cached ->", basename(rds)))
    rm(GD, GM, Y, MK, D); invisible(gc())

  } else if (geno$format == "gd_gm_txt") {
    GD <- read.table(geno$gd_txt, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    GM <- read.table(geno$gm_txt, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    Y  <- read.table(pheno,        header = TRUE, check.names = FALSE)
    stopifnot("SNP" %in% colnames(GM), nrow(GM) == ncol(GD) - 1L)   # GM 行序须与 GD 标记列一致
    X <- as.matrix(GD[, -1, drop = FALSE])                          # 行=样本, 列=标记
    freq <- colMeans(X, na.rm = TRUE) / 2; maf <- pmin(freq, 1 - freq)
    miss <- 1 - colSums(!is.na(X)) / nrow(X)          # 缺失率先算 — 必须在填补之前
    k <- !is.na(maf) & maf >= MAF.min & miss <= MISS.max
    log(sprintf("markers total=%d passing QC=%d (MAF>=%g & miss<=%g)",
                ncol(X), sum(k), MAF.min, MISS.max))
    X <- X[, k, drop = FALSE]; GM <- GM[k, , drop = FALSE]
    naim <- which(is.na(X), arr.ind = TRUE)            # 位点=列 -> colMeans + naim[,2]
    if (nrow(naim) > 0) {
      colm <- colMeans(X, na.rm = TRUE)
      X[naim] <- colm[naim[, 2]]
      X <- round(X)                                    # 保持 0/1/2 整数剂量
    }
    log(sprintf("imputed %d genotype cells (%.2f%%)", nrow(naim),
                100 * nrow(naim) / max(1L, nrow(X) * ncol(X))))
    # cbind 会内部重复传 check.names -> 用 data.frame() 直构, 与 rdata 分支一致
    GD <- data.frame(Taxa = as.character(GD[[1]]), X, check.names = FALSE, stringsAsFactors = FALSE)
    MK <- NULL
    saveRDS(list(GD = GD, GM = GM, Y = Y), rds, compress = FALSE)
    log(paste("input cached ->", basename(rds)))

  } else stop('geno$format 仅支持 "rdata_basecall_hmp" | "gd_gm_txt"; 标准 HapMap 文本请直接用最小脚本 GAPIT(Y=, G=hmp_df) 或先转为 gd_gm_txt')
} else {
  log(paste("using cached input:", basename(rds)))
}

obj <- readRDS(rds)
GD <- obj$GD; GM <- obj$GM; Y <- obj$Y

# ---------- 2. 表型: NaN/NA -> 按列中位数填补 (逐性状记录数量) ----------
Yn <- as.data.frame(Y, stringsAsFactors = FALSE)
tax <- as.character(Yn[[1]]); Yv <- Yn[, -1, drop = FALSE]
for (k in seq_len(ncol(Yv))) {
  v <- suppressWarnings(as.numeric(Yv[[k]])); v[is.nan(v)] <- NA_real_
  med <- median(v, na.rm = TRUE)
  log(sprintf("trait %s: obs=%d miss=%d median=%.4g", names(Yv)[k], sum(!is.na(v)), sum(is.na(v)), med))
  v[is.na(v)] <- med; Yv[[k]] <- v
}
Y <- data.frame(Taxa = tax, Yv, check.names = FALSE)

# ---------- 3. taxa 对齐 (GD 行=样本; 强断言; 小数据集自适应阈值) ----------
gt_taxa <- as.character(GD[[1]])
common  <- intersect(gt_taxa, as.character(Y$Taxa))
frac <- length(common) / max(1L, min(length(gt_taxa), nrow(Y)))
log(sprintf("taxa: genotype=%d phenotype=%d common=%d (%.0f%%)",
            length(gt_taxa), nrow(Y), length(common), 100*frac))
if (length(common) < MIN_COMMON || frac < 0.5)
  stop(sprintf("taxa overlap insufficient (common=%d, fraction=%.2f) — 样本命名体系不同请提供映射表; 确属小数据集可降 MIN_COMMON",
               length(common), frac))
GD <- GD[match(common, gt_taxa), , drop = FALSE]
Y  <- Y[match(common, as.character(Y$Taxa)), , drop = FALSE]
mk <- colnames(GD)[-1]; GM <- GM[match(mk, GM$SNP), , drop = FALSE]
stopifnot(identical(as.character(GD[[1]]), as.character(Y[[1]])), !anyNA(GM$SNP), all(GM$SNP == mk))
log(sprintf("aligned: %d taxa x %d markers x %d traits", length(common), ncol(GD)-1L, ncol(Y)-1L))

# ---------- 4. GAPIT 运行 (输出写到 outdir; GAPIT 以 cwd 为输出目录) ----------
setwd(outdir)
log(sprintf("GAPIT %s start; output dir: %s", model, getwd()))
GS <- GAPIT(Y = Y, GD = GD, GM = GM,
            model = model, PCA.total = pca.total,
            kinship.algorithm = kinship.algorithm,
            SNP.MAF = 0,            # QC 已在外部完成, 避免二次过滤
            cutOff = cutOff, file.output = TRUE, memo = memo)
saveRDS(GS, paste0("GAPIT_", model, "_results.rds"))
log(sprintf("DONE in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
log(paste("result files:", paste(list.files(".", pattern = paste0("(?i)", model)), collapse = " ; ")))
