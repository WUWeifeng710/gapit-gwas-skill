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
# genotype_format 选择:
#   "rdata_basecall_hmp" : Rdata 内为变异行式 HapMap 布局 data.frame
#                          (样本名在数据第1行, 基因型为单字母碱基 A/C/G/T + '-' 缺失, 列为 factor)
#                          —— 2026-09 全量实测路径 (385万位点×320样本)
#   "gd_gm_txt"          : 已是数值型 GD 文本(行=样本, 首列=taxa, 列=SNP, 值0/1/2)
#                          + GM 文本(SNP Chromosome Position, 与 GD 列同序)
#   "hapmap_txt"         : 标准 HapMap 文本(SNP Chr Pos Allele A/G + 0|1|2/ACGT 混合码)
#                          —— 交给 GAPIT 原生 G= 转换; 非标准编码请用第一个分支
# =====================================================================
options(warn = 1)
t0 <- Sys.time()
suppressPackageStartupMessages(library(GAPIT))

## ============================== CONFIG ==============================
wd    <- "/path/to/workdir"                       # 工作目录 (绝对路径)
model <- "BLINK"          # GLM / MLM / CMLM / SUPER / MLMM / FarmCPU / BLINK ...
pca.total <- 3            # BLINK/FarmCPU: PCA 协变量; MLM 演示常用 1
kinship.algorithm <- "VanRaden"                   # MLM/CMLM 等用; BLINK 计算但不用
MAF.min  <- 0.05          # 位点最小 MAF (GAPIT 前 QC)
MISS.max <- 0.20          # 位点最大缺失率
CH       <- 250000L       # rdata 分支的位点分块大小
cutOff   <- 0.01          # 显著性 = cutOff / n_markers (GAPIT 约定)
pheno    <- file.path(wd, "phenotype.txt")        # 首列 Taxa, 行=样本, 余列=数值性状
outdir   <- file.path(wd, paste0(tolower(model), "_results"))
memo     <- paste0(toupper(model), "_proj")
rds      <- file.path(wd, "gapit_input", "gapit_input.rds")   # 输入缓存(重跑免转换)

geno <- list(
  format = "rdata_basecall_hmp",
  rdata  = file.path(wd, "genotype.Rdata"), obj = "myG",
  meta_n = 11L,                                   # 前11列元信息; 样本列自 12 起
  rs = "V1", alleles = "V2", chr = "V3", pos = "V4",   # 元数据列名(按实际改)
  gd_txt = NULL, gm_txt = NULL,                   # format="gd_gm_txt" 时给出
  hapmap = file.path(wd, "genotype.hmp.txt")      # format="hapmap_txt" 时使用
)
## ====================================================================

dir.create(dirname(rds), showWarnings = FALSE)
dir.create(outdir, showWarnings = FALSE)
log <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
BASES <- c("A","C","G","T")

# ---------- 1. 基因型 -> GD(行=样本,首列Taxa,列=标记) + GM ----------
if (!file.exists(rds)) {
  if (geno$format == "rdata_basecall_hmp") {
    log("Loading genotype Rdata ...")
    e <- new.env(); load(geno$rdata, envir = e); x <- get(geno$obj, envir = e)
    sn <- as.character(unlist(x[1, ])); sn <- sn[(geno$meta_n + 1L):length(sn)]   # 第1行样本名
    d <- x[-1, c(geno$rs, geno$alleles, geno$chr, geno$pos), drop = FALSE]        # 去标签行
    G <- x[-1, (geno$meta_n + 1L):ncol(x), drop = FALSE]
    rm(x, e); invisible(gc())
    stopifnot(nrow(d) == nrow(G), length(sn) == ncol(G))
    for (cc in seq_len(ncol(G))) G[[cc]] <- as.character(G[[cc]])   # factor->character(水平序陷阱)

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
      M <- match(gm, BASES, 0L); dim(M) <- dim(gm)          # match() 丢维度, 必须恢复!
      rm(gm)
      E  <- matrix(ac[i:j], nrow = j-i+1L, ncol = ncol(M))  # 列优先广播 (byrow=TRUE 是错的)
      Rm <- matrix(rc[i:j], nrow = j-i+1L, ncol = ncol(M))
      dos <- 2L * (M == E); dos[M == Rm] <- 0L
      dos[!(M == Rm | M == E)] <- NA_integer_               # '-' / 第三碱基 -> 缺失
      freq <- rowMeans(dos, na.rm = TRUE) / 2               # 位点在行 -> row 方向统计!
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

    naim <- which(is.na(D), arr.ind = TRUE)                 # 位点均值填补 -> 取整 0/1/2
    if (nrow(naim) > 0) { rowm <- rowMeans(D, na.rm = TRUE); D[naim] <- rowm[naim[, 1]]; D <- round(D) }
    log(sprintf("imputed %d genotype cells (%.2f%%)", nrow(naim), 100*nrow(naim)/(nrow(D)*ncol(D))))

    rownames(D) <- MK$SNP
    GD <- data.frame(Taxa = sn, t(D), check.names = FALSE, stringsAsFactors = FALSE)
    GM <- MK[, c("SNP", "Chromosome", "Position")]
    Y  <- read.table(pheno, header = TRUE, check.names = FALSE)
    saveRDS(list(GD = GD, GM = GM, Y = Y, MK = MK), rds, compress = FALSE)
    log("input cached -> gapit_input.rds")
    rm(GD, GM, Y, MK, D); invisible(gc())

  } else if (geno$format == "gd_gm_txt") {
    GD <- read.table(geno$gd_txt, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    GM <- read.table(geno$gm_txt, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    Y  <- read.table(pheno,        header = TRUE, check.names = FALSE)
    # 外部 MAF/缺失 QC: 剂量矩阵为 样本x标记 -> 位点统计沿列方向
    X <- as.matrix(GD[, -1, drop = FALSE])
    miss <- 1 - colSums(!is.na(X)) / nrow(X)        # !! 必须先算缺失率, 填补后 NA 消失
    freq <- colMeans(X, na.rm = TRUE) / 2; maf <- pmin(freq, 1 - freq)
    k <- !is.na(maf) & maf >= MAF.min & miss <= MISS.max
    log(sprintf("markers total=%d passing QC=%d", ncol(X), sum(k)))
    X <- X[, k, drop = FALSE]; GM <- GM[k, , drop = FALSE]
    if (anyNA(X)) {                                 # 位点均值填补(列方向) -> 取整 0/1/2
      idx  <- which(is.na(X), arr.ind = TRUE)
      colm <- colMeans(X, na.rm = TRUE)
      X[idx] <- colm[idx[, 2]]; X <- round(X)
      log(sprintf("genotype NA site-mean imputed: %d cells", nrow(idx)))
    }
    GD <- cbind(GD[, 1, drop = FALSE], X); MK <- NULL
    saveRDS(list(GD = GD, GM = GM, Y = Y), rds, compress = FALSE)

  } else if (geno$format == "hapmap_txt") {
    G <- read.table(geno$hapmap, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    Y <- read.table(pheno, header = TRUE, check.names = FALSE)
    stopifnot(ncol(G) > 12)                         # rs alleles chrom pos ... + samples
    saveRDS(list(hapmap = TRUE, G = G, Y = Y), rds, compress = FALSE)
    rm(G); invisible(gc())
    log("hapmap route: cached raw HapMap; passing G= to GAPIT directly (standard coding only, no external QC)")
  } else stop("unknown geno$format")
}

obj <- readRDS(rds)
is_hmp <- isTRUE(obj$hapmap)
if (!is_hmp) { GD <- obj$GD; GM <- obj$GM }
Y <- obj$Y

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

# ---------- 3. taxa 对齐 (GD 行=样本; 强断言) ----------
if (!is_hmp) {
  gt_taxa <- as.character(GD[[1]])
  common  <- intersect(gt_taxa, as.character(Y$Taxa))
  log(sprintf("taxa: genotype=%d phenotype=%d common=%d", length(gt_taxa), nrow(Y), length(common)))
  if (length(common) < max(20, floor(0.5 * length(gt_taxa))))   # 小数据集(全重叠<50)不误杀
    stop("taxa overlap too small — 检查样本命名/映射表")
  GD <- GD[match(common, gt_taxa), , drop = FALSE]
  Y  <- Y[match(common, as.character(Y$Taxa)), , drop = FALSE]
  mk <- colnames(GD)[-1]; GM <- GM[match(mk, GM$SNP), , drop = FALSE]
  stopifnot(identical(as.character(GD[[1]]), as.character(Y[[1]])), !anyNA(GM$SNP), all(GM$SNP == mk))
  log(sprintf("aligned: %d taxa x %d markers x %d traits", length(common), ncol(GD)-1L, ncol(Y)-1L))
} else {
  G <- obj$G
  gt_taxa <- colnames(G)[-(1:11)]
  common  <- intersect(gt_taxa, as.character(Y$Taxa))
  log(sprintf("taxa: genotype=%d phenotype=%d common=%d (GAPIT aligns G/Y internally)",
              length(gt_taxa), nrow(Y), length(common)))
  if (length(common) < max(20, floor(0.5 * length(gt_taxa))))
    stop("taxa overlap too small — 检查样本命名/映射表")
}

# ---------- 4. GAPIT 运行 (输出写到 outdir; GAPIT 以 cwd 为输出目录) ----------
setwd(outdir)
log(sprintf("GAPIT %s start; output dir: %s", model, getwd()))
if (is_hmp) {
  GS <- GAPIT(Y = Y, G = G,                       # HapMap 原生物件, GAPIT 内部转 GD/GM
              model = model, PCA.total = pca.total,
              kinship.algorithm = kinship.algorithm,
              SNP.MAF = 0,                        # hapmap 路线无外部 QC; 需要过滤请在此设 SNP.MAF
              cutOff = cutOff, file.output = TRUE, memo = memo)
} else {
  GS <- GAPIT(Y = Y, GD = GD, GM = GM,
              model = model, PCA.total = pca.total,
              kinship.algorithm = kinship.algorithm,
              SNP.MAF = 0,            # QC 已在外部完成, 避免二次过滤
              cutOff = cutOff, file.output = TRUE, memo = memo)
}
saveRDS(GS, paste0("GAPIT_", model, "_results.rds"))
log(sprintf("DONE in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
log(paste("result files:", paste(list.files(".", pattern = paste0("(?i)", model)), collapse = " ; ")))
