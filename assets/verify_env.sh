#!/usr/bin/env bash
# verify_env.sh — gapit-gwas 技能 Step 0 一键检查
# 用法: bash verify_env.sh   (退出码 0 = 环境就绪; 非 0 = 需按 SKILL.md Step 0 修复)
set -u
MGR=$(command -v mamba || command -v conda) || { echo "[FAIL] no conda/mamba found"; exit 2; }
echo "[ok] conda mgr: $MGR"
if ! $MGR env list | awk '$1=="gapit"{f=1} END{exit !f}'; then
  echo "[NEED-CREATE] env 'gapit' not found — 按 SKILL.md 0.2 创建"; exit 3
fi
echo "[ok] env 'gapit' exists"
$MGR run -n gapit Rscript -e '
ok <- TRUE
cat("R:", R.version.string, "\n")
deps <- c("ape","bigmemory","EMMREML","genetics","gplots","htmltools","magrittr",
          "lme4","MASS","methods","multtest","plotly","RcppArmadillo",
          "scatterplot3d","snowfall","snpStats","BiocManager")
miss <- deps[!vapply(deps, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) { cat("[MISSING]", paste(miss, collapse=", "), "\n"); ok <- FALSE }
v <- tryCatch(as.character(packageVersion("GAPIT")), error = function(e) NA)
if (is.na(v)) { cat("[MISSING] GAPIT\n"); ok <- FALSE } else cat("GAPIT:", v, "\n")
if (ok) cat("ENV-READY\n") else quit(status = 1)' || { echo "[FAIL] env incomplete — 按 SKILL.md 0.3/0.4 修复"; exit 4; }
# GAPIT 可加载 + 导出检查 (触发 .gapit_require_or_install 全部隐式依赖)
$MGR run -n gapit Rscript -e 'suppressPackageStartupMessages(library(GAPIT)); stopifnot("GAPIT" %in% getNamespaceExports("GAPIT") | length(getNamespaceExports("GAPIT"))>0); cat("LOAD-OK exports:", length(getNamespaceExports("GAPIT")), "\n")' || { echo "[FAIL] library(GAPIT) load error"; exit 5; }
echo "[PASS] gapit environment ready"
