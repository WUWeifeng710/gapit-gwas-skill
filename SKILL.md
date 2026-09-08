---
name: gapit-gwas
description: End-to-end GAPIT GWAS in an isolated conda environment - auto-checks mamba/conda and GAPIT install (incl. the tricky multtest/snpStats/EMMREML deps), probes and converts genotype (HapMap text/GD+GM/Rdata/VCF) and phenotype inputs with QC, confirms model and parameters with the user, runs GLM/MLM/CMLM/FarmCPU/BLINK with safe background execution, then summarizes hits and writes an interpretation report. Use this skill whenever a GWAS/GAPIT analysis is requested, so that environment setup and input-format pitfalls never need to be rediscovered.
---

# GAPIT GWAS (Environment to Report, End-to-End)

## Overview

Runs a complete **GAPIT v4.x** GWAS in the dedicated mamba environment `gapit`, from environment verification to a final interpretation report. Idempotent: every step verifies what exists instead of reinstalling.

Steps:
- **Step 0** – Verify/repair the `gapit` conda environment (one command: `bash ${SKILL_DIR}/assets/verify_env.sh`).
- **Step 1** – Probe inputs (metadata only, NEVER dump large data), detect format, prepare GAPIT-ready GD+GM+Y with QC via the tested template `assets/gapit_gwas_run.R`.
- **Step 2** – Confirm analysis parameters with the user. **Never guess silently.**
- **Step 3** – Run GAPIT; long runs detached with a log, monitored via progress markers.
- **Step 4** – Summarize per-trait hits (`assets/summarize.R`), then write the mandatory report.

`${SKILL_DIR}` = `~/.dsh/skills/34.gapit-gwas`. All analysis outputs go to `${out_dir}` (default `./<model>_results/`); intermediates/cache to `${wd}/gapit_input/`. GAPIT writes results **into the current working directory** — the template handles `setwd`.

---

## Step 0 – Environment verify / repair

### 0.1 Quick check (do this first, every time)
```bash
bash ~/.dsh/skills/34.gapit-gwas/assets/verify_env.sh; echo "exit=$?"
```
Exit 0 → go to Step 1. Interpretation of failures:
| exit | meaning | action |
|---|---|---|
| 2 | no conda/mamba | **stop, inform user** (do not self-install a package manager) |
| 3 | env `gapit` missing | create per 0.2 |
| 4 | deps missing/broken (incl. GAPIT not installed) | repair listed packages per 0.3/0.4; if the R report shows `[MISSING] GAPIT`, install GAPIT per 0.5 |
| 5 | `library(GAPIT)` fails at load time | a load-time dep broken (0.3) — GAPIT is installed but cannot load |

If mamba errors mention lockfiles / "no writable cache directory", that is a **sandbox file-policy denial** on the mamba root — request escalated permissions once with a one-sentence justification; do not detour.

### 0.2 Create environment (known-good spec, tested 2026-09)
```bash
mamba create -n gapit -c conda-forge -c bioconda -y \
  "r-base=4.3" "bioconductor-multtest=2.56.0=r43*" \
  r-ape r-bigmemory r-genetics r-gplots r-htmltools r-magrittr r-lme4 \
  r-plotly r-rcpparmadillo r-scatterplot3d r-snowfall r-remotes \
  gcc_linux-64 gxx_linux-64 gfortran_linux-64 make
mamba run -n gapit Rscript -e 'install.packages("EMMREML", repos="https://cloud.r-project.org")'
mamba install -n gapit -c conda-forge -c bioconda -y \
  "bioconductor-snpstats=1.52.0=r43*" r-biocmanager
```
**Why this exact shape**: GAPIT needs R>=4.3; bioconda binaries for `multtest`/`snpstats` top out at r43 builds (multtest was removed from newer Bioconductor — binary r43 avoids the legacy Biobase source-tarball chain); `EMMREML` has no conda build (Fortran source → env compilers); `r-rcpparmadillo/lme4/plotly` as binaries avoid long compiles.

### 0.3 Repair individual gaps
- `multtest` / `snpStats` / `Biobase` → mamba binaries (bioconda, r43 pin) — **never** let R auto-install them.
- `EMMREML` → CRAN source inside env (needs the compilers installed above).
- any `r-*` → `mamba install -n gapit -c conda-forge r-<pkg>`; a package that exists but won't load (ABI mismatch) must be reinstalled via mamba, not from CRAN.
- Verify loadability of **all 17 deps** (verify_env.sh does exactly this list).

### 0.4 Drift fallbacks
- Check current builds: `mamba search --override-channels -c conda-forge -c bioconda bioconductor-multtest`. If r44/r45 binaries appear upstream, you may move to a newer r-base pin accordingly.
- If multtest is unavailable as a binary: create env at r4.3 anyway and `R CMD INSTALL` the Bioconductor 3.19 source tarballs of **Biobase then multtest** (GAPIT README's own recommendation).
- Solver hangs on `defaults` repodata → add `--override-channels -c conda-forge -c bioconda` to every command.

### 0.5 Install GAPIT itself
```bash
cd ${wd} && unzip -o GAPIT-master.zip -d _gapit_src     # local zip preferred
mamba run -n gapit R CMD INSTALL --no-multiarch _gapit_src/GAPIT-master
```
No local copy? Try `git clone --depth 1 https://github.com/jiabowang/GAPIT.git`, then `curl -4 .../codeload.github.com/jiabowang/GAPIT/tar.gz/refs/heads/master`, else **ask the user to download the zip** (server git→github.com frequently times out; raw.githubusercontent often still works).
`R CMD INSTALL` failing during lazy loading with a BiocManager/auto-install trace = a load-time dep is missing → back to 0.3. GAPIT's `GAPIT.0000.R` calls `.gapit_require_or_install()` at load time for gplots, genetics, ape, compiler, grid, bigmemory, EMMREML, scatterplot3d, lme4, multtest, **snpStats** — conda R's staged-install makes that auto-install path fail, so pre-installing everything is mandatory, not optional.

---

## Step 1 – Input inspection & preparation

### 1.1 Required inputs (ask user for paths if unclear)
1. **Genotype** — HapMap text / numeric GD+GM text pair / Rdata(rds) object / VCF.
2. **Phenotype** — tab/CSV: first column = Taxa, one row per sample, remaining columns numeric traits (NaN/NA allowed).
3. Working dir + output dir name (default `<model>_results`).

### 1.2 Metadata-only probes (hard rule: never print bulk genotype data)
```bash
ls -lh ${pheno}; head -3 ${pheno}; wc -l ${pheno}
free -g; nproc; df -h ${wd} | tail -1
```
```r
## Rdata/rds genotype probe — layout, coding census, taxa match. Adapt obj name.
e <- new.env(); nm <- load(geno_file, envir = e); x <- get(nm[1], envir = e)
cat("obj:", nm, "class:", class(x), "dim:", dim(x), "\n")
print(head(as.character(unlist(x[1, ][1:min(16, ncol(x))])), 16))  # header row inside data?
cat("col classes:", paste(sapply(x[1:5, ], class), collapse=","), "\n")
h <- as.vector(as.matrix(x[sample(2:nrow(x), 3000), sample(1:ncol(x), 40)]))
print(sort(table(h), decreasing = TRUE))            # call alphabet: ACGT? 0/1/2? 0|1? B/D? '-'?
y <- read.table(pheno, header = TRUE, check.names = FALSE)
cat("phenotype taxa:", nrow(y), "traits:", ncol(y)-1, "\n")
```
Decide `geno$format` in the run script:
| probe result | route |
|---|---|
| Rdata data.frame of variant-rows, header in row 1, factor sample cols, calls A/C/G/T + `-` | `"rdata_basecall_hmp"` (**fully tested route**, incl. `G/-/A` style alleles) |
| numeric GD text (col1=taxa) + GM text (SNP/Chr/Pos) | `"gd_gm_txt"` |
| standard HapMap text (`A/G` alleles, calls incl. het `B`) | `"hapmap_txt"` → GAPIT-native `G=` conversion |
| VCF | convert first (bcftools/mglr in-env) to gd_gm_txt; log MAF/missing during conversion |

### 1.3 Preparation + execution — use the tested template
```bash
cp ~/.dsh/skills/34.gapit-gwas/assets/gapit_gwas_run.R ${wd}/gapit_gwas_run.R
# edit ONLY the CONFIG block: wd, model, pca.total, kinship.algorithm,
# MAF.min/MISS.max, cutOff, pheno path, geno list (format + paths/columns), outdir, memo
```
The template encapsulates (do not re-implement ad hoc):
- allele parsing: split on `/`, keep sites with **exactly two distinct ACGT tokens** (so `A/C` and `G/-/A` kept; tri-allelic & indel-containing dropped), counts logged;
- dosage: ref=0 / alt=2 (alt = second distinct base); `-`/third-base → NA → site-mean impute → round to 0/1/2, % imputed logged;
- QC: MAF>=0.05 & site-missing<=0.20 defaults, per-chunk kept counts logged;
- unique IDs `rs:pos` + `make.unique`;
- phenotype NaN → per-trait median with obs/miss counts logged;
- strict taxa alignment: intersect → reorder GD **rows** and Y rows identically, GM reordered by `match` to GD marker columns, triple `stopifnot` guard;
- `setwd(outdir)` then `GAPIT(...)` with `SNP.MAF=0` (external QC done), then `saveRDS(GS)`;
- input cache `gapit_input/gapit_input.rds`: re-runs with different model/params skip conversion entirely.

**Coding pitfalls baked into the template — re-verify if you must edit it:**
1. GAPIT **4.x**: `GAPIT(Y=, GD=, GM=)` — `GM` is the genetic **map only**; v3-style `GM=hapmap` fails with "Kinship has to be provided or estimated from genotype!!!".
2. `GD` = data.frame, taxa column first; **do not** `as.matrix()` the whole mixed df (character coercion → `apply(GD,2,sum)` "invalid 'type' (character)").
3. R `match()` drops matrix dims → restore `dim(M) <- dim(gm)`.
4. Per-site alt/ref broadcasting: column-major `matrix(v, nrow=S, ncol=N)`; **`byrow=TRUE` is wrong** (silently corrupts dosage).
5. Per-site stats along the SITE axis of the working orientation: dosage matrix is sites×taxa inside the conversion → `rowMeans/rowSums`; after `t(D)` into GD, sites are columns.
6. Factor genotype columns → convert to character before `as.matrix` (level-order trap).
7. Sample names may live in the **header row stored as data** (row 1), and rs columns may lack coordinates (duplicate IDs).
8. Result-object layout drifts across versions — probe `names(GS)`; `GS$GWAS` is a data.frame in 4.1.0.

---

## Step 2 – Parameter confirmation (ask the user)

Present a compact table and get explicit confirmation (defaults recommended):
| Param | Default | Notes |
|---|---|---|
| `model` | **must ask** | GLM/MLM/CMLM/SUPER/MLMM/FarmCPU/BLINK; BLINK&FarmCPU → PCA covariates; MLM family → kinship |
| traits | all numeric cols | show per-trait n_obs / n_missing from probe before confirming |
| `PCA.total` | 3 (BLINK/FarmCPU), 1–3 (MLM) | |
| `kinship.algorithm` | VanRaden | ignored by BLINK association (still computed) |
| QC | MAF ≥ 0.05, missing ≤ 0.20 | looser missing for GBS/RAD; state implied marker retention if known |
| `cutOff` | 0.01 | implied Bonferroni = 0.01 / n_markers — report it |
| `out_dir` | `./<model>_results` | |

---

## Step 3 – Run GAPIT (detached for long runs)

**Runtime budget** (measured: 320 taxa × 3.53 M markers × 5 traits, BLINK, 64 cores / 251 GB): 181 min total; per-trait association 25–40 min (first trait slowest), multi-trait Circle/High-resolution Manhattan tail adds ~30 min AFTER all per-trait CSVs are saved. Scale ~linearly in markers × traits. Peak RAM ≈ several × markers×taxa×8 B (measured ~125 GB): check `free -g`; if headroom < 2× raise MAF or LD-prune first.
**Rule**: expected wall time > 15 min → run detached:
```bash
cd ${wd} && mkdir -p gapit_input && nohup mamba run -n gapit Rscript gapit_gwas_run.R > gapit_input/run.log 2>&1 &
echo "pid $!"   # wrapper exits; the R child persists
```
Monitor: `pgrep -f "exec/R --no-echo"` + `ps -o etime,%cpu,rss -p <Rpid>` + `tail gapit_input/run.log`. Progress markers in order: `markers passing QC` → `PC created`/`Kinship created` → model banner → per-trait `has been analyzed successfully` → `GAPIT has done all analysis!!!` → script's own `DONE in X min` + `result files:`. A silent stretch of 10–30 min at the end = plotting tail, NOT a hang — check newest mtimes in `outdir` before touching anything; never kill mid-trait. If all `GWAS_Results` CSVs exist, association is done regardless of remaining figures.

---

## Step 4 – Summarize & report

```bash
mamba run -n gapit Rscript ~/.dsh/skills/34.gapit-gwas/assets/summarize.R ${out_dir} 0.01 | tee ${out_dir}/summary_console.txt
```
Produces: per-trait marker count / Bonferroni threshold / significant count / top-5 (SNP, Chr, Pos, P, MAF, Effect, H&B P), `hits_summary.csv`, cross-trait shared significant loci, and the GAPIT-internal `Filter_GWAS_results.csv` echo.

Then write **`${out_dir}/GAPIT_<model>_report.md`** with the six mandatory sections:
1. **Tools & versions** — mamba env name, R, GAPIT, key deps, this skill.
2. **Parameters & rationale** — model choice, covariates, QC (with retained/dropped counts from the log), imputation stats (genotype cells %, per-trait phenotype miss), cutOff & implied Bonferroni.
3. **Commands** — verify/create env commands actually used, full CONFIG of `gapit_gwas_run.R`, launch & monitor commands, runtime.
4. **Results interpretation** — per-trait significance table, effect directions, cross-trait concordant loci (strongest evidence), QQ/NYC sanity notes; state "0 significant, top suggestive P = …" explicitly where applicable.
5. **Biological meaning** — locus positions in the reference build used, candidate intervals (±50–100 kb), trait-design caveats (time-series/saturation effects, structure control adequacy).
6. **Conclusions & next steps** — candidate-gene annotation (GFF proximity), second-method cross-check (e.g., FarmCPU vs BLINK), multi-environment BLUP phenotype rerun, LD-pruning for speed, artifact inventory incl. cache file sizes for cleanup.

---

## Pre-delivery checklist
- [ ] `verify_env.sh` exit 0 (incl. `library(GAPIT)` load)
- [ ] probe logged: dims, call alphabet, taxa overlap (no silent drops)
- [ ] QC & imputation counts logged and quoted in the report
- [ ] alignment `stopifnot` passed
- [ ] >15 min run launched detached; final `DONE in` line observed, not assumed
- [ ] `summarize.R` output captured; `GAPIT_<model>_report.md` written
- [ ] user told: output dir contents, cache location/size, cleanup options
