# gapit-gwas-skill

End-to-end **GAPIT v4.x GWAS** skill for AI agents (dsh / WorkBuddy / Claude-Code-style skill loaders) — from a bare Linux server to a finished interpretation report, in a dedicated, self-verifying mamba environment.

It encodes the hard-won lessons of a real production run (320 taxa × 3.53 M markers × 5 traits, BLINK, 64 cores / 251 GB RAM): the tricky Bioconductor dependency pins, the GAPIT 4.x API pitfalls, and the silent data-corruption traps that cost hours to debug — so they never need to be rediscovered.


## What it does

| Step | Action                                                                                                                             | Asset                     |
| ---- | ---------------------------------------------------------------------------------------------------------------------------------- | ------------------------- |
| 0    | Verify / repair the `gapit` conda environment (17 deps + GAPIT load test, one command, meaningful exit codes)                      | `assets/verify_env.sh`    |
| 1    | Probe genotype & phenotype inputs (metadata only — never dumps bulk data), detect format, build GAPIT-ready GD+GM+Y with logged QC | `assets/gapit_gwas_run.R` |
| 2    | Confirm model & parameters with the user (never guesses silently)                                                                  | —                         |
| 3    | Run GLM / MLM / CMLM / SUPER / MLMM / FarmCPU / BLINK; long runs detached with a log and progress markers                          | `assets/gapit_gwas_run.R` |
| 4    | Summarize per-trait hits (Bonferroni, top-5, cross-trait shared loci) and write a mandatory 6-section interpretation report        | `assets/summarize.R`      |

Every step is **idempotent** — it verifies what exists instead of reinstalling or reconverting.

## Repository layout

```
gapit-gwas-skill/
├── SKILL.md                  # The skill: full playbook the agent follows
└── assets/
    ├── verify_env.sh         # Step 0: one-command env check (exit 0 = ready)
    ├── gapit_gwas_run.R      # Steps 1+3: tested conversion + QC + GAPIT template
    └── summarize.R           # Step 4: per-trait hit summary -> hits_summary.csv
```


## Installation (as an agent skill)

This is first and foremost an **agent skill**: a `SKILL.md` playbook plus tested `assets/` that an LLM agent loads and follows. The easiest install — just hand the repo URL to your agent and ask:

> "Install this skill for me: `https://github.com/<your-account>/gapit-gwas-skill`"

The agent knows where its own skills directory lives (e.g. `~/.dsh/skills/` for dsh, `~/.workbuddy/skills/` for WorkBuddy, `~/.claude/skills/` for Claude Code–style loaders) and will clone it there.

Two things to know after installing:

1. **`${SKILL_DIR}` is hard-coded** in `SKILL.md` as `~/.dsh/skills/gapit-gwas` (three call sites: Step 0.1, Step 1.3, Step 4). If the agent installs under a different name/path, tell it to point those at the actual location — it can do the edit itself.
2. The YAML frontmatter (`name`, `description`) is what lets the agent **auto-trigger** on any GWAS/GAPIT request. Keep it intact.

No server-side dependencies are needed at install time — the skill verifies and repairs the `gapit` mamba environment on first use (Step 0).


## Using it with an agent

Once installed, just describe your analysis in natural language — the description in the frontmatter routes GWAS requests to this skill:

> "Run a GAPIT GWAS on my data using BLINK model. Genotype: `/data/genotype.Rdata` (object `myG`, variant-rows HapMap layout, sample names in row 1). Phenotype: `/data/phenotype.txt`. Work directory: `/data/gwas_run`."

The agent will then, in order:

1. **Step 0** — run `assets/verify_env.sh`; repair or create the `gapit` mamba environment if needed (it stops and asks rather than self-installing a package manager);
2. **Step 1** — probe your inputs (metadata only, never dumps bulk genotype data), detect the format, and prepare the conversion template;
3. **Step 2** — come back to you with a parameter table (model, traits, `PCA.total`, kinship algorithm, QC thresholds, `cutOff`, output dir) and **wait for your explicit confirmation** — it never guesses silently;
4. **Step 3** — launch GAPIT detached with a log, monitor progress markers, and refuse to kill mid-trait runs (the 10–30 min silent plotting tail is not a hang);
5. **Step 4** — summarize per-trait hits and write the mandatory 6-section `GAPIT_<model>_report.md`, then tell you where the outputs and caches are and how to clean up.

What you should have ready: paths to genotype + phenotype files, a rough idea of the genotype format (the agent probes to confirm), and enough free RAM/disk (see Requirements).

## Requirements

- Linux with **mamba** (or conda) — the skill stops and asks if neither exists; it never self-installs a package manager.
- Sufficient RAM: peak usage ≈ several × (markers × taxa × 8 bytes). The reference run peaked at ~125 GB; check `free -g` and raise MAF or LD-prune first if headroom < 2×.
- The GAPIT source (a local `GAPIT-master.zip`, or a shallow clone of [jiabowang/GAPIT](https://github.com/jiabowang/GAPIT)).


## Quick start (manual use, no agent)

The assets also work standalone — handy for debugging or for running on a server without an agent:

```bash
# 0. Verify the environment (exit 0 = ready; 2/3/4/5 tell you exactly what to fix)
bash assets/verify_env.sh; echo "exit=$?"

#    If exit=3, create the known-good environment (tested 2026-09):
mamba create -n gapit -c conda-forge -c bioconda -y \
  "r-base=4.3" "bioconductor-multtest=2.56.0=r43*" \
  r-ape r-bigmemory r-genetics r-gplots r-htmltools r-magrittr r-lme4 \
  r-plotly r-rcpparmadillo r-scatterplot3d r-snowfall r-remotes \
  gcc_linux-64 gxx_linux-64 gfortran_linux-64 make
mamba run -n gapit Rscript -e 'install.packages("EMMREML", repos="https://cloud.r-project.org")'
mamba install -n gapit -c conda-forge -c bioconda -y \
  "bioconductor-snpstats=1.52.0=r43*" r-biocmanager
mamba run -n gapit R CMD INSTALL --no-multiarch /path/to/GAPIT-master

# 1+3. Copy the template, edit ONLY the CONFIG block, launch detached
cp assets/gapit_gwas_run.R /path/to/workdir/
$EDITOR /path/to/workdir/gapit_gwas_run.R
cd /path/to/workdir && mkdir -p gapit_input
nohup mamba run -n gapit Rscript gapit_gwas_run.R > gapit_input/run.log 2>&1 &

#    Monitor (a silent 10-30 min stretch at the end = plotting tail, NOT a hang)
tail -f gapit_input/run.log

# 4. Summarize
mamba run -n gapit Rscript assets/summarize.R /path/to/workdir/gapit_results 0.01
```


## Supported input formats

Set `geno$format` in the CONFIG block:

| `format`             | Input                                                                                                                       | Route                                                                                                                             |
| -------------------- | --------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------- |
| `rdata_basecall_hmp` | Rdata data.frame, variant-rows HapMap layout (sample names in data row 1, calls = A/C/G/T + `-`, factor columns)            | **Fully tested** (3.85 M sites × 320 samples). Handles `G/-/A`-style allele strings; tri-allelic & indel sites dropped and logged |
| `gd_gm_txt`          | Numeric GD text (rows = taxa, first col = Taxa, values 0/1/2) + GM text (SNP/Chromosome/Position, same order as GD columns) | External MAF/missingness QC, then per-site mean imputation                                                                        |
| `hapmap_txt`         | Standard HapMap text                                                                                                        | Passed to GAPIT natively via `G=` (no external QC — set `SNP.MAF` in the GAPIT call if filtering is needed)                       |
| VCF                  | Convert to GD+GM first (bcftools / mglr in-env), then use `gd_gm_txt`                                                       | —                                                                                                                                 |

The phenotype file: first column `Taxa`, one row per sample, remaining columns numeric traits (`NaN`/`NA` allowed — imputed per-trait to the median, counts logged).

## Key design decisions (why it's built this way)

- **r43 binary pins for `multtest`/`snpStats`**: multtest was removed from newer Bioconductor; bioconda r43 binaries avoid the legacy Biobase source-tarball chain. Never let R auto-install these — GAPIT's load-time `.gapit_require_or_install()` fails under conda R's staged install, so pre-installing everything is mandatory.
- **External QC, `SNP.MAF = 0` inside GAPIT**: filtering happens once, in the template, with kept/dropped counts logged — no silent double filtering.
- **Chunked conversion with an input cache**: `gapit_input/gapit_input.rds` lets you re-run with a different model/covariates in minutes instead of re-converting millions of sites. *Changing QC thresholds or input files requires deleting this cache first.*
- **Hard guards over silent coercion**: triple `stopifnot` on taxa/marker alignment, dimension restoration after `match()`, column-major per-site broadcasting — the known silent-corruption traps are asserted, not assumed.

## Runtime expectations

Measured on 320 taxa × 3.53 M markers × 5 traits, BLINK, 64 cores: **181 min total** (per-trait association 25–40 min, first trait slowest; multi-trait Manhattan/QQ plotting adds ~30 min after all CSVs are saved). Scales ~linearly in markers × traits. Anything > 15 min should run detached via `nohup`.

## Output

- `<model>_results/` — GAPIT CSVs, Manhattan/QQ plots, `GAPIT_<model>_results.rds`
- `hits_summary.csv` + `summary_console.txt` — significant hits and cross-trait shared loci
- `GAPIT_<model>_report.md` — the mandatory report: tools & versions, parameters & rationale (with QC/imputation counts), exact commands, results interpretation, biological meaning, conclusions & next steps

## Changelog


### 2026-09-08 — review & bugfix pass

- README: added agent installation & usage sections (install = hand the repo URL to your agent; `${SKILL_DIR}` caveat; natural-language trigger example); corrected the citation to GAPIT Version 4 (Wang & Zhang 2026, *Mol. Biol. Evol.*, msag107).
- Fixed broken `hapmap_txt` route (never cached its object → guaranteed `readRDS` crash; missing `geno$hapmap` config field; wrong GM slice). Now caches the raw HapMap object and calls `GAPIT(Y=, G=)` natively.
- Fixed `gd_gm_txt` QC: missingness was computed **after** imputation (always 0 → missing-rate filter dead) and imputation used per-sample means instead of per-site means. Now: QC first (column-wise), then per-site mean imputation.
- Fixed taxa-overlap guard that rejected small datasets (< 50 taxa) even at 100 % overlap.
- `summarize.R`: cross-trait shared-locus table no longer hard-codes the `P.value` column name.
- `SKILL.md`: added the missing Step 3 heading, corrected the exit-code table (GAPIT-not-installed exits 4, not 5), documented the input-cache staleness rule.

## License

Use freely within your lab / team. GAPIT itself is © its authors — see [jiabowang/GAPIT](https://github.com/jiabowang/GAPIT). If you publish an analysis that used this skill, report the program version and cite GAPIT Version 4:

> Wang, J., Zhang, Z. (2026). GAPIT Version 4: Integration of GWAS into Genomic Prediction. *Molecular Biology and Evolution*, msag107. <https://doi.org/10.1093/molbev/msag107>
