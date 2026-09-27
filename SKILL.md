---
name: gapit-gwas-skill
description: End-to-end GAPIT GWAS in an isolated conda environment - auto-checks mamba/conda and GAPIT install (incl. the tricky multtest/snpStats/EMMREML deps), probes and converts genotype (HapMap text/GD+GM/Rdata/VCF) and phenotype inputs with QC, confirms model and parameters with the user, runs GLM/MLM/CMLM/FarmCPU/BLINK with safe background execution, then summarizes hits and writes an interpretation report. Use this skill whenever a GWAS/GAPIT analysis is requested, so that environment setup and input-format pitfalls never need to be rediscovered.
---