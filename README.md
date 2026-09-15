# OA × Sarcopenia Dual-Track Transcriptomic Atlas

Integrative multi-omics analysis of shared molecular programs between **osteoarthritis (OA)** and **sarcopenia**, built entirely on public data (GEO bulk/single-cell transcriptomes, GWAS summary statistics, eQTL resources).

## Overview

This repository accompanies the manuscript:

> *Integrative multi-omics and genetic analysis reveals shared molecular axes and inter-tissue mirror dysregulation between osteoarthritis and sarcopenia* (submission in preparation)

**Key design features**

- **Direction-aware dual-track framework**: 99 shared genes split into 27 concordant and 72 *mirror* genes (72.7%; one-tailed binomial P = 3.5×10⁻⁶; genome-wide signed-rank ρ = −0.058, P = 5.7×10⁻⁹)
- **10 locked biomarkers**: NELL1, STEAP1, GADD45A, FBLN5, HEMK1, GDE1, RNF14, BNIP3, EGFR, CBR3 — external cartilage validation AUC = 0.818 (GSE57218); seven-classifier robustness panel (external AUC 0.788–0.887); nomogram with bootstrap-corrected C-index = 0.97
- **Triangulated evidence**: immune deconvolution (xCell / ssGSEA / CIBERSORT, independently cross-validated on the Stanford CIBERSORTx platform, 92.0% directional concordance), single-cell/snRNA localization, gene-level Mendelian randomization + colocalization (32 tests; Steiger directionality 32/32), druggable-genome overlap (5/10)

## Repository structure

```
.
├── README.md
├── LICENSE
├── scripts/            # analysis pipeline, numbered in execution order
│   ├── 00–11 ...      # environment setup, GEO download/QC, modules A–H
│                      # (to be added from the analysis workstation before archiving)
│   ├── 12_moduleI_nomogram_单基因AUC_分类器稳健性_20260913.R
│   ├── 13_moduleJ_TFmiRNA调控网络_20260914.R
│   ├── 13b_moduleJ_补丁_TF库重试与TRRUST兜底_20260914.R
│   ├── 14c_moduleK_CellChat_分组标签与presto修正_20260914.R
│   └── 15_sessionInfo_S10工具版本表_20260915.R
└── docs/              # supplementary checklists (TRIPOD+AI / STROBE-MR)
```

## Data availability

All 21 datasets are publicly available from GEO (accession numbers and group definitions: Table 1 of the manuscript; accessed August 2026). GWAS summary statistics, eQTL resources (eQTLGen, GTEx) and the druggable-genome list are public; sources and versions are documented in Methods 2.1.3–2.1.4. No proprietary or restricted-access data are used.

## Reproducing the analysis

1. **Environment**: R 4.6.0 (Windows or Linux). Package versions are frozen in supplementary Table S10 — regenerate with `scripts/15_sessionInfo_S10工具版本表_20260915.R`.
2. **Pipeline order**: scripts are numbered in execution order (`00_setup` → `01_GEO_download_check` → modules A–K). Each module writes to `results/moduleXX/` and logs to a dated log file; hard self-checks (sample seals, gene coverage, direction consistency) halt the run loudly rather than emitting silent empty output.
3. **Seeds**: every stochastic step uses an explicitly fixed seed (documented per module; e.g., classifier seeds 20260914–20260921).
4. **License-restricted components**: the CIBERSORT source code and LM22 signature matrix are *not* redistributed here (Stanford license). Obtain them from the original authors; the independent CIBERSORTx validation was run on the official Stanford online platform.

## Honest-reporting conventions

Selection leakage in internal CV AUC is disclosed in-text (feature selection precedes CV); internal metrics are reference-only, evidential weight rests on external validation. The bootstrap-corrected calibration slope (0.394) is reported as overfitting disclosure. Aging-proxy datasets (GSE167186) are never interpreted as sarcopenia evidence. Zero/negative results are reported rather than dropped.

## License

Code: MIT (see LICENSE). Data remain subject to their original GEO depositor terms.

## Citation

【To be completed upon publication / 投稿后回填】
