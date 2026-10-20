# OA × Sarcopenia Dual-Track Transcriptomic Atlas

Integrative multi-omics analysis of shared molecular programs between **osteoarthritis (OA)** and **sarcopenia**, built entirely on public data (GEO bulk/single-cell transcriptomes, GWAS summary statistics, eQTL resources).

## Overview

This repository accompanies the manuscript:

> *Integrated multi-omics and genetic analyses identify consistent direction-level mirror dysregulation between osteoarthritis and sarcopenia: development and cartilage-tier external validation of an osteoarthritis diagnostic biomarker panel* (prepared for submission to Clinical Rheumatology; Online Resource 1 archived under `docs/`)

**Key design features**

- **Direction-aware dual-track framework**: 99 shared genes split into 27 concordant and 72 *mirror* genes (72.7%; one-tailed binomial P = 3.5×10⁻⁶; genome-wide signed-rank ρ = −0.058, P = 5.7×10⁻⁹)
- **10 locked biomarkers**: NELL1, STEAP1, GADD45A, FBLN5, HEMK1, GDE1, RNF14, BNIP3, EGFR, CBR3 — external cartilage validation AUC = 0.818 (GSE57218); seven-classifier robustness panel (external AUC 0.788–0.887); nomogram with bootstrap-corrected C-index = 0.97
- **Triangulated evidence**: immune deconvolution (xCell / ssGSEA / CIBERSORT, independently cross-validated on the Stanford CIBERSORTx platform, 92.0% directional concordance), single-cell/snRNA localization, gene-level Mendelian randomization + colocalization (32 tests; Steiger directionality 32/32), druggable-genome overlap (5/10)

## Repository structure

```
.
├── README.md
├── LICENSE
├── scripts/            # analysis pipeline, numbered in execution order (final frozen
│                       #  versions; superseded iterations are kept in the local archive)
│   ├── setup_R_environment.R, 00–05  # environment, GEO download/QC, modules A–D
│   │                                 # (preprocessing, dual-track construction, enrichment, ML)
│   ├── 06–09            # modules E–G (external validation, immune infiltration, single-cell)
│   ├── 10–13e           # module H (MR, coloc, SuSiE diagnostics and LD-sensitivity stages)
│   ├── 11–12e           # module I (druggable genome, nomogram, classifier robustness)
│   ├── 13–14c           # modules J–K (TF–miRNA network, CellChat)
│   ├── 15/15b           # sessionInfo → supplementary Table S10 (15b = frozen 2026-10-08 run)
│   ├── 16–21            # ClinRheum-style main figures 1–5 (final via versions)
│   ├── 22–31            # audit, prediction rebuild, figure-input preparation and panel assembly
│   ├── 22_sfig1–38      # supplementary figure redraws and composite layout
│   ├── 23–26_moduleB/I  # LODO robustness, donor-overlap merge correction, coefficient extraction
│   ├── 39–44            # archive-path audit, Table 1 assembly, GSE168505 prediction rebuild,
│   │                    # MAF/SuSiE convention audit
│   └── fix14v4          # CIBERSORTx mixture preparation (Stanford platform validation)
├── tables/             # frozen supplementary tables (S11, S12, S13-Appendix, S16, S19, S20,
│                       #  S20-Appendices), module A status, CIBERSORT two-implementation
│                       #  comparison, group definitions
├── reference/          # local KEGG snapshot (2026-08-23) used for enrichment reproducibility
└── docs/               # Online Resource (ESM_1: S-Fig 1–13 legends, Tables S1–S20 notes,
                        #  TRIPOD+AI 2024 and STROBE-MR 2021 checklists, Appendix D)
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
