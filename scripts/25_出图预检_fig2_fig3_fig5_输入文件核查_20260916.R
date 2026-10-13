## =============================================================================
## 25_出图预检_fig2_fig3_fig5_输入文件核查_20260916.R
## 目的：跑 18/19/21 号图脚本之前，一次性核查全部输入文件是否在本地存在；
##       缺失的自动在 results/ 与项目根下搜索相似文件，给出改名或改路径建议。
## 运行：放在项目根目录，source("25_出图预检_fig2_fig3_fig5_输入文件核查_20260916.R")
## =============================================================================

.script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- grep("--file=", a, value = TRUE)
  if (length(f) >= 1) return(dirname(normalizePath(sub("--file=", "", f[1]))))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable() &&
      nzchar(rstudioapi::getActiveDocumentContext()$path))
    return(dirname(normalizePath(rstudioapi::getActiveDocumentContext()$path)))
  getwd()
}
.d <- .script_dir()
for (i in 1:8) {
  if (any(dir.exists(file.path(.d, c("processed", "results", "data"))))) break
  if (dirname(.d) == .d) break
  .d <- dirname(.d)
}
setwd(.d); message("[setwd] 项目根目录：", .d)

## ---- 各图所需文件清单 -------------------------------------------------------
need <- list(
  "Fig2 (18号)" = c(
    "results/moduleA/fig2a_pca_scores_long.csv",
    "results/moduleA/expr_before_batch.rds",
    "results/moduleA/expr_after_batch.rds",
    "data/sample_group.csv",
    "results/moduleB/layer1_DE_GSE55235.csv",
    "results/moduleB/fig2c_heatmap_logFC_long.csv",
    "results/moduleC/moduleC_gene_tracks.csv",
    "results/moduleB/CORE_concordant_genes.csv",
    "results/moduleB/CORE_mirror_genes.csv",
    "results/moduleB/layer2_shared_DE_genes.csv",
    "results/moduleB/layer2_sensitivity_sarcopenia_only.csv",
    "results/moduleC/enrich_ALL_tracks_significant.csv"
  ),
  "Fig3 (19号)" = c(
    "results/moduleI/moduleI_classifier_robustness.csv",
    "results/moduleI/moduleI_single_gene_AUC.csv",
    "results/moduleD/feature_votes.csv"
  ),
  "Fig3 可选（缺失自动走备路，不报错）" = c(
    "results/moduleD/shap_object.rds",
    "results/moduleI/nomogram_fit.rds"
  ),
  "Fig5 (21号)" = c(
    "results/MR/mr_summary_32tests.csv",
    "results/MR/coloc_abf_summary_24combos.csv",
    "results/MR/coloc_susie_vs_abf.csv",
    "results/MR/region_sumstats_FBLN5_alm.csv",
    "results/MR/lead_snp_info.csv",
    "fix11_LD_1000G_EUR/ld_FBLN5_alm.rds"
  )
)

## ---- 相似文件搜索（在 results/、项目根、moduleI_data/ 下按关键词模糊匹配）----
find_similar <- function(basename_pat) {
  key <- gsub("\\.(csv|rds|tsv|txt)$", "", basename(basename_pat))
  key <- gsub("_32tests|_24combos|fig2a_|fig2c_|layer1_DE_|layer2_", "", key)
  roots <- c("results", ".", "moduleI_data", "moduleE_data", "fix11_LD_1000G_EUR")
  hits <- character(0)
  for (r in roots) {
    if (!dir.exists(r)) next
    fs <- list.files(r, pattern = key, recursive = TRUE, full.names = TRUE,
                     ignore.case = TRUE)
    hits <- c(hits, fs[file.exists(fs)])
  }
  unique(hits)[1:min(3, length(unique(hits)))]
}

cat("\n================ 出图输入预检 ================\n")
all_ok <- TRUE
for (grp in names(need)) {
  cat("\n── ", grp, " ──\n", sep = "")
  for (f in need[[grp]]) {
    if (file.exists(f)) {
      cat("  ✅ ", f, "\n", sep = "")
    } else {
      all_ok <- FALSE
      cat("  ❌ 缺：", f, "\n", sep = "")
      sim <- find_similar(f)
      sim <- sim[!is.na(sim)]
      if (length(sim) > 0) {
        cat("     疑似实际文件（可考虑改名对齐，或把图脚本 CONFIG 路径改成它）：\n")
        for (sf in sim) cat("       → ", sf, "\n", sep = "")
      } else {
        cat("     未搜到相似文件——可能需回对应模块脚本补产出。\n")
      }
    }
  }
}
cat("\n============================================\n")
if (all_ok) {
  cat("全部齐备：可直接依次 source 18/19/21 号脚本。\n")
} else {
  cat("存在缺失：请按上面「疑似实际文件」提示处理后，把本输出截图发我。\n")
}

## ---- 附：results 目录结构速览（方便人工对照）--------------------------------
cat("\n[附] results/ 一级结构：\n")
print(list.dirs("results", recursive = FALSE))
