## ============================================================================
## 39_补充表源文件归档路径核对_20261008.R
## 用途：核对补充材料映射总表 B 登记的全部源文件在本机的存在性（重点 moduleH/I/J/K）
## 运行：在项目根目录（ENG论文）下 source("39_补充表源文件归档路径核对_20261008.R")
## 输出：控制台逐项 ✅/❌ 报告；❌ 项自动全盘递归搜索提示实际位置
## ============================================================================
message("== 脚本版本：39_20261008_归档路径核对版 ==")

## ---- 0. 根目录自检 ---------------------------------------------------------
root <- getwd()
message("工作目录：", root)
if (!dir.exists(file.path(root, "results"))) {
  stop("当前目录下没有 results/ 文件夹——请先 setwd() 到项目根目录（ENG论文）再运行。")
}

## ---- 1. 映射总表 B 登记的显式路径 ------------------------------------------
targets <- c(
  ## S6/S8/S17 —— moduleH
  "results/moduleH/moduleH_MR_all_methods.csv",
  "results/moduleH/moduleH_known_causal_overlap.csv",
  "results/moduleH/moduleH_MR_summary.csv",
  "results/moduleH/moduleH_coloc_summary.csv",
  "results/moduleH/fix11_susie/fix11_susie_summary.csv",
  "results/moduleH/fix11_susie/fix11e_susie_summary_convergence.csv",
  "results/moduleH/fix11e_susie_summary_convergence.csv",   # 备选位置（映射表写法有歧义，两个候选都试）
  ## S11/S12/S13/S18 —— moduleI
  "results/moduleI/moduleI_single_gene_AUC.csv",
  "results/moduleI/moduleI_classifier_robustness.csv",
  "results/moduleI/moduleI_nomogram_validate.txt",
  "results/moduleI/moduleI_druggable_overlap.csv",
  "results/moduleI/moduleI_dsigdb_enrichment.csv",
  ## S14 —— moduleJ
  "results/moduleJ/moduleJ_TFmiRNA_edges.csv",
  "results/moduleJ/moduleJ_TFmiRNA_hub.csv",
  ## S15 —— moduleK
  "results/moduleK/moduleK_biomarker_LR_pairs.csv",
  "results/moduleK/moduleK_GSE169454_compare_summary.csv",
  "results/moduleK/moduleK_GSE152805_compare_summary.csv",
  "results/moduleK/moduleK_GSE167186_compare_summary.csv",
  ## S10（已生成，顺带复核）
  "results/moduleS10/S10_tools_versions.csv",
  "results/moduleS10/S10_sessionInfo_full.txt",
  ## 其他登记路径（顺带复核）
  "results/moduleB/layer2_shared_DE_genes.csv",
  "results/moduleB/CORE_concordant_genes.csv",
  "results/moduleB/CORE_mirror_genes.csv",
  "results/moduleD_v2/moduleD_v2_biomarker_table.csv",
  "moduleA_status_2026-08-22_final.csv"
)

## 封版/订正表（文件名搜索，位置不限——它们可能在项目根目录或其他文件夹）
sealed_names <- c(
  "表S11_单基因诊断效能_方向固定订正版_20260916.csv",
  "表S12_七分类器稳健性汇总_封版_20261008.csv",
  "表S13附_主模型完整系数表_20261003.csv",
  "表S16_21数据集总表_最终工作版_20261003.csv",
  "表S19_LODO全队列矩阵_20261003.csv",
  "表S20_供体重叠合并校正汇总_20261003.csv",
  "表S20附_供体重叠合并校正_OA侧RRA全表_20261003.csv",
  "表S20附_研究级留一_20261003.csv",
  "fix15_双版对照_CIBERSORTx官方_vs_dicepro_20260903.csv"
)

## ---- 2. 逐项核对（❌ 项自动递归搜索实际位置） ------------------------------
all_files <- list.files(root, recursive = TRUE, full.names = TRUE)
find_elsewhere <- function(basename_) {
  hits <- all_files[basename(all_files) == basename_]
  if (length(hits) == 0) return(NULL)
  sub(paste0("^", gsub("([./\\\\])", "\\\\\\1", root), "[/\\\\]?"), "", hits)
}

n_ok <- 0L; n_bad <- 0L
message("\n---- A. 显式路径逐项核对 ----")
for (t in targets) {
  if (file.exists(file.path(root, t))) {
    message("✅ ", t); n_ok <- n_ok + 1L
  } else {
    alt <- find_elsewhere(basename(t))
    if (is.null(alt)) {
      message("❌ ", t, "  << 全盘未找到该文件名")
    } else {
      message("⚠️ ", t, "  << 路径不符，实际位于：", paste(alt, collapse = " | "))
    }
    n_bad <- n_bad + 1L
  }
}

message("\n---- B. 封版/订正表（按文件名全盘搜索）----")
for (nm in sealed_names) {
  alt <- find_elsewhere(nm)
  if (is.null(alt)) {
    message("❌ ", nm, "  << 未找到"); n_bad <- n_bad + 1L
  } else {
    message("✅ ", nm, "  @ ", paste(alt, collapse = " | ")); n_ok <- n_ok + 1L
  }
}

message("\n---- C. 汇总 ----")
message("✅ 在位：", n_ok, " 项；⚠️/❌ 待处理：", n_bad, " 项")
if (n_bad == 0L) message("全部核对通过——可回报「6.1 归档路径核对闭环」。") else message("请把上方 ⚠️/❌ 行贴回，逐项定位。")
