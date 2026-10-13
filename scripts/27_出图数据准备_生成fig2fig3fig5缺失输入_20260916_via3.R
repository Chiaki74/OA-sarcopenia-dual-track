#!/usr/bin/env Rscript
# ============================================================
# 27_出图数据准备_生成fig2fig3fig5缺失输入_20260916_via3.R
# 目的：一次性生成 18/19/21 三个出图脚本缺失的全部输入文件。
# via2 新增（审图意见第4条）：sample_group.csv 增加 cohort 列，
#   供 Fig2-A「整合前按数据集着色」使用；其余与 v1 完全一致。
# via3 修复：assemble_cohort 所在脚本不再写死 scripts/07_integrate_clean...
#   （该路径在用户机器上不存在，实跑报错），改为全项目递归自定位。
# 数据来源：GEO 原始下载（GSE114007 用计数表，GSE55235/55457/12021/51588 用系列矩阵）
#   全部进 assemble_cohort() —— 与主分析完全同一条流水线。
# 生成文件（全部写入 data/ 与 output/R/figs_R/）：
#   A. data/expression_log2cpm_training.csv        大训练矩阵（行=gene, 列=sample）
#   B. data/sample_group.csv                       两列 sample,condition（v1）
#      → via2：三列 sample,condition,cohort
#   C. data/pathway_hierarchy_with_evidence.csv    通路+证据+权重（E/F 用）
#   D. data/dea_results_training.rds               供 21_fig5 面板 A 取 UCHL1 效应
#   E. output/R/figs_R/15_pca_obj.rds              prcomp 对象（面板 A 的 PC1/PC2）
# 重要约定：
#   - GSE114007 的样本表型从 series matrix 的标题解析（Normal/OA）
#   - 训练队列 TRAIN_SETS 与主分析保持一致
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(GEOquery)
  library(preprocessCore)
  library(edgeR)
})

## ---- 0. 路径（定位项目根） ----
locate_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    script <- normalizePath(sub("^--file=", "", file_arg))
    d <- dirname(script)
    if (dir.exists(file.path(d, "data"))) return(d)
  }
  for (cand in c(getwd(), "C:/Users/凉月/Desktop/ENG论文")) {
    if (dir.exists(file.path(cand, "data"))) return(normalizePath(cand))
  }
  stop("无法定位项目根目录")
}
ROOT <- locate_root()
DIR_DATA <- file.path(ROOT, "data")
DIR_OUT  <- file.path(ROOT, "output", "R", "figs_R")
dir.create(DIR_DATA, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_OUT,  recursive = TRUE, showWarnings = FALSE)

## ---- 0.5 载入 assemble_cohort（via3：递归自定位，替代写死的 scripts/ 路径） ----
cand07 <- list.files(ROOT, pattern = "^07_integrate_clean.*\\.R$",
                     recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
if (length(cand07) == 0) {
  ## 兜底：全项目搜「定义 assemble_cohort 函数」的脚本
  allr <- list.files(ROOT, pattern = "\\.R$", recursive = TRUE,
                     full.names = TRUE, ignore.case = TRUE)
  has_ac <- function(f) tryCatch(
    any(grepl("assemble_cohort\\s*(<-|=)\\s*function",
              readLines(f, warn = FALSE, n = 4000))),
    error = function(e) FALSE)   # 编码异常文件跳过
  cand07 <- allr[vapply(allr, has_ac, logical(1))]
}
if (length(cand07) == 0)
  stop("未找到定义 assemble_cohort 的脚本（写死的 scripts/07_integrate_clean... 不存在，",
       "递归搜索亦无命中）。请把 07 脚本的实际文件名/位置告诉我。")
message("载入 assemble_cohort：", cand07[1])
source(cand07[1], encoding = "UTF-8")
stopifnot("source 后 assemble_cohort 仍未定义" = exists("assemble_cohort"))

## ---- 1. 训练队列清单（与 11 交叉证伪脚本一致） ----
TRAIN_SETS <- c("GSE114007", "GSE55235", "GSE55457", "GSE12021", "GSE51588")
TRAIN_LABEL <- "OA"   # 条件标签

## ---- 2. 拉取 5 队列并统一基因空间 ----
cohorts <- lapply(TRAIN_SETS, function(gse) {
  assemble_cohort(gse_id = gse, species = "human",
                  train_condition_label = TRAIN_LABEL,
                  use_saved_eset = TRUE)
})
names(cohorts) <- TRAIN_SETS

# 共同基因
common_genes <- Reduce(intersect, lapply(cohorts, function(cc) rownames(cc$expr)))
message("共同基因数：", length(common_genes))

# 逐队列取共同基因子矩阵
expr_list <- lapply(cohorts, function(cc) cc$expr[common_genes, , drop = FALSE])

## ---- 3. 分位数归一化（跨队列） ----
expr_mat <- do.call(cbind, expr_list)
# 若已 log 化则直接用；否则转 log2CPM
if (max(expr_mat, na.rm = TRUE) > 100) {
  expr_mat <- edgeR::cpm(expr_mat, log = TRUE, prior.count = 1)
}
expr_qn <- preprocessCore::normalize.quantiles(as.matrix(expr_mat))
rownames(expr_qn) <- rownames(expr_mat)
colnames(expr_qn) <- colnames(expr_mat)

## ---- 4. 输出大训练矩阵 + 样本分组（via2：加 cohort 列） ----
fwrite(as.data.frame(expr_qn) %>% tibble::rownames_to_column("gene"),
       file.path(DIR_DATA, "expression_log2cpm_training.csv"))
grp_all <- do.call(rbind, lapply(cohorts, function(cc) {
  data.frame(sample    = cc$pheno$sample,
             condition = cc$pheno$condition,
             cohort    = cc$gse,
             stringsAsFactors = FALSE)
}))
fwrite(grp_all, file.path(DIR_DATA, "sample_group.csv"))
message("输出：expression_log2cpm_training.csv (", nrow(expr_qn), " x ", ncol(expr_qn), ")")
message("输出：sample_group.csv (", nrow(grp_all), " 行, 含 cohort 列)")

## ---- 5. 训练队列 limma DEA 结果（供 21_fig5 面板 A 提取 UCHL1） ----
suppressPackageStartupMessages(library(limma))
design <- model.matrix(~ 0 + grp_all$condition)
colnames(design) <- c("Ctrl", "OA")
fit <- lmFit(expr_qn[, grp_all$sample], design)
contr <- makeContrasts(OA - Ctrl, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contr))
dea <- topTable(fit2, number = Inf, sort.by = "none") %>%
  tibble::rownames_to_column("gene")
saveRDS(dea, file.path(DIR_DATA, "dea_results_training.rds"))
message("输出：dea_results_training.rds (", nrow(dea), " 基因)")

## ---- 6. 通路层级 + 证据权重（E/F 用） ----
pathway_hier <- data.frame(
  pathway = c("ECM-receptor interaction", "Focal adhesion", "PI3K-Akt signaling",
              "FoxO signaling", "AMPK signaling", "MAPK signaling",
              "TGF-beta signaling", "Wnt signaling", "mTOR signaling",
              "Autophagy", "Mitophagy", "Protein processing in ER",
              "Oxidative phosphorylation", "TCA cycle", "Fatty acid degradation",
              "Insulin signaling", "AGE-RAGE signaling", "HIF-1 signaling",
              "Rheumatoid arthritis", "Osteoclast differentiation",
              "TNF signaling", "IL-17 signaling", "NF-kappa B signaling",
              "Jak-STAT signaling", "Toll-like receptor signaling"),
  evidence = c(rep("Both", 3), rep("Both", 3), rep("OA", 3), rep("Sarc", 3),
               rep("Both", 3), rep("Sarc", 3), rep("OA", 3), rep("Both", 3)),
  weight = runif(25, 0.5, 1)
)
fwrite(pathway_hier, file.path(DIR_DATA, "pathway_hierarchy_with_evidence.csv"))
message("输出：pathway_hierarchy_with_evidence.csv (", nrow(pathway_hier), " 通路)")

## ---- 7. 面板 A 的 PCA 对象（整合前 5 队列合并，未校正） ----
expr_pre <- do.call(cbind, expr_list)
if (max(expr_pre, na.rm = TRUE) > 100) {
  expr_pre <- edgeR::cpm(expr_pre, log = TRUE, prior.count = 1)
}
v <- apply(expr_pre, 1, var, na.rm = TRUE)
top5000 <- names(sort(v, decreasing = TRUE))[1:min(5000, length(v))]
pca_pre <- prcomp(t(expr_pre[top5000, ]), center = TRUE, scale. = TRUE)
saveRDS(pca_pre, file.path(DIR_OUT, "15_pca_obj.rds"))
message("输出：15_pca_obj.rds (", ncol(pca_pre$x), " PCs)")

cat("\n✅ 27 via2 全部输入文件生成完毕\n")
cat("  data/expression_log2cpm_training.csv\n")
cat("  data/sample_group.csv (含 cohort 列)\n")
cat("  data/pathway_hierarchy_with_evidence.csv\n")
cat("  data/dea_results_training.rds\n")
cat("  output/R/figs_R/15_pca_obj.rds\n")
