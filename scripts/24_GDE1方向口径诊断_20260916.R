## =============================================================================
## 24_GDE1方向口径诊断_20260916.R
## 目的：模块 E 方向表称 GDE1 在 GSE57218 为「OA 上调（均值口径，与发现层一致）」，
##       而 2026-09-16 复核的方向固定 AUC=0.359（秩次口径，多数 OA 样本低于对照）。
##       同一队列、同一样本集（33 OA / 7 对照）两种统计口径结论相反。
##       本脚本打印全部证据，确认是否为「均值–秩次背离」（少数高表达 OA 样本拉高均值）。
## 运行：source("24_GDE1方向口径诊断_20260916.R")
## =============================================================================

## ---- 根目录定位（与 23 号脚本同口径：processed/results/data 任一存在即可）----
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
if (!dir.exists(file.path(.d, "processed", "GSE57218")))
  stop("在 ", .d, " 下找不到 processed/GSE57218——请把本脚本放到项目根目录再 source")
setwd(.d); message("[setwd] 项目根目录：", .d)

suppressPackageStartupMessages(library(pROC))

## ---- 读取 GSE57218 表达与分组（与 22/23 号脚本同口径）----
expr <- readRDS("processed/GSE57218/GSE57218_expr_gene.rds")
if (!"GDE1" %in% rownames(expr) && "GDE1" %in% colnames(expr)) expr <- t(expr)
pheno <- read.csv("processed/GSE57218/GSE57218_pheno.csv", stringsAsFactors = FALSE)

cond_col <- intersect(c("condition", "condition_raw", "group", "label"), colnames(pheno))[1]
samp_col <- intersect(c("sample", "gsm", "geo_accession", "id"), colnames(pheno))[1]
lab_raw  <- tolower(trimws(as.character(pheno[[cond_col]])))
label    <- ifelse(grepl("preserv", lab_raw), 2L,
            ifelse(grepl("^(oa|case|1|disease)|osteoarthritis", lab_raw), 1L, 0L))
keep <- label != 2L                       # OA_preserved 剔出主验证（33/7 口径）
label <- label[keep]; samp <- as.character(pheno[[samp_col]][keep])

g <- as.numeric(expr["GDE1", samp])
x_case <- g[label == 1]; x_ctrl <- g[label == 0]

cat("\n================ GDE1 @ GSE57218（33 OA / 7 对照）================\n")
cat("样本数：OA =", length(x_case), "，对照 =", length(x_ctrl), "\n\n")
cat("对照值：", paste(round(sort(x_ctrl), 3), collapse = ", "), "\n")
cat("OA 值（升序）：", paste(round(sort(x_case), 3), collapse = ", "), "\n\n")

cat("均值：OA =", round(mean(x_case), 3), "，对照 =", round(mean(x_ctrl), 3),
    " → 均值口径方向：", ifelse(mean(x_case) > mean(x_ctrl), "up（模块E口径）", "down"), "\n")
cat("中位数：OA =", round(median(x_case), 3), "，对照 =", round(median(x_ctrl), 3),
    " → 中位数口径方向：", ifelse(median(x_case) > median(x_ctrl), "up", "down"), "\n\n")

wt <- suppressWarnings(wilcox.test(x_case, x_ctrl))
cat("Wilcoxon 秩和检验 P =", format.pval(wt$p.value, digits = 3), "\n")

r_up   <- as.numeric(pROC::auc(pROC::roc(label, g, levels = c(0, 1), direction = "<", quiet = TRUE)))
r_down <- as.numeric(pROC::auc(pROC::roc(label, g, levels = c(0, 1), direction = ">", quiet = TRUE)))
cat("AUC（发现层 up 固定, direction='<'）：", round(r_up, 3),
    " ← 22 号复核值应为 0.359\n")
cat("AUC（反向, direction='>'）：", round(r_down, 3), "（应 = 1 - 0.359 = 0.641，即旧表 auto 值）\n")

n_pairs_below <- sum(outer(x_case, x_ctrl, `<`)); n_pairs <- length(x_case) * length(x_ctrl)
cat("\n逐对比较：OA < 对照 的对子占 ", n_pairs_below, "/", n_pairs,
    "（", round(100 * n_pairs_below / n_pairs, 1), "%）\n")
n_above_ctrl_max <- sum(x_case > max(x_ctrl))
cat("高于对照最大值的 OA 样本数：", n_above_ctrl_max, "/", length(x_case),
    "（若为少数 → 均值被少数高值样本拉高，坐实均值–秩次背离）\n")

cat("\n判定：")
if (mean(x_case) > mean(x_ctrl) && r_up < 0.5) {
  cat("均值–秩次背离成立。模块 E「up」与复核「反向 0.359」均为真，",
      "差异源于统计口径（均值 vs 秩次），正文按保守口径（秩次/AUC）报告并加调和句。\n")
} else if (mean(x_case) < mean(x_ctrl) && r_up < 0.5) {
  cat("两种口径一致均为 down → 模块 E 表该行为旧跑残留，需重跑模块 E 或人工订正。\n")
} else {
  cat("结果与预期不符，请把本输出完整发回。\n")
}
