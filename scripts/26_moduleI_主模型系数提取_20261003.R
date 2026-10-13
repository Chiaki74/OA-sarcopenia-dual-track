# ============================================================================
# 26_moduleI_主模型系数提取_20261003.R
# TRIPOD+AI 条目 22（模型完整呈现）：提取主模型（10 标志物 logistic 回归）
# 的完整方程——10 个系数＋截距、OR 与 95%CI、Wald P，排版为表 S13 附属。
# ----------------------------------------------------------------------------
# 原理：模块 I 当前训练矩阵（results/moduleI/train_expr_z.csv +
#       train_group.csv，2026-09-16 版，即 GSE51588 供体去重后 117 样本、
#       逐队列 z-score 合并）即最终模型训练数据；glm 极大似然拟合确定性
#       复现主模型系数。表观 AUC 打印供 sanity check（应略高于内部 CV
#       0.959，与外部 GSE57218 0.818 属不同口径，不可直接比较）。
# 输入：results/moduleI/train_expr_z.csv、results/moduleI/train_group.csv
# 输出：results/moduleI/表S13附_主模型完整系数表_20261003.csv
# 运行（项目根目录）：
#   source("26_moduleI_主模型系数提取_20261003.R")
# ============================================================================

BIOM <- c("NELL1","STEAP1","GADD45A","FBLN5","HEMK1","GDE1","RNF14","BNIP3","EGFR","CBR3")

locate_root <- function() {
  d <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  for (i in 1:8) {
    if (dir.exists(file.path(d, "results", "moduleI"))) return(d)
    parent <- dirname(d); if (parent == d) break; d <- parent
  }
  stop("未找到 results/moduleI/——请在项目根目录运行")
}
ROOT <- locate_root()
MIDIR <- file.path(ROOT, "results", "moduleI")

## ---- 读取训练矩阵（自动判别基因在行还是在列） -------------------------------
ex <- read.csv(file.path(MIDIR, "train_expr_z.csv"), check.names = FALSE)
first_col <- ex[[1]]
if (sum(BIOM %in% first_col) >= 8) {          # 基因在行：转置为 样本×基因
  rownames(ex) <- make.unique(as.character(first_col))
  ex[[1]] <- NULL
  X <- as.data.frame(t(as.matrix(ex)))
} else if (sum(BIOM %in% colnames(ex)) >= 8) { # 基因在列
  if (is.character(first_col) && !anyNA(match(first_col, first_col))) {
    rownames(ex) <- make.unique(as.character(first_col)); ex[[1]] <- NULL
  }
  X <- ex
} else stop("无法在 train_expr_z.csv 中定位 10 个标志物——请把该文件传给对话方检查结构")
if (!all(BIOM %in% colnames(X))) stop("标志物缺失：", paste(setdiff(BIOM, colnames(X)), collapse = ", "))
X <- X[, BIOM, drop = FALSE]
message("训练矩阵：", nrow(X), " 样本 × ", ncol(X), " 标志物")

## ---- 读取分组（自动识别样本列与分组列） -------------------------------------
gr <- read.csv(file.path(MIDIR, "train_group.csv"), check.names = FALSE)
sample_col <- names(gr)[vapply(gr, function(col) mean(as.character(col) %in% rownames(X)) > 0.8, logical(1))][1]
if (is.na(sample_col)) stop("train_group.csv 中找不到与训练矩阵匹配的样本列")
rest <- setdiff(names(gr), sample_col)
group_col <- rest[vapply(rest, function(cn) length(unique(gr[[cn]])) == 2, logical(1))][1]
if (is.na(group_col)) stop("未识别出二分类分组列，现有列：", paste(rest, collapse = ", "))
y_raw <- as.character(gr[[group_col]][match(rownames(X), as.character(gr[[sample_col]]))])
if (any(is.na(y_raw))) stop("分组与样本对齐失败")
lv <- unique(y_raw)
case_lv <- lv[grepl("OA|case|Case|病例|1", lv)][1]
if (is.na(case_lv)) stop("无法判定病例水平，分组取值：", paste(lv, collapse = ", "))
y <- as.integer(y_raw == case_lv)
message("分组列：", group_col, "；病例水平：", case_lv, "（n=", sum(y), " 病例 / ", sum(1 - y), " 对照）")

## ---- 拟合主模型并提取系数 ---------------------------------------------------
df <- data.frame(case = y, X, check.names = FALSE)
fit <- glm(case ~ ., family = binomial, data = df)
co <- summary(fit)$coefficients
ci <- suppressWarnings(confint.default(fit))   # Wald 95%CI，确定性
tab <- data.frame(
  变量   = rownames(co),
  系数β  = round(co[, "Estimate"], 4),
  标准误 = round(co[, "Std. Error"], 4),
  OR    = round(exp(co[, "Estimate"]), 4),
  `OR 95%CI` = sprintf("%.4f–%.4f", exp(ci[, 1]), exp(ci[, 2])),
  Wald_z = round(co[, "z value"], 3),
  P值    = signif(co[, "Pr(>|z|)"], 3),
  row.names = NULL, check.names = FALSE)
tab$方向 <- ifelse(tab$系数β > 0, "高表达→OA 风险↑", "高表达→OA 风险↓")
tab$方向[tab$变量 == "(Intercept)"] <- "—"

## ---- 表观 AUC（sanity check） ----------------------------------------------
p_hat <- as.numeric(predict(fit, type = "response"))
# Mann-Whitney 形式 AUC（含 ties 校正）
auc <- (sum(outer(p_hat[y == 1], p_hat[y == 0], `>`)) +
        0.5 * sum(outer(p_hat[y == 1], p_hat[y == 0], `==`))) / (sum(y) * sum(1 - y))
message(sprintf("表观 AUC = %.3f（应≥内部 CV 0.959 附近；与外部 GSE57218 0.818 不同口径）", auc))

## ---- 输出 -------------------------------------------------------------------
out <- file.path(MIDIR, "表S13附_主模型完整系数表_20261003.csv")
write.csv(tab, out, row.names = FALSE, fileEncoding = "UTF-8")
cat("\n主模型完整方程：logit(P) =", round(co["(Intercept)", "Estimate"], 4))
for (g in BIOM) cat(" + (", round(co[g, "Estimate"], 4), ")×", g, sep = "")
cat("\n\n")
print(tab)
message("完成。输出：", out)
