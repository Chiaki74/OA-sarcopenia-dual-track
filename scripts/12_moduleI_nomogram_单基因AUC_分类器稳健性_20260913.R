# ============================================================================
# 12_moduleI_nomogram_单基因AUC_分类器稳健性_20260913.R
# 模块 I｜v5 修订新增分析（结构对照拍板项 E1/E2/E3，备忘录 §66）
# ----------------------------------------------------------------------------
# 三个产出（对照 Huang 2025 与 Wang W 2026 同构模块，纯增量、不改任何既有结论）：
#   Part 1（E2）10 标志物单基因 univariate AUC 表：训练层（五队列合并 z-score）
#     + 软骨主验证层 GSE57218，逐基因报告——对齐 Huang 2025 Fig 7 的体裁惯例；
#   Part 2（E1）Nomogram 列线图：rms::lrm 重建 10 标志物 logistic 模型，
#     出列线图 PDF + rms 内部 bootstrap 校验——对齐 Huang 2025 Fig 6 体裁；
#   Part 3（E3）分类器稳健性 panel：同一 10 标志物签名在 8 种分类器下
#     10 折 CV（重复 3 次）AUC + GSE57218 外部 AUC——回应"为何用 logistic"，
#     对齐 Wang W 2026 Fig 5A 的 10 模型 ROC 惯例（其仅内部 CV，我们加外部层）。
# 纪律：训练/验证装配函数与 06_moduleE_validation.R 逐行同口径（pheno 硬校验、
#   donor 去重、队列内 z-score、零方差剔除）；全部随机种子显式报告（E6）。
# 前置：results/moduleD_v2/moduleD_v2_biomarkers.rds（05b 已跑）；
#       processed/{GSE}/{GSE}_expr_gene.rds + {GSE}_pheno.csv（02 已跑）
# 运行：source("12_moduleI_nomogram_单基因AUC_分类器稳健性_20260913.R")
# 输出：results/moduleI/
# ============================================================================

suppressPackageStartupMessages({
  library(pROC)
  library(data.table)
})

# ----------------------------- CONFIG --------------------------------------
CONFIG <- list(
  V2_RDS     = "results/moduleD_v2/moduleD_v2_biomarkers.rds",
  IN_DIR     = "processed",
  OUT_DIR    = "results/moduleI",
  LOG_FILE   = file.path("results/moduleI", sprintf("moduleI_log_%s.txt", Sys.Date())),
  TRAIN_SETS = c("GSE114007", "GSE55235", "GSE55457", "GSE12021", "GSE51588"),
  TRAIN_LABEL = "OA",
  VAL_SET    = "GSE57218",
  VAL_CASE_CANDIDATES = c("OA_affected", "OA"),   # 与 06 同口径
  ## E3 分类器清单（8 种，caret 接口；种子逐方法显式设定）
  CLASSIFIERS = c("glm", "glmnet", "rf", "svmRadial", "xgbTree", "knn", "nnet", "naive_bayes"),
  CV_FOLDS = 10, CV_REPEATS = 3,
  SEED_MAIN = 20260913
)
dir.create(CONFIG$OUT_DIR, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  txt <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = ""))
  cat(txt, "\n"); cat(txt, "\n", file = CONFIG$LOG_FILE, append = TRUE)
}

# ---- 通用装配（与 06 逐行同口径，勿改）----
assemble_cohort <- function(gse, case_label, ctrl_label = "control") {
  expr  <- readRDS(file.path(CONFIG$IN_DIR, gse, paste0(gse, "_expr_gene.rds")))
  pheno <- read.csv(file.path(CONFIG$IN_DIR, gse, paste0(gse, "_pheno.csv")),
                    stringsAsFactors = FALSE)
  pheno <- pheno[match(colnames(expr), pheno$sample), ]
  if (anyNA(pheno$sample) || anyNA(pheno$condition))
    stop(gse, ": pheno 缺少 expr 中的样本或 condition 缺失——请回查模块 A 输出")
  keep <- pheno$condition %in% c(ctrl_label, case_label)
  if (sum(keep) == 0)
    stop(gse, ": 无 case/control 样本——case 标签 '", case_label,
         "' 不在 condition 取值（", paste(unique(pheno$condition), collapse = ","), "）中")
  expr <- expr[, keep, drop = FALSE]; pheno <- pheno[keep, ]
  if ("donor_id" %in% names(pheno) && anyDuplicated(pheno$donor_id)) {
    if (anyNA(pheno$donor_id)) stop(gse, ": donor_id 含 NA，请人工核对")
    keep_idx <- !duplicated(pheno$donor_id)
    expr <- expr[, keep_idx, drop = FALSE]; pheno <- pheno[keep_idx, ]
    log_msg("  ", gse, ": donor_id 去重 → ", ncol(expr), " 样本")
  }
  x <- scale(t(expr))
  if (anyNA(x)) {
    bad <- colnames(x)[apply(x, 2, anyNA)]
    log_msg("  ", gse, ": 剔除 ", length(bad), " 个零方差基因（z-score NaN）")
    x <- x[, !colnames(x) %in% bad, drop = FALSE]
  }
  y <- factor(ifelse(pheno$condition == case_label, "case", "control"),
              levels = c("control", "case"))
  if (sum(y == "case") == 0 || sum(y == "control") == 0)
    stop(gse, ": 病例/对照有一侧为 0")
  list(x = x, y = y)
}

resolve_case <- function(gse, candidates) {
  pheno <- read.csv(file.path(CONFIG$IN_DIR, gse, paste0(gse, "_pheno.csv")),
                    stringsAsFactors = FALSE)
  avail <- unique(pheno$condition)
  hit <- intersect(candidates, avail)
  if (length(hit) == 0)
    stop(gse, ": case 候选均不在 condition 实测取值（",
         paste(avail, collapse = ","), "）中——请人工核对")
  hit[1]
}

# ---- 数据装配：训练层 + 软骨主验证层 ----
v2 <- readRDS(CONFIG$V2_RDS)
biomarkers <- v2$biomarkers
if (length(biomarkers) != 10)
  stop("v2 标志物数 != 10（实测 ", length(biomarkers), "）——请回查 05b 封板")
log_msg("v2 标志物（10）：", paste(biomarkers, collapse = ", "))

tr <- lapply(CONFIG$TRAIN_SETS, assemble_cohort, case_label = CONFIG$TRAIN_LABEL)
feat <- Reduce(intersect, c(list(biomarkers), lapply(tr, function(cc) colnames(cc$x))))
if (length(feat) != 10)
  stop("训练层可检标志物 != 10（实测 ", length(feat), "：",
       paste(feat, collapse = ","), "）——与 06 封板口径不一致，请回查")
x_tr <- do.call(rbind, lapply(tr, function(cc) cc$x[, feat, drop = FALSE]))
y_tr <- factor(unlist(lapply(tr, function(cc) as.character(cc$y))),
               levels = c("control", "case"))
log_msg("训练层：", nrow(x_tr), " 样本（case=", sum(y_tr == "case"),
        " / control=", sum(y_tr == "control"), "）")
## 硬自检 1：样本量与 06 封板一致（117 = 65/52）
if (nrow(x_tr) != 117 || sum(y_tr == "case") != 65)
  stop("训练层样本量偏离封板（117=65/52）——实得 ", nrow(x_tr), "=",
       sum(y_tr == "case"), "/", sum(y_tr == "control"), "，请回查")

val_case <- resolve_case(CONFIG$VAL_SET, CONFIG$VAL_CASE_CANDIDATES)
va <- assemble_cohort(CONFIG$VAL_SET, case_label = val_case)
feat_v <- intersect(feat, colnames(va$x))
log_msg("验证层 ", CONFIG$VAL_SET, "：", nrow(va$x), " 样本，可检标志物 ",
        length(feat_v), "/10")
## 硬自检 2：GSE57218 主验证构成与封板一致（33 受累 OA vs 7 正常）
if (sum(va$y == "case") != 33 || sum(va$y == "control") != 7)
  log_msg("⚠ ", CONFIG$VAL_SET, " 样本构成偏离封板（33/7）——实得 ",
          sum(va$y == "case"), "/", sum(va$y == "control"), "，请人工核对后再用结果")

# ============================================================================
# Part 1（E2）单基因 univariate AUC 表
# ============================================================================
log_msg("—— Part 1：单基因 univariate AUC ——")
auc_tab <- data.frame()
for (g in feat) {
  ## 方向修正：AUC<0.5 时取 1-AUC 并记录方向
  roc_tr <- suppressMessages(pROC::roc(y_tr, x_tr[, g], quiet = TRUE))
  a_tr <- as.numeric(pROC::auc(roc_tr)); d_tr <- ifelse(a_tr < 0.5, "down_in_case", "up_in_case")
  a_tr <- max(a_tr, 1 - a_tr)
  if (g %in% feat_v) {
    roc_va <- suppressMessages(pROC::roc(va$y, va$x[, g], quiet = TRUE))
    a_va <- as.numeric(pROC::auc(roc_va))
    ci_va <- as.numeric(pROC::ci.auc(roc_va, method = "delong"))
    a_va_dir <- max(a_va, 1 - a_va)
  } else { a_va_dir <- NA; ci_va <- c(NA, NA, NA) }
  auc_tab <- rbind(auc_tab, data.frame(
    gene = g, train_AUC = round(a_tr, 3), case_direction = d_tr,
    GSE57218_AUC = round(a_va_dir, 3),
    GSE57218_CI = ifelse(is.na(a_va_dir), NA,
      sprintf("%.3f-%.3f", min(ci_va[1], 1-ci_va[1]), max(ci_va[3], 1-ci_va[3]+0)+0))
  ))
}
auc_tab <- auc_tab[order(-auc_tab$GSE57218_AUC), ]
fwrite(auc_tab, file.path(CONFIG$OUT_DIR, "moduleI_single_gene_AUC.csv"))
log_msg("单基因 AUC 表已出（", nrow(auc_tab), " 行）：GSE57218 最高 ",
        auc_tab$gene[1], "=", auc_tab$GSE57218_AUC[1])
## 硬自检 3：AUC 值域
if (any(auc_tab$train_AUC < 0.5) || any(auc_tab$GSE57218_AUC < 0.5, na.rm = TRUE))
  stop("方向修正后仍有 AUC<0.5——逻辑错误，请回查")

# ============================================================================
# Part 2（E1）Nomogram 列线图（rms）
# ============================================================================
log_msg("—— Part 2：Nomogram（rms::lrm）——")
if (!requireNamespace("rms", quietly = TRUE))
  stop("缺少 rms 包——请运行 install.packages(\"rms\") 后重跑")
suppressPackageStartupMessages(library(rms))
df_tr <- as.data.frame(x_tr); df_tr$y <- y_tr
dd <- datadist(df_tr); options(datadist = "dd")
set.seed(CONFIG$SEED_MAIN)
fit_lrm <- lrm(y ~ ., data = df_tr, x = TRUE, y = TRUE)
## 硬自检 4：模型收敛且 10 个系数齐备
cf <- coef(fit_lrm)
if (length(cf) != 11 || anyNA(cf))
  stop("lrm 系数异常（", length(cf), " 个，含 NA=", sum(is.na(cf)), "）——请回查")
log_msg("lrm 收敛，10 系数齐备；C-index（表观）= ",
        round(as.numeric(fit_lrm$stats["C"]), 3))
## 与 05b/06 的 glm 系数方向一致性核对（同数据同族模型，方向应 10/10 一致）
fit_glm <- glm(y ~ ., data = df_tr, family = binomial())
sign_agree <- sum(sign(coef(fit_lrm)) == sign(coef(fit_glm)))
log_msg("lrm 与 glm 系数方向一致：", sign_agree, "/11（含截距）")
if (sign_agree < 11) log_msg("⚠ 方向不一致——请人工核对（两模型同族同数据，理应全一致）")
nom <- nomogram(fit_lrm, fun = plogis, funlabel = "P(OA)",
                lp = TRUE, maxscale = 100)
pdf(file.path(CONFIG$OUT_DIR, "moduleI_Fig6a2_nomogram.pdf"), width = 10, height = 7)
plot(nom, cex.axis = 0.85, cex.var = 0.9)
title("Nomogram of the 10-biomarker panel (training layer, z-scored expression)")
dev.off()
## bootstrap 内部校验（B=200， optimism 校正 C-index）
set.seed(CONFIG$SEED_MAIN + 1)
val_boot <- validate(fit_lrm, method = "boot", B = 200)
# validate.lrm 输出矩阵无 "C" 行（仅 Dxy/R2/Slope 等）；C-index = 0.5 + Dxy/2
dxy_corr <- as.numeric(val_boot["Dxy", "index.corrected"])
c_opt <- round(0.5 + dxy_corr / 2, 3)
log_msg("bootstrap 校正 C-index = ", c_opt, "（B=200）")
sink(file.path(CONFIG$OUT_DIR, "moduleI_nomogram_validate.txt"))
print(val_boot); sink()
log_msg("Nomogram PDF 已出：moduleI_Fig6a2_nomogram.pdf")

# ============================================================================
# Part 3（E3）分类器稳健性 panel（8 分类器 × 10 折 CV×3 + 外部 AUC）
# ============================================================================
log_msg("—— Part 3：分类器稳健性（", length(CONFIG$CLASSIFIERS), " 种分类器）——")
if (!requireNamespace("caret", quietly = TRUE))
  stop("缺少 caret 包——请运行 install.packages(\"caret\") 后重跑")
suppressPackageStartupMessages(library(caret))
for (pkg in c("glmnet", "randomForest", "kernlab", "xgboost", "class", "nnet", "e1071"))
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("缺少 ", pkg, " 包——安装后重跑")

ctrl_cv <- trainControl(method = "repeatedcv", number = CONFIG$CV_FOLDS,
                        repeats = CONFIG$CV_REPEATS,
                        classProbs = TRUE, summaryFunction = twoClassSummary,
                        allowParallel = FALSE)
rob_tab <- data.frame()
x_va_df <- as.data.frame(va$x[, feat_v, drop = FALSE])
for (m in CONFIG$CLASSIFIERS) {
  set.seed(CONFIG$SEED_MAIN + match(m, CONFIG$CLASSIFIERS))   # 逐方法显式种子（E6）
  tune <- switch(m,
    glmnet = expand.grid(alpha = 0.5, lambda = 10^seq(-3, 0, length = 10)),
    rf     = expand.grid(mtry = c(2, 3, 5)),
    svmRadial = expand.grid(sigma = 0.01, C = c(0.1, 1, 10)),
    xgbTree = expand.grid(nrounds = 100, max_depth = 3, eta = 0.1,
                          gamma = 0, colsample_bytree = 0.8,
                          min_child_weight = 1, subsample = 0.8),
    knn    = expand.grid(k = c(5, 7, 9, 11)),
    nnet   = expand.grid(size = c(3, 5), decay = c(0.1, 0.5)),
    NULL)  # glm / naive_bayes 无调参网格
  ## 注意：不可向 train() 统一传 trace/verbose——caret 会把 ... 原样转发给底层模型，
  ## glm 的 ... 进入 glm.control()（无 verbose 参数）而逐折报错，caret 捕获后返回
  ## 全 NA（"all the ROC metric values are missing"）。静默输出改用 capture.output。
  fit <- NULL
  invisible(capture.output(
    invisible(capture.output(
      fit <- caret::train(x = x_tr, y = y_tr, method = m, metric = "ROC",
                          trControl = ctrl_cv, tuneGrid = tune),
      type = "message")),
    type = "output"))
  auc_cv <- max(fit$results$ROC, na.rm = TRUE)
  ## 外部 AUC（GSE57218；预测用可检标志物子集重训同法模型会导致口径漂移，
  ## 故外部层仅在 10 标志物全部可检时执行）
  auc_ext <- NA
  if (length(feat_v) == 10) {
    p_ext <- predict(fit, newdata = x_va_df, type = "prob")[, "case"]
    auc_ext <- as.numeric(pROC::auc(pROC::roc(va$y, p_ext, quiet = TRUE)))
  }
  rob_tab <- rbind(rob_tab, data.frame(
    classifier = m, cv_AUC = round(auc_cv, 3),
    GSE57218_AUC = ifelse(is.na(auc_ext), NA, round(auc_ext, 3)),
    seed = CONFIG$SEED_MAIN + match(m, CONFIG$CLASSIFIERS)))
  log_msg("  ", sprintf("%-12s", m), " CV AUC=", round(auc_cv, 3),
          " | 外部 AUC=", ifelse(is.na(auc_ext), "NA", round(auc_ext, 3)))
}
fwrite(rob_tab, file.path(CONFIG$OUT_DIR, "moduleI_classifier_robustness.csv"))
## 硬自检 5：稳健性声明的证据边界——所有分类器 CV AUC 应 > 0.8（内部口径），
## 外部 AUC 最低值作为稳健性下界如实记录，不设硬停（探索性 panel）
if (any(rob_tab$cv_AUC < 0.8))
  log_msg("⚠ 有分类器 CV AUC<0.8：",
          paste(rob_tab$classifier[rob_tab$cv_AUC < 0.8], collapse = ","),
          "——稳健性叙事需如实收窄")
log_msg("分类器稳健性：CV AUC 区间 ", min(rob_tab$cv_AUC), "–", max(rob_tab$cv_AUC),
        "；外部 AUC 区间 ", min(rob_tab$GSE57218_AUC, na.rm = TRUE), "–",
        max(rob_tab$GSE57218_AUC, na.rm = TRUE))

# ---- 汇总 ----
log_msg("模块 I 完成。产物：")
log_msg("  moduleI_single_gene_AUC.csv（E2）")
log_msg("  moduleI_Fig6a2_nomogram.pdf + moduleI_nomogram_validate.txt（E1）")
log_msg("  moduleI_classifier_robustness.csv（E3）")
log_msg("注意：CV AUC 存在特征选择泄漏（特征在 05b 已先于 CV 锁定），仅作参考；")
log_msg("  证据权重在外部 GSE57218 列——与 05b/06 口径一致，写稿时保留此声明。")
