# ===== 补丁 20260913：从 validate 处续跑（会话内对象仍在，无需重跑前段）=====
## bootstrap 内部校验（B=200）——修正：validate.lrm 无 "C" 行，用 Dxy 换算
set.seed(CONFIG$SEED_MAIN + 1)
val_boot <- validate(fit_lrm, method = "boot", B = 200)
dxy_corr <- as.numeric(val_boot["Dxy", "index.corrected"])
c_opt <- round(0.5 + dxy_corr / 2, 3)
log_msg("bootstrap 校正 C-index = ", c_opt, "（B=200，由 Dxy=", round(dxy_corr,3), " 换算）")
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
  set.seed(CONFIG$SEED_MAIN + match(m, CONFIG$CLASSIFIERS))
  tune <- switch(m,
    glmnet = expand.grid(alpha = 0.5, lambda = 10^seq(-3, 0, length = 10)),
    rf     = expand.grid(mtry = c(2, 3, 5)),
    svmRadial = expand.grid(sigma = 0.01, C = c(0.1, 1, 10)),
    xgbTree = expand.grid(nrounds = 100, max_depth = 3, eta = 0.1,
                          gamma = 0, colsample_bytree = 0.8,
                          min_child_weight = 1, subsample = 0.8),
    knn    = expand.grid(k = c(5, 7, 9, 11)),
    nnet   = expand.grid(size = c(3, 5), decay = c(0.1, 0.5)),
    NULL)
  fit <- suppressWarnings(suppressMessages(
    caret::train(x = x_tr, y = y_tr, method = m, metric = "ROC",
                 trControl = ctrl_cv, tuneGrid = tune, trace = FALSE,
                 verbose = FALSE)))
  auc_cv <- max(fit$results$ROC, na.rm = TRUE)
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
if (any(rob_tab$cv_AUC < 0.8))
  log_msg("⚠ 有分类器 CV AUC<0.8：",
          paste(rob_tab$classifier[rob_tab$cv_AUC < 0.8], collapse = ","),
          "——稳健性叙事需如实收窄")
log_msg("分类器稳健性：CV AUC 区间 ", min(rob_tab$cv_AUC), "–", max(rob_tab$cv_AUC),
        "；外部 AUC 区间 ", min(rob_tab$GSE57218_AUC, na.rm = TRUE), "–",
        max(rob_tab$GSE57218_AUC, na.rm = TRUE))
log_msg("模块 I 全部完成。")
