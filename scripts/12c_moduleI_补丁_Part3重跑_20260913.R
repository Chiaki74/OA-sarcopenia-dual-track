# ===== 补丁 12c 20260913：仅重跑 Part 3（会话内对象仍在，validate 已跑通）=====
# 根因：上一版向 caret::train() 统一传了 trace=FALSE / verbose=FALSE，
#   caret 将 ... 原样转发给底层模型函数；glm 的 ... 进入 control=list(...)
#   → do.call(glm.control, ...)，而 glm.control 没有 verbose 参数 → 每个重抽样
#   折都报 "unused argument" 错误，caret 捕获后该折返回 NA →
#   "Something is wrong; all the ROC metric values are missing"。
#   （rf/knn 等同样不接受这两个参数，会连环失败。）
# 修复：不向 train() 传通用 trace/verbose；打印输出改由 capture.output 吞掉。
# 前置：本会话已存在 CONFIG / log_msg / x_tr / y_tr / va / feat_v / CV_FOLDS 等。

log_msg("—— Part 3（12c 补丁版）：分类器稳健性（", length(CONFIG$CLASSIFIERS), " 种分类器）——")
suppressPackageStartupMessages(library(caret))
for (pkg in c("glmnet", "randomForest", "kernlab", "xgboost", "class", "nnet", "naivebayes"))
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("缺少 ", pkg, " 包——install.packages(\"", pkg, "\") 后重跑本补丁")

ctrl_cv <- trainControl(method = "repeatedcv", number = CONFIG$CV_FOLDS,
                        repeats = CONFIG$CV_REPEATS,
                        classProbs = TRUE, summaryFunction = twoClassSummary,
                        allowParallel = FALSE)

## —— 先单发诊断 glm（不吞任何输出），确认管道通畅 ——
diag_fit <- caret::train(x = x_tr, y = y_tr, method = "glm", metric = "ROC",
                         trControl = ctrl_cv)
diag_auc <- max(diag_fit$results$ROC, na.rm = TRUE)
if (is.na(diag_auc))
  stop("诊断 glm 仍全 NA——请把上方完整控制台输出（含真实报错）发回")
log_msg("诊断 glm CV ROC = ", round(diag_auc, 3), " ✓ 管道通畅，开始 8 分类器全量")

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
  fit <- NULL
  invisible(capture.output(
    invisible(capture.output(
      fit <- caret::train(x = x_tr, y = y_tr, method = m, metric = "ROC",
                          trControl = ctrl_cv, tuneGrid = tune),
      type = "message")),
    type = "output"))
  auc_cv <- max(fit$results$ROC, na.rm = TRUE)
  ## 外部 AUC（GSE57218；仅 10 标志物全部可检时执行，口径与 05b/06 一致）
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
log_msg("模块 I 全部完成。产物：moduleI_single_gene_AUC.csv / moduleI_Fig6a2_nomogram.pdf",
        " / moduleI_nomogram_validate.txt / moduleI_classifier_robustness.csv")
log_msg("注意：CV AUC 存在特征选择泄漏（特征先于 CV 锁定），仅作参考；",
        "证据权重在外部 GSE57218 列——写稿时保留此声明。")
