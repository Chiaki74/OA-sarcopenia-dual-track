# ===== 补丁 12e 20260913：xgbTree 诊断 + 补跑（会话内对象仍在）=====
# 12d 中 xgbTree 报 "Stopping"（caret+xgboost 版本兼容问题的典型笼统报错）。
# 本补丁：1) 不吞输出单跑，暴露真实报错；2) 依次尝试修复参数组合；
# 3) 成功则把 xgbTree 行并入 moduleI_classifier_robustness.csv 并重写。

log_msg("—— xgbTree 诊断（12e）——")
suppressPackageStartupMessages(library(caret)); library(xgboost)
log_msg("caret ", as.character(packageVersion("caret")),
        " | xgboost ", as.character(packageVersion("xgboost")))

ctrl_cv <- trainControl(method = "repeatedcv", number = CONFIG$CV_FOLDS,
                        repeats = CONFIG$CV_REPEATS,
                        classProbs = TRUE, summaryFunction = twoClassSummary,
                        allowParallel = FALSE)
tune_xgb <- expand.grid(nrounds = 100, max_depth = 3, eta = 0.1, gamma = 0,
                        colsample_bytree = 0.8, min_child_weight = 1, subsample = 0.8)

## 第 1 步：裸跑暴露真实报错（不吞任何输出）
log_msg("第 1 步：裸跑 xgbTree（真实报错见下，若无报错即已成功）")
set.seed(CONFIG$SEED_MAIN + match("xgbTree", CONFIG$CLASSIFIERS))
fit_xgb <- tryCatch(
  caret::train(x = x_tr, y = y_tr, method = "xgbTree", metric = "ROC",
               trControl = ctrl_cv, tuneGrid = tune_xgb),
  error = function(e) { log_msg("裸跑报错：", conditionMessage(e)); NULL })

## 第 2 步：若失败，尝试 nthread=1 + verbose=0（xgboost 2.x 兼容修复）
if (is.null(fit_xgb)) {
  log_msg("第 2 步：尝试 nthread=1 + verbose=0 修复")
  set.seed(CONFIG$SEED_MAIN + match("xgbTree", CONFIG$CLASSIFIERS))
  fit_xgb <- tryCatch(
    caret::train(x = x_tr, y = y_tr, method = "xgbTree", metric = "ROC",
                 trControl = ctrl_cv, tuneGrid = tune_xgb,
                 nthread = 1, verbose = 0),
    error = function(e) { log_msg("修复尝试报错：", conditionMessage(e)); NULL })
}

if (is.null(fit_xgb)) {
  log_msg("xgbTree 在本环境不可用——以 7 分类器定稿（glm/glmnet/rf/svmRadial/knn/nnet/naive_bayes），",
          "写稿时按 7 种算法口径陈述，不回避、不粉饰。")
} else {
  auc_cv <- max(fit_xgb$results$ROC, na.rm = TRUE)
  p_ext <- predict(fit_xgb, newdata = as.data.frame(va$x[, feat_v, drop = FALSE]),
                   type = "prob")[, "case"]
  auc_ext <- as.numeric(pROC::auc(pROC::roc(va$y, p_ext, quiet = TRUE)))
  log_msg("✓ xgbTree 补跑成功： CV AUC=", round(auc_cv, 3),
          " | 外部 AUC=", round(auc_ext, 3))
  rob_tab <- read.csv(file.path(CONFIG$OUT_DIR, "moduleI_classifier_robustness.csv"))
  rob_tab <- rob_tab[rob_tab$classifier != "xgbTree", ]
  rob_tab <- rbind(rob_tab, data.frame(
    classifier = "xgbTree", cv_AUC = round(auc_cv, 3),
    GSE57218_AUC = round(auc_ext, 3),
    seed = CONFIG$SEED_MAIN + match("xgbTree", CONFIG$CLASSIFIERS)))
  rob_tab <- rob_tab[match(CONFIG$CLASSIFIERS, rob_tab$classifier), ]
  fwrite(rob_tab, file.path(CONFIG$OUT_DIR, "moduleI_classifier_robustness.csv"))
  log_msg("CSV 已重写（8 分类器）。CV AUC 区间 ", min(rob_tab$cv_AUC), "–",
          max(rob_tab$cv_AUC), "；外部 AUC 区间 ", min(rob_tab$GSE57218_AUC),
          "–", max(rob_tab$GSE57218_AUC))
}
log_msg("12e 完成。")
