# ============================================================================
# 05b_moduleD_v2_multicohort.R
# 模块 D v2｜多队列训练 + 平台预筛 + 分层验证（2026-08-23，备忘录 §八-27）
# ----------------------------------------------------------------------------
# 回炉背景（v1 结果触发计划书 §4.5「验证导向」预案）：
#   v1 单队列（GSE114007, n=38）训练，外部血液 GSE48556 AUC=0.642（0.524–0.755）；
#   方向核对显示血液可测性集中于「下调轴」（IER2/BNIP3/GADD45A 一致，
#   FAM162A 平台缺失、GDE1 翻转）——5 特征中 2 个在血液失效。
# v2 设计（方案 A，研究者授权）：
#   1) 多队列发现层：5 个 OA 训练集逐队列 z-score 后合并（≈120 样本），
#      提升签名稳健性；GSE51588 双部位按 donor_id 去重（每供体留 1 样本，
#      防伪重复）；
#   2) 平台预筛：特征空间 = 99 双轨基因 ∩ 各训练队列可检 ∩ GPL6947（GSE48556
#      血液平台）可检——先验平台信息，不构成验证集泄漏；
#   3) 四法交集重选标志物（LASSO + SVM-RFE + RF + XGBoost），口径同 v1；
#   4) GSE48556 保持纯验证（不碰标签）；分层报告：血液=探索性验证，
#      软骨验证集（GSE57218/168505/206848）留模块 E 作主验证；
#   5) 评估三件套（ROC/校准/DCA）+ 方向一致性核对，口径同 v1；
#      内部 CV AUC 存在特征选择泄漏，仅作参考，证据权重在外部验证。
# ----------------------------------------------------------------------------
# 前置：results/moduleC/moduleC_gene_tracks.csv；5 个训练集与 GSE48556
#       已经 02_moduleA_preprocess.R 预处理
# 运行：source("05b_moduleD_v2_multicohort.R")；结果输出 results/moduleD_v2/
# ============================================================================

suppressPackageStartupMessages({
  library(glmnet)
  library(e1071)
  library(randomForest)
  library(xgboost)
  library(pROC)
  library(caret)
  library(data.table)
})
## 审查S1：caret 的 svmLinear 后端是 kernlab::ksvm（Suggests 依赖，可能缺失），前置检查
if (!requireNamespace("kernlab", quietly = TRUE))
  stop("缺少 kernlab 包（SVM-RFE 后端）——请运行 install.packages(\"kernlab\") 后重跑")
set.seed(20260823)

# ----------------------------- CONFIG --------------------------------------
CONFIG <- list(
  TRACK_FILE = "results/moduleC/moduleC_gene_tracks.csv",
  IN_DIR  = "processed",
  OUT_DIR = "results/moduleD_v2",
  LOG_FILE = file.path("results/moduleD_v2", sprintf("moduleD_v2_log_%s.txt", Sys.Date())),
  TRAIN_SETS = c("GSE114007", "GSE55235", "GSE55457", "GSE12021", "GSE51588"),
  TRAIN_LABEL = "OA",
  VALIDATION_GSE   = "GSE48556",   # 血液探索性验证（GARP，106 OA vs 33 对照）
  VALIDATION_LABEL = "OA",
  CV_REPEATS = 10, CV_FOLDS = 10,
  N_BOOT = 2000,
  MAX_BIOMARKERS = 10
)
dir.create(CONFIG$OUT_DIR, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  txt <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = ""))
  cat(txt, "\n"); cat(txt, "\n", file = CONFIG$LOG_FILE, append = TRUE)
}

# ----------------------------- 双轨基因表读取与校验 --------------------------
if (!file.exists(CONFIG$TRACK_FILE))
  stop("缺少 ", CONFIG$TRACK_FILE, "——请先运行 04_moduleC_enrichment.R")
track_tab <- read.csv(CONFIG$TRACK_FILE, stringsAsFactors = FALSE)
valid_tracks <- c("concordant", "mirror_OAup_muscleDown", "mirror_OAdown_muscleUp")
if (!all(c("gene", "track", "direction_OA", "direction_muscle", "in_WGCNA_key") %in% names(track_tab)))
  stop("moduleC_gene_tracks.csv 列名口径漂移，请核对模块 C 输出")
if (any(!track_tab$track %in% valid_tracks) ||
    any(!track_tab$direction_OA %in% c("up", "down")) ||
    any(!track_tab$direction_muscle %in% c("up", "down")) ||
    any(is.na(track_tab$gene)) || anyDuplicated(track_tab$gene))
  stop("track/direction/gene 列口径异常，请核对模块 C 输出")
core_genes <- track_tab$gene
log_msg("双轨核心基因：", length(core_genes), "（预期 99）")

# ----------------------------- 单队列装配（z-score + 供体去重）---------------
assemble_cohort <- function(gse, label) {
  expr  <- readRDS(file.path(CONFIG$IN_DIR, gse, paste0(gse, "_expr_gene.rds")))
  pheno <- read.csv(file.path(CONFIG$IN_DIR, gse, paste0(gse, "_pheno.csv")),
                    stringsAsFactors = FALSE)
  pheno <- pheno[match(colnames(expr), pheno$sample), ]
  ## 审查S2：pheno 缺样本/脏 condition 会在 keep 过滤中被静默剔除，必须硬校验
  if (anyNA(pheno$sample) || anyNA(pheno$condition))
    stop(gse, ": pheno 缺少 expr 中的样本或 condition 缺失——请回查模块 A 输出")
  if (!all(pheno$condition %in% c("control", label, "other")))
    stop(gse, ": condition 存在越界取值：",
         paste(setdiff(pheno$condition, c("control", label, "other")), collapse = ","))
  keep <- pheno$condition %in% c("control", label)
  expr <- expr[, keep, drop = FALSE]; pheno <- pheno[keep, ]
  ## GSE51588 等双部位队列：每供体保留 1 个样本，防伪重复（v2 新增）
  if ("donor_id" %in% names(pheno) && anyDuplicated(pheno$donor_id)) {
    ## 审查A2：donor_id 含 NA 时 duplicated 会误删；保留样本留痕
    if (anyNA(pheno$donor_id)) stop(gse, ": donor_id 含 NA，请人工核对")
    before <- ncol(expr)
    keep_idx <- !duplicated(pheno$donor_id)
    kept_samples <- pheno$sample[keep_idx]
    expr <- expr[, keep_idx, drop = FALSE]; pheno <- pheno[keep_idx, ]
    log_msg("  ", gse, ": donor_id 去重 ", before, " → ", ncol(expr),
            " 样本；保留：", paste(kept_samples, collapse = ", "))
  }
  x <- t(expr)                                  # 样本 × 基因（全基因，稍后取交集）
  x <- scale(x)                                 # 队列内逐基因 z-score
  y <- factor(ifelse(pheno$condition == label, "case", "control"),
              levels = c("control", "case"))
  if (sum(y == "case") == 0 || sum(y == "control") == 0)
    stop(gse, ": 病例/对照有一侧为 0——请检查 pheno condition")
  list(x = x, y = y, gse = gse)
}

# ----------------------------- 1. 装配训练集 + 平台预筛 ----------------------
log_msg("—— 装配 ", length(CONFIG$TRAIN_SETS), " 个训练队列 ——")
cohorts <- lapply(CONFIG$TRAIN_SETS, assemble_cohort, label = CONFIG$TRAIN_LABEL)

## 平台预筛：特征须在所有训练队列 + 血液平台均可检（先验信息，无泄漏）
val_expr_file <- file.path(CONFIG$IN_DIR, CONFIG$VALIDATION_GSE,
                           paste0(CONFIG$VALIDATION_GSE, "_expr_gene.rds"))
if (!file.exists(val_expr_file)) stop("缺少验证队列预处理文件：", val_expr_file)
val_genes_all <- rownames(readRDS(val_expr_file))
feat <- Reduce(intersect, c(list(core_genes, val_genes_all),
                            lapply(cohorts, function(cc) colnames(cc$x))))
dropped <- setdiff(core_genes, feat)
log_msg("平台预筛：", length(core_genes), " → ", length(feat), " 特征（剔除 ",
        length(dropped), " 个不可检基因：",
        paste(head(dropped, 10), collapse = ", "), if (length(dropped) > 10) " ..." else "", ")")
if (length(feat) < 20) log_msg("⚠ 特征数 <20，交集过窄，请检查各队列基因注释口径")

## 合并（逐队列 z-score 已做，合并即调和；组成留痕）
x_all <- do.call(rbind, lapply(cohorts, function(cc) cc$x[, feat, drop = FALSE]))
y_all <- factor(unlist(lapply(cohorts, function(cc) as.character(cc$y))),
                levels = c("control", "case"))
batch <- factor(rep(CONFIG$TRAIN_SETS, sapply(cohorts, function(cc) nrow(cc$x))))
if (anyDuplicated(rownames(x_all))) stop("合并后样本名重复——请检查各队列 GSM 编号")
## 审查A1：队列内零方差基因经 scale 产生 NaN 列，必须硬校验
if (anyNA(x_all)) {
  bad <- colnames(x_all)[apply(x_all, 2, anyNA)]
  stop("合并矩阵含 NaN（零方差基因经 z-score 产生）：",
       paste(head(bad, 10), collapse = ", "), "——请从特征空间剔除后重跑")
}
comp <- data.frame(cohort = CONFIG$TRAIN_SETS,
                   n = sapply(cohorts, function(cc) nrow(cc$x)),
                   case = sapply(cohorts, function(cc) sum(cc$y == "case")),
                   control = sapply(cohorts, function(cc) sum(cc$y == "control")))
write.csv(comp, file.path(CONFIG$OUT_DIR, "cohort_composition.csv"), row.names = FALSE)
log_msg("合并训练集：", nrow(x_all), " 样本 × ", ncol(x_all), " 特征（case=",
        sum(y_all == "case"), " / control=", sum(y_all == "control"), "）")
## z-score 调和自检：各队列各基因合并后均值应≈0（抽样打印前 3 个基因的批次均值）
chk <- sapply(split(seq_len(nrow(x_all)), batch), function(idx) mean(x_all[idx, 1]))
log_msg("z-score 调和自检（第 1 特征各队列均值）：",
        paste(sprintf("%s=%.2f", names(chk), chk), collapse = " "))

# ----------------------------- 2. 四法特征选择（口径同 v1）-------------------
run_lasso <- function(x, y) {
  cv_fit <- cv.glmnet(x, y, family = "binomial", alpha = 1,
                      nfolds = CONFIG$CV_FOLDS, type.measure = "auc")
  coefs <- coef(cv_fit, s = "lambda.1se")
  setdiff(rownames(coefs)[as.vector(coefs != 0)], "(Intercept)")
}
run_svm_rfe <- function(x, y) {
  ctrl <- rfeControl(functions = caretFuncs, method = "repeatedcv",
                     number = CONFIG$CV_FOLDS, repeats = CONFIG$CV_REPEATS)
  rfe_fit <- rfe(x, y, sizes = seq_len(min(20, ncol(x))), rfeControl = ctrl,
                 method = "svmLinear")
  predictors(rfe_fit)
}
run_rf_importance <- function(x, y) {
  fit <- randomForest(x, y, ntree = 1000, importance = TRUE)
  names(sort(importance(fit, type = 2)[, 1], decreasing = TRUE))
}
run_xgb_importance <- function(x, y) {
  dtrain <- xgb.DMatrix(x, label = as.integer(y) - 1)
  fit <- xgb.train(params = list(objective = "binary:logistic", eval_metric = "auc"),
                   data = dtrain, nrounds = 200, verbose = 0)
  xgb.importance(colnames(x), model = fit)$Feature
}

log_msg("运行四种特征选择器 …")
sel_lasso <- run_lasso(x_all, y_all);          log_msg("LASSO 选出：", length(sel_lasso))
sel_rfe   <- run_svm_rfe(x_all, y_all);        log_msg("SVM-RFE 选出：", length(sel_rfe))
sel_rf    <- run_rf_importance(x_all, y_all)
sel_xgb   <- run_xgb_importance(x_all, y_all)

freq <- table(c(sel_lasso, sel_rfe, head(sel_rf, 10), head(sel_xgb, 10)))
biomarkers <- intersect(sel_lasso, sel_rfe)
if (length(biomarkers) < 2) {
  biomarkers <- names(freq[freq >= 2])
  log_msg("⚠ LASSO∩RFE 不足 2 个，改用频次≥2 规则")
}
if (length(biomarkers) > CONFIG$MAX_BIOMARKERS)
  biomarkers <- names(sort(freq[biomarkers], decreasing = TRUE))[seq_len(CONFIG$MAX_BIOMARKERS)]
freq_df <- data.frame(gene = names(freq), votes = as.integer(freq),
                      selected = names(freq) %in% biomarkers)
write.csv(freq_df, file.path(CONFIG$OUT_DIR, "feature_selection_votes.csv"), row.names = FALSE)
log_msg("✅ v2 候选标志物（", length(biomarkers), " 个）：", paste(biomarkers, collapse = ", "))
if (length(biomarkers) < 2) stop("候选标志物不足 2 个——请回查特征空间与训练集")

## 标志物轨归属注释表
anno_idx <- match(biomarkers, track_tab$gene)
biomarker_table <- data.frame(
  gene = biomarkers,
  track = track_tab$track[anno_idx],
  direction_OA = track_tab$direction_OA[anno_idx],
  direction_muscle = track_tab$direction_muscle[anno_idx],
  in_WGCNA_key = track_tab$in_WGCNA_key[anno_idx],
  votes = as.integer(freq[biomarkers]),
  stringsAsFactors = FALSE
)
write.csv(biomarker_table, file.path(CONFIG$OUT_DIR, "moduleD_v2_biomarker_table.csv"),
          row.names = FALSE)
log_msg("标志物轨分布：", paste(sprintf("%s=%d", names(table(biomarker_table$track)),
                                        as.integer(table(biomarker_table$track))), collapse = " / "))

# ----------------------------- 3. 评估三件套 ---------------------------------
fit_final <- function(x, y, feats)
  suppressWarnings(glm(y ~ ., data = data.frame(y = y, x[, feats, drop = FALSE]),
                       family = binomial))
predict_prob <- function(fit, x, feats)
  as.numeric(predict(fit, newdata = data.frame(x[, feats, drop = FALSE]), type = "response"))

report_performance <- function(prob, y, tag) {
  ## 审查A4：glm 完全分离时 predict 可能返回 NA，roc() 报错信息不直观，前置拦截
  if (anyNA(prob)) stop(tag, ": 预测概率含 NA（疑似 glm 完全分离）——请回查标志物与数据")
  ## 显式固定方向（control < case），防 pROC auto 静默翻转掩盖外部验证失败
  roc_obj <- roc(y, prob, levels = c("control", "case"), direction = "<", quiet = TRUE)
  ci_obj <- ci.auc(roc_obj, method = "bootstrap", boot.n = CONFIG$N_BOOT)
  brks <- unique(quantile(prob, probs = seq(0, 1, 0.1)))
  if (length(brks) < 4) brks <- unique(seq(min(prob), max(prob), length.out = 6))
  if (length(brks) >= 2) {
    cal_df <- data.frame(bin = cut(prob, breaks = brks, include.lowest = TRUE),
                         obs = as.integer(y) - 1, pred = prob)
    cal <- aggregate(cbind(obs, pred) ~ bin, data = cal_df, mean)
    write.csv(cal, file.path(CONFIG$OUT_DIR, paste0("calibration_", tag, ".csv")), row.names = FALSE)
  } else log_msg("  ⚠ ", tag, " 预测值退化，校准曲线跳过")
  thresholds <- seq(0.01, 0.8, by = 0.01)
  dca <- do.call(rbind, lapply(thresholds, function(pt) {
    pred_pos <- prob >= pt
    data.frame(threshold = pt,
               net_benefit = sum(pred_pos & y == "case") / length(y) -
                 sum(pred_pos & y == "control") / length(y) * (pt / (1 - pt)))
  }))
  write.csv(dca, file.path(CONFIG$OUT_DIR, paste0("dca_", tag, ".csv")), row.names = FALSE)
  log_msg(tag, ": AUC = ", round(as.numeric(auc(roc_obj)), 3),
          " (bootstrap 95%CI ", round(ci_obj[1], 3), "–", round(ci_obj[3], 3), ")")
  attr(roc_obj, "ci_bootstrap") <- ci_obj
  roc_obj
}

## 内部验证：重复 10×10 CV out-of-fold（注意：特征选择在 CV 外，AUC 偏乐观，仅参考）
cv_predict <- function(x, y, feats) {
  ctrl <- trainControl(method = "repeatedcv", number = CONFIG$CV_FOLDS,
                       repeats = CONFIG$CV_REPEATS, savePredictions = "final",
                       classProbs = TRUE)
  fit <- train(x[, feats, drop = FALSE], y, method = "glm",
               trControl = ctrl, family = binomial)
  oof <- fit$pred[fit$pred$Resample != "", ]
  aggregate(case ~ rowIndex, data = oof, mean)
}

oof <- cv_predict(x_all, y_all, biomarkers)
y_oof <- y_all[match(oof$rowIndex, seq_along(y_all))]
roc_internal <- report_performance(oof$case, y_oof, "internal_cv_multicohort")

# ----------------------------- 4. 外部验证（GSE48556 血液，纯验证）-----------
val <- assemble_cohort(CONFIG$VALIDATION_GSE, CONFIG$VALIDATION_LABEL)
## 平台预筛已保证标志物在血液中全部可检；仍做一次防御性交集
## 审查A5：用显式标志变量判断外部验证是否执行（防同会话重跑时 exists() 残留误报）
common <- intersect(biomarkers, colnames(val$x))
ext_done <- length(common) >= 2
if (length(common) != length(biomarkers))
  log_msg("⚠ 血液队列缺 ", length(biomarkers) - length(common), " 个标志物（平台预筛后不应发生，请回查）")
if (ext_done) {
  fit <- fit_final(x_all, y_all, common)
  prob_val <- predict_prob(fit, val$x[, common, drop = FALSE], common)
  roc_external <- report_performance(prob_val, val$y,
                                     paste0("external_blood_", CONFIG$VALIDATION_GSE))
  saveRDS(list(roc = roc_external, fit = fit),
          file.path(CONFIG$OUT_DIR, "roc_external_blood.rds"))

  ## 方向一致性核对（血液实测 vs 软骨发现层）
  dir_check <- do.call(rbind, lapply(common, function(g) {
    obs <- ifelse(mean(val$x[val$y == "case", g]) > mean(val$x[val$y == "control", g]),
                  "up", "down")
    idx <- match(g, track_tab$gene)
    data.frame(gene = g, track = track_tab$track[idx],
               direction_OA_discovery = track_tab$direction_OA[idx],
               direction_blood_observed = obs,
               consistent = obs == track_tab$direction_OA[idx],
               stringsAsFactors = FALSE)
  }))
  write.csv(dir_check, file.path(CONFIG$OUT_DIR,
            paste0("direction_check_", CONFIG$VALIDATION_GSE, ".csv")), row.names = FALSE)
  log_msg("方向一致性（血液 vs 软骨发现层）：", sum(dir_check$consistent), "/",
          nrow(dir_check), "（", round(100 * mean(dir_check$consistent), 1), "%）")
}

saveRDS(list(biomarkers = biomarkers, biomarker_table = biomarker_table,
             votes = freq_df, features_space = feat,
             cohort_composition = comp, roc_internal = roc_internal),
        file.path(CONFIG$OUT_DIR, "moduleD_v2_biomarkers.rds"))

cat("\n模块 D v2 完成（多队列训练 + 平台预筛）。产物：\n",
    "  ", CONFIG$OUT_DIR, "/moduleD_v2_biomarker_table.csv —— 标志物轨归属/方向/置信度总表\n",
    "  ", CONFIG$OUT_DIR, "/cohort_composition.csv —— 五队列组成留痕\n",
    "  ", CONFIG$OUT_DIR, "/feature_selection_votes.csv —— 四法投票留痕\n",
    "  ", CONFIG$OUT_DIR, "/dca_*.csv / calibration_*.csv —— 临床效用与校准\n",
    if (ext_done)
      paste0("  ", CONFIG$OUT_DIR, "/direction_check_", CONFIG$VALIDATION_GSE,
             ".csv —— 血液方向一致性核对\n") else "",
    "  ", CONFIG$OUT_DIR, "/moduleD_v2_biomarkers.rds —— v2 锁定标志物（模块 E/H 输入）\n",
    "下一步：血液 AUC 与 v1（0.642）对比；软骨主验证在模块 E 进行。\n",
    "注意：内部 CV AUC 因特征选择在 CV 外而偏乐观，正文仅作参考、权重在外部验证。\n",
    sep = "")
