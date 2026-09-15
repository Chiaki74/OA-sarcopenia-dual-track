# ============================================================================
# 05_moduleD_ML.R
# 模块 D｜机器学习标志物筛选与诊断模型（计划书 v1.4a §3.3.4）
# ----------------------------------------------------------------------------
# 设计（对齐 TRIPOD+AI 报告规范，Collins et al. 2024 BMJ）：
#   1) 输入特征 = 双轨核心基因全量 99 个（备忘录 §八-24/26：同向轨 27 + 镜像轨
#      72 并行进入下游；特征带轨标签与方向注释；不用全转录组，避免维度灾难
#      与过拟合——小样本 bulk 队列的关键防线）；
#   2) 三种互补特征选择器取交集：LASSO（线性稀疏）+ SVM-RFE（非线性边界）
#      + RF/XGBoost 重要性（非线性交互）→ 候选标志物（预期 3–6 个）；
#      入选标志物标注轨归属（concordant/mirror↑↓/mirror↓↑）与 WGCNA 置信度；
#   3) 模型在发现队列内部用重复交叉验证评估，再在独立外部队列
#      （默认 GSE48556 外周血）做外部验证（TRIPOD 要求的最高证据等级）；
#      外部验证同时报告方向一致性（备忘录 §八-24 双轨框架要求）；
#   4) 报告三件套：ROC/AUC（区分度）+ 校准曲线（校准度）+ DCA（临床效用），
#      AUC 置信区间用 2000 次 bootstrap；
#   5) 全部随机种子固定；中间结果落盘可复现。
# ----------------------------------------------------------------------------
# 前置：results/moduleC/moduleC_gene_tracks.csv（模块 C 产物，双轨 99 基因总表）；
#       processed/<GSE>/ 表达矩阵与 pheno（condition 已人工核对）
# 运行：source("05_moduleD_ML.R")；结果输出 results/moduleD/
# ============================================================================

suppressPackageStartupMessages({
  library(glmnet)       # LASSO
  library(e1071)        # SVM（RFE 用）
  library(randomForest)
  library(xgboost)
  library(pROC)
  library(caret)        # 重复 CV 框架 + 校准
  library(data.table)
  library(ggplot2)
})
set.seed(20260820)

# ----------------------------- CONFIG --------------------------------------
CONFIG <- list(
  TRACK_FILE = "results/moduleC/moduleC_gene_tracks.csv",  # 双轨 99 基因总表（2026-08-23 适配）
  IN_DIR  = "processed",
  OUT_DIR = "results/moduleD",
  LOG_FILE = file.path("results/moduleD", sprintf("moduleD_log_%s.txt", Sys.Date())),  # 审查A5：日志入 OUT_DIR
  ZSCORE_FEATURES = TRUE,   # 审查A4：各队列内部逐基因 z-score，提升跨队列/跨平台可迁移性
  ## 发现队列（训练+内部验证）：默认取模块 A 处理后样本量最大的 OA 数据集
  DISCOVERY_GSE   = "GSE114007",
  DISCOVERY_LABEL = "OA",
  ## 外部验证队列：GSE48556 外周血 PBMC——OA 队列（GARP 研究，106 OA vs 33 对照，
  ## 全女性、同胞对家系相关；GEO 已核实，2026-08-22 更正此前的 sarcopenia 误标）。
  ## 用途：验证共享标志物在 OA 人群血液中的可检测性——计划书 §3.3.5 模块 E 的
  ## "血液转化卖点"。标签必须是 OA；若错配为 sarcopenia，外部验证必然失败。
  VALIDATION_GSE   = "GSE48556",
  VALIDATION_LABEL = "OA",
  CV_REPEATS = 10, CV_FOLDS = 10,
  N_BOOT = 2000,
  MAX_BIOMARKERS = 10     # 交集过大时按出现频次+重要性截断
)

log_msg <- function(...) {
  txt <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = ""))
  cat(txt, "\n"); cat(txt, "\n", file = CONFIG$LOG_FILE, append = TRUE)
}
dir.create(CONFIG$OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ----------------------------- 数据装配 -------------------------------------
assemble_xy <- function(gse, core_genes, disease_label) {
  expr  <- readRDS(file.path(CONFIG$IN_DIR, gse, paste0(gse, "_expr_gene.rds")))
  pheno <- read.csv(file.path(CONFIG$IN_DIR, gse, paste0(gse, "_pheno.csv")),
                    stringsAsFactors = FALSE)
  pheno <- pheno[match(colnames(expr), pheno$sample), ]
  keep <- pheno$condition %in% c("control", disease_label)
  expr <- expr[, keep, drop = FALSE]; pheno <- pheno[keep, ]
  genes_use <- intersect(core_genes, rownames(expr))
  missing <- setdiff(core_genes, genes_use)
  if (length(missing) > 0)
    log_msg(gse, ": ", length(missing), " 个核心基因在该队列缺失：",
            paste(head(missing, 10), collapse = ", "))
  x <- t(expr[genes_use, , drop = FALSE])          # 样本 × 基因
  if (CONFIG$ZSCORE_FEATURES)
    x <- scale(x)   # 队列内逐基因标准化：跨组织（软骨→血液）迁移的关键（审查A4）
  y <- factor(ifelse(pheno$condition == disease_label, "case", "control"),
              levels = c("control", "case"))
  if (sum(y == "case") == 0 || sum(y == "control") == 0)
    stop(gse, ": 标签 '", disease_label, "' 下病例/对照有一侧为 0——",
         "请检查 DISCOVERY/VALIDATION 标签与该队列 pheno 的 condition 是否匹配")
  list(x = x, y = y, genes = genes_use)
}

# ---- 双轨基因总表读取与口径校验（2026-08-23 适配，备忘录 §八-26）----
if (!file.exists(CONFIG$TRACK_FILE))
  stop("缺少 ", CONFIG$TRACK_FILE, "——请先运行 04_moduleC_enrichment.R")
track_tab <- read.csv(CONFIG$TRACK_FILE, stringsAsFactors = FALSE)
valid_tracks <- c("concordant", "mirror_OAup_muscleDown", "mirror_OAdown_muscleUp")
if (!all(c("gene", "track", "direction_OA", "direction_muscle", "in_WGCNA_key") %in% names(track_tab)))
  stop("moduleC_gene_tracks.csv 列名口径漂移，请核对模块 C 输出")
if (any(!track_tab$track %in% valid_tracks))
  stop("存在未知轨标签：", paste(setdiff(track_tab$track, valid_tracks), collapse = ","))
if (any(!track_tab$direction_OA %in% c("up", "down")) ||
    any(!track_tab$direction_muscle %in% c("up", "down")) ||
    any(is.na(track_tab$direction_OA)) || any(is.na(track_tab$direction_muscle)))
  stop("direction 列存在非 up/down 或 NA 取值，请核对模块 C 输出（审查A3）")
if (any(is.na(track_tab$gene)) || anyDuplicated(track_tab$gene))
  stop("基因列存在 NA 或重复，请核对模块 C 输出")
core_genes <- track_tab$gene
log_msg("核心基因数：", length(core_genes), "（双轨框架：同向轨 ",
        sum(track_tab$track == "concordant"), " + 镜像轨 ",
        sum(track_tab$track != "concordant"), "；预期 27+72=99）")
disc <- assemble_xy(CONFIG$DISCOVERY_GSE, core_genes, CONFIG$DISCOVERY_LABEL)
log_msg("发现队列 ", CONFIG$DISCOVERY_GSE, "：", nrow(disc$x), " 样本 × ",
        ncol(disc$x), " 特征（case=", sum(disc$y == "case"), "）")

# ----------------------------- 1. LASSO -------------------------------------
run_lasso <- function(x, y) {
  cv_fit <- cv.glmnet(x, y, family = "binomial", alpha = 1,
                      nfolds = CONFIG$CV_FOLDS, type.measure = "auc")
  coefs <- coef(cv_fit, s = "lambda.1se")   # 1se 规则：更稀疏、更保守
  sel <- rownames(coefs)[as.vector(coefs != 0)]
  setdiff(sel, "(Intercept)")
}

# ----------------------------- 2. SVM-RFE -----------------------------------
run_svm_rfe <- function(x, y) {
  ## caret::rfe 封装递归特征消除（内部以线性 SVM 权重排序），重复 CV 选最优特征数
  ctrl <- rfeControl(functions = caretFuncs, method = "repeatedcv",
                     number = CONFIG$CV_FOLDS, repeats = CONFIG$CV_REPEATS)
  rfe_fit <- rfe(x, y, sizes = seq_len(min(20, ncol(x))), rfeControl = ctrl,
                 method = "svmLinear")
  predictors(rfe_fit)
}

# ----------------------------- 3. RF / XGBoost ------------------------------
run_rf_importance <- function(x, y) {
  fit <- randomForest(x, y, ntree = 1000, importance = TRUE)
  imp <- importance(fit, type = 2)[, 1]     # MeanDecreaseGini
  names(sort(imp, decreasing = TRUE))
}
run_xgb_importance <- function(x, y) {
  dtrain <- xgb.DMatrix(x, label = as.integer(y) - 1)
  fit <- xgb.train(params = list(objective = "binary:logistic", eval_metric = "auc"),
                   data = dtrain, nrounds = 200, verbose = 0)
  imp <- xgb.importance(colnames(x), model = fit)
  imp$Feature
}

# ----------------------------- 特征交集 -------------------------------------
log_msg("运行三种特征选择器 …")
sel_lasso <- run_lasso(disc$x, disc$y);          log_msg("LASSO 选出：", length(sel_lasso))
sel_rfe   <- run_svm_rfe(disc$x, disc$y);        log_msg("SVM-RFE 选出：", length(sel_rfe))
sel_rf    <- run_rf_importance(disc$x, disc$y)   # 全排序
sel_xgb   <- run_xgb_importance(disc$x, disc$y)  # 重要性排序

## 交集策略：LASSO ∩ RFE 为主；若为空/过少，并入 RF+XGBoost top-10 中
## 在 ≥2 种方法中出现的基因（频次表全部留痕）
freq <- table(c(sel_lasso, sel_rfe, head(sel_rf, 10), head(sel_xgb, 10)))
biomarkers <- intersect(sel_lasso, sel_rfe)
if (length(biomarkers) < 2) {
  biomarkers <- names(freq[freq >= 2])
  log_msg("⚠ LASSO∩RFE 不足 2 个，改用频次≥2 规则")
}
if (length(biomarkers) > CONFIG$MAX_BIOMARKERS) {
  biomarkers <- names(sort(freq[biomarkers], decreasing = TRUE))[seq_len(CONFIG$MAX_BIOMARKERS)]
}
freq_df <- data.frame(gene = names(freq), votes = as.integer(freq),
                      selected = names(freq) %in% biomarkers)
write.csv(freq_df, file.path(CONFIG$OUT_DIR, "feature_selection_votes.csv"), row.names = FALSE)
log_msg("✅ 候选标志物（", length(biomarkers), " 个）：", paste(biomarkers, collapse = ", "))

# ---- 标志物轨归属注释表（双轨框架核心产物，备忘录 §八-24/26）----
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
write.csv(biomarker_table, file.path(CONFIG$OUT_DIR, "moduleD_biomarker_table.csv"),
          row.names = FALSE)
log_msg("标志物轨分布：", paste(sprintf("%s=%d", names(table(biomarker_table$track)),
                                        as.integer(table(biomarker_table$track))), collapse = " / "))

# ----------------------------- 模型评估三件套 -------------------------------
## 逻辑回归作为最终可解释模型（小样本下比复杂模型更稳，TRIPOD 推荐）
fit_final <- function(x, y, feats) {
  df <- data.frame(y = y, x[, feats, drop = FALSE])
  suppressWarnings(glm(y ~ ., data = df, family = binomial))
}
predict_prob <- function(fit, x, feats) {
  df <- data.frame(x[, feats, drop = FALSE])
  as.numeric(predict(fit, newdata = df, type = "response"))
}

report_performance <- function(prob, y, tag) {
  ## 审查S1：显式固定方向（control < case），禁用 auto——
  ## 否则外部验证翻车（AUC<0.5）会被 pROC 静默翻转成 >0.5，掩盖失败信号
  roc_obj <- roc(y, prob, levels = c("control", "case"), direction = "<", quiet = TRUE)
  ## AUC 置信区间：显式 bootstrap（roc(ci=TRUE) 默认 DeLong，boot.n 会被静默吞掉）
  ci_obj <- ci.auc(roc_obj, method = "bootstrap", boot.n = CONFIG$N_BOOT)
  ## 校准曲线（10 分位；预测值有并列时分位数断点不唯一，去重后不足则降级等宽分箱；
  ## 审查A1：退化预测（常数 prob）下断点仍可能不唯一，此时跳过校准输出而非崩溃）
  brks <- unique(quantile(prob, probs = seq(0, 1, 0.1)))
  if (length(brks) < 4) brks <- unique(seq(min(prob), max(prob), length.out = 6))
  if (length(brks) >= 2) {
    cal_df <- data.frame(bin = cut(prob, breaks = brks, include.lowest = TRUE),
                         obs = as.integer(y) - 1, pred = prob)
    cal <- aggregate(cbind(obs, pred) ~ bin, data = cal_df, mean)
    write.csv(cal, file.path(CONFIG$OUT_DIR, paste0("calibration_", tag, ".csv")), row.names = FALSE)
  } else {
    cal <- NULL
    log_msg("  ⚠ ", tag, " 预测值退化（近常数），校准曲线跳过")
  }
  ## DCA（净获益曲线，阈值 0–0.8）
  thresholds <- seq(0.01, 0.8, by = 0.01)
  dca <- do.call(rbind, lapply(thresholds, function(pt) {
    pred_pos <- prob >= pt
    tp <- sum(pred_pos & y == "case"); fp <- sum(pred_pos & y == "control")
    n <- length(y)
    data.frame(threshold = pt,
               net_benefit = tp / n - fp / n * (pt / (1 - pt)))
  }))
  write.csv(dca, file.path(CONFIG$OUT_DIR, paste0("dca_", tag, ".csv")), row.names = FALSE)
  log_msg(tag, ": AUC = ", round(as.numeric(auc(roc_obj)), 3),
          " (bootstrap 95%CI ", round(ci_obj[1], 3), "–", round(ci_obj[3], 3), ")")
  attr(roc_obj, "ci_bootstrap") <- ci_obj
  roc_obj
}

## 内部验证：重复 10×10 CV 的 out-of-fold 预测
cv_predict <- function(x, y, feats) {
  ctrl <- trainControl(method = "repeatedcv", number = CONFIG$CV_FOLDS,
                       repeats = CONFIG$CV_REPEATS, savePredictions = "final",
                       classProbs = TRUE)
  fit <- train(x[, feats, drop = FALSE], y, method = "glm",
               trControl = ctrl, family = binomial)
  oof <- fit$pred[fit$pred$Resample != "", ]
  aggregate(case ~ rowIndex, data = oof, mean)  # 每样本跨重复平均
}

# ----------------------------- 主流程 ---------------------------------------
if (length(biomarkers) < 2) stop("候选标志物不足 2 个，模块 D 无法建模——请回查模块 B 核心集")

## 内部验证
oof <- cv_predict(disc$x, disc$y, biomarkers)
y_oof <- disc$y[match(oof$rowIndex, seq_along(disc$y))]
roc_internal <- report_performance(oof$case, y_oof, "internal_cv")

## 外部验证（GSE48556 外周血）
val_file <- file.path(CONFIG$IN_DIR, CONFIG$VALIDATION_GSE,
                      paste0(CONFIG$VALIDATION_GSE, "_expr_gene.rds"))
if (file.exists(val_file)) {
  val <- assemble_xy(CONFIG$VALIDATION_GSE, biomarkers, CONFIG$VALIDATION_LABEL)
  common <- intersect(biomarkers, val$genes)
  if (length(common) >= max(2, floor(length(biomarkers) * 0.7))) {
    fit <- fit_final(disc$x, disc$y, common)
    prob_val <- predict_prob(fit, val$x, common)
    roc_external <- report_performance(prob_val, val$y, paste0("external_", CONFIG$VALIDATION_GSE))
    saveRDS(list(roc = roc_external, fit = fit),
            file.path(CONFIG$OUT_DIR, "roc_external.rds"))

    ## 方向一致性核对（双轨框架要求，备忘录 §八-24）：
    ## 外部血液队列中各标志物的实测方向 vs 发现层 direction_OA。
    ## 注意：血液与软骨方向不一致不构成模型失败（组织特异性正是镜像叙事的一部分），
    ## 仅作如实报告；全反或不相关才需回查标签/批次。
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
            nrow(dir_check), " 一致（",
            round(100 * mean(dir_check$consistent), 1), "%）——详见 direction_check CSV")
  } else {
    log_msg("⚠ 外部队列仅覆盖 ", length(common), "/", length(biomarkers),
            " 个标志物，外部验证不执行——考虑换用全基因模型或降低标志物数")
  }
} else {
  log_msg("⚠ 外部验证队列 ", CONFIG$VALIDATION_GSE, " 尚未预处理，跳过外部验证")
}

saveRDS(list(biomarkers = biomarkers, biomarker_table = biomarker_table,
             votes = freq_df, roc_internal = roc_internal),
        file.path(CONFIG$OUT_DIR, "moduleD_biomarkers.rds"))

## 审查A2：产物清单按外部验证是否实际执行动态组织，防误导
ext_done <- exists("roc_external")
cat("\n模块 D 完成（双轨框架）。产物：\n",
    "  ", CONFIG$OUT_DIR, "/moduleD_biomarker_table.csv —— 标志物轨归属/方向/WGCNA 置信度总表\n",
    "  ", CONFIG$OUT_DIR, "/feature_selection_votes.csv —— 四法投票留痕\n",
    "  ", CONFIG$OUT_DIR, "/dca_*.csv / calibration_*.csv —— 临床效用与校准\n",
    if (ext_done)
      paste0("  ", CONFIG$OUT_DIR, "/direction_check_", CONFIG$VALIDATION_GSE,
             ".csv —— 血液方向一致性核对\n") else "",
    "  ", CONFIG$OUT_DIR, "/moduleD_biomarkers.rds —— 锁定标志物（模块 E/H 的输入）\n",
    if (ext_done)
      paste0("下一步：外部验证 AUC 若 <0.7，回查 direction_check_", CONFIG$VALIDATION_GSE,
             " 与批次问题；\n方向不一致不必然是失败——血液 vs 软骨的组织特异性正是镜像叙事的一部分。\n")
    else
      "下一步：外部验证未执行（队列缺失或标志物覆盖不足），请先解决后重跑。\n",
    sep = "")
