# ============================================================================
# 05c_blood_exploratory_subsignature.R
# 路径 a｜血液探索性子签名（2026-08-23，备忘录 §八-28）
# ----------------------------------------------------------------------------
# ★ 探索性分析声明（hypothesis-generating）：
#   本脚本全部产物标注 EXPLORATORY——GSE48556 是唯一公开 OA 血液队列，
#   子签名的筛选与评估在同一队列内完成，虽有嵌套 CV 防泄漏，仍无独立
#   血液验证队列，结论只能作假设生成，正文须写入 Limitations。
# 设计：
#   1) 特征空间 = 84 个平台可检双轨基因（moduleD_v2 rds 的 features_space）；
#   2) 嵌套重复 10×10 CV：每折内仅用训练折做 limma DE（名义 p<0.05，
#      上限 20 个按 p 排序截断）→ glm 拟合 → 测试折预测；
#      特征选择在 CV 内部重做，OOF AUC 无选择泄漏；
#   3) 对照报告：全数据 DE 子签名（标注"表观、偏乐观"）与 v2 十基因
#      面板在血液中的表现回顾（引自模块 D 结果，不重算）；
#   4) 方向对照：血液 DE 方向 vs 软骨 direction_OA，检验"下调轴可测性"
#      模式是否稳定。
# ----------------------------------------------------------------------------
# 前置：results/moduleC/moduleC_gene_tracks.csv；
#       results/moduleD_v2/moduleD_v2_biomarkers.rds；
#       processed/GSE48556/ 已预处理
# 运行：source("05c_blood_exploratory_subsignature.R")；输出 results/moduleD_v2/blood_exploratory/
# ============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(pROC)
  library(data.table)
})
set.seed(20260823)

CONFIG <- list(
  TRACK_FILE = "results/moduleC/moduleC_gene_tracks.csv",
  V2_RDS     = "results/moduleD_v2/moduleD_v2_biomarkers.rds",
  IN_DIR  = "processed",
  BLOOD_GSE = "GSE48556",
  BLOOD_LABEL = "OA",
  OUT_DIR = "results/moduleD_v2/blood_exploratory",
  LOG_FILE = file.path("results/moduleD_v2/blood_exploratory",
                       sprintf("blood_exploratory_log_%s.txt", Sys.Date())),
  DE_P_CUT = 0.05,        # 名义 p（探索性，不做 FDR——假设生成口径）
  MAX_FEATS_PER_FOLD = 20,
  CV_REPEATS = 10, CV_FOLDS = 10,
  N_BOOT = 2000
)
dir.create(CONFIG$OUT_DIR, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  txt <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = ""))
  cat(txt, "\n"); cat(txt, "\n", file = CONFIG$LOG_FILE, append = TRUE)
}
log_msg("★ 探索性分析（EXPLORATORY / hypothesis-generating）——备忘录 §八-28 路径 a")

# ---- 读取双轨表与 v2 特征空间 ----
track_tab <- read.csv(CONFIG$TRACK_FILE, stringsAsFactors = FALSE)
v2 <- readRDS(CONFIG$V2_RDS)
feat_space <- v2$features_space
if (is.null(feat_space) || length(feat_space) < 10)
  stop("v2 rds 中 features_space 缺失或过短——请确认 05b 已成功运行")
log_msg("特征空间：", length(feat_space), " 个平台可检双轨基因（来自 moduleD_v2）")

# ---- 装配血液队列（z-score；无 family_id，走非配对 fallback，备忘录 §八-16）----
expr  <- readRDS(file.path(CONFIG$IN_DIR, CONFIG$BLOOD_GSE,
                           paste0(CONFIG$BLOOD_GSE, "_expr_gene.rds")))
pheno <- read.csv(file.path(CONFIG$IN_DIR, CONFIG$BLOOD_GSE,
                            paste0(CONFIG$BLOOD_GSE, "_pheno.csv")), stringsAsFactors = FALSE)
pheno <- pheno[match(colnames(expr), pheno$sample), ]
if (anyNA(pheno$condition)) stop("pheno 缺少 expr 中的样本——请回查模块 A 输出")
keep <- pheno$condition %in% c("control", CONFIG$BLOOD_LABEL)
expr <- expr[, keep, drop = FALSE]; pheno <- pheno[keep, ]
genes_use <- intersect(feat_space, rownames(expr))
log_msg("血液队列：", ncol(expr), " 样本（case=", sum(pheno$condition == CONFIG$BLOOD_LABEL),
        "）；可检特征 ", length(genes_use), "/", length(feat_space))
## 审查 S2：交集过小会让 100 折全部空跑
if (length(genes_use) < 10)
  stop("血液中可检双轨基因仅 ", length(genes_use), " 个（<10）——交集过窄，请回查平台注释")
x <- scale(t(expr[genes_use, , drop = FALSE]))   # 样本 × 基因，队列内逐基因 z-score
if (anyNA(x)) stop("z-score 产生 NaN（零方差基因）——请回查")
y <- factor(ifelse(pheno$condition == CONFIG$BLOOD_LABEL, "case", "control"),
            levels = c("control", "case"))

# ---- 折内 DE 初筛函数（limma，非配对） ----
de_select <- function(x_tr, y_tr) {
  design <- model.matrix(~ y_tr)
  fit <- eBayes(lmFit(t(x_tr), design))   # limma 需要 基因×样本
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
  tt <- tt[tt$P.Value < CONFIG$DE_P_CUT, , drop = FALSE]
  if (nrow(tt) == 0) return(character(0))
  head(rownames(tt), CONFIG$MAX_FEATS_PER_FOLD)
}

# ---- 嵌套重复 CV：每折内重选特征 + 重训 glm ----
log_msg("运行嵌套 ", CONFIG$CV_REPEATS, "×", CONFIG$CV_FOLDS, " CV（折内 DE 重选，防泄漏）…")
oof_list <- list(); fold_genes <- list(); k <- 0
n_notconverged <- 0   # 审查 A4：glm 不收敛/分离折计数留痕
for (rep in seq_len(CONFIG$CV_REPEATS)) {
  folds <- caret::createFolds(y, k = CONFIG$CV_FOLDS, returnTrain = FALSE)
  for (fi in seq_along(folds)) {
    k <- k + 1
    te <- folds[[fi]]; tr <- setdiff(seq_len(nrow(x)), te)
    gs <- de_select(x[tr, , drop = FALSE], y[tr])
    fold_genes[[k]] <- gs
    if (length(gs) < 2) next
    fit <- suppressWarnings(glm(y ~ ., family = binomial,
                                data = data.frame(y = y[tr], x[tr, gs, drop = FALSE])))
    if (!fit$converged) n_notconverged <- n_notconverged + 1
    pr <- as.numeric(predict(fit, newdata = data.frame(x[te, gs, drop = FALSE]),
                             type = "response"))
    if (!anyNA(pr)) oof_list[[length(oof_list) + 1]] <- data.frame(idx = te, prob = pr)
  }
}
if (n_notconverged > 0)
  log_msg("  ⚠ ", n_notconverged, "/", k, " 折 glm 不收敛（疑似分离），OOF 口径已留痕")
## 审查 S2：全部折空跑时显式报错，而非 aggregate(NULL) 误导性崩溃
if (length(oof_list) == 0)
  stop("所有折折内 DE 入选基因均 <2——血液中无有效信号，请放宽 DE_P_CUT 或回查数据")
oof_all <- do.call(rbind, oof_list)
oof <- aggregate(prob ~ idx, data = oof_all, mean)
y_oof <- y[oof$idx]
roc_cv <- roc(y_oof, oof$prob, levels = c("control", "case"), direction = "<", quiet = TRUE)
ci_cv <- ci.auc(roc_cv, method = "bootstrap", boot.n = CONFIG$N_BOOT)
log_msg("★ 嵌套 CV OOF AUC = ", round(as.numeric(auc(roc_cv)), 3),
        " (bootstrap 95%CI ", round(ci_cv[1], 3), "–", round(ci_cv[3], 3),
        ")——特征选择无泄漏的探索性估计；", length(oof_list), " 折有效",
        "（口径声明：z-score 为全队列无监督预处理；CI 未校正重复 CV 折间相关性，偏窄——审查 S1/A5）")

## 折内入选频次（稳定性指标）
fg <- table(unlist(fold_genes))
stab <- data.frame(gene = names(fg), folds_selected = as.integer(fg),
                   stability = round(as.integer(fg) / k, 3))
idx <- match(stab$gene, track_tab$gene)
stab$track <- track_tab$track[idx]
stab$direction_OA <- track_tab$direction_OA[idx]
stab <- stab[order(-stab$folds_selected), ]
write.csv(stab, file.path(CONFIG$OUT_DIR, "EXPLORATORY_blood_fold_selection_stability.csv"),
          row.names = FALSE)
log_msg("折内入选稳定性 top5：",
        paste(sprintf("%s(%.0f%%)", head(stab$gene, 5), 100 * head(stab$stability, 5)),
              collapse = ", "))

# ---- 全数据 DE 子签名（表观结果，偏乐观，仅供展示对照） ----
design <- model.matrix(~ y)
fit_full <- eBayes(lmFit(t(x), design))
tt_full <- topTable(fit_full, coef = 2, number = Inf, sort.by = "P")
blood_de <- tt_full[tt_full$P.Value < CONFIG$DE_P_CUT, ]
blood_de$gene <- rownames(blood_de)
idx <- match(blood_de$gene, track_tab$gene)
if (anyNA(idx)) stop("血液 DE 基因不在双轨表中——口径漂移，请回查（审查 A1）")   # 审查 A1
blood_de$track <- track_tab$track[idx]
blood_de$direction_OA <- track_tab$direction_OA[idx]
## 审查 A3：x 已 z-score，logFC 实为标准化均值差，改名防误导
blood_de$blood_z_diff <- blood_de$logFC; blood_de$logFC <- NULL
blood_de$blood_direction <- ifelse(blood_de$blood_z_diff > 0, "up", "down")
blood_de$consistent_with_cartilage <- blood_de$blood_direction == blood_de$direction_OA
write.csv(blood_de, file.path(CONFIG$OUT_DIR, "EXPLORATORY_blood_DE_dualtrack_genes.csv"),
          row.names = FALSE)
## 审查 A2：空表防 NaN%
if (nrow(blood_de) > 0) {
  log_msg("全数据血液 DE（名义 p<0.05，表观/偏乐观）：", nrow(blood_de), "/", length(genes_use),
          " 个；与软骨方向一致 ", sum(blood_de$consistent_with_cartilage), " 个（",
          round(100 * mean(blood_de$consistent_with_cartilage), 1), "%）")
} else log_msg("全数据血液 DE：0 个名义显著——血液信号弱，如实报告")
## 下调轴可测性模式检验
dn <- blood_de[blood_de$direction_OA == "down", ]
up <- blood_de[blood_de$direction_OA == "up", ]
if (nrow(dn) > 0 && nrow(up) > 0)
  log_msg("下调轴模式：软骨↓基因血液一致率 ", round(100 * mean(dn$consistent_with_cartilage), 1),
          "% vs 软骨↑基因 ", round(100 * mean(up$consistent_with_cartilage), 1), "%")

## 与 v2 十基因面板交集
overlap <- intersect(v2$biomarkers, blood_de$gene)
log_msg("v2 十基因中在血液名义 DE 的：", length(overlap), " 个（",
        paste(overlap, collapse = ", "), "）")

cat("\n路径 a 完成（★探索性，hypothesis-generating）。产物：\n",
    "  ", CONFIG$OUT_DIR, "/EXPLORATORY_blood_fold_selection_stability.csv —— 嵌套 CV 折内入选稳定性\n",
    "  ", CONFIG$OUT_DIR, "/EXPLORATORY_blood_DE_dualtrack_genes.csv —— 血液 DE × 双轨方向对照\n",
    "正文口径：嵌套 CV OOF AUC 可报（标注探索性）；全数据 DE 结果仅作展示、须声明偏乐观；\n",
    "无第二独立血液队列，全部结论 hypothesis-generating，写入 Limitations。\n", sep = "")
