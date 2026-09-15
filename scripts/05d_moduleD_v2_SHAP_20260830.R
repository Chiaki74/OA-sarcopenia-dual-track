# ============================================================================
# 05d_moduleD_v2_SHAP_20260830.R
# 模块 D v2 补丁｜最终 logistic 模型的 SHAP 可解释性分析（2026-08-30，备忘录 §50-1）
# ----------------------------------------------------------------------------
# 裁定背景：写作前清点时研究者拍板"SHAP 做"（同类 ML 标志物论文标配，
#   防审稿人点名模型可解释性）。
# 技术口径：
#   1) 模块 D v2 最终模型 = 10 标志物 logistic 回归（glm, binomial），
#      在 link scale 上是线性模型 → SHAP 有精确解析解：
#        phi_j(i) = beta_j * (x_ij - mean(x_j))    （baseline = 训练集均值）
#      满足效率性：sum_j phi_j(i) + intercept = linear predictor(i)
#      不依赖 kernelshap / shapviz / treeshap 等外部包，规避依赖风险。
#   2) 训练矩阵不在 moduleD_v2_biomarkers.rds 中 → 本脚本逐字复用 05b 的
#      装配链（逐队列 z-score + GSE51588 供体去重 + 平台预筛）重建 117 样本
#      训练矩阵，并用 roc_external_blood.rds 内封存 fit 的系数做逐位一致性
#      自检（装配若漂移，系数必然不同，直接 stop）。
#   3) 标志物取自 moduleD_v2_biomarkers.rds（v2 锁定口径），不重跑四法选择。
# 前置：05b_moduleD_v2_multicohort.R 已实跑（results/moduleD_v2/ 产物在位）；
#       5 个训练集 + GSE48556 已经 02_moduleA_preprocess.R 预处理。
# 运行：source("05d_moduleD_v2_SHAP_20260830.R")；结果输出 results/moduleD_v2/
# 产物：moduleD_v2_SHAP_values.csv     —— 117 样本 × 10 标志物全量 SHAP 值
#       moduleD_v2_SHAP_importance.csv —— mean|SHAP| 排序 + beta 方向 + 轨归属
#       moduleD_v2_SHAP_beeswarm.pdf   —— 蜜蜂图 + 重要性条形图（Fig 5 新 panel）
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
})
set.seed(20260830)

# ----------------------------- CONFIG --------------------------------------
CONFIG <- list(
  TRACK_FILE = "results/moduleC/moduleC_gene_tracks.csv",
  IN_DIR  = "processed",
  OUT_DIR = "results/moduleD_v2",
  LOG_FILE = file.path("results/moduleD_v2", sprintf("moduleD_v2_SHAP_log_%s.txt", Sys.Date())),
  TRAIN_SETS = c("GSE114007", "GSE55235", "GSE55457", "GSE12021", "GSE51588"),
  TRAIN_LABEL = "OA",
  VALIDATION_GSE = "GSE48556"   # 仅用于平台预筛（与 05b 完全同口径）
)
dir.create(CONFIG$OUT_DIR, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  txt <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = ""))
  cat(txt, "\n"); cat(txt, "\n", file = CONFIG$LOG_FILE, append = TRUE)
}

# ----------------------------- 0. 锁定产物读取 -------------------------------
bm_file <- file.path(CONFIG$OUT_DIR, "moduleD_v2_biomarkers.rds")
if (!file.exists(bm_file)) stop("缺少 ", bm_file, "——请先运行 05b_moduleD_v2_multicohort.R")
bm_pack <- readRDS(bm_file)
biomarkers <- bm_pack$biomarkers
if (length(biomarkers) < 2) stop("v2 锁定标志物不足 2 个——请回查 05b 产物")
log_msg("v2 锁定标志物（", length(biomarkers), " 个）：", paste(biomarkers, collapse = ", "))

ext_file <- file.path(CONFIG$OUT_DIR, "roc_external_blood.rds")
if (!file.exists(ext_file)) stop("缺少 ", ext_file, "——05b 外部验证未执行，无法做系数自检")

# ----------------------------- 1. 重建训练矩阵（复用 05b 装配链）------------
## —— 以下装配逻辑与 05b 逐行一致，任何改动必须双改并留痕 ——
track_tab <- read.csv(CONFIG$TRACK_FILE, stringsAsFactors = FALSE)
core_genes <- track_tab$gene

assemble_cohort <- function(gse, label) {
  expr  <- readRDS(file.path(CONFIG$IN_DIR, gse, paste0(gse, "_expr_gene.rds")))
  pheno <- read.csv(file.path(CONFIG$IN_DIR, gse, paste0(gse, "_pheno.csv")),
                    stringsAsFactors = FALSE)
  pheno <- pheno[match(colnames(expr), pheno$sample), ]
  if (anyNA(pheno$sample) || anyNA(pheno$condition))
    stop(gse, ": pheno 缺少 expr 中的样本或 condition 缺失——请回查模块 A 输出")
  if (!all(pheno$condition %in% c("control", label, "other")))
    stop(gse, ": condition 存在越界取值")
  keep <- pheno$condition %in% c("control", label)
  expr <- expr[, keep, drop = FALSE]; pheno <- pheno[keep, ]
  if ("donor_id" %in% names(pheno) && anyDuplicated(pheno$donor_id)) {
    if (anyNA(pheno$donor_id)) stop(gse, ": donor_id 含 NA，请人工核对")
    keep_idx <- !duplicated(pheno$donor_id)
    expr <- expr[, keep_idx, drop = FALSE]; pheno <- pheno[keep_idx, ]
    log_msg("  ", gse, ": donor_id 去重 -> ", ncol(expr), " 样本")
  }
  x <- scale(t(expr))
  y <- factor(ifelse(pheno$condition == label, "case", "control"),
              levels = c("control", "case"))
  list(x = x, y = y, gse = gse)
}

log_msg("—— 重建 ", length(CONFIG$TRAIN_SETS), " 个训练队列（05b 同口径）——")
cohorts <- lapply(CONFIG$TRAIN_SETS, assemble_cohort, label = CONFIG$TRAIN_LABEL)

val_genes_all <- rownames(readRDS(file.path(
  CONFIG$IN_DIR, CONFIG$VALIDATION_GSE, paste0(CONFIG$VALIDATION_GSE, "_expr_gene.rds"))))
feat <- Reduce(intersect, c(list(core_genes, val_genes_all),
                            lapply(cohorts, function(cc) colnames(cc$x))))
x_all <- do.call(rbind, lapply(cohorts, function(cc) cc$x[, feat, drop = FALSE]))
y_all <- factor(unlist(lapply(cohorts, function(cc) as.character(cc$y))),
                levels = c("control", "case"))
batch <- factor(rep(CONFIG$TRAIN_SETS, sapply(cohorts, function(cc) nrow(cc$x))))
if (anyNA(x_all)) stop("重建矩阵含 NaN——与 05b 实跑环境不一致，请回查")
if (!all(biomarkers %in% colnames(x_all)))
  stop("锁定标志物不在重建特征空间内：",
       paste(setdiff(biomarkers, colnames(x_all)), collapse = ", "))
log_msg("重建训练集：", nrow(x_all), " 样本 × ", ncol(x_all), " 特征（case=",
        sum(y_all == "case"), " / control=", sum(y_all == "control"), "）")

## 重建保真自检 1：队列组成 vs 05b 存盘
comp_saved <- read.csv(file.path(CONFIG$OUT_DIR, "cohort_composition.csv"),
                       stringsAsFactors = FALSE)
comp_now <- data.frame(cohort = CONFIG$TRAIN_SETS,
                       n = sapply(cohorts, function(cc) nrow(cc$x)),
                       case = sapply(cohorts, function(cc) sum(cc$y == "case")),
                       control = sapply(cohorts, function(cc) sum(cc$y == "control")))
if (!isTRUE(all.equal(comp_saved, comp_now, check.attributes = FALSE)))
  stop("队列组成与 05b 存盘不一致——装配漂移，请回查 processed/ 是否变动")
log_msg("✅ 自检1：队列组成与 05b 存盘逐格一致")

## 重建保真自检 2：重拟合系数 vs roc_external_blood.rds 封存 fit 逐位一致
x_bm <- x_all[, biomarkers, drop = FALSE]
fit_shap <- suppressWarnings(glm(y_all ~ .,
                                 data = data.frame(y_all = y_all, x_bm),
                                 family = binomial))
fit_saved <- readRDS(ext_file)$fit
if (!identical(names(coef(fit_saved)), names(coef(fit_shap))) ||
    !isTRUE(all.equal(as.numeric(coef(fit_saved)), as.numeric(coef(fit_shap)),
                      tolerance = 1e-8)))
  stop("重拟合系数与 05b 封存 fit 不一致（容差 1e-8）——装配漂移，请回查")
log_msg("✅ 自检2：重拟合 11 个系数（含截距）与 05b 封存 fit 逐位一致（容差 1e-8）")

# ----------------------------- 2. 精确线性 SHAP（link scale）----------------
beta <- coef(fit_shap)[-1]                       # 10 个标志物系数
baseline <- colMeans(x_bm)                       # 训练集均值作 baseline
SHAP <- sweep(x_bm, 2L, baseline, `-`)           # x_ij - mean(x_j)
SHAP <- sweep(SHAP, 2L, beta, `*`)               # phi_ij = beta_j * (x_ij - mean_j)

## 效率性硬校验：sum_j phi_ij + intercept == linear predictor_i
lp_pred <- as.numeric(predict(fit_shap, newdata = data.frame(x_bm)))
## as.numeric 去名：rowSums/系数名残留会使 all.equal 的 names 校验误报
lp_shap <- as.numeric(rowSums(SHAP)) + as.numeric(coef(fit_shap)[1])
if (!isTRUE(all.equal(lp_pred, lp_shap, tolerance = 1e-8)))
  stop("SHAP 效率性校验失败（max|diff| = ",
       format(max(abs(lp_pred - lp_shap)), digits = 3), "）——实现有误，请回查")
log_msg("✅ 自检3：SHAP 效率性成立（sum(phi)+intercept == 线性预测值，容差 1e-8；max|diff| = ",
        format(max(abs(lp_pred - lp_shap)), digits = 3), "）")

# ----------------------------- 3. 产物 1：全量 SHAP 值表 ---------------------
shap_out <- data.frame(sample = rownames(SHAP), cohort = batch,
                       label = y_all, SHAP, check.names = FALSE)
write.csv(shap_out, file.path(CONFIG$OUT_DIR, "moduleD_v2_SHAP_values.csv"),
          row.names = FALSE)
log_msg("产物1：moduleD_v2_SHAP_values.csv（", nrow(shap_out), " × ",
        ncol(shap_out) - 3, " 标志物 SHAP 值）")

# ----------------------------- 4. 产物 2：重要性表 ---------------------------
anno_idx <- match(biomarkers, track_tab$gene)
importance <- data.frame(
  gene = biomarkers,
  beta = as.numeric(beta),
  mean_abs_shap = colMeans(abs(SHAP)),
  mean_shap_case = colMeans(SHAP[y_all == "case", , drop = FALSE]),
  mean_shap_control = colMeans(SHAP[y_all == "control", , drop = FALSE]),
  track = track_tab$track[anno_idx],
  direction_OA = track_tab$direction_OA[anno_idx],
  stringsAsFactors = FALSE
)
importance <- importance[order(-importance$mean_abs_shap), ]
importance$rank <- seq_len(nrow(importance))
rownames(importance) <- NULL
write.csv(importance, file.path(CONFIG$OUT_DIR, "moduleD_v2_SHAP_importance.csv"),
          row.names = FALSE)
log_msg("产物2：moduleD_v2_SHAP_importance.csv（mean|SHAP| 排序）：")
for (i in seq_len(nrow(importance)))
  log_msg("  ", importance$rank[i], ". ", importance$gene[i],
          "  mean|SHAP|=", round(importance$mean_abs_shap[i], 4),
          "  beta=", round(importance$beta[i], 4),
          "  轨=", importance$track[i], "/OA ", importance$direction_OA[i])

## beta 方向 vs 发现层方向一致性留痕（beta>0 应≈OA 上调；z-score 口径下非严格等价，仅参考）
dir_align <- sum(sign(importance$beta) == ifelse(importance$direction_OA == "up", 1, -1))
log_msg("beta 方向与发现层 OA 方向一致：", dir_align, "/", nrow(importance),
        "（参考口径：多变量 logistic 系数含共线性调整，不一致不必然矛盾）")

# ----------------------------- 5. 产物 3：蜜蜂图 + 重要性条形图 --------------
long <- do.call(rbind, lapply(biomarkers, function(g) {
  data.frame(gene = g, shap = SHAP[, g], value = x_bm[, g],
             stringsAsFactors = FALSE)
}))
gene_order <- importance$gene[nrow(importance):1]   # 顶部为最重要
long$gene <- factor(long$gene, levels = gene_order)

p_bee <- ggplot(long, aes(x = shap, y = gene, color = value)) +
  geom_jitter(width = 0, height = 0.22, size = 1.1, alpha = 0.65) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  scale_color_gradient2(low = "#3B6FB5", mid = "grey80", high = "#C0392B",
                        midpoint = 0, name = "Feature value\n(z-score)") +
  labs(x = "SHAP value (contribution to log-odds of OA)", y = NULL,
       title = "SHAP summary of the 10-biomarker logistic model",
       subtitle = sprintf("Multi-cohort training set, n = %d; exact linear SHAP on link scale",
                          nrow(x_bm))) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        panel.grid.major.y = element_blank())

p_bar <- ggplot(importance,
                aes(x = mean_abs_shap, y = factor(gene, levels = gene_order),
                    fill = beta > 0)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "#3B6FB5"),
                    labels = c(`TRUE` = "beta > 0 (risk)", `FALSE` = "beta < 0 (protective)"),
                    name = "Coefficient") +
  labs(x = "mean(|SHAP value|)", y = NULL,
       title = "Global feature importance (mean |SHAP|)") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "bottom")

pdf(file.path(CONFIG$OUT_DIR, "moduleD_v2_SHAP_beeswarm.pdf"), width = 7.5, height = 5.5)
print(p_bee); print(p_bar)
dev.off()
log_msg("产物3：moduleD_v2_SHAP_beeswarm.pdf（2 页：蜜蜂图 + 重要性条形图）")

# ----------------------------- 收尾 ------------------------------------------
cat("\n模块 D v2 SHAP 补丁完成。产物：\n",
    "  ", CONFIG$OUT_DIR, "/moduleD_v2_SHAP_values.csv —— 全量 SHAP 值（117×10）\n",
    "  ", CONFIG$OUT_DIR, "/moduleD_v2_SHAP_importance.csv —— 重要性排序 + beta 方向\n",
    "  ", CONFIG$OUT_DIR, "/moduleD_v2_SHAP_beeswarm.pdf —— Fig 5 新 panel 素材\n",
    "写作纪律（备忘录 §50-1）：正文只报 SHAP 排序，不与四法投票/LASSO 系数序强行对齐；\n",
    "  Methods 注明 exact linear SHAP on the logit scale with training-mean baseline。\n",
    sep = "")
