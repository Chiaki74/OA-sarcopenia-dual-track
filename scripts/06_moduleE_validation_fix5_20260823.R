# ============================================================================
# 06_moduleE_validation_fix5_20260823.R
# 模块 E｜外部验证（2026-08-23，计划书 §4.5 + 备忘录 §八-27/28 分层口径）
# ----------------------------------------------------------------------------
# 分层验证框架（血液降级后的新口径）：
#   【主验证 · 软骨对软骨】GSE57218（7 正常 vs 33 OA 受累）/ GSE168505 /
#     GSE206848 —— AUC > 0.7 的验收标准仅适用于本层；
#     逐标志物方向一致性要求：≥2 个独立队列方向一致（§4.5）；
#   【探索性 · 血液】GSE48556 结果引自模块 D（v2：AUC=0.605），本脚本不重算；
#   【探索性 · 肌肉镜像方向】GSE136344/GSE117525/GSE144304（frailty/MetS 代理，
#     非 sarcopenia 确诊队列）——检验 10 个标志物的肌肉方向是否与
#     direction_muscle 一致（双轨框架的组织间外部证据）。
# 模型：在 v2 五队列训练集（逐队列 z-score 合并）上重训 glm（与 05b 同口径），
#       再应用于各验证队列（队列内 z-score）。
# ----------------------------------------------------------------------------
# 前置：results/moduleD_v2/moduleD_v2_biomarkers.rds（05b 已运行）；
#       results/moduleC/moduleC_gene_tracks.csv；
#       5 个训练集 + 3 个软骨验证集 + 3 个肌肉验证集均已预处理
# 运行：source("06_moduleE_validation.R")；输出 results/moduleE/
# ============================================================================

suppressPackageStartupMessages({
  library(pROC)
  library(data.table)
})
set.seed(20260823)

CONFIG <- list(
  TRACK_FILE = "results/moduleC/moduleC_gene_tracks.csv",
  V2_RDS     = "results/moduleD_v2/moduleD_v2_biomarkers.rds",
  IN_DIR  = "processed",
  OUT_DIR = "results/moduleE",
  LOG_FILE = file.path("results/moduleE", sprintf("moduleE_log_%s.txt", Sys.Date())),
  TRAIN_SETS  = c("GSE114007", "GSE55235", "GSE55457", "GSE12021", "GSE51588"),
  TRAIN_LABEL = "OA",
  ## 主验证（软骨）：case 标签候选按优先级排列，运行时实测匹配
  OA_VAL_SETS = c("GSE57218", "GSE168505", "GSE206848"),
  OA_CASE_CANDIDATES = list(
    GSE57218  = c("OA_affected", "OA"),   # RAAK：受累软骨 vs 健康；OA_preserved 不进主分析
    GSE168505 = c("OA"),
    GSE206848 = c("OA")
  ),
  ## 探索性肌肉镜像验证（frailty/MetS 代理，非 sarcopenia 确诊）
  MUSCLE_VAL_SETS = c("GSE136344", "GSE117525", "GSE144304"),
  MUSCLE_CASE_CANDIDATES = list(
    GSE136344 = c("metabolic_syndrome", "frailty", "frail"),
    GSE117525 = c("frailty", "frail"),   # 2026-08-23 实测修正：标签为 frailty 非 frail
    GSE144304 = c("frailty", "frail")    # 同上
  ),
  MUSCLE_CTRL_CANDIDATES = list(
    GSE136344 = c("aging", "control", "young"),  # 实测：aging=12 为 MetS 的老年对照
    GSE117525 = c("control", "healthy", "non_frail", "nonfrail"),
    GSE144304 = c("control", "healthy", "non_frail", "nonfrail")
  ),
  MIN_GENE_COVER = 0.7,   # 验证队列至少覆盖 70% 标志物才执行模型验证
  N_BOOT = 2000
)
dir.create(CONFIG$OUT_DIR, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  txt <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = ""))
  cat(txt, "\n"); cat(txt, "\n", file = CONFIG$LOG_FILE, append = TRUE)
}

# ---- 通用装配（与 05b 同口径：pheno 硬校验 + donor 去重 + 队列内 z-score）----
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
  ## 2026-08-23 修正v2：零方差基因（z-score NaN）剔除而非硬停——
  ## 常数基因本就不携带方向/区分信息（GSE136344/GSE144304 实机触发）
  if (anyNA(x)) {
    bad <- colnames(x)[apply(x, 2, anyNA)]
    log_msg("  ", gse, ": 剔除 ", length(bad), " 个零方差基因（z-score NaN）：",
            paste(head(bad, 5), collapse = ", "), if (length(bad) > 5) " ..." else "")
    x <- x[, !colnames(x) %in% bad, drop = FALSE]
  }
  y <- factor(ifelse(pheno$condition == case_label, "case", "control"),
              levels = c("control", "case"))
  if (sum(y == "case") == 0 || sum(y == "control") == 0)
    stop(gse, ": 病例/对照有一侧为 0")
  list(x = x, y = y)
}

## case/control 标签实测匹配（候选制，防口径漂移静默失败）
resolve_case <- function(gse, candidates, what = "case") {
  pheno <- read.csv(file.path(CONFIG$IN_DIR, gse, paste0(gse, "_pheno.csv")),
                    stringsAsFactors = FALSE)
  avail <- unique(pheno$condition)
  hit <- intersect(candidates, avail)
  if (length(hit) == 0)
    stop(gse, ": ", what, " 候选（", paste(candidates, collapse = "/"),
         "）均不在 condition 实测取值（", paste(avail, collapse = ","), "）中——请人工核对")
  if (length(hit) > 1 || hit[1] != candidates[1])
    log_msg("  ", gse, ": ", what, " 标签实测匹配为 '", hit[1], "'（候选优先级自动选择）")
  hit[1]
}

# ---- 读取 v2 标志物与双轨表 ----
v2 <- readRDS(CONFIG$V2_RDS)
biomarkers <- v2$biomarkers
track_tab <- read.csv(CONFIG$TRACK_FILE, stringsAsFactors = FALSE)
## 审查 S2：标志物必须全部在双轨表中，否则 direction 注释静默变 NA 污染验收表
if (!all(biomarkers %in% track_tab$gene))
  stop("v2 标志物不全在 moduleC_gene_tracks.csv 中：",
       paste(setdiff(biomarkers, track_tab$gene), collapse = ", "), "——请回查上游")
log_msg("v2 标志物（", length(biomarkers), " 个）：", paste(biomarkers, collapse = ", "))

# ---- 1. 重训最终模型（五队列合并，口径同 05b）----
log_msg("—— 重训最终模型（", length(CONFIG$TRAIN_SETS), " 训练队列合并）——")
tr <- lapply(CONFIG$TRAIN_SETS, assemble_cohort, case_label = CONFIG$TRAIN_LABEL)
feat_tr <- Reduce(intersect, c(list(biomarkers), lapply(tr, function(cc) colnames(cc$x))))
if (length(feat_tr) < length(biomarkers))
  log_msg("⚠ ", length(biomarkers) - length(feat_tr), " 个标志物在部分训练队列缺失，用 ",
          length(feat_tr), " 个交集特征重训")
x_all <- do.call(rbind, lapply(tr, function(cc) cc$x[, feat_tr, drop = FALSE]))
y_all <- factor(unlist(lapply(tr, function(cc) as.character(cc$y))),
                levels = c("control", "case"))
log_msg("训练集：", nrow(x_all), " 样本（case=", sum(y_all == "case"), "）")
## 注：最终模型在各验证队列按其可检标志物子集逐队列重训（审查 F1 修复口径），
## 此处不再拟合单一全特征模型

# ---- 2. 主验证：软骨三队列 ----
log_msg("—— 主验证（软骨对软骨，AUC>0.7 标准适用本层）——")
auc_report <- data.frame()
dir_long <- data.frame()
for (gse in CONFIG$OA_VAL_SETS) {
  ## 审查 S1：软骨层 resolve_case 同样降级跳过（与肌肉层一致），防单队列标签漂移终止全部主验证
  case_lab <- tryCatch(resolve_case(gse, CONFIG$OA_CASE_CANDIDATES[[gse]]),
                       error = function(e) { log_msg("⚠ ", gse, "：", conditionMessage(e)); NULL })
  if (is.null(case_lab)) next
  vv <- tryCatch(assemble_cohort(gse, case_lab),
                 error = function(e) { log_msg("⚠ ", gse, " 装配失败：", conditionMessage(e)); NULL })
  if (is.null(vv)) next
  common <- intersect(feat_tr, colnames(vv$x))
  log_msg(gse, "（case='", case_lab, "'）：", nrow(vv$x), " 样本（case=",
          sum(vv$y == "case"), "），标志物覆盖 ", length(common), "/", length(feat_tr))
  ## 逐基因方向核对（无论模型是否执行都做）
  for (g in common) {
    obs <- ifelse(mean(vv$x[vv$y == "case", g]) > mean(vv$x[vv$y == "control", g]), "up", "down")
    idx <- match(g, track_tab$gene)
    dir_long <- rbind(dir_long, data.frame(
      cohort = gse, gene = g, track = track_tab$track[idx],
      direction_OA_discovery = track_tab$direction_OA[idx],
      direction_observed = obs, consistent = obs == track_tab$direction_OA[idx]))
  }
  ## 模型 AUC（覆盖率达标才执行）
  if (length(common) >= max(2, ceiling(length(feat_tr) * CONFIG$MIN_GENE_COVER))) {
    ## 审查 F1（致命）：predict 无法处理 newdata 缺列——覆盖率为 70%~100% 时
    ## 必须按 common 逐队列重训，与"交集特征重训"口径一致
    fit_v <- suppressWarnings(glm(y ~ ., family = binomial,
                                  data = data.frame(y = y_all, x_all[, common, drop = FALSE])))
    if (!fit_v$converged) log_msg("  ⚠ ", gse, " 重训 glm 不收敛（疑似分离），AUC 口径留痕")
    prob <- as.numeric(predict(fit_v, newdata = data.frame(vv$x[, common, drop = FALSE]),
                               type = "response"))
    if (anyNA(prob)) { log_msg("⚠ ", gse, " 预测含 NA（完全分离），跳过 AUC"); next }
    roc_v <- roc(vv$y, prob, levels = c("control", "case"), direction = "<", quiet = TRUE)
    ## 审查建议 4：小样本队列 bootstrap 偶发报错，降级为只报点估计
    ci_v <- tryCatch(ci.auc(roc_v, method = "bootstrap", boot.n = CONFIG$N_BOOT),
                     error = function(e) { log_msg("  ⚠ ", gse, " bootstrap CI 失败，仅报点估计")
                                           c(NA, as.numeric(auc(roc_v)), NA) })
    auc_report <- rbind(auc_report, data.frame(
      cohort = gse, n_case = sum(vv$y == "case"), n_control = sum(vv$y == "control"),
      genes_used = length(common), AUC = round(as.numeric(auc(roc_v)), 3),
      CI_low = round(ci_v[1], 3), CI_high = round(ci_v[3], 3)))
    log_msg("  ", gse, " 模型 AUC = ", round(as.numeric(auc(roc_v)), 3),
            " (", round(ci_v[1], 3), "–", round(ci_v[3], 3), ")",
            ifelse(as.numeric(auc(roc_v)) > 0.7, " ✓ 达标", " ✗ 未达 0.7"))
  } else log_msg("  ", gse, " 标志物覆盖不足 70%，模型验证不执行（仅方向核对）")
}
## 审查 S3：空表守卫——三队列全失败时显式报错，而非 aggregate 误导性崩溃
if (nrow(dir_long) == 0)
  stop("主验证三软骨队列均装配失败，无验收依据——请回查 processed/ 与 condition 口径")
if (nrow(auc_report) > 0) {
  write.csv(auc_report, file.path(CONFIG$OUT_DIR, "moduleE_cartilage_AUC.csv"), row.names = FALSE)
} else {
  log_msg("⚠ 无队列达成覆盖率门槛，moduleE_cartilage_AUC.csv 不生成（仅方向核对产物）")
}
write.csv(dir_long, file.path(CONFIG$OUT_DIR, "moduleE_cartilage_direction_by_gene.csv"),
          row.names = FALSE)

## §4.5 逐标志物验收：≥2 个独立队列方向一致
verdict <- aggregate(consistent ~ gene, data = dir_long,
                     FUN = function(v) c(n = length(v), k = sum(v)))
verdict$n_cohorts <- verdict$consistent[, "n"]
verdict$n_consistent <- verdict$consistent[, "k"]
verdict$consistent <- NULL
verdict$pass_direction <- verdict$n_consistent >= 2
idx <- match(verdict$gene, track_tab$gene)
verdict$track <- track_tab$track[idx]
write.csv(verdict, file.path(CONFIG$OUT_DIR, "moduleE_biomarker_verdict.csv"), row.names = FALSE)
log_msg("§4.5 方向验收（≥2 队列一致）：", sum(verdict$pass_direction), "/",
        nrow(verdict), " 个标志物通过")

# ---- 3. 探索性：肌肉镜像方向验证 ----
log_msg("—— 探索性：肌肉镜像方向验证（frailty/MetS 代理队列，非 sarcopenia 确诊）——")
mir_long <- data.frame()
for (gse in CONFIG$MUSCLE_VAL_SETS) {
  case_lab <- tryCatch(resolve_case(gse, CONFIG$MUSCLE_CASE_CANDIDATES[[gse]]),
                       error = function(e) { log_msg("⚠ ", gse, "：", conditionMessage(e)); NULL })
  if (is.null(case_lab)) next
  ctrl_lab <- tryCatch(resolve_case(gse, CONFIG$MUSCLE_CTRL_CANDIDATES[[gse]], what = "control"),
                       error = function(e) { log_msg("⚠ ", gse, "：", conditionMessage(e)); NULL })
  if (is.null(ctrl_lab)) next
  vv <- tryCatch(assemble_cohort(gse, case_lab, ctrl_lab),
                 error = function(e) { log_msg("⚠ ", gse, " 装配失败：", conditionMessage(e)); NULL })
  if (is.null(vv)) next
  common <- intersect(feat_tr, colnames(vv$x))
  for (g in common) {
    obs <- ifelse(mean(vv$x[vv$y == "case", g]) > mean(vv$x[vv$y == "control", g]), "up", "down")
    idx <- match(g, track_tab$gene)
    mir_long <- rbind(mir_long, data.frame(
      cohort = gse, gene = g, track = track_tab$track[idx],
      direction_muscle_expected = track_tab$direction_muscle[idx],
      direction_observed = obs,
      consistent = obs == track_tab$direction_muscle[idx]))
  }
  sub <- mir_long[mir_long$cohort == gse, ]
  log_msg("  ", gse, "（case='", case_lab, "'，", nrow(vv$x), " 样本）：肌肉方向一致 ",
          sum(sub$consistent), "/", nrow(sub))
}
if (nrow(mir_long) > 0) {
  write.csv(mir_long, file.path(CONFIG$OUT_DIR, "moduleE_muscle_mirror_direction.csv"),
            row.names = FALSE)
  ## 按轨汇总镜像一致性（镜像轨基因的反向模式是否在外部肌肉队列复现）
  for (tr_name in unique(mir_long$track)) {
    sub <- mir_long[mir_long$track == tr_name, ]
    log_msg("  轨 ", tr_name, "：肌肉方向一致率 ",
            round(100 * mean(sub$consistent), 1), "%（", sum(sub$consistent), "/",
            nrow(sub), " 基因×队列）")
  }
}

cat("\n模块 E 完成（分层验证）。产物：\n",
    "  ", CONFIG$OUT_DIR, "/moduleE_cartilage_AUC.csv —— 【主验证】软骨三队列模型 AUC\n",
    "  ", CONFIG$OUT_DIR, "/moduleE_cartilage_direction_by_gene.csv —— 逐基因×队列方向核对\n",
    "  ", CONFIG$OUT_DIR, "/moduleE_biomarker_verdict.csv —— §4.5 逐标志物验收结论\n",
    "  ", CONFIG$OUT_DIR, "/moduleE_muscle_mirror_direction.csv —— 【探索性】肌肉镜像方向验证\n",
    "判读：主验证看软骨 AUC（>0.7）与 verdict 表；肌肉层按轨看一致率（镜像轨应呈反向复现）。\n",
    sep = "")
