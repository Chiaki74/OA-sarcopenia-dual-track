# ============================================================================
# 40_Table1_十标志物三档分层诊断效能_组装_20261010.R
# 版本：v1.0　日期：2026-10-10　编码：UTF-8
# ----------------------------------------------------------------------------
# 任务：组装主文 Table 1「10 标志物三档分层诊断效能汇总」（英文投稿版表头）。
#   逐标志物行 + PANEL 汇总行 + GSE206848 注释行。
# 数据源（全部为已封版/已归档产物，本脚本只读不写）：
#   ① 补充表_封版/表S11_单基因诊断效能_方向固定订正版_20260916.csv
#      → 训练层单基因 AUC、GSE57218 方向固定 AUC + DeLong 95%CI、方向一致性
#   ② results/moduleC/moduleC_gene_tracks.csv
#      → 轨道归属（concordant / mirror）与发现层双病方向
#   ③ results/moduleE/moduleE_cartilage_direction_by_gene.csv
#      → GSE168505 等软骨队列逐基因方向核对
#   ④ results/moduleE/moduleE_biomarker_verdict.csv
#      → 软骨层方向验收（≥2 独立队列一致）
#   ⑤ results/moduleD_v2/direction_check_GSE48556.csv
#      → 血液探索层逐基因方向一致性
#   ⑥ results/moduleE/moduleE_muscle_mirror_direction.csv
#      → 肌肉代理探索层逐基因逐队列方向一致性
#   ⑦ results/moduleE/moduleE_cartilage_AUC.csv（可选）
#      → PANEL 行软骨层 AUC；缺失时用封版常数并打印警告
# 产出：
#   results/moduleE/Table1_10markers_three_tier_20261010.csv（英文表头，投稿用）
#   + 控制台预览 + 封版数字自检（与正文口径不一致时显式 WARNING）
# 口径锁定（与中文终稿 v16 / 英文版一致，脚本内自检用）：
#   GSE57218 方向一致 6/10；软骨层方向验收 7/10；血液 6/10；
#   肌肉代理总体 9/20；GSE136344 交集 0/0；GSE206848 覆盖 0/10 未执行。
# ============================================================================

message("== 脚本版本：40_Table1_组装 v1.0（2026-10-10）==")

## ---- 0. 环境 ----
ROOT <- "C:/Users/凉月/Desktop/ENG论文"
if (dir.exists(ROOT)) setwd(ROOT) else message("⚠ 默认根目录不存在，使用当前工作目录：", getwd())

options(stringsAsFactors = FALSE)
BIOMARKERS <- c("NELL1","STEAP1","GADD45A","FBLN5","HEMK1","GDE1","RNF14","BNIP3","EGFR","CBR3")

find_file <- function(basename, roots = c(".", "results", "补充表_封版")) {
  ## 按文件名递归查找（与脚本 39 的 sealed_names 策略一致），返回首个命中路径或 NA
  hits <- character(0)
  for (r in roots) {
    if (dir.exists(r))
      hits <- c(hits, list.files(r, pattern = paste0("^", gsub("\\.", "\\\\.", basename), "$"),
                                 recursive = TRUE, full.names = TRUE))
  }
  if (length(hits) == 0) return(NA_character_)
  if (length(hits) > 1) message("  ⚠ ", basename, " 命中多处，取第一个：", hits[1])
  hits[1]
}

read_or_na <- function(path, label) {
  if (is.na(path) || !file.exists(path)) {
    message("  ❌ ", label, " 未找到 —— 对应列标记「待补」")
    return(NULL)
  }
  message("  ✅ ", label, " <- ", path)
  read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM")
}

## ---- 1. 读取六个数据源 ----
message("—— 数据源定位 ——")
s11   <- read_or_na(find_file("表S11_单基因诊断效能_方向固定订正版_20260916.csv"), "① 表S11（封版）")
tracks<- read_or_na(find_file("moduleC_gene_tracks.csv"), "② moduleC 轨道表")
dirE  <- read_or_na(find_file("moduleE_cartilage_direction_by_gene.csv"), "③ moduleE 软骨逐基因方向")
verd  <- read_or_na(find_file("moduleE_biomarker_verdict.csv"), "④ moduleE 方向验收")
blood <- read_or_na(find_file("direction_check_GSE48556.csv"), "⑤ 血液方向核对")
mir   <- read_or_na(find_file("moduleE_muscle_mirror_direction.csv"), "⑥ 肌肉代理方向")
aucE  <- read_or_na(find_file("moduleE_cartilage_AUC.csv"), "⑦ moduleE 软骨 panel AUC（可选）")

if (is.null(s11)) stop("表 S11 封版 csv 为硬依赖，未找到——请先确认归档位置（补充表_封版/）")

## ---- 2. 逐标志物行 ----
rows <- lapply(BIOMARKERS, function(g) {
  s <- s11[s11$gene == g, ]
  if (nrow(s) != 1) stop("S11 中 ", g, " 行数异常：", nrow(s))

  ## 轨道与发现层方向
  track <- disc_oa <- disc_muscle <- "待补"
  if (!is.null(tracks)) {
    tr <- tracks[tracks$gene == g, ]
    if (nrow(tr) == 1) {
      track <- as.character(tr$track)
      if ("direction_OA" %in% names(tr))     disc_oa     <- as.character(tr$direction_OA)
      if ("direction_muscle" %in% names(tr)) disc_muscle <- as.character(tr$direction_muscle)
    }
  }

  ## GSE57218 单基因 AUC（方向固定）+ CI
  g57218 <- sprintf("%.3f (%.3f–%.3f)", s$GSE57218_AUC_方向固定, s$CI_low_DeLong, s$CI_high_DeLong)

  ## GSE168505 方向（moduleE 软骨逐基因表）
  g168505 <- "待补"
  if (!is.null(dirE)) {
    d <- dirE[dirE$cohort == "GSE168505" & dirE$gene == g, ]
    if (nrow(d) == 1) g168505 <- ifelse(d$consistent, "consistent", "reversed")
  }

  ## 软骨层方向验收（≥2 队列）
  accept <- "待补"
  if (!is.null(verd)) {
    v <- verd[verd$gene == g, ]
    if (nrow(v) == 1) accept <- sprintf("%d/%d%s", v$n_consistent, v$n_cohorts,
                                        ifelse(v$pass_direction, " ✓", ""))
  }

  ## 血液 GSE48556 方向
  bld <- "待补"
  if (!is.null(blood)) {
    b <- blood[blood$gene == g, ]
    if (nrow(b) == 1) bld <- ifelse(b$consistent, "consistent", "reversed")
  }

  ## 肌肉代理逐队列一致数
  mus <- "待补"
  if (!is.null(mir)) {
    m <- mir[mir$gene == g, ]
    if (nrow(m) > 0) mus <- sprintf("%d/%d", sum(m$consistent), nrow(m))
    else mus <- "0/0"
  }

  data.frame(
    Gene = g,
    Track = track,
    `Discovery direction (OA / muscle)` = paste0(disc_oa, " / ", disc_muscle),
    `Training AUC (5-cohort)` = sprintf("%.3f", s$train_AUC_五队列合并),
    `GSE57218 AUC (95% CI), direction-fixed` = g57218,
    `GSE168505 direction` = g168505,
    `Cartilage direction acceptance` = accept,
    `Blood GSE48556 direction` = bld,
    `Muscle-proxy concordance` = mus,
    check.names = FALSE)
})
tab <- do.call(rbind, rows)

## ---- 3. PANEL 汇总行与 GSE206848 注释行 ----
panel_57218  <- "0.818 (0.614–1.000)"   # 封版口径（DeLong）；⑦存在时以文件复核
panel_168505 <- "0.750"
if (!is.null(aucE)) {
  a5 <- aucE[aucE$cohort == "GSE57218", ]
  if (nrow(a5) == 1) {
    file_val <- sprintf("%.3f (%.3f–%.3f)", a5$AUC, a5$CI_low, a5$CI_high)
    if (file_val != panel_57218)
      message("⚠ ⑦文件值 ", file_val, " 与封版口径 ", panel_57218, " 不一致——以正文封版口径为准并留痕")
  }
  a16 <- aucE[aucE$cohort == "GSE168505", ]
  if (nrow(a16) == 1 && sprintf("%.3f", a16$AUC) != panel_168505)
    message("⚠ ⑦GSE168505 文件值 ", a16$AUC, " 与封版口径 0.750 不一致——以封版为准并留痕")
}
panel_row <- data.frame(
  Gene = "PANEL (10 markers, logistic)",
  Track = "—",
  `Discovery direction (OA / muscle)` = "—",
  `Training AUC (5-cohort)` = "0.959 (internal CV, reference only)",
  `GSE57218 AUC (95% CI), direction-fixed` = panel_57218,
  `GSE168505 direction` = paste0("AUC ", panel_168505, " (n=7, down-weighted)"),
  `Cartilage direction acceptance` = "7/10 pass (≥2 cohorts)",
  `Blood GSE48556 direction` = "AUC 0.605 (0.499–0.713); 6/10 consistent",
  `Muscle-proxy concordance` = "9/20 gene×cohort (45%)",
  check.names = FALSE)
note_row <- data.frame(
  Gene = "GSE206848",
  Track = "—",
  `Discovery direction (OA / muscle)` = "—",
  `Training AUC (5-cohort)` = "—",
  `GSE57218 AUC (95% CI), direction-fixed` = "not evaluated (marker coverage 0/10)",
  `GSE168505 direction` = "—",
  `Cartilage direction acceptance` = "—",
  `Blood GSE48556 direction` = "—",
  `Muscle-proxy concordance` = "—",
  check.names = FALSE)
tab_full <- rbind(tab, panel_row, note_row)

## ---- 4. 封版数字自检 ----
message("—— 封版口径自检 ——")
n_57218_ok <- sum(s11$direction_consistent == "yes")
cat(sprintf("GSE57218 方向一致：%d/10（封版口径 6/10）%s\n", n_57218_ok,
            ifelse(n_57218_ok == 6, " ✓", " ✗ WARNING")))
if (!is.null(verd)) {
  n_accept <- sum(verd$pass_direction)
  cat(sprintf("软骨层方向验收通过：%d/10（封版口径 7/10）%s\n", n_accept,
              ifelse(n_accept == 7, " ✓", " ✗ WARNING")))
}
if (!is.null(blood)) {
  n_blood <- sum(blood$consistent)
  cat(sprintf("血液方向一致：%d/10（封版口径 6/10）%s\n", n_blood,
              ifelse(n_blood == 6, " ✓", " ✗ WARNING")))
}
if (!is.null(mir)) {
  cat(sprintf("肌肉代理总体：%d/%d（封版口径 9/20；GSE136344 交集 0/0 不贡献）%s\n",
              sum(mir$consistent), nrow(mir),
              ifelse(sum(mir$consistent) == 9, " ✓", " ✗ WARNING（若总对子数≠20 请核对 GSE136344 是否被正确排除）")))
}

## ---- 5. 写出与预览 ----
dir.create("results/moduleE", recursive = TRUE, showWarnings = FALSE)
out_csv <- "results/moduleE/Table1_10markers_three_tier_20261010.csv"
write.csv(tab_full, out_csv, row.names = FALSE, fileEncoding = "UTF-8")
cat("\n== 脚本版本：40_Table1_组装 v1.0 完成 ==\n产出：", out_csv, "\n\n")
print(tab_full, row.names = FALSE)
cat("\n注：表体为英文投稿表头；任何「待补」列说明对应数据源未在常规位置找到，",
    "请把该文件路径反馈给我，我出补丁版（41）补列。\n", sep = "")
