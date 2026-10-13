# ============================================================================
# 41_Table1_补丁_血液列修正_20261010.R
# 版本：v1.0　日期：2026-10-10　编码：UTF-8
# ----------------------------------------------------------------------------
# 问题（40 首跑控制台）：血液列取自旧版 results/moduleD/direction_check_GSE48556.csv
#   （v1 单队列训练时代产物，仅 3/10 基因在册），正确文件为 moduleD_v2 同名文件
#   （10 基因版，封版口径 6/10 方向一致）。
# 本补丁：①显式锁定 moduleD_v2 路径；②校验该文件含全部 10 个标志物且一致数=6；
#   ③重写 Table1 csv 的血液列并复打预览。
# 前置：40 已生成 results/moduleE/Table1_10markers_three_tier_20261010.csv
# ============================================================================

message("== 脚本版本：41_Table1_补丁_血液列修正 v1.0（2026-10-10）==")

ROOT <- "C:/Users/凉月/Desktop/ENG论文"
if (dir.exists(ROOT)) setwd(ROOT) else message("⚠ 默认根目录不存在，使用当前工作目录：", getwd())

BIOMARKERS <- c("NELL1","STEAP1","GADD45A","FBLN5","HEMK1","GDE1","RNF14","BNIP3","EGFR","CBR3")
tab_csv  <- "results/moduleE/Table1_10markers_three_tier_20261010.csv"
blood_ok <- "results/moduleD_v2/direction_check_GSE48556.csv"

## ---- 1. 路径锁定与完整性校验 ----
if (!file.exists(blood_ok)) {
  ## 兜底：全盘找同名文件，但只接受含全部 10 标志物的那份
  cands <- list.files(".", pattern = "^direction_check_GSE48556\\.csv$",
                      recursive = TRUE, full.names = TRUE)
  cands <- cands[order(!grepl("moduleD_v2", cands))]  # moduleD_v2 优先
  blood_ok <- NA_character_
  for (f in cands) {
    d <- read.csv(f, check.names = FALSE, fileEncoding = "UTF-8-BOM")
    if (all(BIOMARKERS %in% d$gene)) { blood_ok <- f; break }
  }
  if (is.na(blood_ok)) stop("找不到含全部 10 个标志物的 direction_check_GSE48556.csv——请反馈实际路径")
}
blood <- read.csv(blood_ok, check.names = FALSE, fileEncoding = "UTF-8-BOM")
message("血液方向文件锁定：", blood_ok)
stopifnot(all(BIOMARKERS %in% blood$gene))
n_blood <- sum(blood$consistent[blood$gene %in% BIOMARKERS])
cat(sprintf("血液方向一致：%d/10（封版口径 6/10）%s\n", n_blood,
            ifelse(n_blood == 6, " ✓", " ✗ WARNING——文件内容仍与封版不符，请反馈我核查")))
if (n_blood != 6) stop("自检未过，不写表")

## ---- 2. 重写血液列 ----
tab <- read.csv(tab_csv, check.names = FALSE, fileEncoding = "UTF-8-BOM")
blood_col <- "Blood GSE48556 direction"
stopifnot(blood_col %in% names(tab))
for (g in BIOMARKERS) {
  b <- blood[blood$gene == g, ]
  tab[tab$Gene == g, blood_col] <- ifelse(b$consistent, "consistent", "reversed")
}

## ---- 3. 写回与预览 ----
write.csv(tab, tab_csv, row.names = FALSE, fileEncoding = "UTF-8")
cat("\n== 41 补丁完成，血液列已修正 ==\n产出：", tab_csv, "\n\n")
print(tab[, c("Gene", "Track", "GSE57218 AUC (95% CI), direction-fixed",
              "Cartilage direction acceptance", blood_col, "Muscle-proxy concordance")],
      row.names = FALSE)
