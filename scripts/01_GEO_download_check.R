################################################################################
## OA × 肌少症整合生信研究 — 第 1 周任务③：GEO 批量下载 + 分组核对脚本
## 文件名：01_GEO_download_check.R
## 版本：v1.0（2026-08-18）  配套：数据集信息核对表_预填版_20260813.csv / 部署操作指引 Step 3
##
## 功能：
##   1) 批量下载 21 个 GEO 数据集（18 bulk + 3 单细胞/单核；GSE57218 已替换 GSE129147，
##      后者仅留备用不下载）
##   2) 逐数据集提取样本元数据（pData）存为 CSV，供人工分组复核
##   3) 样本数与核对表实测口径自动比对，输出 download_check_report.csv
##   4) GSE167186 触发"93 vs 89 样本数差异"专项核对提醒（备忘录 §八-7）
##   5) 断点续跑：已完成的数据集自动跳过
##
## 用法（Windows / RStudio，先完成 setup_R_environment.R）：
##   方式一：RStudio 中 Source 本文件；
##   方式二：Rscript 01_GEO_download_check.R
##
## 耗时与磁盘预估：
##   - bulk 芯片/RNA-seq 矩阵：单个几十 MB ~ 500 MB，合计约 2–5 GB，约 0.5–2 小时；
##   - 单细胞/单核原始补充文件（GSE152805、GSE167186-sn）可达数 GB，默认【不下载】，
##     待第 5 周模块 G 启动前将下方 DOWNLOAD_SUPP 改为 TRUE 再补下；
##   - 建议工作盘预留 ≥ 30 GB。
################################################################################

## ============================================================================
## 0. 全局配置
## ============================================================================
DATA_DIR       <- file.path("data", "GEO")          # 数据存放目录（自动创建）
DOWNLOAD_SUPP  <- FALSE                             # TRUE = 补下 GSE152805/GSE167186/GSE169454 的单细胞补充大文件
SUPP_GSES      <- c("GSE152805", "GSE167186", "GSE169454")   # 需要补充大文件的单细胞/单核数据集
REPORT_FILE    <- paste0("download_check_report_", format(Sys.Date(), "%Y%m%d"), ".csv")
LOG_FILE       <- paste0("download_log_", format(Sys.Date(), "%Y%m%d"), ".log")
## 人工复核时请对照《数据集信息核对表_预填版_20260813.csv》（建议放在同一工作目录）

options(timeout = 1800)                             # 大文件下载放宽到 30 分钟
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(GEOquery); library(data.table)
})

msg <- function(...) {
  line <- sprintf(...)
  cat(line, "\n")
  cat(line, "\n", file = LOG_FILE, append = TRUE)
}

## ============================================================================
## 1. 数据集清单（内嵌，口径与核对表 CSV 完全一致，2026-08-18 版）
##    expected_total = 系列矩阵预期样本总数；NA = 需人工复核；
##    n_eset = 预期 eSet 个数（>1 表示多平台/多子系列超系列）
## ============================================================================
manifest <- data.frame(
  gse = c("GSE114007", "GSE55235", "GSE55457", "GSE12021", "GSE51588",
          "GSE169077", "GSE168505", "GSE57218", "GSE206848", "GSE48556",
          "GSE1428", "GSE25941", "GSE111016", "GSE111010", "GSE111006",
          "GSE136344", "GSE117525", "GSE144304",
          "GSE152805", "GSE169454", "GSE167186"),
  side = c(rep("OA", 10), rep("肌肉", 8), rep("单细胞/单核", 3)),
  tissue = c("软骨", "滑膜", "滑膜", "滑膜", "软骨下骨",
             "软骨", "软骨", "软骨", "滑膜", "外周血PBMC",
             "股外侧肌", "股外侧肌", "骨骼肌", "骨骼肌", "骨骼肌",
             "股外侧肌", "股外侧肌", "股外侧肌",
             "OA软骨+滑膜", "软骨细胞", "股外侧肌"),
  platform = c("GPL11154+GPL18573", "GPL96", "GPL96", "GPL96+GPL97", "GPL13497",
               "GPL96", "GPL16791", "GPL6947", "GPL570", "GPL6947",
               "GPL96", "GPL570", "GPL16791", "GPL16791", "GPL16791",
               "GPL5175", "GPL20880", "GPL18573",
               "GPL20301", "GPL16791", "GPL20301+GPL24676"),
  expected_total = c(38, 30, 33, NA, 50,
                     11, 7, 73, 16, 139,
                     22, 36, 40, 39, 40,
                     NA, 259, 80,
                     9, 7, 93),
  n_eset = c(2, 1, 1, 2, 1,
             1, 1, 1, 1, 1,
             1, 1, 1, 1, 1,
             1, 1, 1,
             1, 1, 2),
  role = c("训练集(主)", "训练集", "训练集", "训练集(仅RRA)", "训练集",
           "验证集(方向验证)", "验证集", "验证集(替换GSE129147)", "验证集", "验证集(血液卖点)",
           "训练集", "训练集", "训练集(表型锚点)", "训练/验证", "训练/验证",
           "验证集", "验证集", "验证集",
           "单细胞定位(主)", "单细胞验证", "单核定位(主)"),
  check_note = c("两平台批次，预处理需校正",
                 "仅取正常+OA（另有10例RA）",
                 "仅取正常+OA（另有13例RA）",
                 "GPL96口径9正常/10OA；与GSE55457供体重叠，仅训练侧RRA",
                 "50样本=25供体×2部位；统计按供体",
                 "5个体/池混合样本，仅方向验证",
                 "样本量小(3/4)",
                 "RAAK队列：7健康+33OA受累+33OA保留；主分析取7vs33非配对",
                 "仅取正常+OA（另有2例RA）",
                 "全女性GARP队列，同胞对家系相关，统计需考虑",
                 "全男性；衰老+轻度肌少症",
                 "纯衰老代理表型（勿标肌少症）",
                 "确诊肌少症20/20，表型锚点",
                 "肌少症仅9例；与GSE111016同系列注意批次",
                 "肌少症仅4例，统计效力弱",
                 "三组：年轻/健康老年/老年+代谢综合征，全男性",
                 "总259含随访；基线样本作主分析",
                 "另含26例年轻人",
                 "scRNA-seq，约3.7万细胞",
                 "scRNA-seq 3正常/4OA",
                 "★GEO挂93样本vs论文口径89，须逐样本对清单（备忘录§八-7）"),
  stringsAsFactors = FALSE
)

## GSE129147（备用，默认不下载；其膝内配对设计留作补充材料方向一致性证据）
## 如需下载，手动运行：GEOquery::getGEO("GSE129147", destdir = file.path(DATA_DIR, "GSE129147"))

## ============================================================================
## 2. 单数据集处理函数
## ============================================================================
process_one <- function(row) {
  gse <- row$gse
  gse_dir <- file.path(DATA_DIR, gse)
  done_marker <- file.path(gse_dir, ".done")
  supp_marker <- file.path(gse_dir, ".done_supp")
  need_supp   <- DOWNLOAD_SUPP && gse %in% SUPP_GSES   # 仅单细胞数据集且开关打开时才下大文件

  if (file.exists(done_marker)) {
    msg("  ⏭ %s 已完成（断点续跑跳过；如需重跑请删除 %s）", gse, done_marker)
    return(data.frame(gse = gse, status = "skipped_done", n_samples = NA,
                      expected = row$expected_total, match = "SKIPPED", note = ""))
  }
  dir.create(gse_dir, showWarnings = FALSE, recursive = TRUE)

  ## 2.1 下载系列矩阵（多 eSet 超系列自动逐个保存）
  eset_list <- tryCatch(
    getGEO(gse, GSEMatrix = TRUE, getGPL = FALSE, destdir = gse_dir),
    error = function(e) { msg("  【下载失败】%s：%s", gse, conditionMessage(e)); NULL }
  )
  if (is.null(eset_list))
    return(data.frame(gse = gse, status = "download_failed", n_samples = NA,
                      expected = row$expected_total, match = "FAIL",
                      note = "网络或GEO服务器问题，稍后重跑本脚本即可续传"))
  if (!is.list(eset_list)) eset_list <- list(eset_list)

  ## 2.2 逐 eSet 保存表达矩阵 + 样本元数据
  total_samples <- 0L
  for (i in seq_along(eset_list)) {
    es <- eset_list[[i]]
    gpl <- annotation(es)
    pdata <- pData(es)
    total_samples <- total_samples + nrow(pdata)
    suffix <- if (length(eset_list) > 1) paste0("_", gpl) else ""
    # 手动加 ID 列：不用 keep.rownames（行名为序号时 data.table 不建 rn 列会报错，
    # GSE114007 实机触发；2026-08-22 修复）
    fwrite(cbind(gsm = rownames(pdata), as.data.table(pdata)),
           file.path(gse_dir, paste0(gse, suffix, "_pData.csv")))
    expr <- exprs(es)
    fwrite(cbind(feature_id = rownames(expr), as.data.table(expr)),
           file.path(gse_dir, paste0(gse, suffix, "_expr.csv.gz")))
    msg("    · eSet %d/%d（%s）：%d 样本 × %d 特征", i, length(eset_list), gpl, nrow(pdata), nrow(expr))
  }

  ## 2.3 可选：下载补充文件（单细胞原始数据等大文件）
  supp_note <- ""
  if (need_supp) {
    msg("    · 下载补充文件（大文件，可能耗时较长）...")
    tryCatch({
      getGEOSuppFiles(gse, makeDirectory = FALSE, baseDir = gse_dir)
      writeLines(format(Sys.time()), supp_marker)
    }, error = function(e) msg("    【补充文件下载失败】%s（可稍后重跑，仅补此项）", conditionMessage(e)))
  } else if (gse %in% SUPP_GSES && !file.exists(supp_marker)) {
    supp_note <- "单细胞/单核补充文件未下载（DOWNLOAD_SUPP=FALSE），模块G启动前改TRUE重跑即可只补大文件"
  }

  ## 2.4 样本数比对
  match_flag <- if (is.na(row$expected_total)) "MANUAL_CHECK"
                else if (total_samples == row$expected_total) "OK"
                else "MISMATCH"
  note_parts <- c(supp_note,
                  if (match_flag == "MISMATCH")
                    sprintf("样本数与核对表不符（预期%s，实际%d），请对照 pData CSV 人工复核",
                            row$expected_total, total_samples) else "")
  note <- paste(note_parts[nzchar(note_parts)], collapse = "；")

  writeLines(format(Sys.time()), done_marker)   # 完成标记（断点续跑依据）
  data.frame(gse = gse, status = "downloaded", n_samples = total_samples,
             expected = row$expected_total, match = match_flag, note = note)
}

## ============================================================================
## 3. 主循环
## ============================================================================
msg("============================================")
msg("GEO 批量下载启动：%s（共 %d 个数据集）", format(Sys.time()), nrow(manifest))
msg("数据目录：%s", normalizePath(DATA_DIR))

results <- vector("list", nrow(manifest))
for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, ]
  msg("----------------------------------------")
  msg("[%d/%d] %s（%s / %s / %s）", i, nrow(manifest), row$gse, row$side, row$tissue, row$role)
  msg("  核对备注：%s", row$check_note)
  results[[i]] <- tryCatch(
    process_one(row),
    error = function(e) {
      msg("  【处理失败】%s：%s", row$gse, conditionMessage(e))
      data.frame(gse = row$gse, status = "error", n_samples = NA,
                 expected = row$expected_total, match = "FAIL", note = conditionMessage(e))
    }
  )
  ## ★ GSE167186 专项提醒（备忘录 §八-7，研究者 2026-08-18 指示）
  if (row$gse == "GSE167186" && results[[i]]$status %in% c("downloaded", "skipped_done")) {
    msg("")
    msg("  ★★★ 提醒：GSE167186 样本数核对（93 vs 89）—— 请现在打开")
    msg("      %s 下的 pData CSV，逐样本对清单：", file.path(DATA_DIR, "GSE167186"))
    msg("      确认多出的 4 个记录属于：技术重复 / 质控淘汰样本 / 辅助子系列，")
    msg("      并将无关记录剔除。预期口径：bulk 72（19年轻/29非肌少老年/24肌少老年）")
    msg("      + snRNA-seq 17 供体（143,051 核）。此结果将写入论文 Methods。")
    msg("")
  }
}

## ============================================================================
## 4. 汇总报告
## ============================================================================
report <- rbindlist(results, fill = TRUE)
report <- merge(report, manifest[, c("gse", "side", "tissue", "platform", "role", "check_note")],
                by = "gse", all.x = TRUE)
setcolorder(report, c("gse", "side", "tissue", "platform", "role",
                      "status", "n_samples", "expected", "match", "note", "check_note"))
write.csv(report, REPORT_FILE, row.names = FALSE, fileEncoding = "UTF-8")

n_ok   <- sum(report$match == "OK")
n_bad  <- sum(report$match %in% c("MISMATCH", "FAIL"))
n_skip <- sum(report$match == "SKIPPED")

msg("")
msg("==================== 下载核对总结 ====================")
msg("数据集总数 %d：样本数吻合 %d，需人工复核/失败 %d，跳过（已完成） %d",
    nrow(report), n_ok, n_bad + sum(report$match == "MANUAL_CHECK"), n_skip)
if (n_bad > 0) {
  msg("须处理清单：")
  for (g in report$gse[report$match %in% c("MISMATCH", "FAIL")]) msg("  ✘ %s", g)
}
msg("报告文件：%s", REPORT_FILE)
msg("下一步：逐数据集打开 *_pData.csv 做分组变量人工复核（对照核对表'分组变量'列），")
msg("        复核通过后即可进入模块 A 预处理。")
msg("=====================================================")
