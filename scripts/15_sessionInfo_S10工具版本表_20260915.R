# ============================================================================
# 15_sessionInfo_S10工具版本表_20260915.R
# 生成补充材料表 S10（分析工具与 R 包版本清单）——投稿前随代码库一并归档
# ----------------------------------------------------------------------------
# 产出：results/moduleS10/S10_tools_versions.csv（包/版本/所属模块/用途）
#   + S10_sessionInfo_full.txt（完整 sessionInfo 留档）
# 用法：在完成全部分析的 R 环境中 source() 本脚本；缺装的包会标「未安装」
#   并触发提醒（S10 须如实反映实际运行环境）。
# ============================================================================

suppressPackageStartupMessages(library(data.table))

OUT_DIR <- "results/moduleS10"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                     paste0(..., collapse = "")))

## 包 → 模块/用途映射（与备忘录各模块口径一致；版本以本机实测为准）
PKG_MAP <- list(
  "环境与数据获取" = c("GEOquery", "data.table", "readxl", "curl", "httr"),
  "模块A 预处理/QC" = c("affy", "oligo", "limma", "sva", "preprocessCore", "impute"),
  "模块B 差异表达"  = c("limma", "DESeq2", "edgeR"),
  "模块C 双轨整合"  = c("RobustRankAggreg", "WGCNA", "ggplot2", "pheatmap"),
  "模块D 机器学习"  = c("glmnet", "caret", "randomForest", "e1071", "xgboost",
                    "kernlab", "nnet", "naivebayes", "pROC"),
  "模块E 验证"     = c("pROC", "rms", "dcurves", "ResourceSelection"),
  "模块F 免疫浸润"  = c("xCell", "GSVA", "dicepro"),
  "模块G 单细胞"   = c("Seurat", "SeuratObject", "harmony", "UCell", "Matrix"),
  "模块H MR/共定位" = c("TwoSampleMR", "ieugwasr", "coloc", "MRPRESSO"),
  "模块I 列线图/稳健性" = c("rms", "caret", "pROC"),
  "模块J 调控网络"  = c("enrichR", "igraph"),
  "模块K 细胞通讯"  = c("CellChat", "future", "NMF", "ggalluvial"),
  "通用可视化"     = c("ggplot2", "ggrepel", "patchwork", "circlize", "ComplexHeatmap")
)

rows <- list()
for (mod in names(PKG_MAP)) {
  for (pkg in unique(PKG_MAP[[mod]])) {
    ver <- tryCatch(as.character(packageVersion(pkg)),
                    error = function(e) NA_character_)
    rows[[length(rows) + 1]] <- data.frame(
      package = pkg, version = ifelse(is.na(ver), "未安装", ver),
      module = mod, stringsAsFactors = FALSE)
  }
}
s10 <- unique(rbindlist(rows))
missing_pkgs <- s10$package[s10$version == "未安装"]
if (length(missing_pkgs) > 0)
  log_msg("⚠ 以下包在本环境未安装（S10 如实标注；若某模块实际未用到可人工删行）：",
          paste(missing_pkgs, collapse = ", "))

fwrite(s10, file.path(OUT_DIR, "S10_tools_versions.csv"))
sink(file.path(OUT_DIR, "S10_sessionInfo_full.txt"))
cat("R version: ", R.version.string, "\nOS: ", Sys.info()["sysname"], " ",
    Sys.info()["release"], "\n\n", sep = "")
print(sessionInfo())
sink()

log_msg("S10 已生成：", nrow(s10), " 行（", length(unique(s10$package)),
        " 个包）→ results/moduleS10/S10_tools_versions.csv")
log_msg("完整 sessionInfo 留档 → S10_sessionInfo_full.txt")
log_msg("归档动作：把 S10_tools_versions.csv 收入补充材料；本脚本与版本表一并入库。")
