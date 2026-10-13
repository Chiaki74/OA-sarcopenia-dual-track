## =============================================================================
## 26_出图预检2_实际文件清单_20260916.R
## 目的：列出图脚本缺失文件对应的实际本地文件名，供对齐 CONFIG 用。
## 运行：source("26_出图预检2_实际文件清单_20260916.R")，把全部输出截图发我
## =============================================================================

.script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- grep("--file=", a, value = TRUE)
  if (length(f) >= 1) return(dirname(normalizePath(sub("--file=", "", f[1]))))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable() &&
      nzchar(rstudioapi::getActiveDocumentContext()$path))
    return(dirname(rstudioapi::getActiveDocumentContext()$path))
  getwd()
}
.d <- .script_dir()
for (i in 1:8) {
  if (any(dir.exists(file.path(.d, c("processed", "results", "data"))))) break
  if (dirname(.d) == .d) break
  .d <- dirname(.d)
}
setwd(.d); message("[setwd] 项目根目录：", .d)

show <- function(title, files) {
  cat("\n── ", title, " ──\n", sep = "")
  if (length(files) == 0) cat("  （无）\n") else
    for (f in files) cat("  ", f, "\n", sep = "")
}

## 1. Fig2 缺：PCA / 批次前后矩阵 / 分组表 / 热图长表
show("PCA 相关（全项目搜 pca）",
     list.files(".", pattern = "pca", recursive = TRUE, ignore.case = TRUE))
show("批次前后表达矩阵（搜 before|after|batch）",
     list.files(".", pattern = "before|after|batch", recursive = TRUE, ignore.case = TRUE))
show("分组表（搜 sample_group|group.*csv）",
     list.files(".", pattern = "sample_group|sample.*group", recursive = TRUE, ignore.case = TRUE))
show("热图长表（搜 heatmap|heat）",
     list.files(".", pattern = "heatmap|heat", recursive = TRUE, ignore.case = TRUE))

## 2. Fig3 缺：feature_votes
show("特征投票表（搜 vote）",
     list.files(".", pattern = "vote", recursive = TRUE, ignore.case = TRUE))
show("SHAP 相关（搜 shap，可选备路）",
     list.files(".", pattern = "shap", recursive = TRUE, ignore.case = TRUE))
show("列线图相关（搜 nomogram|nomo，可选备路）",
     list.files(".", pattern = "nomogram|nomo", recursive = TRUE, ignore.case = TRUE))

## 3. Fig5：moduleH 全部文件（共定位/SuSiE/区域汇总/lead SNP 的真实文件名）
show("results/moduleH 全部（递归）",
     list.files("results/moduleH", recursive = TRUE))

## 4. moduleD / moduleD_v2 全部文件（feature_votes 与 SHAP 可能在这里）
show("results/moduleD 全部（递归）",
     list.files("results/moduleD", recursive = TRUE))
show("results/moduleD_v2 全部（递归）",
     list.files("results/moduleD_v2", recursive = TRUE))

## 5. data/ 目录（分组表可能在这里）
if (dir.exists("data")) show("data/ 全部（递归）", list.files("data", recursive = TRUE))

cat("\n把以上全部输出截图发我。\n")
