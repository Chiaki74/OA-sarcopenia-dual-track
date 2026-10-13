# =============================================================================
# 44_口径核查_MAF与SuSiE_20261011.R
# 功能：主图审查遗留两项本地核实的证据搜集（只读，不改任何文件）
#   Part A（P1-3）：扫描项目内全部 .R 脚本，打印含 maf / 5e-8 / r2 / F-stat
#     的行（文件名 + 行号 + 内容），用于判定 MR 仪器筛选是否真实执行 MAF>1%
#   Part B（P2-2）：列出 results/ 下 moduleF 相关文件，并扫描 results/
#     内 <5MB 文本文件中含 converged / susie 的行，用于确认
#     正文「10/10 converged」与 Fig 5C「2 可算 + 3 不可算」的分母关系
# 运行：放项目根目录，source("44_口径核查_MAF与SuSiE_20261011.R")，
#       把控制台输出原样发回（截图或复制文本均可）
# 作者：刘晟宇    日期：2026-10-11
# =============================================================================

message("== 脚本版本：44_口径核查_MAF与SuSiE_20261011 ==")

## ---- 0. 项目根定位 -----------------------------------------------------------
.script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- grep("--file=", a, value = TRUE)
  if (length(f) >= 1) return(dirname(normalizePath(sub("--file=", "", f[1]))))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable() &&
      nzchar(rstudioapi::getActiveDocumentContext()$path))
    return(dirname(normalizePath(rstudioapi::getActiveDocumentContext()$path)))
  getwd()
}
.d <- .script_dir()
for (i in 1:8) {
  if (any(dir.exists(file.path(.d, c("processed", "results", "data"))))) break
  if (dirname(.d) == .d) break
  .d <- dirname(.d)
}
setwd(.d); message("[44] 项目根目录：", .d)

scan_lines <- function(files, pattern, label, max_per_file = 30L) {
  cat("\n========== ", label, " ==========\n")
  n_hit <- 0L
  for (f in files) {
    L <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
    if (length(L) == 0) next
    idx <- grep(pattern, L, ignore.case = TRUE)
    if (length(idx) == 0) next
    n_hit <- n_hit + 1L
    cat("\n--- ", f, " （", length(idx), " 行命中）---\n", sep = "")
    for (i in head(idx, max_per_file))
      cat(sprintf("  L%d: %s\n", i, trimws(L[i])))
    if (length(idx) > max_per_file)
      cat("  ...（其余 ", length(idx) - max_per_file, " 行省略）\n")
  }
  if (n_hit == 0) cat("（无任何命中）\n")
  invisible(n_hit)
}

## ---- Part A：MR 仪器筛选口径（全部 .R 脚本） ---------------------------------
r_files <- list.files(".", pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
r_files <- r_files[!grepl("/renv/|/packrat/", r_files)]
cat("\n[Part A] 共扫描 ", length(r_files), " 个 .R 脚本\n")

scan_lines(r_files, "maf", "A1: 含 MAF 的行（核心证据）")
scan_lines(r_files, "5e-8|5\\*10\\^-8|5e-08|clump|r2|F.stat|F_stat|Fstat",
           "A2: 含 5e-8 / clump / r2 / F-stat 的行（旁证）", max_per_file = 20L)

## ---- Part B：SuSiE 收敛口径 ---------------------------------------------------
cat("\n[Part B] results/ 下 moduleF / MR / coloc 相关文件清单：\n")
res_all <- list.files("results", recursive = TRUE, full.names = TRUE)
focus <- res_all[grepl("moduleF|coloc|susie|MR", res_all, ignore.case = TRUE)]
if (length(focus) == 0) cat("  （results/ 下未见 moduleF/coloc/susie/MR 字样文件）\n")
for (f in focus) cat("  ", f, "\n")

txt_files <- res_all[grepl("\\.(txt|log|csv|md)$", res_all, ignore.case = TRUE)]
sz <- file.size(txt_files); sz[is.na(sz)] <- 0
txt_files <- txt_files[sz < 5e6]   # 跳过大矩阵（如 region_sumstats 合并文件）
scan_lines(txt_files, "converg|susie", "B1: results/ 文本中含 converged / susie 的行")

cat("\n[44] 完成。请把以上全部输出发回（截图或复制文本）。\n")
cat("[44] 若 B 部清单中出现 SuSiE 汇总表（含 PPH4 / converged 列的 csv/txt），",
    "请将该文件一并发回。\n")
