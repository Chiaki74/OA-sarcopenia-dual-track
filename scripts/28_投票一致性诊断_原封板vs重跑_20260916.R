## =============================================================================
## 28_投票一致性诊断_原封板vs重跑_20260916.R
## 目的：Fig3 面板 (a)/(b) 画的是 27 号重跑的逐算法投票表；必须确认其 top10
##       与论文封板的 10 个标志物一致，否则图与文不符。本脚本只读不写。
## 运行：source("28_投票一致性诊断_原封板vs重跑_20260916.R")，全部输出发我
## =============================================================================

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
setwd(.d); message("[setwd] ", .d)

BM <- c("NELL1", "STEAP1", "GADD45A", "FBLN5", "HEMK1",
        "GDE1", "RNF14", "BNIP3", "EGFR", "CBR3")

## ---- 1. 27 号重跑的逐算法投票表 ----
f1 <- "results/moduleD_v2/feature_votes_by_algorithm.csv"
stopifnot(file.exists(f1))
v1 <- read.csv(f1, check.names = FALSE, stringsAsFactors = FALSE)
alg <- setdiff(colnames(v1), "gene")
v1$n_votes <- rowSums(v1[, alg, drop = FALSE])
cat("====== [1] 重跑投票表 ======\n")
cat("票数分布（0-4 票各多少基因）：\n"); print(table(v1$n_votes))
top10 <- v1[order(-v1$n_votes, v1$gene), ][seq_len(10), ]
cat("\n重跑 top10（Fig3 面板 a/b 将展示这些基因）：\n")
print(top10[, c("gene", alg, "n_votes")], row.names = FALSE)
cat("\n10 个封板标志物在重跑中的逐算法表现：\n")
print(v1[v1$gene %in% BM, c("gene", alg, "n_votes")], row.names = FALSE)
cat("\n>>> 判定1：重跑 top10 与封板 10 基因完全一致？",
    setequal(top10$gene, BM), "\n")
cat(">>> 判定2：封板 10 基因是否恰好是票数最高的 10 个（无并列挤出）？")
nb <- sum(v1$n_votes >= min(top10$n_votes))
cat(" 票数>=top10 最低票的基因数 =", nb, "（=10 则无并列问题）\n")

## ---- 2. 原封板投票表（05b 当时输出）----
f0 <- "results/moduleD_v2/feature_selection_votes.csv"
cat("\n====== [2] 原封板 feature_selection_votes.csv ======\n")
if (!file.exists(f0)) {
  cat("未找到", f0, "——请把该文件实际路径发我\n")
} else {
  v0 <- read.csv(f0, check.names = FALSE, stringsAsFactors = FALSE)
  cat("列：", paste(colnames(v0), collapse = ", "), "\n")
  cat("票数分布：\n"); print(table(v0$votes))
  cat("selected 汇总：\n"); print(table(v0$selected))
  cat("\n封板 10 基因在原表中的 votes/selected：\n")
  print(v0[v0$gene %in% BM, ], row.names = FALSE)
  sel_genes <- v0$gene[v0$selected %in% c(TRUE, 1, "TRUE", "1", "true")]
  cat("\n>>> 判定3：原表 selected=TRUE 的基因与封板 10 一致？",
      setequal(sel_genes, BM), "（selected 共", length(sel_genes), "个）\n")
  cat(">>> 判定4：原表 selected 恰好等于 votes==4？",
      setequal(sel_genes, v0$gene[v0$votes == 4]), "\n")
  ## 重跑 vs 原封板的票数对照（只看封板 10 基因）
  mg <- merge(v0[, c("gene", "votes")], v1[, c("gene", "n_votes")], by = "gene")
  mg <- mg[mg$gene %in% BM, ]
  cat("\n封板 10 基因：原 votes vs 重跑 n_votes：\n")
  print(mg, row.names = FALSE)
}
cat("\n====== 诊断结束，请把全部输出发我 ======\n")
