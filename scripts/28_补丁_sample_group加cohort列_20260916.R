## ===========================================================================
## 28_补丁_sample_group加cohort列_20260916.R
## ---------------------------------------------------------------------------
## 目的：给 v1 版 27 号脚本已生成的 data/sample_group.csv 追加 cohort 列
##       （Fig2-A「Before 按数据集着色」的唯一缺口），不重跑整条流水线。
## 背景：27 via2/via3 的全量重建路线依赖 assemble_cohort，其实际定义在
##       05b_moduleD_v2_multicohort.R 内——source 它会触发整条 moduleD
##       流水线（且撞上 caret::train 版本问题），代价大、无必要。
## 原理：v1 的 grp_all 按 TRAIN_SETS 顺序逐队列 rbind（备忘录实载），
##       样本块顺序确定；各队列样本数 38/20/20/19/20，合计 117，
##       与 05b 控制台日志「合并训练集：117 样本」完全一致。
## 防呆：尾块 20 个样本必须与 GSE51588 donor_id 去重后保留集完全一致
##       （05b 控制台日志实锤的 20 个 GSM），不一致立即报错终止，
##       绝不静默错标——此时请把 v1 版 27 脚本发我改精确方案。
## 运行：秒级；只改写 data/sample_group.csv 一个文件（先自动备份原版）。
## ===========================================================================

sg_path <- file.path(getwd(), "data", "sample_group.csv")
stopifnot("找不到 data/sample_group.csv，请确认 setwd 到项目根目录" =
            file.exists(sg_path))
sg <- read.csv(sg_path, check.names = FALSE, stringsAsFactors = FALSE)

## 样本列防御性探测（正常应为 "sample"）
s_col <- if ("sample" %in% names(sg)) "sample" else names(sg)[1]
if (s_col != "sample") {
  names(sg)[names(sg) == s_col] <- "sample"
  message("[补丁] 样本列名非 sample，已按首列处理：", s_col)
}

TRAIN_SETS <- c("GSE114007", "GSE55235", "GSE55457", "GSE12021", "GSE51588")
N_PER      <- c(38L, 20L, 20L, 19L, 20L)   # 合计 117

if (nrow(sg) != sum(N_PER))
  stop("sample_group.csv 行数 = ", nrow(sg), " ≠ 117——块结构假设不成立，",
       "请把 v1 版 27 脚本（或本报错）发我，改精确方案。")

## 硬校验：尾块 20 样本 == GSE51588 保留集（05b 日志实锤）
GSE51588_GSM <- c("GSM1248759", "GSM1248760", "GSM1248761", "GSM1248762",
                  "GSM1248763", "GSM1248769", "GSM1248770", "GSM1248771",
                  "GSM1248772", "GSM1248773", "GSM1248774", "GSM1248775",
                  "GSM1248776", "GSM1248777", "GSM1248778", "GSM1248779",
                  "GSM1248785", "GSM1248786", "GSM1248787", "GSM1248788")
tail20 <- tail(sg$sample, 20)
if (!setequal(tail20, GSE51588_GSM))
  stop("尾块 20 样本与 GSE51588 保留集不一致——块顺序假设不成立，",
       "请把 v1 版 27 脚本（或本报错）发我，改精确方案。\n",
       "  实得尾块前 5 个：", paste(head(tail20, 5), collapse = ", "))
message("[补丁] 尾块硬校验通过：20/20 匹配 GSE51588 保留集")

## 备份原版（仅首次）后写入 cohort 列
bak <- sub("\\.csv$", "_无cohort备份.csv", sg_path)
if (!file.exists(bak)) file.copy(sg_path, bak)
sg$cohort <- rep(TRAIN_SETS, times = N_PER)
write.csv(sg, sg_path, row.names = FALSE, quote = FALSE)

## 各块首尾 GSM 打印，供人工二次确认
starts <- cumsum(c(1L, head(N_PER, -1L)))
for (i in seq_along(TRAIN_SETS)) {
  blk <- sg$sample[starts[i]:(starts[i] + N_PER[i] - 1L)]
  cat(sprintf("%-10s n=%2d  首 = %s  尾 = %s\n",
              TRAIN_SETS[i], N_PER[i], head(blk, 1), tail(blk, 1)))
}
cat("\n✅ cohort 列已写入：", sg_path,
    "\n   原版备份：", bak, "\n")
