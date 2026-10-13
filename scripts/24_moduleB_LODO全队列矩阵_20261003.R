# ============================================================================
# 24_moduleB_LODO全队列矩阵_20261003.R
# 全队列留一法（leave-one-dataset-out, LODO）矩阵：
#   肌肉侧 5 队列 + OA 侧 5 训练队列逐个剔除重跑 RRA，量化每个队列对
#   ①本侧 RRA 显著集合、②双轨 99 候选集/72 镜像轨 的影响力。
#   用于回答：GSE111006（36 对照/4 肌少症）的留一效应在全部队列中是否特殊。
# ----------------------------------------------------------------------------
# 前置：23_moduleB_LODO_GSE111006敏感性_20261003.R 已验证全量复跑与主分析
#       CSV 逐分一致（run_rra 复现忠实）。本脚本复用同一 run_rra。
# 输入：results/moduleB/layer1_DE_all_datasets.rds
#       results/moduleB/layer2_RRA_muscle.csv、layer2_RRA_OA.csv（主分析基准）
# 输出：results/moduleB/LODO/LODO全队列矩阵.csv
#       results/moduleB/LODO/LODO全队列矩阵_摘要.txt（可直接引用进补充材料）
# 运行（项目根目录）：
#   source("24_moduleB_LODO全队列矩阵_20261003.R")
# ============================================================================

suppressPackageStartupMessages(library(RobustRankAggreg))

RRA_CUT            <- 0.05
MUSCLE_NOMINAL_CUT <- 0.05

## ---- 定位项目根目录 --------------------------------------------------------
locate_root <- function() {
  d <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  for (i in 1:8) {
    if (dir.exists(file.path(d, "results", "moduleB"))) return(d)
    parent <- dirname(d)
    if (parent == d) break
    d <- parent
  }
  stop("未找到 results/moduleB/——请在项目根目录运行")
}
ROOT  <- locate_root()
MBDIR <- file.path(ROOT, "results", "moduleB")
OUTD  <- file.path(MBDIR, "LODO")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
message("项目根目录：", ROOT)

## ---- run_rra：与 03_moduleB_shared_genes.R L144-163 逐行同构 ---------------
run_rra <- function(de_list, side_name) {
  de_list <- de_list[!vapply(de_list, is.null, logical(1))]
  if (length(de_list) < 2) stop(side_name, " 侧有效数据集 < 2，无法 RRA")
  up_lists   <- lapply(de_list, function(tt) tt$gene[order(tt$P.Value, -abs(tt$logFC))][tt$logFC[order(tt$P.Value, -abs(tt$logFC))] > 0])
  down_lists <- lapply(de_list, function(tt) tt$gene[order(tt$P.Value, -abs(tt$logFC))][tt$logFC[order(tt$P.Value, -abs(tt$logFC))] < 0])
  universe <- Reduce(union, lapply(de_list, `[[`, "gene"))
  up   <- aggregateRanks(glist = up_lists,   N = length(universe))
  down <- aggregateRanks(glist = down_lists, N = length(universe))
  colnames(up) <- c("gene", "RRA_score_up"); colnames(down) <- c("gene", "RRA_score_down")
  m <- merge(up, down, by = "gene")
  m$direction <- ifelse(m$RRA_score_up < m$RRA_score_down, "up", "down")
  m$RRA_score <- pmin(m$RRA_score_up, m$RRA_score_down)
  m$RRA_adj   <- p.adjust(m$RRA_score, method = "BH")
  m[order(m$RRA_adj), ]
}

## ---- 双轨共享集构建：与 03 脚本 L243-249 同构 -------------------------------
build_shared <- function(rra_oa_tab, rra_mu_tab) {
  s <- merge(rra_oa_tab[rra_oa_tab$RRA_adj < RRA_CUT, c("gene", "direction")],
             rra_mu_tab[rra_mu_tab$RRA_score < MUSCLE_NOMINAL_CUT, c("gene", "direction")],
             by = "gene", suffixes = c("_OA", "_muscle"))
  s$concordant <- s$direction_OA == s$direction_muscle
  s
}

## ---- 读取数据与主分析基准 ---------------------------------------------------
de_all <- readRDS(file.path(MBDIR, "layer1_DE_all_datasets.rds"))
de_oa <- de_all$OA; de_mu <- de_all$muscle
rra_oa_main <- read.csv(file.path(MBDIR, "layer2_RRA_OA.csv"),     check.names = FALSE)
rra_mu_main <- read.csv(file.path(MBDIR, "layer2_RRA_muscle.csv"), check.names = FALSE)

## 主分析基准集合
oa511   <- rra_oa_main$gene[rra_oa_main$RRA_adj < RRA_CUT]
mu_adj  <- rra_mu_main$gene[rra_mu_main$RRA_adj < RRA_CUT]          # 29
mu_nom  <- rra_mu_main$gene[rra_mu_main$RRA_score < MUSCLE_NOMINAL_CUT]
shared_main <- build_shared(rra_oa_main, rra_mu_main)               # 99
cand99  <- shared_main$gene
mirror72 <- shared_main$gene[!shared_main$concordant]
message("基准：OA 511=", length(oa511), "；肌肉 adj 集=", length(mu_adj),
        "；候选=", length(cand99), "（镜像 ", length(mirror72), "）")

## ---- 单侧留一 ---------------------------------------------------------------
lodo_one <- function(side, drop, de_side, rra_main, rra_other_main) {
  rra_lodo <- run_rra(de_side[setdiff(names(de_side), drop)], paste0(side, "-no", drop))
  ## 本侧显著集合保持
  if (side == "OA") {
    main_set <- oa511; lodo_set <- rra_lodo$gene[rra_lodo$RRA_adj < RRA_CUT]
  } else {
    main_set <- mu_adj; lodo_set <- rra_lodo$gene[rra_lodo$RRA_adj < RRA_CUT]
  }
  retained <- length(intersect(main_set, lodo_set))
  jaccard  <- retained / max(length(union(main_set, lodo_set)), 1)
  ## 全基因组秩相关
  common <- intersect(rra_main$gene, rra_lodo$gene)
  o <- rra_main[match(common, rra_main$gene), ]; l <- rra_lodo[match(common, rra_lodo$gene), ]
  rho <- cor(-log10(o$RRA_score), -log10(l$RRA_score), method = "spearman")
  ## 双轨候选重构（本侧 LODO × 对侧主分析）
  shared_l <- if (side == "OA") build_shared(rra_lodo, rra_other_main)
              else              build_shared(rra_other_main, rra_lodo)
  cand_l   <- shared_l$gene
  mir_l    <- shared_l$gene[!shared_l$concordant]
  cand_ret <- length(intersect(cand99, cand_l))
  mir_ret  <- length(intersect(mirror72, mir_l))
  mir_pct  <- 100 * length(mir_l) / max(length(cand_l), 1)
  mir_p    <- if (length(cand_l) > 0)
    binom.test(length(mir_l), length(cand_l), p = 0.5, alternative = "greater")$p.value else NA_real_
  data.frame(
    侧 = side, 剔除队列 = drop, 剩余队列数 = length(de_side) - 1,
    本侧显著集_主 = length(main_set), 本侧显著集_LODO = length(lodo_set),
    本侧显著集_保留 = retained, 本侧显著集_保留率 = round(100 * retained / max(length(main_set), 1), 1),
    本侧显著集_Jaccard = round(jaccard, 3),
    全基因组秩相关rho = round(rho, 4),
    候选数_主 = length(cand99), 候选数_LODO = length(cand_l),
    候选保留 = cand_ret, 候选保留率 = round(100 * cand_ret / length(cand99), 1),
    镜像数_LODO = length(mir_l), 镜像占比_LODO = round(mir_pct, 1),
    镜像二项P_LODO = signif(mir_p, 3),
    镜像轨保留 = mir_ret, 镜像轨保留率 = round(100 * mir_ret / length(mirror72), 1),
    stringsAsFactors = FALSE)
}

res <- list()
for (g in names(de_mu)) res[[paste0("mu_", g)]] <- lodo_one("肌肉", g, de_mu, rra_mu_main, rra_oa_main)
for (g in names(de_oa)) res[[paste0("oa_", g)]] <- lodo_one("OA",  g, de_oa, rra_oa_main, rra_mu_main)
mat <- do.call(rbind, res)
rownames(mat) <- NULL
write.csv(mat, file.path(OUTD, "LODO全队列矩阵.csv"), row.names = FALSE)

## ---- 摘要 -------------------------------------------------------------------
fmt <- function(r) sprintf(
  "  剔 %-10s：本侧显著集保留 %s/%s（%s%%，Jaccard=%s）；候选 %s→%s（保留 %s%%）；镜像轨保留 %s%%；镜像占比 %s%%（P=%s）；ρ=%s",
  r$剔除队列, r$本侧显著集_保留, r$本侧显著集_主, r$本侧显著集_保留率, r$本侧显著集_Jaccard,
  r$候选数_主, r$候选数_LODO, r$候选保留率, r$镜像轨保留率, r$镜像占比_LODO,
  format(r$镜像二项P_LODO, digits = 2), r$全基因组秩相关rho)

mu_rows <- mat[mat$侧 == "肌肉", ]; oa_rows <- mat[mat$侧 == "OA", ]
fmt_rows <- function(df) vapply(seq_len(nrow(df)), function(i) fmt(df[i, ]), character(1))
sm <- c(
  "全队列留一法（LODO）矩阵：各队列对 RRA 显著集与双轨候选集的影响力",
  sprintf("生成时间：%s", Sys.time()),
  sprintf("基准：OA 511 基因；肌肉轨校正显著 %d；双轨候选 %d（镜像 %d，%.1f%%）",
          length(mu_adj), length(cand99), length(mirror72), 100 * length(mirror72) / length(cand99)),
  "",
  "【肌肉侧留一】", fmt_rows(mu_rows),
  "",
  "【OA 侧留一】", fmt_rows(oa_rows),
  "",
  sprintf("肌肉侧候选保留率范围：%.1f–%.1f%%（GSE111006=%.1f%%）；镜像轨保留率范围：%.1f–%.1f%%（GSE111006=%.1f%%）",
          min(mu_rows$候选保留率), max(mu_rows$候选保留率), mu_rows$候选保留率[mu_rows$剔除队列 == "GSE111006"],
          min(mu_rows$镜像轨保留率), max(mu_rows$镜像轨保留率), mu_rows$镜像轨保留率[mu_rows$剔除队列 == "GSE111006"]),
  sprintf("全部留一中镜像占比范围：%.1f–%.1f%%，二项检验最大 P=%s（均<0.05 则镜像主导结论全留一稳健）",
          min(mat$镜像占比_LODO), max(mat$镜像占比_LODO), format(max(mat$镜像二项P_LODO), digits = 2))
)
writeLines(sm, file.path(OUTD, "LODO全队列矩阵_摘要.txt"))
cat(paste0(sm, collapse = "\n"), "\n")
message("完成。输出：", file.path(OUTD, "LODO全队列矩阵.csv"))
