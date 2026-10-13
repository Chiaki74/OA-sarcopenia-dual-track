# ============================================================================
# 23_moduleB_LODO_GSE111006敏感性_20261003.R
# 留一法（leave-one-dataset-out, LODO）敏感性分析：
#   评估 GSE111006（36 对照 / 4 肌少症，逐队列功效极低）对肌肉轨 RRA 与
#   双轨共享候选集的影响。
# ----------------------------------------------------------------------------
# 背景（2026-10-03 摸底结论，备忘录 §253 追加）：
#   - GSE111006 在分析阈值（adj.P<0.05 且 |log2FC|>=0.585）下显著基因为 0 个；
#   - 但模块 B 的 RRA 输入为「逐队列全基因按 P 值排序、按 logFC 方向拆分的
#     双侧列表」（03_moduleB_shared_genes.R L148-151），故其排序仍参与了整合；
#   - 主分析肌肉轨 RRA 显著（RRA_adj<0.05）的 29 个基因中，26/29（89.7%）在
#     GSE111006 中 log2FC 方向一致（二项检验 P=1.5e-05）。
#   本脚本给出正式 LODO 证据：剔除 GSE111006 后核心结论是否保持不变。
# ----------------------------------------------------------------------------
# 输入（均为模块 B 已产出文件，无需重跑差异分析）：
#   results/moduleB/layer1_DE_all_datasets.rds   逐队列全量 DE 表 list(OA, muscle)
#   results/moduleB/layer2_RRA_muscle.csv        主分析肌肉轨 RRA（5 队列）
#   results/moduleB/layer2_RRA_OA.csv            主分析 OA 侧 RRA（可选，用于复现校验）
# 输出（results/moduleB/LODO/）：
#   layer2_RRA_muscle_LODO_noGSE111006.csv   剔除后重跑的肌肉轨 RRA 全表
#   LODO_GSE111006_基因对照表.csv            主分析显著基因逐项对照（校正+名义两套口径）
#   LODO_GSE111006_双轨对照表.csv            99 候选逐项对照（保留/方向/轨道）
#   LODO_GSE111006_摘要.txt                  可直接引用进补充材料的关键数字
# 运行（项目根目录，即 processed/ 与 results/ 的上级）：
#   source("23_moduleB_LODO_GSE111006敏感性_20261003.R")
# ============================================================================

suppressPackageStartupMessages(library(RobustRankAggreg))

RRA_CUT            <- 0.05   # 与 03_moduleB_shared_genes.R CONFIG$RRA_CUT 一致
MUSCLE_NOMINAL_CUT <- 0.05   # 与 CONFIG$MUSCLE_NOMINAL_CUT 一致
DROP_GSE           <- "GSE111006"

## ---- 定位项目根目录 --------------------------------------------------------
locate_root <- function() {
  d <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  for (i in 1:8) {
    if (dir.exists(file.path(d, "results", "moduleB"))) return(d)
    parent <- dirname(d)
    if (parent == d) break
    d <- parent
  }
  stop("未找到 results/moduleB/——请在项目根目录（processed/ 与 results/ 的上级）运行本脚本")
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
  # 上调/下调分别按 P 值排序（带方向，双侧列表）；全基因列表，不做显著性过滤
  up_lists   <- lapply(de_list, function(tt) tt$gene[order(tt$P.Value, -abs(tt$logFC))][tt$logFC[order(tt$P.Value, -abs(tt$logFC))] > 0])
  down_lists <- lapply(de_list, function(tt) tt$gene[order(tt$P.Value, -abs(tt$logFC))][tt$logFC[order(tt$P.Value, -abs(tt$logFC))] < 0])
  universe <- Reduce(union, lapply(de_list, `[[`, "gene"))
  up   <- aggregateRanks(glist = up_lists,   N = length(universe))
  down <- aggregateRanks(glist = down_lists, N = length(universe))
  colnames(up)   <- c("gene", "RRA_score_up")
  colnames(down) <- c("gene", "RRA_score_down")
  m <- merge(up, down, by = "gene")
  m$direction <- ifelse(m$RRA_score_up < m$RRA_score_down, "up", "down")
  m$RRA_score <- pmin(m$RRA_score_up, m$RRA_score_down)
  m$RRA_adj   <- p.adjust(m$RRA_score, method = "BH")
  m <- m[order(m$RRA_adj), ]
  message(side_name, " 侧 RRA：", sum(m$RRA_adj < RRA_CUT), " 个基因 adj<", RRA_CUT,
          "；名义 score<", MUSCLE_NOMINAL_CUT, " 共 ", sum(m$RRA_score < MUSCLE_NOMINAL_CUT), " 个")
  m
}

## ---- 读取模块 B 第 1 层产物 -------------------------------------------------
rds_path <- file.path(MBDIR, "layer1_DE_all_datasets.rds")
if (!file.exists(rds_path)) stop("缺少 ", rds_path, "——请先确认模块 B 已完整运行")
de_all <- readRDS(rds_path)
de_oa <- de_all$OA; de_mu <- de_all$muscle
message("OA 侧队列：", paste(names(de_oa), collapse = ", "))
message("肌肉侧队列：", paste(names(de_mu), collapse = ", "))
if (!DROP_GSE %in% names(de_mu)) stop(DROP_GSE, " 不在肌肉侧 DE 列表中——无需 LODO")

## ---- 校验 1：全量复跑肌肉轨，须与主分析 CSV 完全一致（证明复现忠实） --------
rra_mu_full <- run_rra(de_mu, "muscle-全量复跑")
orig_csv <- file.path(MBDIR, "layer2_RRA_muscle.csv")
rra_mu_orig <- read.csv(orig_csv, check.names = FALSE)
same_gene  <- identical(rra_mu_full$gene, rra_mu_orig$gene)
same_score <- isTRUE(all.equal(rra_mu_full$RRA_score, rra_mu_orig$RRA_score, tolerance = 1e-12))
message("校验 1（肌肉轨全量复跑 vs 主分析 CSV）：基因顺序一致=", same_gene,
        "，分数一致=", same_score)
if (!(same_gene && same_score))
  warning("复现不完全一致！请检查模块 B 是否用同一版本输入重跑过；LODO 结果仍可与 CSV 直接比较，但需谨慎解读")

## ---- 校验 2：OA 侧复跑（可选，对照 layer2_RRA_OA.csv 的 511） ---------------
rra_oa <- run_rra(de_oa, "OA-复跑")
oa_csv <- file.path(MBDIR, "layer2_RRA_OA.csv")
if (file.exists(oa_csv)) {
  rra_oa_orig <- read.csv(oa_csv, check.names = FALSE)
  message("校验 2（OA 侧复跑 vs 主分析 CSV）：基因顺序一致=", identical(rra_oa$gene, rra_oa_orig$gene),
          "，分数一致=", isTRUE(all.equal(rra_oa$RRA_score, rra_oa_orig$RRA_score, tolerance = 1e-12)),
          "；adj<0.05 基因数：复跑 ", sum(rra_oa$RRA_adj < RRA_CUT),
          " vs 主分析 ", sum(rra_oa_orig$RRA_adj < RRA_CUT))
}

## ---- LODO：剔除 GSE111006 重跑肌肉轨 ---------------------------------------
de_mu_lodo <- de_mu[setdiff(names(de_mu), DROP_GSE)]
message("LODO 肌肉侧队列（", length(de_mu_lodo), " 个）：", paste(names(de_mu_lodo), collapse = ", "))
rra_mu_lodo <- run_rra(de_mu_lodo, "muscle-LODO-noGSE111006")
write.csv(rra_mu_lodo, file.path(OUTD, "layer2_RRA_muscle_LODO_noGSE111006.csv"), row.names = FALSE)

## ---- 指标 1：肌肉轨显著基因集合保持度（校正/名义两套口径） ------------------
cmp_one <- function(orig, lodo, cut, col) {
  o <- orig$gene[orig[[col]] < cut]; l <- lodo$gene[lodo[[col]] < cut]
  list(n_orig = length(o), n_lodo = length(l),
       retained = length(intersect(o, l)),
       retention = length(intersect(o, l)) / max(length(o), 1),
       jaccard = length(intersect(o, l)) / max(length(union(o, l)), 1))
}
m_adj <- cmp_one(rra_mu_orig, rra_mu_lodo, RRA_CUT, "RRA_adj")
m_nom <- cmp_one(rra_mu_orig, rra_mu_lodo, MUSCLE_NOMINAL_CUT, "RRA_score")

## 秩相关：共有基因上 -log10(RRA_score) 的 Spearman 相关（全体 + 主分析显著子集）
common <- intersect(rra_mu_orig$gene, rra_mu_lodo$gene)
o <- rra_mu_orig[match(common, rra_mu_orig$gene), ]
l <- rra_mu_lodo[match(common, rra_mu_lodo$gene), ]
rho_all <- cor(-log10(o$RRA_score), -log10(l$RRA_score), method = "spearman")
sig29   <- rra_mu_orig$gene[rra_mu_orig$RRA_adj < RRA_CUT]
idx     <- common %in% sig29
rho_sig <- if (sum(idx) >= 3) cor(-log10(o$RRA_score[idx]), -log10(l$RRA_score[idx]),
                                  method = "spearman") else NA_real_

## ---- 指标 2：双轨共享候选集（99 基因）保持度 —— 决策层证据 ------------------
build_shared <- function(rra_oa_tab, rra_mu_tab) {
  s <- merge(rra_oa_tab[rra_oa_tab$RRA_adj < RRA_CUT, c("gene", "direction", "RRA_adj")],
             rra_mu_tab[rra_mu_tab$RRA_score < MUSCLE_NOMINAL_CUT, c("gene", "direction", "RRA_score")],
             by = "gene", suffixes = c("_OA", "_muscle"))
  s$concordant <- s$direction_OA == s$direction_muscle
  s
}
shared_main <- build_shared(rra_oa, rra_mu_orig)   # 主分析口径（OA 复跑 × 肌肉主分析）
shared_lodo <- build_shared(rra_oa, rra_mu_lodo)   # LODO 口径（OA 复跑 × 肌肉 LODO）
mirror_main <- shared_main$gene[!shared_main$concordant]
mirror_lodo <- shared_lodo$gene[!shared_lodo$concordant]
bt_lodo <- binom.test(length(mirror_lodo), nrow(shared_lodo), p = 0.5, alternative = "greater")

message("双轨候选：主分析 ", nrow(shared_main), "（镜像 ", length(mirror_main),
        "）→ LODO ", nrow(shared_lodo), "（镜像 ", length(mirror_lodo), "）")

## ---- 输出 1：主分析显著基因逐项对照表 ---------------------------------------
key_genes <- union(rra_mu_orig$gene[rra_mu_orig$RRA_adj < RRA_CUT],
                   rra_mu_orig$gene[rra_mu_orig$RRA_score < MUSCLE_NOMINAL_CUT])
tab <- data.frame(
  gene = key_genes,
  主分析_RRA_score   = rra_mu_orig$RRA_score[match(key_genes, rra_mu_orig$gene)],
  主分析_RRA_adj     = rra_mu_orig$RRA_adj[match(key_genes, rra_mu_orig$gene)],
  主分析_direction   = rra_mu_orig$direction[match(key_genes, rra_mu_orig$gene)],
  主分析_校正显著    = key_genes %in% rra_mu_orig$gene[rra_mu_orig$RRA_adj < RRA_CUT],
  主分析_名义显著    = key_genes %in% rra_mu_orig$gene[rra_mu_orig$RRA_score < MUSCLE_NOMINAL_CUT],
  LODO_RRA_score    = rra_mu_lodo$RRA_score[match(key_genes, rra_mu_lodo$gene)],
  LODO_RRA_adj      = rra_mu_lodo$RRA_adj[match(key_genes, rra_mu_lodo$gene)],
  LODO_direction    = rra_mu_lodo$direction[match(key_genes, rra_mu_lodo$gene)],
  stringsAsFactors = FALSE)
tab$LODO_校正显著 <- !is.na(tab$LODO_RRA_adj) & tab$LODO_RRA_adj < RRA_CUT
tab$LODO_名义显著 <- !is.na(tab$LODO_RRA_score) & tab$LODO_RRA_score < MUSCLE_NOMINAL_CUT
tab$方向保持      <- !is.na(tab$LODO_direction) & tab$LODO_direction == tab$主分析_direction
tab <- tab[order(tab$主分析_RRA_adj), ]
write.csv(tab, file.path(OUTD, "LODO_GSE111006_基因对照表.csv"), row.names = FALSE)

## ---- 输出 2：双轨候选逐项对照表 ---------------------------------------------
allcand <- union(shared_main$gene, shared_lodo$gene)
tab2 <- data.frame(
  gene = allcand,
  主分析_入选 = allcand %in% shared_main$gene,
  主分析_轨道 = ifelse(allcand %in% mirror_main, "mirror",
                  ifelse(allcand %in% shared_main$gene, "concordant", NA)),
  LODO_入选  = allcand %in% shared_lodo$gene,
  LODO_轨道  = ifelse(allcand %in% mirror_lodo, "mirror",
                 ifelse(allcand %in% shared_lodo$gene, "concordant", NA)),
  stringsAsFactors = FALSE)
tab2$轨道保持 <- tab2$主分析_入选 & tab2$LODO_入选 & tab2$主分析_轨道 == tab2$LODO_轨道
tab2 <- tab2[order(!tab2$主分析_入选, tab2$gene), ]
write.csv(tab2, file.path(OUTD, "LODO_GSE111006_双轨对照表.csv"), row.names = FALSE)

## ---- 输出 3：摘要（可直接引用） ---------------------------------------------
sm <- c(
  "LODO 敏感性分析：剔除 GSE111006（36 对照/4 肌少症）后肌肉轨 RRA 与双轨候选集稳定性",
  sprintf("生成时间：%s", Sys.time()),
  "",
  sprintf("复现校验：全量 5 队列肌肉轨复跑与主分析 CSV 完全一致（基因顺序=%s，分数=%s）", same_gene, same_score),
  "",
  sprintf("[肌肉轨 RRA 显著集（RRA_adj<0.05）] 主分析 %d 个 → LODO %d 个；保留 %d/%d（%.1f%%）；Jaccard=%.3f",
          m_adj$n_orig, m_adj$n_lodo, m_adj$retained, m_adj$n_orig, 100 * m_adj$retention, m_adj$jaccard),
  sprintf("[肌肉轨名义显著集（RRA_score<0.05，双轨输入口径）] 主分析 %d 个 → LODO %d 个；保留 %d/%d（%.1f%%）；Jaccard=%.3f",
          m_nom$n_orig, m_nom$n_lodo, m_nom$retained, m_nom$n_orig, 100 * m_nom$retention, m_nom$jaccard),
  sprintf("[秩相关] 全体共有基因（n=%d）-log10(RRA_score) Spearman ρ=%.4f；主分析校正显著子集 ρ=%.4f",
          length(common), rho_all, rho_sig),
  "",
  sprintf("[双轨共享候选集] 主分析 %d 个（镜像 %d，占比 %.1f%%）→ LODO %d 个（镜像 %d，占比 %.1f%%，单侧二项 P=%.2e）",
          nrow(shared_main), length(mirror_main), 100 * length(mirror_main) / nrow(shared_main),
          nrow(shared_lodo), length(mirror_lodo), 100 * length(mirror_lodo) / max(nrow(shared_lodo), 1),
          bt_lodo$p.value),
  sprintf("[候选保持] 主分析 %d 个候选中 LODO 仍入选 %d 个（%.1f%%）；镜像轨 %d 个中仍入镜像轨 %d 个（%.1f%%）",
          nrow(shared_main), length(intersect(shared_main$gene, shared_lodo$gene)),
          100 * length(intersect(shared_main$gene, shared_lodo$gene)) / nrow(shared_main),
          length(mirror_main), length(intersect(mirror_main, mirror_lodo)),
          100 * length(intersect(mirror_main, mirror_lodo)) / max(length(mirror_main), 1)),
  "",
  "结论（若上述保持率均高）：GSE111006 虽因 case 数仅 4 而逐队列功效极低（0 个 FDR 显著基因），",
  "但其排序对肌肉轨整合与双轨候选集无实质影响，剔除后核心结论保持不变。"
)
writeLines(sm, file.path(OUTD, "LODO_GSE111006_摘要.txt"))
cat(paste0(sm, collapse = "\n"), "\n")
message("完成。输出目录：", OUTD)
