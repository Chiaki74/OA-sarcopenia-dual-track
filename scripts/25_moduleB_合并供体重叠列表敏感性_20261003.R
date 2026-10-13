# ============================================================================
# 25_moduleB_合并供体重叠列表敏感性_20261003.R
# 供体重叠校正敏感性分析：
#   GSE55457 与 GSE12021 存在供体重叠（正文 L51 已注明），RRA 零假设要求
#   各列表独立；主分析把二者作为两个独立列表，构成双重计数。全队列 LODO
#   （24 号脚本）显示：剔除其中任一者，OA 侧 511 仅剩 ~2%、候选 99→1，
#   证实 511 的相当部分显著性建立在该重复证据上。
#   本脚本将二者按 Stouffer 法合并为单个「研究级」列表（OA 侧变为 4 个
#   相互独立的列表），重跑 RRA 与双轨共享判定，检验核心结论在消除双重
#   计数后是否保持；并在合并结构下做研究级留一。
# ----------------------------------------------------------------------------
# 合并方法：对每个基因，Z = Σt_i / √k（t 为逐队列 limma t 统计量，k=1或2），
#   P = 2·Φ(-|Z|)，方向 = sign(Z)；仅出现于一个队列的基因按 k=1 处理。
#   合并表与 run_rra 输入格式兼容（gene / P.Value / logFC=Z 三列起效）。
# 输入：results/moduleB/layer1_DE_all_datasets.rds
#       results/moduleB/layer2_RRA_OA.csv、layer2_RRA_muscle.csv（主分析基准）
# 输出：results/moduleB/LODO_merged/
#   layer2_RRA_OA_mergedPair.csv        合并后 OA 侧 RRA 全表（4 列表）
#   合并列表_OA显著集对照.csv           511 与合并后显著集的逐项对照
#   合并列表_双轨候选对照.csv           99 候选逐项对照
#   合并列表_研究级LODO.csv             合并结构下 4 个研究级列表的留一矩阵
#   合并列表_摘要.txt                   可直接引用进补充材料的关键数字
# 运行（项目根目录）：
#   source("25_moduleB_合并供体重叠列表敏感性_20261003.R")
# ============================================================================

suppressPackageStartupMessages(library(RobustRankAggreg))

RRA_CUT            <- 0.05
MUSCLE_NOMINAL_CUT <- 0.05
PAIR <- c("GSE55457", "GSE12021")

locate_root <- function() {
  d <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  for (i in 1:8) {
    if (dir.exists(file.path(d, "results", "moduleB"))) return(d)
    parent <- dirname(d); if (parent == d) break; d <- parent
  }
  stop("未找到 results/moduleB/——请在项目根目录运行")
}
ROOT  <- locate_root()
MBDIR <- file.path(ROOT, "results", "moduleB")
OUTD  <- file.path(MBDIR, "LODO_merged")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

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

build_shared <- function(rra_oa_tab, rra_mu_tab) {
  s <- merge(rra_oa_tab[rra_oa_tab$RRA_adj < RRA_CUT, c("gene", "direction")],
             rra_mu_tab[rra_mu_tab$RRA_score < MUSCLE_NOMINAL_CUT, c("gene", "direction")],
             by = "gene", suffixes = c("_OA", "_muscle"))
  s$concordant <- s$direction_OA == s$direction_muscle
  s
}

## ---- 读取 -------------------------------------------------------------------
de_all <- readRDS(file.path(MBDIR, "layer1_DE_all_datasets.rds"))
de_oa <- de_all$OA; de_mu <- de_all$muscle
rra_oa_main <- read.csv(file.path(MBDIR, "layer2_RRA_OA.csv"),     check.names = FALSE)
rra_mu_main <- read.csv(file.path(MBDIR, "layer2_RRA_muscle.csv"), check.names = FALSE)
oa511 <- rra_oa_main$gene[rra_oa_main$RRA_adj < RRA_CUT]
shared_main <- build_shared(rra_oa_main, rra_mu_main)
cand99   <- shared_main$gene
mirror72 <- shared_main$gene[!shared_main$concordant]
message("基准：OA 511；候选 ", length(cand99), "（镜像 ", length(mirror72), "）")
stopifnot(all(PAIR %in% names(de_oa)))

## ---- Stouffer 合并 GSE55457 + GSE12021 ---------------------------------------
merge_pair <- function(tt1, tt2) {
  mg <- merge(tt1[, c("gene", "t")], tt2[, c("gene", "t")],
              by = "gene", all = TRUE, suffixes = c("_1", "_2"))
  t1 <- mg$t_1; t2 <- mg$t_2
  k  <- as.integer(!is.na(t1)) + as.integer(!is.na(t2))
  t1[is.na(t1)] <- 0; t2[is.na(t2)] <- 0
  Z <- (t1 + t2) / sqrt(k)
  P <- 2 * pnorm(-abs(Z))
  data.frame(gene = mg$gene, logFC = Z, P.Value = P, stringsAsFactors = FALSE)
}
tt_pair <- merge_pair(de_oa[[PAIR[1]]], de_oa[[PAIR[2]]])
message("合并列表基因数：", nrow(tt_pair),
        "（", PAIR[1], " ", nrow(de_oa[[PAIR[1]]]), "；", PAIR[2], " ", nrow(de_oa[[PAIR[2]]]), "）")

## ---- 合并结构 OA 侧 RRA（4 个独立列表） --------------------------------------
de_oa_m <- de_oa[setdiff(names(de_oa), PAIR)]
de_oa_m[["GSE55457+GSE12021"]] <- tt_pair
rra_oa_m <- run_rra(de_oa_m, "OA-合并供体重叠对")
write.csv(rra_oa_m, file.path(OUTD, "layer2_RRA_OA_mergedPair.csv"), row.names = FALSE)
oa_m_sig <- rra_oa_m$gene[rra_oa_m$RRA_adj < RRA_CUT]

## ---- 双轨候选重构（合并 OA × 主分析肌肉轨） -----------------------------------
shared_m <- build_shared(rra_oa_m, rra_mu_main)
cand_m   <- shared_m$gene
mirror_m <- shared_m$gene[!shared_m$concordant]
bt_m <- binom.test(length(mirror_m), length(cand_m), p = 0.5, alternative = "greater")

## ---- 研究级留一（合并结构，4 列表逐个剔除） ------------------------------------
lodo_m <- do.call(rbind, lapply(names(de_oa_m), function(g) {
  r <- run_rra(de_oa_m[setdiff(names(de_oa_m), g)], paste0("OA合并-no", g))
  sig <- r$gene[r$RRA_adj < RRA_CUT]
  sh  <- build_shared(r, rra_mu_main)
  mir <- sh$gene[!sh$concordant]
  p   <- if (nrow(sh) > 0) binom.test(length(mir), nrow(sh), p = 0.5, alternative = "greater")$p.value else NA_real_
  data.frame(剔除 = g, OA显著集 = length(sig),
             OA显著集_与合并结构交集 = length(intersect(oa_m_sig, sig)),
             候选数 = nrow(sh), 镜像数 = length(mir),
             镜像占比 = round(100 * length(mir) / max(nrow(sh), 1), 1),
             镜像二项P = signif(p, 3),
             候选与合并结构交集 = length(intersect(cand_m, sh$gene)),
             镜像与合并结构交集 = length(intersect(mirror_m, mir)),
             stringsAsFactors = FALSE)
}))
write.csv(lodo_m, file.path(OUTD, "合并列表_研究级LODO.csv"), row.names = FALSE)

## ---- 逐项对照表 ---------------------------------------------------------------
allg <- union(oa511, oa_m_sig)
cmp <- data.frame(
  gene = allg,
  主分析_RRA_adj  = rra_oa_main$RRA_adj[match(allg, rra_oa_main$gene)],
  主分析_511入选  = allg %in% oa511,
  合并_RRA_adj    = rra_oa_m$RRA_adj[match(allg, rra_oa_m$gene)],
  合并_入选       = allg %in% oa_m_sig,
  stringsAsFactors = FALSE)
write.csv(cmp, file.path(OUTD, "合并列表_OA显著集对照.csv"), row.names = FALSE)

allc <- union(cand99, cand_m)
cmp2 <- data.frame(
  gene = allc,
  主分析_入选 = allc %in% cand99,
  主分析_轨道 = ifelse(allc %in% mirror72, "mirror", ifelse(allc %in% cand99, "concordant", NA)),
  合并_入选   = allc %in% cand_m,
  合并_轨道   = ifelse(allc %in% mirror_m, "mirror", ifelse(allc %in% cand_m, "concordant", NA)),
  stringsAsFactors = FALSE)
write.csv(cmp2, file.path(OUTD, "合并列表_双轨候选对照.csv"), row.names = FALSE)

## ---- 摘要 ---------------------------------------------------------------------
sm <- c(
  "供体重叠校正敏感性：GSE55457+GSE12021 合并为单一研究级列表（Stouffer 法）",
  sprintf("生成时间：%s", Sys.time()),
  "",
  sprintf("[OA 侧显著集] 主分析（5 列表，含重复证据）%d 个 → 合并结构（4 独立列表）%d 个；",
          length(oa511), length(oa_m_sig)),
  sprintf("  主分析 511 中保留 %d（%.1f%%）；Jaccard=%.3f",
          length(intersect(oa511, oa_m_sig)), 100 * length(intersect(oa511, oa_m_sig)) / length(oa511),
          length(intersect(oa511, oa_m_sig)) / length(union(oa511, oa_m_sig))),
  "",
  sprintf("[双轨候选集] 主分析 %d（镜像 %d，%.1f%%）→ 合并结构 %d（镜像 %d，%.1f%%，单侧二项 P=%.2e）",
          length(cand99), length(mirror72), 100 * length(mirror72) / length(cand99),
          length(cand_m), length(mirror_m), 100 * length(mirror_m) / max(length(cand_m), 1), bt_m$p.value),
  sprintf("  99 候选中合并结构仍入选 %d（%.1f%%）；72 镜像中仍为镜像 %d（%.1f%%）",
          length(intersect(cand99, cand_m)), 100 * length(intersect(cand99, cand_m)) / length(cand99),
          length(intersect(mirror72, mirror_m)), 100 * length(intersect(mirror72, mirror_m)) / length(mirror72)),
  "",
  "[合并结构研究级留一]",
  apply(lodo_m, 1, function(r) sprintf(
    "  剔 %-22s：OA显著集 %s（与合并结构交集 %s）；候选 %s（镜像 %s，%s%%，P=%s）；候选交集 %s；镜像交集 %s",
    r["剔除"], r["OA显著集"], r["OA显著集_与合并结构交集"], r["候选数"], r["镜像数"],
    r["镜像占比"], format(as.numeric(r["镜像二项P"]), digits = 2),
    r["候选与合并结构交集"], r["镜像与合并结构交集"])),
  "",
  "判读要点：若合并结构下镜像占比仍显著>50% 且研究级留一不再出现塌缩，",
  "则镜像主导结论在消除供体重叠双重计数后成立；成员清单变动需在补充材料如实报告。"
)
writeLines(sm, file.path(OUTD, "合并列表_摘要.txt"))
cat(paste0(sm, collapse = "\n"), "\n")
message("完成。输出目录：", OUTD)
