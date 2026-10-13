## ============================================================================
## 34_moduleF_免疫浸润图整批重绘_20261003.R
## 目的：S-Fig 5/6/7 整批重绘 + S-Fig 11a（Fig7a verdict 热图）重设计。
## 修复旧图系统性问题：①标题模板乱码；②Y 轴标题乱码；③control/case 配色
## 与全文规范相反（统一为 control 蓝 / case 红，同 S-Fig 1）；④补 low_power
## 标注（GSE111006/GSE111010）；⑤箱线图面板内补 adj.P 标注；⑥signed 热图
## 色标改为对称发散色、修正刻度；⑦Fig7a 行标签过密不可读 → 仅展示
## verdict != none 的 97 行、按层级分面。
## 输入：results/moduleF/ 下 8 个 CSV（scores_long ×3、group_compare ×3、
##       crosscohort_direction_summary、biomarker_celltype_spearman）
## 输出：figures/moduleF_redraw_20261003/ 下 28 个 PDF（+PNG）
## 运行：Rscript 34_moduleF_免疫浸润图整批重绘_20261003.R
## ============================================================================

suppressPackageStartupMessages({library(ggplot2); library(dplyr)})

MF  <- file.path("results", "moduleF")
OUT <- file.path("figures", "moduleF_redraw_20261003")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- 统一主题与配色（control 蓝 / case 红，与 S-Fig 1 一致） -----------------
COL_CTRL <- "#7A9EB1"; COL_CASE <- "#C05746"
SAGE <- "#7A9E7E"; TERRA <- "#C07763"
CASE_LVLS <- c("OA", "sarcopenia", "aging")          # case 侧取值
METHODS <- c("ssGSEA" = "ssgsea", "xCell" = "xcell", "CIBERSORT" = "cibersort")
TIERS_EXP <- c("Muscle_sarcopenia", "Blood_exploratory", "Aging_proxy")
LOW_POWER <- c("GSE111006", "GSE111010")
TOP_N <- 6

theme_fig <- function() theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        strip.text = element_text(face = "bold", size = 8),
        panel.grid.minor = element_blank())

save_both <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  message("  ✓ ", name)
}

## ---- 读数 -------------------------------------------------------------------
read1 <- function(n) read.csv(file.path(MF, n), check.names = FALSE)
scores <- setNames(lapply(METHODS, function(m) read1(sprintf("moduleF_%s_scores_long.csv", m))), names(METHODS))
cmp    <- setNames(lapply(METHODS, function(m) read1(sprintf("moduleF_%s_group_compare.csv", m))), names(METHODS))
cc     <- read1("moduleF_crosscohort_direction_summary.csv")
sp     <- read1("moduleF_biomarker_celltype_spearman.csv")

for (m in names(METHODS)) {
  scores[[m]]$grp <- ifelse(scores[[m]]$condition %in% CASE_LVLS, "case", "control")
  scores[[m]]$grp <- factor(scores[[m]]$grp, levels = c("control", "case"))
  cmp[[m]]$gse_lab <- ifelse(cmp[[m]]$gse %in% LOW_POWER,
                             paste0(cmp[[m]]$gse, "\n(low power)"), cmp[[m]]$gse)
}

## ---- S-Fig 5：组间箱线图（3 方法 × 3 探索层 + OA_training 存档） -------------
message("== S-Fig 5 箱线图 ==")
for (m in names(METHODS)) {
  for (tr in c(TIERS_EXP, "OA_training")) {
    cm <- cmp[[m]][cmp[[m]]$tier == tr, ]
    if (nrow(cm) == 0) next
    top_cells <- cm$cell_type[order(cm$adjP)][seq_len(min(TOP_N, nrow(cm)))] |> unique()
    sc <- scores[[m]][scores[[m]]$tier == tr & scores[[m]]$cell_type %in% top_cells, ]
    if (nrow(sc) == 0) next
    lab <- unique(cm[, c("gse", "gse_lab")])
    sc$gse_lab <- lab$gse_lab[match(sc$gse, lab$gse)]
    cm_k <- cm[cm$cell_type %in% top_cells, ]
    cm_k$adjP_lab <- sprintf("adj.P=%.2g", cm_k$adjP)
    sc$cell_type <- factor(sc$cell_type, levels = top_cells)
    cm_k$cell_type <- factor(cm_k$cell_type, levels = top_cells)
    gse_order <- unique(sc$gse_lab[order(match(sc$gse, lab$gse))])
    sc$gse_lab <- factor(sc$gse_lab, levels = gse_order)
    cm_k$gse_lab <- factor(lab$gse_lab[match(cm_k$gse, lab$gse)], levels = gse_order)

    tag <- ifelse(tr %in% TIERS_EXP, "S-Fig 5", "archive")
    p <- ggplot(sc, aes(grp, score, fill = grp)) +
      geom_boxplot(outlier.size = 0.5, alpha = 0.85, width = .62) +
      geom_text(data = cm_k, aes(x = 1.5, y = Inf, label = adjP_lab),
                inherit.aes = FALSE, vjust = 1.6, size = 2.4, colour = "grey30") +
      facet_grid(cell_type ~ gse_lab, scales = "free_y") +
      scale_fill_manual(values = c(control = COL_CTRL, case = COL_CASE), name = NULL) +
      labs(title = sprintf("%s — %s", m, tr),
           subtitle = sprintf("%s; top %d cell types by within-cohort adj.P (BH); low-power flagged",
                              tag, length(top_cells)),
           x = NULL, y = sprintf("%s score", m)) +
      theme_fig() +
      theme(legend.position = "bottom",
            axis.text.x = element_text(angle = 0, size = 7.5),
            strip.text.y = element_text(size = 7.5, angle = 0),
            panel.spacing = unit(0.4, "lines"))
    n_gse <- length(unique(sc$gse)); n_ct <- length(top_cells)
    save_both(p, sprintf("SFig5_%s_boxplot_%s_20261003", tolower(METHODS[m]), tr),
              w = 1.6 * n_gse + 2.2, h = 1.35 * n_ct + 1.4)
  }
}

## ---- S-Fig 7：signed 热图（逐方法，10 队列全列） ------------------------------
message("== S-Fig 7 signed 热图 ==")
panel_letter <- c(ssGSEA = "a", xCell = "b", CIBERSORT = "c")   # 与图注 (a)ssGSEA/(b)xCell/(c)CIBERSORT 一致
for (m in names(METHODS)) {
  cm <- cmp[[m]]
  cm$signed <- sign(cm$delta) * pmin(-log10(cm$adjP), 16)
  lim <- max(abs(cm$signed), na.rm = TRUE)
  p <- ggplot(cm, aes(gse, cell_type, fill = signed)) +
    geom_tile(colour = "grey90", linewidth = 0.2) +
    scale_fill_gradient2(low = "#4C72B0", mid = "white", high = COL_CASE,
                         midpoint = 0, limits = c(-lim, lim),
                         name = "signed\n-log10(adj.P)") +
    labs(title = sprintf("S-Fig 7%s. %s signed -log10(adj.P) across cohorts",
                         panel_letter[m], m),
         subtitle = "red = higher in case, blue = lower in case; BH within cohort × method",
         x = NULL, y = NULL) +
    theme_fig() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 7.5))
  save_both(p, sprintf("SFig7%s_%s_signed_heatmap_20261003", panel_letter[m], tolower(METHODS[m])),
            w = 6.85, h = max(3.5, 0.16 * length(unique(cm$cell_type)) + 1.6))
}

## ---- S-Fig 6：标志物 × 细胞类型相关热图（3 方法 × 4 层） ----------------------
message("== S-Fig 6 相关热图 ==")
GENES10 <- c("NELL1","STEAP1","GADD45A","FBLN5","HEMK1","GDE1","RNF14","BNIP3","EGFR","CBR3")
for (m in names(METHODS)) {
  for (tr in c(TIERS_EXP, "OA_training")) {
    d <- sp[sp$method == m & sp$tier == tr, ]
    if (nrow(d) == 0) next
    med <- d %>% group_by(gene, cell_type) %>%
      summarise(med_rho = median(rho, na.rm = TRUE),
                n_sig = sum(adjP < 0.05, na.rm = TRUE),
                n_test = sum(!is.na(adjP)), .groups = "drop")
    med$gene <- factor(med$gene, levels = GENES10)
    p <- ggplot(med, aes(cell_type, gene, fill = med_rho)) +
      geom_tile(colour = "grey90", linewidth = 0.2) +
      geom_text(aes(label = ifelse(n_sig > 0, paste0(n_sig, "/", n_test), "")),
                size = 2.2, colour = "grey20") +
      scale_fill_gradient2(low = "#4C72B0", mid = "white", high = COL_CASE,
                           midpoint = 0, limits = c(-1, 1), name = "median\nrho") +
      labs(title = sprintf("%s — %s", m, tr),
           subtitle = "biomarker × cell-type Spearman; median rho across cohorts; label = n sig / n tested",
           x = NULL, y = NULL) +
      theme_fig() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
    save_both(p, sprintf("SFig6_%s_corr_heatmap_%s_20261003", tolower(METHODS[m]), tr),
              w = max(5, 0.16 * length(unique(med$cell_type)) + 2.4), h = 4.2)
  }
}

## ---- S-Fig 11a 重设计：verdict 汇总热图（仅非 none 行） ------------------------
message("== S-Fig 11a verdict 热图（重设计）==")
cc2 <- cc[cc$verdict != "none", ]
cc2$verdict <- factor(cc2$verdict,
  levels = c("replicated_up_in_case", "replicated_down_in_case", "conflicting", "exploratory"))
cc2$sig_lab <- sprintf("%d↑/%d↓", cc2$n_sig_up, cc2$n_sig_down)
## 缩短分面条标签，避免竖排条内裁切（全称见正文/图注）
tier_lab <- c(OA_training = "OA_training", Muscle_sarcopenia = "Muscle_sarc.",
              Blood_exploratory = "Blood_expl.", Aging_proxy = "Aging_proxy")
cc2$tier <- factor(tier_lab[as.character(cc2$tier)],
                   levels = c("OA_training", "Muscle_sarc.", "Blood_expl.", "Aging_proxy"))
p <- ggplot(cc2, aes(method, cell_type, fill = verdict)) +
  geom_tile(colour = "grey85", linewidth = 0.2) +
  geom_text(aes(label = sig_lab), size = 2.2, colour = "grey15") +
  facet_grid(tier ~ ., scales = "free_y", space = "free_y") +   # 按层级行数比例分配高度
  scale_fill_manual(values = c(replicated_up_in_case = COL_CASE,
                               replicated_down_in_case = "#4C72B0",
                               conflicting = "#B25D92",
                               exploratory = "#D9CFBE"),
                    labels = c("replicated up", "replicated down", "conflicting", "exploratory"),
                    name = NULL) +
  labs(title = "S-Fig 11a (redraw). Cross-cohort direction-consistency verdicts",
       subtitle = sprintf("%d non-none entries shown (of 468); replicated signals (16) all in OA_training; label = n sig up / n sig down",
                          nrow(cc2)),
       x = NULL, y = NULL) +
  theme_fig() +
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 8.5),
        strip.text.y = element_text(size = 7.5, face = "bold"),
        plot.margin = margin(6, 26, 6, 6),
        legend.position = "bottom")
save_both(p, "SFig11a_verdict_heatmap_redraw_20261003", w = 9.2,
          h = max(6, 0.13 * nrow(cc2) + 2))

message("\n== 全部完成，产物在 ", OUT, " ==")
message("自检锚点：replicated 应 16 条且全在 OA_training；CIBERSORT 2 条 = Mast cells resting(up) + T cells CD4 memory resting(down)")
