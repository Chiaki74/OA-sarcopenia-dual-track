## ============================================================================
## 35_sfig13_转化拓展重绘_20261003.R
## 目的：S-Fig 13 两件重绘——
##   (a) Fig10a 实图损坏（标题称 10/10 命中、条形却全标 not matched），按
##       moduleI_druggable_overlap.csv 重绘 markers_only 层比对结果；
##   (b) Fig10b 标题与底部 "none passed FDR" 注释被裁断，一并重绘。
## 口径：markers_only 层 = 仅 10 个标志物；命中 = 可成药并集（Finan ∪ DGIdb）。
##       图注预期：5/10 命中（NELL1、STEAP1、FBLN5、EGFR、CBR3），跑后请核对。
## 输入：results/moduleI/moduleI_druggable_overlap.csv
##       results/moduleI/moduleI_dsigdb_enrichment.csv
## 输出：figures/SFig13a_druggable_overlap_20261003.pdf/.png
##       figures/SFig13b_DSigDB_candidates_20261003.pdf/.png
## 运行：Rscript 35_sfig13_转化拓展重绘_20261003.R
## ============================================================================

suppressPackageStartupMessages(library(ggplot2))

MARKERS10 <- c("NELL1","STEAP1","GADD45A","FBLN5","HEMK1",
               "GDE1","RNF14","BNIP3","EGFR","CBR3")
MI <- file.path("results", "moduleI")
SAGE  <- "#7A9E7E"; TERRA <- "#C07763"; TAUPE <- "#8a7f76"

theme_fig <- function() theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        panel.grid.minor = element_blank())

dir.create("figures", showWarnings = FALSE)

## ---- (a) 可成药并集比对 -----------------------------------------------------
f_ov <- file.path(MI, "moduleI_druggable_overlap.csv")
stopifnot(file.exists(f_ov))
ov <- read.csv(f_ov, check.names = FALSE)
ov <- ov[ov$gene %in% MARKERS10, ]                 # markers_only 层
stopifnot(nrow(ov) == 10)
ov$hit <- ov$in_universe %in% TRUE
n_hit <- sum(ov$hit)
message(sprintf("(a) markers_only 命中 %d/10：%s",
                n_hit, paste(ov$gene[ov$hit], collapse = ", ")))
message("    图注预期 5/10 = NELL1/STEAP1/FBLN5/EGFR/CBR3 —— 不符请回传控制台输出")

ov$gene <- factor(ov$gene, levels = rev(MARKERS10))
pa <- ggplot(ov, aes(1, gene, fill = hit)) +
  geom_col(width = .7, orientation = "y") +
  geom_text(data = ov[ov$hit, ],  aes(0.97, gene, label = "druggable"),
            hjust = 1, size = 2.8, colour = "white") +
  geom_text(data = ov[!ov$hit, ], aes(0.03, gene, label = "not matched"),
            hjust = 0, size = 2.8, colour = "grey35") +
  scale_fill_manual(values = c(`TRUE` = SAGE, `FALSE` = "grey80"),
                    labels = c("not matched", "in druggable universe"), name = NULL) +
  scale_x_continuous(limits = c(0, 1.12), expand = c(0, 0)) +
  labs(title = "S-Fig 13a. Druggable-genome overlap (markers-only)",
       subtitle = sprintf("%d/10 biomarkers in Finan & DGIdb druggable union (n = 5,928)", n_hit),
       x = NULL, y = NULL) +
  theme_fig() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom")
ggsave(file.path("figures", "SFig13a_druggable_overlap_20261003.pdf"), pa,
       width = 5.2, height = 3.6)
ggsave(file.path("figures", "SFig13a_druggable_overlap_20261003.png"), pa,
       width = 5.2, height = 3.6, dpi = 300)

## ---- (b) DSigDB 名义显著候选 ------------------------------------------------
f_ds <- file.path(MI, "moduleI_dsigdb_enrichment.csv")
stopifnot(file.exists(f_ds))
ds <- read.csv(f_ds, check.names = FALSE)
n_tested <- nrow(ds); n_fdr <- sum(ds$p_adjust_BH < 0.05, na.rm = TRUE)
ds <- head(ds[order(ds$p_value), ], 20)
## 截断超长药物名（旧版长名把面板挤出画布）
ds$drug_lab <- ifelse(nchar(ds$drug) > 42, paste0(substr(ds$drug, 1, 42), "…"), ds$drug)
ds$drug_lab <- factor(ds$drug_lab, levels = rev(ds$drug_lab))

pb <- ggplot(ds, aes(drug_lab, -log10(p_value), fill = overlap)) +
  geom_col(width = .75) +
  scale_fill_gradient(low = "grey85", high = SAGE, name = "overlap n") +
  coord_flip() +
  labs(title = "S-Fig 13b. DSigDB drug perturbation candidates",
       subtitle = sprintf("exploratory ORA; top 20 by nominal P; no term passed FDR (BH, %d tested)", n_tested),
       x = NULL, y = "-log10(nominal P)") +
  theme_fig() +
  theme(plot.margin = margin(6, 14, 6, 6),
        axis.text.y = element_text(size = 7))
ggsave(file.path("figures", "SFig13b_DSigDB_candidates_20261003.pdf"), pb,
       width = 6.85, height = 5.4)
ggsave(file.path("figures", "SFig13b_DSigDB_candidates_20261003.png"), pb,
       width = 6.85, height = 5.4, dpi = 300)
message(sprintf("✓ (b) 完成：检验条目 %d，FDR<0.05 条目 %d（图注口径：无条目通过 FDR）",
                n_tested, n_fdr))
message("    图注另称 4,026 条目中 561 条进入检验 —— 若 n_tested != 561 请回传控制台输出")
