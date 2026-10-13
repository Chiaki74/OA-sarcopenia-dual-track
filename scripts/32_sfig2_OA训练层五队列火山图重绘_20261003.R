## ============================================================================
## 32_sfig2_OA训练层五队列火山图重绘_20261003.R
## 目的：S-Fig 2 规范版重绘——仅 5 个 OA 训练队列（GSE55235/GSE55457/GSE51588/
##       GSE12021/GSE114007），替换误收 10 队列且残留 "Fig 2b" 标题的旧图。
## 口径：与正文 Methods 及 S-Fig 2 图注一致——
##       显著 = adj.P.Val < 0.05 且 |log2FC| >= 0.585（约 1.5 倍）；
##       分面计数即图注 "445–3,423" 的分析口径数字（跑后请核对落在该区间）。
## 输入：layer1_DE_GSE<号>.csv（列：gene, logFC, AveExpr, t, P.Value, adj.P.Val）
##       自动在 results/moduleB/、results/、当前目录 三处定位。
## 输出：figures/SFig2_volcano_OA_training_5cohorts_20261003.pdf（174mm）+ .png
## 运行：Rscript 32_sfig2_OA训练层五队列火山图重绘_20261003.R
## ============================================================================

suppressPackageStartupMessages(library(ggplot2))

## ---- 配置 ------------------------------------------------------------------
COHORTS <- c("GSE55235", "GSE55457", "GSE51588", "GSE12021", "GSE114007")  # 面板字母顺序即图注 (a)-(e)
LFC_CUT <- 0.585
FDR_CUT <- 0.05
COL_UP   <- "#C05746"   # 与主图 Fig 2b 一致：红 = 上调
COL_DOWN <- "#6C9A8B"   # 绿 = 下调
COL_NS   <- "grey75"
OUT_PDF  <- file.path("figures", "SFig2_volcano_OA_training_5cohorts_20261003.pdf")

find_de <- function(gse) {
  cand <- c(file.path("results", "moduleB", sprintf("layer1_DE_%s.csv", gse)),
            file.path("results", sprintf("layer1_DE_%s.csv", gse)),
            sprintf("layer1_DE_%s.csv", gse))
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0) stop("未找到 ", gse, " 的 layer1_DE 表（已查 results/moduleB、results、当前目录）")
  hit[1]
}

## ---- 读数 + 计数 ------------------------------------------------------------
panels <- lapply(COHORTS, function(gse) {
  d <- read.csv(find_de(gse), check.names = FALSE)
  d <- d[!grepl("^__", d$gene), ]                      # 去 RSEM 对齐伪行
  d <- d[!is.na(d$adj.P.Val) & !is.na(d$logFC), ]
  d$cls <- ifelse(d$adj.P.Val < FDR_CUT & d$logFC >=  LFC_CUT, "Up",
           ifelse(d$adj.P.Val < FDR_CUT & d$logFC <= -LFC_CUT, "Down", "NS"))
  n_up <- sum(d$cls == "Up"); n_dn <- sum(d$cls == "Down")
  message(sprintf("%s: up %d | down %d | 合计显著 %d", gse, n_up, n_dn, n_up + n_dn))
  d$panel <- sprintf("%s\nup %d | down %d", gse, n_up, n_dn)
  d
})
de <- do.call(rbind, panels)
de$panel <- factor(de$panel, levels = unique(de$panel))   # 固定面板顺序
de$cls   <- factor(de$cls, levels = c("Up", "Down", "NS"))

## 面板字母条：在 strip 前缀 (a)-(e)
lv <- levels(de$panel)
levels(de$panel) <- paste0("(", letters[seq_along(lv)], ") ", lv)

## ---- 绘图 -------------------------------------------------------------------
p <- ggplot(de, aes(logFC, -log10(P.Value), colour = cls)) +
  geom_point(size = 0.5, alpha = 0.6, stroke = 0) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed",
             linewidth = 0.3, colour = "grey50") +
  scale_colour_manual(values = c(Up = COL_UP, Down = COL_DOWN, NS = COL_NS),
                      breaks = c("Up", "Down"),
                      labels = c(sprintf("Up (adj.P<0.05, log2FC>%s)", LFC_CUT),
                                 sprintf("Down (adj.P<0.05, log2FC<-%s)", LFC_CUT)),
                      name = NULL) +
  facet_wrap(~ panel, ncol = 3, scales = "free") +
  labs(title = "S-Fig 2. Per-cohort volcano plots — OA training layer",
       subtitle = "Five OA training cohorts; significance: adj.P < 0.05 & |log2FC| ≥ 0.585 (~1.5-fold)",
       x = "log2 fold change", y = "-log10(P)") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 8.5),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        panel.grid.minor = element_blank())

dir.create("figures", showWarnings = FALSE)
ggsave(OUT_PDF, p, width = 6.85, height = 4.9)          # 174 mm 版式
ggsave(sub("\\.pdf$", ".png", OUT_PDF), p, width = 6.85, height = 4.9, dpi = 300)
message("✓ S-Fig 2 重绘完成 -> ", OUT_PDF)
message("!! 请核对：5 个队列显著合计数应落在 445–3,423 区间（图注口径），且与正文一致。")
