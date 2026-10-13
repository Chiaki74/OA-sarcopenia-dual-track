## ============================================================================
## 38_fig8b_pseudobulk重绘_20261006.R
## 用途：重绘 figures/Fig8b_pseudobulk_oLT_vs_MT.pdf —— 原版标题过长被裁断
##   （"GSE152805 = within-donor oLT vs MT (paired n=..." 在右缘截断），
##   且原图例把正向条统标 "higher in MT"，对 GSE169454/GSE167186 行不成立。
## 改动：
##   1) 标题缩短为一行；三对比的方向定义移入副标题（两行，不裁断）
##   2) 图例改为中性 "Delta > 0" / "Delta < 0"，方向语义由副标题按对比分别定义：
##        GSE152805: Delta = MT - oLT  （paired n=3, descriptive only）
##        GSE169454: Delta = OA - normal（unpaired n=7, exploratory）
##        GSE167186: Delta = Old - Young（aging proxy, BH 校正后无显著）
##      ※ GSE152805 的 median_A=MT 系由原图图例（正向=higher in MT）与 CSV
##        delta=median_A-median_B 反推；若与 moduleG 脚本实际设定相反，
##        只需把下方 DIR_152805 翻转即可
##   3) 行序、行标签（FEATURE | GSE | contrast）与配色语义与原图一致
## 输入：results/moduleG/G2/moduleG_G2_pseudobulk_stats.csv（自动多路径查找）
## 输出：figures/Fig8b_pseudobulk_oLT_vs_MT_v2_20261006.pdf / .png（300dpi）
## ============================================================================

PROJECT <- "C:/Users/凉月/Desktop/ENG论文"
suppressPackageStartupMessages(library(ggplot2))

TERRA <- "#C07763"; SAGE <- "#7A9E7E"; INK <- "#3a3532"

## ---- 输入（多路径容错） ------------------------------------------------------
cand <- c(file.path(PROJECT, "results/moduleG/G2/moduleG_G2_pseudobulk_stats.csv"),
          file.path(PROJECT, "results/moduleG/moduleG_G2_pseudobulk_stats.csv"))
hit <- cand[file.exists(cand)]
if (length(hit) == 0) {
  hit <- list.files(file.path(PROJECT, "results"),
                    pattern = "moduleG_G2_pseudobulk_stats\\.csv$",
                    recursive = TRUE, full.names = TRUE)
}
stopifnot("找不到 moduleG_G2_pseudobulk_stats.csv" = length(hit) >= 1)
d <- read.csv(hit[1], stringsAsFactors = FALSE)
message("读取: ", hit[1], "（", nrow(d), " 行）")

## ---- 仅 10 标志物行（与原图一致，track_score 3 行不入图） ---------------------
m <- d[d$feature_type == "marker", ]
stopifnot(nrow(m) == 30)
m$row_lab <- sprintf("%s | %s | %s", m$feature, m$gse, m$comparison)
m <- m[order(m$log2FC_or_delta, decreasing = TRUE), ]
m$row_lab <- factor(m$row_lab, levels = rev(m$row_lab))   # 顶部 = 最大正值
m$dir <- ifelse(m$log2FC_or_delta >= 0, "pos", "neg")

## ---- 方向语义开关（如与 moduleG 设定相反，翻转此值即可） ----------------------
DIR_152805 <- "Delta = MT - oLT"   # 反推自原图图例 + CSV；翻转则改为 "Delta = oLT - MT"

sub <- sprintf(paste0(
  "GSE152805: %s (paired n=3, descriptive);  GSE169454: Delta = OA - normal (n=7)\n",
  "GSE167186: Delta = Old - Young (aging proxy; BH none < 0.05);  ",
  "Delta = donor-level pseudobulk median difference"), DIR_152805)

p <- ggplot(m, aes(log2FC_or_delta, row_lab, fill = dir)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.4) +
  scale_fill_manual(values = c(pos = TERRA, neg = SAGE),
                    labels = c(pos = "Delta > 0", neg = "Delta < 0"), name = NULL) +
  labs(title = "Fig 8b. Pseudobulk marker contrasts",
       subtitle = sub,
       x = "log2FC / delta", y = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 7, colour = "grey30"),
        axis.text.y = element_text(size = 6.5, colour = INK))

OUT <- file.path(PROJECT, "figures", "Fig8b_pseudobulk_oLT_vs_MT_v2_20261006")
ggsave(paste0(OUT, ".pdf"), p, width = 8.5, height = 7)
ggsave(paste0(OUT, ".png"), p, width = 8.5, height = 7, dpi = 300)
message("写出: ", OUT, ".pdf / .png")
message("自检锚点：行数 30；顶部应为 STEAP1|GSE169454 (+1.88)，底部应为 GADD45A|GSE169454 (-2.32)")
