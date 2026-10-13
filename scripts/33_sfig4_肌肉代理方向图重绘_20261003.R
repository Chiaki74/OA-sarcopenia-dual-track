## ============================================================================
## 33_sfig4_肌肉代理方向图重绘_20261003.R
## 目的：S-Fig 4（旧 Fig6f）规范版重绘——修复旧版标题右侧被裁断问题，
##       标题改为 S-Fig 4 编号，副标题双口径动态计算，加宽画布。
## 口径：与 plot_fig6_v4 完全一致——consistent 列已按轨道预期编码
##       （镜像轨 expected 已翻转，核实于 2026-09-03），本脚本不改任何数值。
## 输入：results/moduleE/moduleE_muscle_mirror_direction.csv
## 输出：figures/SFig4_muscle_mirror_direction_20261003.pdf + .png
## 运行：Rscript 33_sfig4_肌肉代理方向图重绘_20261003.R
## ============================================================================

suppressPackageStartupMessages(library(ggplot2))

IN  <- file.path("results", "moduleE", "moduleE_muscle_mirror_direction.csv")
OUT <- file.path("figures", "SFig4_muscle_mirror_direction_20261003.pdf")
SAGE  <- "#7A9E7E"   # matches expectation
TERRA <- "#C07763"   # opposite to expectation

stopifnot(file.exists(IN))
d <- read.csv(IN, check.names = FALSE)

## 列识别（同 plot_fig6_v4 逻辑）
guess_col <- function(df, cands) {
  hit <- intersect(cands, colnames(df))
  if (length(hit) == 0) return(NA_character_) else return(hit[1])
}
gcol <- guess_col(d, c("gene", "Gene"))
dcol <- guess_col(d, c("same_direction", "direction_consistent", "consistent", "same_dir"))
stopifnot(!is.na(gcol), !is.na(dcol))

d$consistent <- as.logical(d[[dcol]])
d <- d[order(d$consistent, d[[gcol]]), ]
lab <- as.character(d[[gcol]])
extra <- intersect(c("cohort", "gse", "dataset", "comparison"), colnames(d))[1]
if (anyDuplicated(lab)) lab <- if (!is.na(extra)) paste(lab, "|", d[[extra]]) else make.unique(lab)
d$gene <- factor(lab, levels = unique(lab))

tcol <- intersect(c("track", "Track"), colnames(d))[1]
if (!is.na(tcol)) {
  d$track_grp <- ifelse(grepl("mirror", d[[tcol]], ignore.case = TRUE),
                        "Mirror track", "Concordant track")
} else {
  d$track_grp <- "all"
}

n_ok   <- sum(d$consistent); n_all <- nrow(d)
## 基因级统计：CSV 的 gene 列可能已预拼接 "GENE|cohort"，统计前剥离 cohort 后缀
gene_base      <- sub("\\s*\\|.*$", "", lab)
n_gene_any <- length(unique(gene_base[d$consistent]))
n_gene_all <- length(unique(gene_base))

p <- ggplot(d, aes(gene, 1, fill = consistent)) +
  geom_col(width = .75) +
  facet_grid(track_grp ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c(`TRUE` = SAGE, `FALSE` = TERRA),
                    labels = c("opposite to expectation", "matches expectation"),
                    name = "direction") +
  coord_flip() +
  labs(title = "S-Fig 4. Muscle proxy: direction vs track-specific expectation",
       subtitle = sprintf("%d/%d gene-cohort combinations match expectation; %d/%d biomarkers match in at least one cohort",
                          n_ok, n_all, n_gene_any, n_gene_all),
       x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.y = element_text(size = 8),
        strip.text.y = element_text(size = 8, face = "bold"),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        plot.margin = margin(6, 16, 6, 6))   # 右侧留白，防标题裁断

dir.create("figures", showWarnings = FALSE)
ggsave(OUT, p, width = 8.2, height = 2 + 0.28 * n_all)
ggsave(sub("\\.pdf$", ".png", OUT), p, width = 8.2, height = 2 + 0.28 * n_all, dpi = 300)
message(sprintf("✓ S-Fig 4 重绘完成 -> %s（组合级 %d/%d；基因级 %d/%d；应分别为 9/20 与 7/10）",
                OUT, n_ok, n_all, n_gene_any, n_gene_all))
