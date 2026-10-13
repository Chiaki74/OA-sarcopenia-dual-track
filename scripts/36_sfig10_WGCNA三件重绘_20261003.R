## ============================================================================
## 36_sfig10_WGCNA三件重绘_20261003.R
## 目的：S-Fig 10 三件规范版重绘——
##   ① 修复 plot_fig3_v2 的 bug：read_module_file() 未带入 GS 列，
##      d0$GS 为 NULL → n_hc 静默为 0（副标题误显示 "0 high-confidence genes"）。
##      真值已核验：OA 2,119 / muscle 304（layer3 CSV 直算）。
##   ② 标题改 S-Fig 10a/b/c；肌肉侧 21 个模块 x 轴标签改 45° 防重叠。
## 输入：results/moduleB/layer3_WGCNA_modules_{OA,muscle}.csv（列 gene/module/kME_own/GS）
##       results/moduleB/CORE_concordant_genes.csv、CORE_mirror_genes.csv（c 图用）
## 输出：figures/SFig10a_WGCNA_OA_modules_20261003.pdf
##       figures/SFig10b_WGCNA_muscle_modules_20261003.pdf
##       figures/SFig10c_candidate_module_distribution_20261003.pdf（+png）
## 运行：Rscript 36_sfig10_WGCNA三件重绘_20261003.R
## ============================================================================

suppressPackageStartupMessages(library(ggplot2))

MB  <- file.path("results", "moduleB")
KEY <- list(OA = c("blue","brown","green","magenta","pink","red","tan","yellow"),
            muscle = c("brown","red","salmon","yellow"))   # moduleB 日志封板口径
KME_CUT <- 0.7; GS_CUT <- 0.3
SAGE <- "#7A9E7E"; TERRA <- "#C07763"; TAUPE <- "#8a7f76"; INK <- "#3a3532"

theme_fig <- function() theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        panel.grid.minor = element_blank())

dir.create("figures", showWarnings = FALSE)

## 修复点：读入保留 GS 列
read_mod <- function(tag) {
  f <- file.path(MB, paste0("layer3_WGCNA_modules_", tag, ".csv"))
  stopifnot(file.exists(f))
  d <- read.csv(f, check.names = FALSE)
  stopifnot(all(c("gene","module","kME_own","GS") %in% colnames(d)))
  d
}

## ---- a/b：模块规模条形图 ----------------------------------------------------
for (tag in c("OA", "muscle")) {
  d0 <- read_mod(tag)
  key <- KEY[[tag]]
  n_hc <- sum(d0$module %in% key & abs(d0$kME_own) >= KME_CUT & d0$GS >= GS_CUT)
  d <- d0[d0$module != "grey", ]
  tab <- as.data.frame(table(d$module), stringsAsFactors = FALSE)
  colnames(tab) <- c("module", "n"); tab <- tab[order(tab$n, decreasing = TRUE), ]
  tab$module <- factor(tab$module, levels = tab$module)
  tab$is_key <- tab$module %in% key

  mods <- levels(tab$module)
  ok <- mods %in% grDevices::colours()
  cols <- ifelse(mods %in% key, mods, "#DDD6C8")
  cols[mods %in% key & !ok] <- TAUPE
  cols <- setNames(cols, mods)

  suffix <- ifelse(tag == "OA", "a", "b")
  p <- ggplot(tab, aes(module, n, fill = module, alpha = is_key)) +
    geom_col(width = .72) +
    geom_text(aes(label = n), angle = 90, hjust = -0.15, size = 2.6, colour = INK) +
    scale_fill_manual(values = cols, guide = "none") +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = .45), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, .12))) +
    labs(title = sprintf("S-Fig 10%s. WGCNA modules — %s (confidence annotation layer)",
                         suffix, tag),
         subtitle = sprintf("%d key modules (solid) with %d high-confidence genes (|kME|>=%.1f & GS>=%.1f); %d non-grey modules total",
                            length(key), n_hc, KME_CUT, GS_CUT, nrow(tab)),
         x = NULL, y = "genes per module") +
    theme_fig() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5))
  out <- file.path("figures", sprintf("SFig10%s_WGCNA_%s_modules_20261003.pdf", suffix, tag))
  ggsave(out, p, width = 6.85, height = 4.2)
  ggsave(sub("\\.pdf$", ".png", out), p, width = 6.85, height = 4.2, dpi = 300)
  message(sprintf("✓ S-Fig 10%s 完成：%s 高置信基因 %d（应%s）",
                  suffix, tag, n_hc, ifelse(tag=="OA", "为 2,119", "为 304")))
}

## ---- c：双轨候选模块分布气泡图 ----------------------------------------------
read_genes <- function(label) {
  f <- file.path(MB, paste0("CORE_", label, "_genes.csv"))
  if (!file.exists(f)) return(character(0))
  d <- read.csv(f, check.names = FALSE)
  unique(d[[intersect(c("gene","Gene","x"), colnames(d))[1]]])
}
g1 <- read_genes("concordant"); g2 <- read_genes("mirror")
cand <- c(setNames(rep("Concordant", length(g1)), g1),
          setNames(rep("Mirror", length(g2)), g2))
stopifnot(length(cand) > 0)

hit_list <- list()
for (tag in c("OA", "muscle")) {
  d <- read_mod(tag)
  d <- d[d$gene %in% names(cand) & d$module != "grey", ]
  if (nrow(d) > 0) hit_list[[tag]] <- data.frame(dataset = tag, gene = d$gene,
                                                 module = d$module, track = cand[d$gene])
}
hh <- do.call(rbind, hit_list)
cnt <- as.data.frame(table(hh$dataset, hh$module, hh$track), stringsAsFactors = FALSE)
colnames(cnt) <- c("dataset", "module", "track", "n"); cnt <- cnt[cnt$n > 0, ]

pc <- ggplot(cnt, aes(track, module)) +
  geom_point(aes(size = n, colour = track), alpha = .9) +
  facet_wrap(~ dataset, scales = "free_y") +
  scale_colour_manual(values = c(Concordant = SAGE, Mirror = TERRA), guide = "none") +
  scale_size_continuous(range = c(2, 8), name = "genes") +
  labs(title = "S-Fig 10c. Candidate genes across WGCNA modules (confidence annotation)",
       subtitle = "Concordant vs mirror track membership per module; grey excluded",
       x = NULL, y = "module") +
  theme_fig()
out3 <- file.path("figures", "SFig10c_candidate_module_distribution_20261003.pdf")
ggsave(out3, pc, width = 8, height = 5)
ggsave(sub("\\.pdf$", ".png", out3), pc, width = 8, height = 5, dpi = 300)
message("✓ S-Fig 10c 完成（候选命中模块 ", nrow(cnt), " 组）")
