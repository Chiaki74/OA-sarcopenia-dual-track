# =============================================================================
# 22_sfig1_数据预处理QC_ClinRheum风格_20261003.R
# S-Fig 1：18 个 bulk 转录组数据集"标准化后"逐数据集 QC（表达密度 + PCA）
# -----------------------------------------------------------------------------
# 背景：
#   原 S-Fig 1（Fig2a_post_PCA.pdf）仅覆盖 6 队列、且为 2026-09-07 旧脚本产物，
#   标题里烧录了 "Fig 2a"。本脚本基于 moduleA 落盘的
#     processed/{GSE}/{GSE}_expr_gene.rds   （标准化后基因表达矩阵，行=基因 列=样本）
#     processed/{GSE}/{GSE}_pheno.csv       （列含 sample / condition / ... ）
#   以 ClinRheum 共享风格重绘全部 18 个 bulk 数据集的标准化后 QC，作为正式 S-Fig 1。
#
#   重要说明：moduleA 只落盘"标准化后"矩阵（expr_pre 标准化前未保存），
#     故本图 = "标准化后" QC（密度 + PCA）；"批次校正前后"对照由主图 Fig 2a 承担，
#     二者互补、不重复。本图与主图 Fig 2a（合并 OA 矩阵 before/after PCA）口径不同。
#
# 运行（两步）：
#   第 1 步 冒烟（默认只跑 2 个数据集，先验证结构 + 看预览）：
#       source("22_sfig1_数据预处理QC_ClinRheum风格_20261003.R")
#   第 2 步 全量（确认冒烟无误后，取消冒烟再跑一次）：
#       SMOKE_GSES <- NULL
#       source("22_sfig1_数据预处理QC_ClinRheum风格_20261003.R")
#
# 输出（FIG_DIR = results/supp_fig/，自动创建）：
#   SFig1A_postnorm_density.pdf / .tiff   —— 18 面板标准化后表达密度
#   SFig1B_postnorm_PCA.pdf     / .tiff   —— 18 面板标准化后 PCA（top 5000 高变基因）
#   SFig1_preprocessing_QC.pdf            —— 2 页合并版（第1页=密度 A，第2页=PCA B）
#   _preview_SFig1A.png / _preview_SFig1B.png  —— 预览 PNG（供快速查看）
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

## ===========================================================================
## 0. 内联 ClinRheum 主题（自包含，不依赖外部主题文件，避免路径问题）
##    色值与 16_style_theme_共享主题_ClinRheum风格 via3 完全一致
## ===========================================================================
CR_COL <- list(ctrl = "#5090F0",   # 对照蓝
               case = "#F07070",   # 疾病红（OA）
               muscle = "#F0A050", # 肌少橙
               ns   = "#C8C8C8")   # 其他/灰

theme_cr <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "sans") %+replace%
    theme(
      axis.line         = element_line(color = "black", linewidth = 0.4),
      axis.ticks        = element_line(color = "black", linewidth = 0.3),
      axis.text         = element_text(color = "black", size = rel(0.90)),
      axis.title        = element_text(color = "black", size = rel(1.00)),
      plot.title        = element_text(face = "bold", size = rel(1.10), hjust = 0),
      plot.subtitle     = element_text(size = rel(0.90), hjust = 0),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      legend.text       = element_text(size = rel(0.90)),
      legend.title      = element_text(size = rel(0.95)),
      strip.background  = element_blank(),
      strip.text        = element_text(face = "bold", size = rel(0.92)),
      panel.border      = element_blank(),
      plot.margin       = margin(4, 6, 4, 4)
    )
}

save_cr <- function(p, file, width_mm = 174, height_mm = 210, dpi = 600) {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  ggsave(paste0(file, ".pdf"),  p, width = width_mm, height = height_mm,
         units = "mm", device = cairo_pdf)
  ggsave(paste0(file, ".tiff"), p, width = width_mm, height = height_mm,
         units = "mm", dpi = dpi, compression = "lzw")
  message("[save_cr] ", basename(file), ".pdf/.tiff (", width_mm, "x", height_mm, " mm)")
}

## 防御性列名探测（用户本地列名以 moduleA 实际落盘为准）
pick_col <- function(df, candidates, label = "") {
  hit <- intersect(candidates, names(df))
  if (length(hit) > 1L) hit <- hit[1]
  if (length(hit) == 0L) {
    rx <- grep(paste(candidates, collapse = "|"), names(df), ignore.case = TRUE, value = TRUE)
    hit <- rx[1]
  }
  if (length(hit) == 0L || is.na(hit))
    stop("[", label, "] 列名探测失败；实际列：", paste(names(df), collapse = ", "))
  hit[1]
}

## ===========================================================================
## 1. 根目录定位 + 输入输出路径
## ===========================================================================
locate_root <- function(markers = c("processed", "data", "results", "GEO")) {
  d <- getwd()
  for (i in 1:8) {
    if (any(dir.exists(file.path(d, markers)))) return(d)
    p <- dirname(d); if (p == d) break; d <- p
  }
  warning("未自动定位根目录，使用当前工作目录：", getwd())
  getwd()
}
ROOT     <- locate_root()
PROC_DIR <- file.path(ROOT, "processed")
FIG_DIR  <- file.path(ROOT, "results", "supp_fig")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
message("[S-Fig1] ROOT = ", ROOT)
message("[S-Fig1] processed 目录存在：", dir.exists(PROC_DIR), "  （", PROC_DIR, "）")

## ===========================================================================
## 2. 18 个 bulk 数据集 + 轨道分组（facet 顺序按此）
## ===========================================================================
TRACK <- list(
  "OA train" = c("GSE55235", "GSE55457", "GSE12021", "GSE51588", "GSE114007"),
  "OA valid" = c("GSE169077", "GSE168505", "GSE57218", "GSE206848"),
  "OA blood" = "GSE48556",
  "Muscle"   = c("GSE1428", "GSE25941", "GSE111016", "GSE111010", "GSE111006",
                 "GSE136344", "GSE117525", "GSE144304")
)
track_of <- function(g) names(TRACK)[vapply(TRACK, function(v) g %in% v, logical(1))][1]
ALL_GSE  <- unlist(TRACK, use.names = FALSE)

## 冒烟开关：默认只跑 2 个（1 个 affy OA + 1 个肌肉 RNA-seq），确认无误后设 NULL 跑全量
if (!exists("SMOKE_GSES")) SMOKE_GSES <- c("GSE55235", "GSE111006")
GSES <- if (is.null(SMOKE_GSES)) ALL_GSE else intersect(ALL_GSE, SMOKE_GSES)
message("[S-Fig1] 运行模式：", if (is.null(SMOKE_GSES)) "★全量 18★" else paste0("冒烟 ", length(GSES), "（", paste(GSES, collapse = ", "), "）"))

## ===========================================================================
## 3. condition → Case / Control / Other（与主图 Fig 2a 的 Case/Control 口径一致）
##    moduleA standardize_condition 的可能取值全部覆盖
## ===========================================================================
map_group <- function(cond) {
  c0 <- tolower(as.character(cond))
  out <- rep("Other", length(c0))
  out[c0 %in% c("control", "young")] <- "Control"
  out[c0 %in% c("oa", "oa_preserved", "oa_affected",
                "sarcopenia", "pre_sarcopenia",
                "aging", "frailty", "metabolic_syndrome")] <- "Case"
  out[is.na(c0) | c0 == "other"] <- "Other"
  factor(out, levels = c("Control", "Case", "Other"))
}
GRP_COL <- c("Control" = CR_COL$ctrl, "Case" = CR_COL$case, "Other" = CR_COL$ns)

## ===========================================================================
## 4. 读单个数据集 → 返回 density 长表 + PCA 长表
## ===========================================================================
read_one <- function(gse) {
  f_expr  <- file.path(PROC_DIR, gse, paste0(gse, "_expr_gene.rds"))
  f_pheno <- file.path(PROC_DIR, gse, paste0(gse, "_pheno.csv"))
  if (!file.exists(f_expr))  { message("  [skip] ", gse, " 缺文件：", basename(f_expr));  return(NULL) }
  if (!file.exists(f_pheno)) { message("  [skip] ", gse, " 缺文件：", basename(f_pheno)); return(NULL) }

  m <- readRDS(f_expr)
  m <- as.matrix(m); storage.mode(m) <- "numeric"
  ph <- read.csv(f_pheno, stringsAsFactors = FALSE, check.names = FALSE)
  scol <- pick_col(ph, c("sample", "sample_id", "gsm", "geo_accession", "title"), label = paste0(gse, "-sample"))
  ccol <- pick_col(ph, c("condition", "group", "cond", "phenotype"),            label = paste0(gse, "-condition"))

  grp <- map_group(ph[[ccol]]); names(grp) <- ph[[scol]]

  common <- intersect(colnames(m), ph[[scol]])
  if (length(common) < 3L) {
    message("  [skip] ", gse, " 样本匹配 <3：colnames(expr) 与 pheno$", scol, " 对不上（示例 expr 列名：",
            paste(head(colnames(m), 3), collapse = "/"), "；pheno 样本：", paste(head(ph[[scol]], 3), collapse = "/"), "）")
    return(NULL)
  }
  m <- m[, common, drop = FALSE]; grp <- grp[common]

  message(sprintf("  [%s] %d 基因 × %d 样本 | 分组 %s",
                  gse, nrow(m), ncol(m),
                  paste0(names(table(grp)), "=", as.integer(table(grp)), collapse = " ")))

  ## ---- density：每样本一条曲线（预计算 256 点；样本>120 时抽 120） ----
  smp <- colnames(m)
  if (length(smp) > 120L) { set.seed(1); smp <- sample(smp, 120L) }
  dens <- do.call(rbind, lapply(smp, function(s) {
    d <- density(stats::na.omit(m[, s]), n = 256)
    data.frame(gse = gse, sample = s, group = grp[[s]], x = d$x, y = d$y,
               stringsAsFactors = FALSE)
  }))

  ## ---- PCA：top 5000 高变基因，center=TRUE, scale.=FALSE（与主图 Fig 2a 一致） ----
  v <- apply(m, 1, var, na.rm = TRUE); v[is.na(v)] <- 0
  m2 <- m[order(v, decreasing = TRUE)[1:min(5000, nrow(m))], , drop = FALSE]
  pc <- prcomp(t(m2), center = TRUE, scale. = FALSE)
  ve <- round(100 * summary(pc)$importance[2, 1:2], 1)
  pca <- data.frame(gse = gse, sample = colnames(m2), group = grp[colnames(m2)],
                    PC1 = pc$x[, 1], PC2 = pc$x[, 2], ve1 = ve[1], ve2 = ve[2],
                    stringsAsFactors = FALSE)

  list(dens = dens, pca = pca, n = ncol(m))
}

## ===========================================================================
## 5. 主流程：读取全部 → 合并 → 出图
## ===========================================================================
message("[S-Fig1] 开始读取 ", length(GSES), " 个数据集 …")
res <- lapply(GSES, read_one)
names(res) <- GSES
res <- res[!vapply(res, is.null, logical(1))]
if (length(res) == 0L) stop("没有可用数据集——请确认 processed/{GSE}/{GSE}_expr_gene.rds 与 _pheno.csv 存在。")
message("[S-Fig1] 成功读取 ", length(res), " 个数据集：", paste(names(res), collapse = ", "))

dens_all <- do.call(rbind, lapply(res, `[[`, "dens"))
pca_all  <- do.call(rbind, lapply(res, `[[`, "pca"))
nsamp    <- vapply(res, `[[`, integer(1), "n")

## facet 标签（含轨道 + 样本量 / PC 解释方差），按 TRACK 顺序排
trk <- vapply(names(res), track_of, character(1))
lab_d <- setNames(sprintf("%s\n%s · n=%d", names(res), trk, nsamp), names(res))
lab_p <- setNames(sprintf("%s\n%s · PC1 %.1f%% PC2 %.1f%%", names(res), trk,
                          vapply(res, function(r) r$pca$ve1[1], numeric(1)),
                          vapply(res, function(r) r$pca$ve2[1], numeric(1))), names(res))
ord <- names(res)  # res 已按 GSES（=TRACK）顺序
dens_all$lab <- factor(lab_d[dens_all$gse], levels = lab_d[ord])
pca_all$lab  <- factor(lab_p[pca_all$gse],  levels = lab_p[ord])

ncol_f <- if (length(res) <= 2L) length(res) else 3L
nrow_f <- ceiling(length(res) / ncol_f)
h_mm   <- max(70, min(234, 40 + nrow_f * 30))   # 按行数自适应高度，≤刊限 234mm

## ---- Panel A：标准化后表达密度 ---------------------------------------------
p_dens <- ggplot(dens_all, aes(x = x, y = y, group = sample, color = group)) +
  geom_line(linewidth = 0.25, alpha = 0.55, na.rm = TRUE) +
  facet_wrap(~ lab, ncol = ncol_f, scales = "free") +
  scale_color_manual(values = GRP_COL, name = "Group", drop = FALSE) +
  labs(title = "A   Expression density after normalization",
       x = "Normalized expression (log2)", y = "Density") +
  theme_cr() +
  theme(legend.position = "top",
        axis.text.x = element_text(size = rel(0.75)),
        axis.text.y = element_text(size = rel(0.75)))

## ---- Panel B：标准化后 PCA（top 5000 高变基因，点=样本，无样本名标签） -------
p_pca <- ggplot(pca_all, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 0.9, alpha = 0.85, na.rm = TRUE) +
  facet_wrap(~ lab, ncol = ncol_f, scales = "free") +
  scale_color_manual(values = GRP_COL, name = "Group", drop = FALSE) +
  labs(title = "B   PCA of post-normalization expression (top 5,000 variable genes)",
       x = "PC1", y = "PC2") +
  theme_cr() +
  theme(legend.position = "top",
        axis.text.x = element_text(size = rel(0.75)),
        axis.text.y = element_text(size = rel(0.75)))

## ---- 导出：单面板 PDF+TIFF、2 页合并 PDF、预览 PNG ---------------------------
save_cr(p_dens, file.path(FIG_DIR, "SFig1A_postnorm_density"), height_mm = h_mm)
save_cr(p_pca,  file.path(FIG_DIR, "SFig1B_postnorm_PCA"),     height_mm = h_mm)

cairo_pdf(file.path(FIG_DIR, "SFig1_preprocessing_QC.pdf"),
          width = 174 / 25.4, height = h_mm / 25.4)
print(p_dens); print(p_pca); dev.off()
message("[save_cr] SFig1_preprocessing_QC.pdf（2 页：A 密度 + B PCA）")

ggsave(file.path(FIG_DIR, "_preview_SFig1A.png"), p_dens, width = 174, height = h_mm, units = "mm", dpi = 150)
ggsave(file.path(FIG_DIR, "_preview_SFig1B.png"), p_pca,  width = 174, height = h_mm, units = "mm", dpi = 150)

message("============================================================")
if (is.null(SMOKE_GSES)) {
  message("[S-Fig1] ★全量完成★ 输出目录：", FIG_DIR)
  message("  请查看 SFig1_preprocessing_QC.pdf / _preview_SFig1*.png")
} else {
  message("[S-Fig1] 冒烟完成（", length(res), " 个数据集）。请查看 ", FIG_DIR, " 下的 _preview_SFig1*.png")
  message("  确认无误后跑全量 18：")
  message("    SMOKE_GSES <- NULL")
  message("    source(\"22_sfig1_数据预处理QC_ClinRheum风格_20261003.R\")")
}
message("============================================================")
