## ============================================================================
## 37_补充图拼版_20261006.R
## 用途：把多文件补充图拼为单页 PDF（174 mm 版式，600 dpi）
##   S-Fig 5  = 9 个 boxplot（行=方法 ssGSEA/xCell/CIBERSORT，列=层级 肌肉/血液/aging）→ (a)-(i)
##   S-Fig 6  = 9 个相关热图，同版式
##   S-Fig 7  = 3 个 signed 热图竖向堆叠 → (a) ssGSEA (b) xCell (c) CIBERSORT
##   S-Fig 9  = 3 个 CellChat diffInteraction + 面板标题 + 红/蓝边说明
##   S-Fig 12 = 8 个 Fig8 文件 + 18 个 UCell 文件 → (a)-(e) 分节
## 输入（自动检查，缺件会列出清单后停止）：
##   figures/moduleF_redraw_20261003/SFig5_*_20261003.pdf  ×9   （脚本 34 产物）
##   figures/moduleF_redraw_20261003/SFig6_*_20261003.pdf  ×9
##   figures/moduleF_redraw_20261003/SFig7{a,b,c}_*_20261003.pdf ×3
##   results/moduleK/moduleK_GSE{169454,152805,167186}_diffInteraction.pdf ×3
##   figures/Fig8_umap_clusters_GSE*.pdf ×3、Fig8_umap_markers_GSE*.pdf ×3、
##   figures/Fig8a_marker_localization_dotplot.pdf、figures/Fig8b_pseudobulk_oLT_vs_MT.pdf
##   results/moduleG/G2/<GSE>/<GSE>_UCell_*_{FeaturePlot,VlnPlot_*}.pdf ×18
## 输出：figures/assembled_20261006/SFig{5,6,7,9,12}_assembled_20261006.pdf
##        + 同名 _preview.png（150 dpi，便于快速目检/回传）
## 依赖：magick、pdftools、png（pdftools 内置 poppler，Windows 无需装 Ghostscript）
## ============================================================================

PROJECT <- "C:/Users/凉月/Desktop/ENG论文"   # 项目根目录
message("======================================================")
message("== 脚本版本：37b_20261006_SFig9注释修正版 ==")
message("== 若本行未出现，说明跑的是旧文件！ ==")
message("======================================================")
FIGDIR  <- file.path(PROJECT, "figures")
REDRAW  <- file.path(FIGDIR, "moduleF_redraw_20261003")
G2DIR   <- file.path(PROJECT, "results", "moduleG", "G2")
KDIR    <- file.path(PROJECT, "results", "moduleK")
OUTDIR  <- file.path(FIGDIR, "assembled_20261006")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

for (pkg in c("magick", "pdftools", "png")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(magick); library(pdftools); library(png)

DPI      <- 600                       # 刊规组合彩图 ≥600 dpi
MM2PX    <- DPI / 25.4
PAGE_W   <- round(174 * MM2PX)        # 4110 px
MARGIN   <- round(3 * MM2PX)          # 页边距 3 mm
CONTENTW <- PAGE_W - 2 * MARGIN
GAP      <- round(2 * MM2PX)          # 面板间隙 2 mm
LETTER_SZ <- round(10 / 72 * DPI)     # 面板字母 ≈10 pt
BAR_SZ    <- round(9 / 72 * DPI)      # 分节标题 ≈9 pt

## ---- 工具函数 ---------------------------------------------------------------
## 读 PDF 首页 → 白底压平 → 缩放到目标像素宽
read_panel <- function(path, target_w_px) {
  bmp <- pdftools::pdf_render_page(path, page = 1, dpi = 600, antialias = TRUE)
  tmp <- tempfile(fileext = ".png")
  png::writePNG(bmp, tmp)
  img <- magick::image_read(tmp)
  unlink(tmp)
  img <- magick::image_flatten(magick::image_background(img, "white"))
  magick::image_resize(img, paste0(target_w_px, "x"))
}

## 面板字母：以置顶横条形式加在面板上方（不遮挡面板内标题）
add_letter <- function(img, letter) {
  bar <- label_bar(paste0("(", letter, ")"), w_px = magick::image_info(img)$width,
                   size = LETTER_SZ)
  magick::image_append(c(bar, img), stack = TRUE)
}

## 文字横条（分节标题 / 注释）
label_bar <- function(text, w_px = CONTENTW, size = BAR_SZ, weight = 700,
                      gravity = "west", height_px = NULL) {
  h <- if (is.null(height_px)) round(size * 1.7) else height_px
  img <- magick::image_blank(w_px, h, color = "white")
  magick::image_annotate(img, text, gravity = gravity,
                         location = paste0("+", round(0.5 * MM2PX), "+0"),
                         size = size, font = "sans", weight = weight, color = "black")
}

## 横向拼一行（面板间加白缝）
hstack <- function(imgs, gap = GAP) {
  imgs <- lapply(imgs, function(im) magick::image_border(im, "white",
        paste0(gap, "x", round(gap / 2))))
  magick::image_append(do.call(c, imgs), stack = FALSE)
}

## 竖向堆叠并贴到 174mm 白画布
to_page <- function(imgs, gap = GAP) {
  body <- magick::image_append(do.call(c, lapply(imgs, function(im)
    magick::image_border(im, "white", paste0("0x", round(gap / 2))))), stack = TRUE)
  info <- magick::image_info(body)
  canvas <- magick::image_blank(PAGE_W, info$height + 2 * MARGIN, color = "white")
  magick::image_composite(canvas, body, offset = paste0("+", MARGIN, "+", MARGIN))
}

## 保存 PDF（600 dpi）+ 150dpi 预览 PNG
save_page <- function(canvas, name) {
  pdff <- file.path(OUTDIR, paste0(name, ".pdf"))
  pngf <- file.path(OUTDIR, paste0(name, "_preview.png"))
  magick::image_write(canvas, pdff, format = "pdf", density = DPI)
  prev <- magick::image_resize(canvas, paste0(round(174 / 25.4 * 150), "x"))
  magick::image_write(prev, pngf, format = "png", density = 150)
  message("  写出: ", basename(pdff), "  (", round(magick::image_info(canvas)$width / MM2PX), " x ",
          round(magick::image_info(canvas)$height / MM2PX), " mm)")
}

## ---- 输入文件清单与存在性检查 ------------------------------------------------
METHODS <- c("ssgsea", "xcell", "cibersort")
TIERS3  <- c("Muscle_sarcopenia", "Blood_exploratory", "Aging_proxy")

f_sfig5 <- as.vector(t(outer(METHODS, TIERS3, function(m, t)
  file.path(REDRAW, sprintf("SFig5_%s_boxplot_%s_20261003.pdf", m, t)))))
f_sfig6 <- as.vector(t(outer(METHODS, TIERS3, function(m, t)
  file.path(REDRAW, sprintf("SFig6_%s_corr_heatmap_%s_20261003.pdf", m, t)))))
f_sfig7 <- file.path(REDRAW, sprintf("SFig7%s_%s_signed_heatmap_20261003.pdf",
                                     c("a", "b", "c"), METHODS))
f_sfig9 <- file.path(KDIR, sprintf("moduleK_GSE%s_diffInteraction.pdf",
                                   c("169454", "152805", "167186")))
GSE3 <- c("GSE152805", "GSE169454", "GSE167186")
f_fig8b_v2 <- file.path(FIGDIR, "Fig8b_pseudobulk_oLT_vs_MT_v2_20261006.pdf")
f_fig8 <- c(file.path(FIGDIR, sprintf("Fig8_umap_clusters_%s.pdf", GSE3)),
            file.path(FIGDIR, sprintf("Fig8_umap_markers_%s.pdf", GSE3)),
            file.path(FIGDIR, "Fig8a_marker_localization_dotplot.pdf"),
            if (file.exists(f_fig8b_v2)) f_fig8b_v2 else
              file.path(FIGDIR, "Fig8b_pseudobulk_oLT_vs_MT.pdf"))
SIGS <- c("concordant", "mirror_OAup_muscleDown", "mirror_OAdown_muscleUp")

## VlnPlot 后缀按数据集不同（tissue / condition），用通配解析真实文件名
resolve_vln <- function() {
  unlist(lapply(GSE3, function(g) unlist(lapply(SIGS, function(s) {
    hit <- list.files(file.path(G2DIR, g),
                      pattern = sprintf("^%s_UCell_%s_VlnPlot_.*\\.pdf$", g, s),
                      full.names = TRUE)
    if (length(hit) != 1) return(NA_character_)
    hit
  }))))
}
f_ucell <- unlist(lapply(GSE3, function(g) unlist(lapply(SIGS, function(s)
  c(file.path(G2DIR, g, sprintf("%s_UCell_%s_FeaturePlot.pdf", g, s)))))))
f_vln <- resolve_vln()

all_files <- c(f_sfig5, f_sfig6, f_sfig7, f_sfig9, f_fig8, f_ucell, f_vln)
missing <- all_files[!file.exists(all_files) | is.na(all_files)]
if (length(missing) > 0) {
  message("!! 缺少以下 ", length(missing), " 个输入文件，请先确认路径/文件名：")
  for (x in missing) message("   - ", x)
  stop("输入文件不齐，拼版中止。")
}
message("== 输入检查通过：S-Fig5 ×9、S-Fig6 ×9、S-Fig7 ×3、S-Fig9 ×3、Fig8 ×8、UCell ×18 ==")

## ---- S-Fig 5：3×3 箱线图网格（行=方法，列=层级；a–i 行优先） ------------------
message("== 拼 S-Fig 5 ==")
cell_w <- round((CONTENTW - 2 * GAP) / 3)
letters9 <- letters[1:9]
rows5 <- lapply(seq_along(METHODS), function(i) {
  panels <- lapply(1:3, function(j) {
    k <- (i - 1) * 3 + j
    add_letter(read_panel(f_sfig5[k], cell_w), letters9[k])
  })
  hstack(panels)
})
save_page(to_page(rows5), "SFig5_assembled_20261006")

## ---- S-Fig 6：3×3 相关热图网格（同版式） --------------------------------------
message("== 拼 S-Fig 6 ==")
rows6 <- lapply(seq_along(METHODS), function(i) {
  panels <- lapply(1:3, function(j) {
    k <- (i - 1) * 3 + j
    add_letter(read_panel(f_sfig6[k], cell_w), letters9[k])
  })
  hstack(panels)
})
save_page(to_page(rows6), "SFig6_assembled_20261006")

## ---- S-Fig 7：3 张 signed 热图竖向堆叠 ---------------------------------------
message("== 拼 S-Fig 7 ==")
rows7 <- lapply(1:3, function(i) add_letter(read_panel(f_sfig7[i], CONTENTW), letters[i]))
save_page(to_page(rows7), "SFig7_assembled_20261006")

## ---- S-Fig 9：CellChat ×3 ----------------------------------------------------
message("== 拼 S-Fig 9 ==")
cap9 <- c("(a) GSE169454 — OA vs Normal (chondrocyte scRNA-seq, main tier)",
          "(b) GSE152805 — MT vs oLT (within-donor, descriptive only)",
          "(c) GSE167186 — Old vs Young (aging proxy)")
half9 <- round((CONTENTW - GAP) / 2)
p9a <- read_panel(f_sfig9[1], round(CONTENTW * 0.72))
p9b <- read_panel(f_sfig9[2], half9)
p9c <- read_panel(f_sfig9[3], half9)
col_b <- magick::image_append(c(label_bar(cap9[2], w_px = half9, size = round(BAR_SZ * 0.9)),
                                p9b), stack = TRUE)
col_c <- magick::image_append(c(label_bar(cap9[3], w_px = half9, size = round(BAR_SZ * 0.9)),
                                p9c), stack = TRUE)
rows9 <- list(
  label_bar(cap9[1]),
  p9a,
  hstack(list(col_b, col_c)),
  label_bar(paste0("Note: red = stronger, blue = weaker (2nd group vs 1st, CellChat differential network); ",
                   "edge width = strength."),
            size = round(BAR_SZ * 0.85), weight = 400)
)
save_page(to_page(rows9), "SFig9_assembled_20261006")

## ---- S-Fig 12：8 Fig8 + 18 UCell，分节 (a)-(e) --------------------------------
message("== 拼 S-Fig 12（面板多，需时稍长）==")
third_w <- round((CONTENTW - 2 * GAP) / 3)
half_w  <- round((CONTENTW - GAP) / 2)
sig_lab <- c("concordant", "mirror: OA-up / muscle-down", "mirror: OA-down / muscle-up")
gse_lab <- c(GSE152805 = "GSE152805 — OA cartilage scRNA-seq (VlnPlot by tissue)",
             GSE169454 = "GSE169454 — chondrocyte scRNA-seq (VlnPlot by condition)",
             GSE167186 = "GSE167186 — muscle snRNA-seq, aging proxy (VlnPlot by condition)")

rows12 <- list(label_bar("(a) Harmony-integrated UMAP clusters"))
rows12 <- c(rows12, list(hstack(lapply(f_fig8[1:3], read_panel, target_w_px = third_w))))
rows12 <- c(rows12, list(label_bar("(b) FeaturePlot of the 10 biomarkers")))
rows12 <- c(rows12, list(hstack(lapply(f_fig8[4:6], read_panel, target_w_px = third_w))))
rows12 <- c(rows12, list(label_bar("(c) Biomarker × cell-type localization dotplot")))
rows12 <- c(rows12, list(read_panel(f_fig8[7], CONTENTW)))
rows12 <- c(rows12, list(label_bar("(d) Donor-level pseudobulk comparison")))
rows12 <- c(rows12, list(read_panel(f_fig8[8], CONTENTW)))
rows12 <- c(rows12, list(label_bar("(e) UCell module scores of the three gene-set tracks")))

for (gi in seq_along(GSE3)) {
  g <- GSE3[gi]
  rows12 <- c(rows12, list(label_bar(gse_lab[g], size = round(BAR_SZ * 0.95))))
  hdr <- hstack(lapply(sig_lab, label_bar, w_px = third_w,
                       size = round(BAR_SZ * 0.8), weight = 400, gravity = "center"))
  rows12 <- c(rows12, list(hdr))
  fp <- lapply(SIGS, function(s)
    read_panel(file.path(G2DIR, g, sprintf("%s_UCell_%s_FeaturePlot.pdf", g, s)), third_w))
  rows12 <- c(rows12, list(hstack(fp)))
  vp <- lapply(SIGS, function(s) {
    hit <- list.files(file.path(G2DIR, g),
                      pattern = sprintf("^%s_UCell_%s_VlnPlot_.*\\.pdf$", g, s),
                      full.names = TRUE)
    read_panel(hit[1], third_w)
  })
  rows12 <- c(rows12, list(hstack(vp)))
}
save_page(to_page(rows12), "SFig12_assembled_20261006")

message("\n== 全部拼版完成，产物在 ", OUTDIR, " ==")
message("自检锚点：S-Fig5/6 各 9 面板 (a)-(i)；S-Fig7 = a ssGSEA / b xCell / c CIBERSORT；",
        "S-Fig9 = a GSE169454 / b GSE152805 / c GSE167186；S-Fig12 = (a)-(e) 分节、UCell 18 件。")
