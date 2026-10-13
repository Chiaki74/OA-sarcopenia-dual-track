## ===========================================================================
## 31_拼图_Fig5_小图拼接_20260916_via3.R
## via3（2026-09-23）：投稿合规升版——DPI 300→600（刊规组合彩图 ≥600dpi），
##       全部像素常数同步 ×2，版式物理尺寸不变。
##       前置：须先用 16_style_theme via3 重跑 21_fig5 via11 生成 600dpi 面板。
## ---------------------------------------------------------------------------
## 目的：把 21_fig5 via6 输出的 4 个小图拼回整版 Fig5（化零为整）
## via2 相对 via1：配合 21 via6 的统一边距（左 14pt、上 6pt）与
##       plot.title.position="plot"，切换 Fig2 29 via7 同款「图内同行盖章」：
##       废弃留白条，字母盖在每面板左上角标题同行留白区，与标题同基线。
## 版式：沿用 21 原 patchwork 设计——(A|B)/(C|D) 两行，宽比 1.7:1、高比 1.6:1：
##       上行 A_MR森林图(110×123) + B_共定位PPH4热图(64×123)，
##       下行 C_SuSiE一致性散点(87×77) + D_LocusZoom(87×77)，整版宽 174mm。
## 策略：magick 位图拼接——不重算、不开任何图形设备。
## 输入：figures/ClinRheum风格/Fig5_panel{A..D}_*.tiff（21 via6 save_cr 产物）
## 输出：figures/ClinRheum风格/Fig5_拼接整版_174mm.tiff（600dpi LZW，投稿用）
##       figures/ClinRheum风格/Fig5_拼接整版_preview.png（快速目检用）
## 依赖：magick（缺则自动安装）；16_style_theme via2（仅用其 locate_root）
## ===========================================================================

## --- 0. 环境 ----------------------------------------------------------------
if (!requireNamespace("magick", quietly = TRUE)) install.packages("magick")
library(magick)

## 主题仅需 locate_root；按实际主题文件名调整（via 递增时同步改这里）
.theme <- file.path(getwd(), "16_style_theme_共享主题_ClinRheum风格_20260915_via3.R")
if (!file.exists(.theme)) {
  .cand <- list.files(getwd(), pattern = "^16_style_theme.*\\.R$", full.names = TRUE)
  if (length(.cand) == 0) stop("未找到 16_style_theme 主题脚本，请将其放在项目根目录")
  .theme <- sort(.cand, decreasing = TRUE)[1]   # 取文件名最大（通常=最高 via）
}
source(.theme, encoding = "UTF-8")

ROOT     <- locate_root()
FIG_DIR  <- file.path(ROOT, "figures", "ClinRheum风格")
DPI      <- 600
MM2PX    <- DPI / 25.4
PAGE_W   <- round(174 * MM2PX)          # 整版宽 174mm -> 2055px

## --- 1. 读入 4 个小图（缺失即列清单报错） ------------------------------------
tags <- c(A = "A_MR森林图",
          B = "B_共定位PPH4热图",
          C = "C_SuSiE一致性散点",
          D = "D_LocusZoom")
paths <- file.path(FIG_DIR, paste0("Fig5_panel", tags, ".tiff"))
missing_f <- paths[!file.exists(paths)]
if (length(missing_f) > 0) {
  stop("[拼图] 缺少 ", length(missing_f), " 个小图，请先跑 21_fig5 via11：\n  ",
       paste(missing_f, collapse = "\n  "), call. = FALSE)
}
panels <- lapply(paths, magick::image_read)
names(panels) <- names(tags)   # paste0 不保名，必须显式命名，否则 panels$A = NULL

## --- 2. 字母盖章（via2：图内同行盖章——21 via6 已统一左 14pt 留白， --------------
##       标题从 x≈58px 起笔，其左侧 x8–33px 标题同行 y 带内无图形元素）
LAB_SIZE <- 72                               # 字号 72px@600dpi（=36px@300dpi，物理字号不变）
LAB_X    <- 16                               # 字母左缘 x（=8px@300dpi）
LAB_Y    <- 36                               # 字母 y（=18px@300dpi）
stamp_on <- function(row_img, letters, xs, y = LAB_Y) {
  for (i in seq_along(letters)) {
    row_img <- magick::image_annotate(row_img, letters[i], gravity = "northwest",
                                      location = paste0("+", xs[i], "+", y),
                                      size = LAB_SIZE, font = "sans", weight = 700,
                                      color = "black")
  }
  row_img
}

## --- 3. 行内横拼（列间白缝）→ 统一拉宽到 174mm → 行间白缝竖拼 ----------------
GAP_H <- 52                                          # 列间白缝 px（=26px@300dpi）
GAP_V <- 32                                          # 行间白缝 px（=16px@300dpi）

hcat <- function(imgs) {                             # 横向拼接，中间夹白缝
  out <- imgs[[1]]
  for (img in imgs[-1]) {
    sp <- magick::image_blank(GAP_H, magick::image_info(out)$height, color = "white")
    out <- magick::image_append(c(out, sp, img))
  }
  magick::image_resize(out, paste0(PAGE_W, "x"))     # 统一整版宽，等比缩放
}

## 双列行第二列左缘 x 坐标：按 hcat 实际缩放比折算（A/B、C/D 原始宽度 + 列缝）
.col_x <- function(imgL, imgR) {
  lw <- magick::image_info(imgL)$width; rw <- magick::image_info(imgR)$width
  round((lw + GAP_H) * PAGE_W / (lw + GAP_H + rw))
}
xB <- .col_x(panels$A, panels$B)
xD <- .col_x(panels$C, panels$D)

row1 <- hcat(list(panels$A, panels$B))                        # A 110 + B 64 mm
row2 <- hcat(list(panels$C, panels$D))                        # C 87 + D 87 mm

row1 <- stamp_on(row1, c("A", "B"), c(LAB_X, xB + LAB_X))
row2 <- stamp_on(row2, c("C", "D"), c(LAB_X, xD + LAB_X))

vcat <- function(rows) {                             # 竖向拼接，行间夹窄白缝
  out <- rows[[1]]
  for (r in rows[-1]) {
    sp <- magick::image_blank(PAGE_W, GAP_V, color = "white")
    out <- magick::image_append(c(out, sp, r), stack = TRUE)
  }
  out
}
fig5 <- vcat(list(row1, row2))

## --- 4. 导出 ------------------------------------------------------------------
fin <- magick::image_info(fig5)
out_tiff <- file.path(FIG_DIR, "Fig5_拼接整版_174mm.tiff")
out_png  <- file.path(FIG_DIR, "Fig5_拼接整版_preview.png")
magick::image_write(fig5, out_tiff, format = "tiff", compression = "LZW")
magick::image_write(fig5, out_png,  format = "png")

message("\n========== Fig5 拼接完成 ==========")
message("整版尺寸：", fin$width, " × ", fin$height, " px",
        "（", round(fin$width / MM2PX, 1), " × ", round(fin$height / MM2PX, 1), " mm @600dpi）")
message("投稿图：", out_tiff)
message("预览图：", out_png)
