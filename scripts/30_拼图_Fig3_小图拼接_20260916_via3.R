## ===========================================================================
## 30_拼图_Fig3_小图拼接_20260916_via3.R
## via3（2026-09-23）：投稿合规升版——DPI 300→600（刊规组合彩图 ≥600dpi），
##       全部像素常数同步 ×2，版式物理尺寸不变。
##       前置：须先用 16_style_theme via3 重跑 19_fig3 via12 生成 600dpi 面板。
## ---------------------------------------------------------------------------
## via2（2026-09-21，配合 19 via7 版式重构）：
##   三排版式 ABC / DE / F（列线图整宽底排，原 E/F 内容对调）；
##   字母仍走细留白条盖章（19 面板无内部标题，沿用 via1 策略）；
##   行循环通用化（rows_def 驱动），不再逐行硬编码。
## 目的：把 19_fig3 via7 输出的 6 个小图拼回整版 Fig3（化零为整）
## 版式：第 1 排 A/B/C（58×70mm ×3）；第 2 排 D/E（85×70mm ×2）；
##       第 3 排 F 列线图整宽（174×76mm）。整版 174×约 229mm ≤ 234mm。
## 字母策略：细留白条（LAB_H=48px）+ 粗体 36px 字母，面板左缘 +8px，
##       与 Fig2 定版外观一致，不压任何图形元素。
## 策略：magick 位图拼接——不重算、不开任何图形设备，
##       从根上规避 RStudio "Plots window may be too small" 设备报错。
## 输入：figures/ClinRheum风格/Fig3_panel{A..F}_*.tiff（19 via7 save_cr 产物）
## 输出：figures/ClinRheum风格/Fig3_拼接整版_174mm.tiff（600dpi LZW，投稿用）
##       figures/ClinRheum风格/Fig3_拼接整版_preview.png（快速目检用）
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

## --- 1. 读入 6 个小图（缺失即列清单报错；via7⑦：E/F 内容已对调） -------------
tags <- c(A = "A_四算法投票",
          B = "B_四算法Venn",
          C = "C_SHAP全局重要性",
          D = "D_SHAP蜜蜂图",
          E = "E_七分类器稳健性",          # via7⑦：原 F 升 E（85x70mm）
          F = "F_列线图")                  # via7⑦：原 E 升 F（整宽 174x76mm）
paths <- file.path(FIG_DIR, paste0("Fig3_panel", tags, ".tiff"))
missing_f <- paths[!file.exists(paths)]
if (length(missing_f) > 0) {
  stop("[拼图] 缺少 ", length(missing_f), " 个小图，请先跑 19_fig3 via12：\n  ",
       paste(missing_f, collapse = "\n  "), call. = FALSE)
}
panels <- lapply(paths, magick::image_read)
names(panels) <- names(tags)   # paste0 不保名，必须显式命名

## --- 2. 字母盖章（细留白条：字母在每排面板正上方左缘，不压图） ----------------
LAB_H <- 96                                # 留白条高 96px@600dpi（=48px@300dpi ≈4.1mm）
stamp_row <- function(row_img, letters, xs) {
  w <- magick::image_info(row_img)$width
  strip <- magick::image_blank(w, LAB_H, color = "white")
  for (i in seq_along(letters)) {
    strip <- magick::image_annotate(strip, letters[i], gravity = "northwest",
                                    location = paste0("+", xs[i], "+4"),
                                    size = 72, font = "sans", weight = 700,
                                    color = "black")
  }
  magick::image_append(c(strip, row_img), stack = TRUE)
}

## --- 3. 行内横拼（列间白缝）→ 统一拉宽到 174mm → 行间白缝竖拼 ----------------
GAP_H <- 52                                          # 列间白缝 px（=26px@300dpi）
GAP_V <- 24                                          # 行间白缝 px（=12px@300dpi）

hcat <- function(imgs) {                             # 横向拼接，中间夹白缝
  out <- imgs[[1]]
  for (img in imgs[-1]) {
    sp <- magick::image_blank(GAP_H, magick::image_info(out)$height, color = "white")
    out <- magick::image_append(c(out, sp, img))
  }
  magick::image_resize(out, paste0(PAGE_W, "x"))     # 统一整版宽，等比缩放
}

## 通用：n 列行各列左缘 x 坐标（按 hcat 实际缩放比折算；单面板行返回 c(0)）
.cols_x <- function(imgs) {
  ws <- sapply(imgs, function(im) magick::image_info(im)$width)
  total <- sum(ws) + GAP_H * (length(imgs) - 1)
  scale <- PAGE_W / total
  round(c(0, cumsum(ws[-length(ws)] + GAP_H)) * scale)
}

## via2：三排版式定义（成员 = panels 名字；字母按阅读顺序逐排盖章）
rows_def <- list(
  list(members = c("A", "B", "C"), letters = c("A", "B", "C")),
  list(members = c("D", "E"),      letters = c("D", "E")),
  list(members = c("F"),           letters = c("F"))
)
rows <- lapply(rows_def, function(rd) {
  imgs <- panels[rd$members]
  row_img <- hcat(imgs)
  stamp_row(row_img, rd$letters, .cols_x(imgs) + 16)
})

vcat <- function(rows) {                             # 竖向拼接，行间夹窄白缝
  out <- rows[[1]]
  for (r in rows[-1]) {
    sp <- magick::image_blank(PAGE_W, GAP_V, color = "white")
    out <- magick::image_append(c(out, sp, r), stack = TRUE)
  }
  out
}
fig3 <- vcat(rows)

## --- 4. 导出 ------------------------------------------------------------------
fin <- magick::image_info(fig3)
out_tiff <- file.path(FIG_DIR, "Fig3_拼接整版_174mm.tiff")
out_png  <- file.path(FIG_DIR, "Fig3_拼接整版_preview.png")
magick::image_write(fig3, out_tiff, format = "tiff", compression = "LZW")
magick::image_write(fig3, out_png,  format = "png")

message("\n========== Fig3 拼接完成（via2 三排版式） ==========")
message("整版尺寸：", fin$width, " × ", fin$height, " px",
        "（", round(fin$width / MM2PX, 1), " × ", round(fin$height / MM2PX, 1), " mm @600dpi）")
message("投稿图：", out_tiff)
message("预览图：", out_png)
