## ===========================================================================
## 29_拼图_Fig2_小图拼接_20260916_via8.R
## via8（2026-09-23）：投稿合规升版——DPI 300→600（刊规组合彩图 ≥600dpi），
##       全部像素常数同步 ×2（字母/缝隙/坐标），版式物理尺寸不变。
##       前置：须先用 16_style_theme via3 重跑 18_fig2 via20 生成 600dpi 面板。
## ---------------------------------------------------------------------------
## 目的：把 18_fig2 via17 输出的 6 个小图拼回整版 Fig2（化零为整）
## via7 相对 via6：图号-标题同行化（用户实图反馈「字母与标题不平行、不对齐」）：
##       废弃 gutter 留白条方案；配合 18 via17 统一的 plot.margin（左 14pt、
##       上 6pt，F 顶 8pt），标题全部从同一 x 起笔，左上角自然空出
##       ≈58px×50px 纯留白区；字母直接盖在该留白区，与标题同行同基线。
## 策略：magick 位图拼接——不重算、不开任何图形设备，
##       从根上规避 RStudio "Plots window may be too small" 设备报错。
## 布局：与原 patchwork 设计一致 "AA / BC / DE / FF"（A、F 通栏，BC、DE 双列），
##       整版宽 174mm，行间/列间留白缝，每面板左上角加粗体标签 A–F。
## 输入：figures/ClinRheum风格/Fig2_panel{A..F}_*.tiff（save_cr 双格式产物，600dpi LZW）
## 输出：figures/ClinRheum风格/Fig2_拼接整版_174mm.tiff（600dpi LZW，投稿用）
##       figures/ClinRheum风格/Fig2_拼接整版_preview.png（快速目检用）
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

## --- 1. 读入 6 个小图（缺失即列清单报错） ------------------------------------
tags <- c(A = "A_批次校正前后PCA",
          B = "B_代表队列火山图",
          C = "C_双轨镜像热图",
          D = "D_镜像占比条形图",
          E = "E_分轨富集气泡图",
          F = "F_富集compareCluster点图")
paths <- file.path(FIG_DIR, paste0("Fig2_panel", tags, ".tiff"))
missing_f <- paths[!file.exists(paths)]
if (length(missing_f) > 0) {
  stop("[拼图] 缺少 ", length(missing_f), " 个小图，请先跑 18_fig2 via20：\n  ",
       paste(missing_f, collapse = "\n  "), call. = FALSE)
}
panels <- lapply(paths, magick::image_read)
names(panels) <- names(tags)   # paste0 不保名，必须显式命名，否则 panels$A = NULL

## --- 2. 字母盖章（via7：盖在每行图左上角标题同行留白区，与标题同基线） ----------
## 18 via17 已统一所有面板左 margin=14pt（≈58px）：标题从 x≈58 起笔，
## 其左侧 x8–33px、标题同行的 y 带内无任何图形元素（y 轴标签/轴标题均在
## 面板纵中部），字母盖章安全无碰撞。
LAB_SIZE <- 72                               # 字号 72px@600dpi（=36px@300dpi，物理字号不变）
LAB_X    <- 16                               # 字母左缘 x（=8px@300dpi）
LAB_Y    <- 36                               # 第 1–3 行字母 y（=18px@300dpi）
LAB_Y_F  <- 52                               # 第 4 行 F（=26px@300dpi）
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
GAP_V <- 32                # via8：×2 适配 600dpi（=16px@300dpi）

hcat <- function(imgs) {                             # 横向拼接，中间夹白缝
  out <- imgs[[1]]
  for (img in imgs[-1]) {
    sp <- magick::image_blank(GAP_H, magick::image_info(out)$height, color = "white")
    out <- magick::image_append(c(out, sp, img))
  }
  magick::image_resize(out, paste0(PAGE_W, "x"))     # 统一整版宽，等比缩放
}

row1 <- magick::image_resize(panels$A, paste0(PAGE_W, "x"))
row2 <- hcat(list(panels$B, panels$C))
row3 <- hcat(list(panels$D, panels$E))
row4 <- magick::image_resize(panels$F, paste0(PAGE_W, "x"))

## 双列行第二列左缘 x 坐标：按 hcat 实际缩放比折算（B/C、D/E 原始宽度 + 列缝）
.col_x <- function(imgL, imgR) {
  lw <- magick::image_info(imgL)$width; rw <- magick::image_info(imgR)$width
  round((lw + GAP_H) * PAGE_W / (lw + GAP_H + rw))
}
xC <- .col_x(panels$B, panels$C)
xE <- .col_x(panels$D, panels$E)

row1 <- stamp_on(row1, "A", LAB_X)
row2 <- stamp_on(row2, c("B", "C"), c(LAB_X, xC + LAB_X))
row3 <- stamp_on(row3, c("D", "E"), c(LAB_X, xE + LAB_X))
row4 <- stamp_on(row4, "F", LAB_X, y = LAB_Y_F)

vcat <- function(rows) {                             # 竖向拼接，行间夹窄白缝
  out <- rows[[1]]
  for (r in rows[-1]) {
    sp <- magick::image_blank(PAGE_W, GAP_V, color = "white")
    out <- magick::image_append(c(out, sp, r), stack = TRUE)
  }
  out
}
fig2 <- vcat(list(row1, row2, row3, row4))

## --- 4. 导出 ------------------------------------------------------------------
fin <- magick::image_info(fig2)
out_tiff <- file.path(FIG_DIR, "Fig2_拼接整版_174mm.tiff")
out_png  <- file.path(FIG_DIR, "Fig2_拼接整版_preview.png")
magick::image_write(fig2, out_tiff, format = "tiff", compression = "LZW")
magick::image_write(fig2, out_png,  format = "png")

message("\n========== Fig2 拼接完成 ==========")
message("整版尺寸：", fin$width, " × ", fin$height, " px",
        "（", round(fin$width / MM2PX, 1), " × ", round(fin$height / MM2PX, 1), " mm @600dpi）")
message("投稿图：", out_tiff)
message("预览图：", out_png)
