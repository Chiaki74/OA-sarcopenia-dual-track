# =============================================================================
# 16_style_theme_共享主题_ClinRheum风格_20260915_via3.R
# via3（2026-09-23）：save_cr 默认 TIFF 分辨率 300→600 dpi —— Clin Rheumatol 官方
#   要求组合彩图 ≥600 dpi（submission-guidelines / Artwork）；PDF 矢量版不受影响。
#   其余主题、色板、字体不变。旧版 via2 保留勿删。
# Clin Rheumatology 范文图风共享主题（所有出图脚本 source 本文件，保证全篇一致）
# 风格基准（像素级采样+画像解析）：
#   [R1] Clin Rheumatol 2026;45 RA锰代谢×ML诊断标志物 DOI 10.1007/s10067-026-08081-3
#   [R2] Clin Rheumatol 2025;44(12) M-CSF1×OA MR+共定位 DOI 10.1007/s10067-025-07765-6
#   [R3] Clin Rheumatol 2026;45(5) 膝OA scRNA          DOI 10.1007/s10067-026-08030-0
# 铁律：除具体内容外，样式与范文完全一致 —— 白底无网格、无衬线字体、
#        面板标签=左上角大写粗体字母、统计量入图、颜色语义化且跨图一致
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)   # 面板组合 + 左上角大写粗体字母标签（范文同款）
})

## ---------------------------------------------------------------------------
## 1. 全篇统一色板（CR_COL）—— 颜色必须语义化，任何图不得自定义新色
##    色值来自范文图版像素采样（量化 ±16），与范文肉眼不可辨
## ---------------------------------------------------------------------------
CR_COL <- list(
  ## 分组双色（R1 范文箱线图采样：Control=蓝 / Case=红，全篇统一）
  ctrl   = "#5090F0",   # 对照蓝
  case   = "#F07070",   # 疾病红（OA）
  muscle = "#F0A050",   # 肌少症橙（第二疾病，同族暖色）

  ## 三轨方向着色（双轨/三轨分析专用，全篇固定）
  concordant = "#5090F0",   # 同向轨 = 蓝
  mirror     = "#F07070",   # 镜像轨 = 红（核心卖点色）
  specific   = "#BFBFBF",   # 疾病特异/非显著 = 灰

  ## 火山图（R1 范文：柔和红蓝 + 灰）
  volc_up   = "#F0B0B0",
  volc_down = "#90D0F0",
  volc_ns   = "#C8C8C8",

  ## 共定位/locuszoom（R2 范文 Fig6 采样：navy→红 编码 LD r²，紫色菱形=lead SNP）
  ld_low  = "#101060",
  ld_high = "#F01010",
  lead    = "#7A00C0",

  ## 热图渐变（R2 范文 PPH4 热图：青绿→白→橙红）
  heat_low  = "#2A9D8F",
  heat_mid  = "#FFFFFF",
  heat_high = "#E76F51",

  ## 流程图分区底色（R2 范文 Fig1：米黄=数据 / 粉=分析路径 / 蓝=方法 / 青=收束）
  zone_data     = "#FDF3D8",
  zone_analysis = "#FDE8E4",
  zone_method   = "#E8F0FA",
  zone_result   = "#E4F2F5",
  box_border    = "#5A5A5A",
  arrow         = "#404040"
)

## ---------------------------------------------------------------------------
## 2. 共享主题 theme_cr() —— 白底、无网格、黑色轴线、无衬线字体
## ---------------------------------------------------------------------------
theme_cr <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "sans") %+replace%
    theme(
      axis.line        = element_line(color = "black", linewidth = 0.4),
      axis.ticks       = element_line(color = "black", linewidth = 0.3),
      axis.text        = element_text(color = "black", size = rel(0.90)),
      axis.title       = element_text(color = "black", size = rel(1.00)),
      plot.title       = element_text(face = "bold", size = rel(1.05), hjust = 0),
      plot.subtitle    = element_text(size = rel(0.90), hjust = 0),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      legend.text       = element_text(size = rel(0.85)),
      legend.title      = element_text(size = rel(0.90)),
      strip.background  = element_blank(),
      strip.text        = element_text(face = "bold", size = rel(0.95)),
      panel.border      = element_blank(),
      plot.margin       = margin(4, 6, 4, 4)
    )
}

## 面板标签主题：patchwork 加 tag 后套用 —— 左上角大写粗体字母（R1/R3 范文同款）
CR_TAG <- theme(
  plot.tag          = element_text(face = "bold", size = 14, family = "sans"),
  plot.tag.position = c(0.005, 0.995)
)

## 显著性星号（R1 范文：图内直接打星）
cr_stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***",
                             ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

## ---------------------------------------------------------------------------
## 3. 统一导出 save_cr() —— 双格式 PDF(矢量)+TIFF(600dpi,LZW)，宽度按刊规
##    单栏 84 mm / 双栏 174 mm；高 ≤234 mm（Clin Rheumatol 图版规范）
## ---------------------------------------------------------------------------
save_cr <- function(p, file, width_mm = 174, height_mm = 120, dpi = 600) {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  ggsave(paste0(file, ".pdf"),  p, width = width_mm, height = height_mm,
         units = "mm", device = cairo_pdf)
  ggsave(paste0(file, ".tiff"), p, width = width_mm, height = height_mm,
         units = "mm", dpi = dpi, compression = "lzw")
  message("[save_cr] ", file, ".pdf / .tiff (", width_mm, "x", height_mm, " mm)")
  invisible(file)
}

## ---------------------------------------------------------------------------
## 4. 防御性列名探测（用户本地中间结果列名不完全已知，按候选名正则匹配）
## ---------------------------------------------------------------------------
pick_col <- function(df, candidates, required = TRUE, label = "") {
  hit <- intersect(candidates, names(df))
  if (length(hit) > 1L) hit <- hit[1]   # 多列精确命中时按候选优先级取第一
  if (length(hit) == 0) {
    rx <- grep(paste(candidates, collapse = "|"), names(df),
               ignore.case = TRUE, value = TRUE)
    hit <- rx[1]
  }
  if (length(hit) == 0 || is.na(hit)) {
    if (required) stop("列名探测失败 [", label, "]：候选 ", paste(candidates, collapse = "/"),
                       "；实际列：", paste(names(df), collapse = ", "))
    return(NULL)
  }
  hit[1]
}

## 根目录定位（与既往模块脚本同款逻辑：向上找含 data/ 或 results/ 的目录）
locate_root <- function(markers = c("data", "results", "01_raw", "GEO")) {
  d <- getwd()
  for (i in 1:8) {
    if (any(dir.exists(file.path(d, markers)))) return(d)
    parent <- dirname(d)
    if (parent == d) break
    d <- parent
  }
  warning("未自动定位根目录，使用当前工作目录：", getwd())
  getwd()
}

message("[16_style_theme] Clin Rheumatology 风格主题已加载（CR_COL / theme_cr / CR_TAG / save_cr / pick_col / locate_root）")
