# =============================================================================
# 17_fig1_流程图_ClinRheum风格_20260915_via9.R
# via9（2026-10-11，P1-3 销项：MR 仪器筛选行对齐 Methods 与主管线实际口径）：
#   证据（44 号核查脚本全库扫描）：moduleH_TwoSampleMR_coloc_pipeline.R 及
#   fix2–fix11e 中 MAF 仅用于 z_to_beta_se 换算 / eaf 近似 / coloc 数据列，
#   仪器筛选无 MAF>1% 过滤；实际为 iv_p_threshold=5e-8、clump_r2=0.001
#   （kb=10000，分级放宽 0.01/0.1）、check_f_stat 弱工具检验。
#   故图内 "MAF > 1%" 删除、"P < 5e-8" 补上（ASCII 写法与 p_mirror 同款，
#   保字体兼容）；几何零改动。
# via8（2026-10-11，via7 实跑预览两处溢出线修复）：
#   ① 机器学习框 4 行正文超出 15mm 框高（末行被底边横穿）→ 折回 3 行：
#     nomogram 并入第 3 行，"classifier robustness" 字样撤出图 1
#     （七分类器稳健性见 Fig 3E，流程图层面无信息损失）
#   ② 三层验证框红色 note 被虚线底边横穿 → 框底边 67→64.5mm
#     （仍在粉色区内，距区界 2.5mm）；B4→C1 直箭起点 y 66.5→64、
#     红弧 curved1 起点 y 67→64.5 同步下移
# via7（2026-10-11，主图审查 P0-3 / P1-4 / P1-5 / P2-3 文字口径修正）：
#   ① P0-3 上游调控框数据库名对齐 Methods：TRRUST / RegNetwork + miRWalk
#     → TRRUST / ChEA_2022 (TF) + miRTarBase (miRNA), via Enrichr（正文 4 行→5 行，框高富余）
#   ② P1-4 MR 框 "6 nominally significant loci" → "signals"
#     （6 个为基因-结局检验，共定位已排除共享因果位点，"loci" 具误导性）
#   ③ P2-3a 结局 GWAS 框 "OA-related outcomes" → "OA- & sarcopenia-related outcomes"
#     （4 个结局中 3 个为肌少症性状 ALM/握力/步速）
#   ④ P2-3b 机器学习框 "4-algorithm consensus" 展开为 Methods 口径
#     "LASSO∩SVM-RFE primary + RF/XGBoost top-10 voting"（正文 3→4 行，框高富余）
#   ⑤ P1-5 三层验证框补 muscle-proxy tier 一行（GSE117525/144304/136344，exploratory）；
#     框底边 69→67 mm 容纳新行，B4→C1 直箭起点 y 68.5→66.5、红弧 curved1 起点 y 69→67 同步
#   ⑥ P1-3（MR 仪器筛选 MAF > 1% 行）暂不改动——待 MR 主脚本核实后另行处理
# via6（2026-10-02，措辞与正文口径一致化 + 头注/自引用版本号修正）：
#   ① 机器学习框 4-algorithm intersection 改 consensus——严格四算法交集仅
#     3 个基因（Fig3B Venn 中心），10 个标志物为「LASSO∩SVM-RFE 为主、
#     RF/XGBoost 前10投票为辅」交叉策略产物，与正文口径一致；
#   ② 修正 via5 遗留：头注文件名 via2→via6、自引用路径 via4→via6。
# via4（2026-09-23）：主题引用改为自动探测最高 via —— 配合 16_style_theme via3（TIFF 600dpi，
#   刊规组合彩图要求），本脚本内容逻辑零改动，重跑即得 600dpi 面板/整版。
# 功能：Fig 1 研究流程图（单面板、无字母标签）——21 个公共数据集整合分析总体设计
# 风格锚点：Clin Rheumatol 2025;44 M-CSF1×OA MR 范文 Fig.1（img_p03_x36）像素级复刻
# via2 修订（2026-09-22，对照范文 + 高分期刊流程图惯例调研结论）：
#   ① 顶区框题溢出碰撞修复：4 个数据框全部改两行标题（题\n(n = X)），
#     正文起始 y 按标题行数自适应（build_boxes 内 n_title 逻辑）
#   ② 主箭头戳穿区标题修复：左/右/底区标题改为居中、缩至 3.1/3.0mm，
#     数据层入区箭头改走两侧干净通道（x=10 / x=96），不再穿过任何文字
#   ③ 红色收束弯箭头范文同款改造：加粗(1.1)+缩短+左右对称（B4→C2、R3→C2），
#     不再长距离横扫、不再与黑色直箭交叉
#   ④ 全部箭头加粗（直箭 0.85、箭头 2.2mm），逼近范文视觉权重
#   ⑤ 顶区框压缩（26→21mm 高），整版画布 200→190mm
#   ⑥ 惯例调研（7 篇 IF≥7 流程图 0/7 使用 "Fig. X" 编号徽章）→ 不加徽章；
#     顶部保持 "(21 datasets)" 口径（用户决定 2026-09-22）
#   ⑦ "-- aging proxy --" → "— aging proxy —"
# 实现：纯 ggplot2（geom_rect + geom_segment/geom_curve(arrow=)），不依赖 DiagrammeR 等
# 输入：无外部数据文件——全部数字为稿件冻结口径，集中在 CONFIG 块声明（如需改动只改这里）
# 输出：figures/ClinRheum风格/Fig1_研究流程图.pdf / .tiff（174 × 190 mm，双栏整版）
# 作者：R 绘图工程师（subagent）  日期：2026-09-15  via2/via3：2026-09-22  种子：set.seed(20260915)
# via3 修订（2026-09-22，via2 实跑预览两处微调）：
#   ⑧ B1 共享基因框正文 4→3 行（binomial P 并入 mirror 行）——via2 末行压框底边
#   ⑨ 底区标题改 geom_label 区底色底板：两条红弧必经中央标题带（几何不可绕行），
#     底板让弧从标题底下穿过，文字完整可读（范文红箭头穿框底同款处理）
# =============================================================================

## ---- 0. 定位并加载共享主题（脚本与主题同目录；找不到则回退 getwd）---------------
## via5（2026-09-30）：修复底区两条红色收束弧与居中区标题互相遮盖——
##   底区标题「Causal inference & translational convergence」由居中单行
##   改为左对齐双行（x=8, hjust=0），让出中央走廊；B4→C1 黑色直箭
##   由 x=44 移至 x=56（避开左置标题）；红弧路径不变，收束点仍落 C2 顶边，
##   弧与标题零重叠。

this_file <- tryCatch(
  normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])),
  error = function(e) NA_character_)
if (is.na(this_file) || !file.exists(this_file))
  this_file <- file.path(getwd(), "17_fig1_流程图_ClinRheum风格_20260915_via9.R")  # via9
SCRIPT_DIR  <- dirname(this_file)
theme_cand <- c(list.files(SCRIPT_DIR, pattern = "^16_style_theme.*\\.R$", full.names = TRUE),
                list.files(getwd(), pattern = "^16_style_theme.*\\.R$", full.names = TRUE))
if (length(theme_cand) == 0)
  stop("找不到共享主题 16_style_theme 脚本；请将本脚本与主题文件放在同一目录后重跑。")
theme_path <- sort(theme_cand, decreasing = TRUE)[1]   # 文件名最大 = 最高 via（via3 起 TIFF=600dpi）
source(theme_path)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})
set.seed(20260915)

## ---- 1. CONFIG：冻结口径数字 + 输出设置 -----------------------------------------
## 以下数字与 v6 主稿「图注与表注 → 主文图表（方案 A）」及结果节完全一致；
## 投稿前若任何中间结果更新，仅需修改本块并重跑。
CONFIG <- list(
  n_bulk        = 18,                    # GEO bulk 转录组数
  n_sc          = 3,                     # sc/snRNA-seq 数据集数
  n_gwas        = 4,                     # 结局 GWAS 数
  n_shared      = 99,                    # 共享候选基因
  n_mirror      = 72,                    # 镜像轨基因
  n_conc        = 27,                    # 同向轨基因
  mirror_pct    = 72.7,                  # 镜像占比 %
  p_mirror      = "3.5e-6",              # 单侧二项检验 P（ASCII 写法，保证各设备字体兼容）
  n_terms       = 35,                    # 分轨富集显著条目（FDR<0.05）
  n_biomarkers  = 10,                    # 四算法交集标志物
  auc_ext       = "0.818 (0.614-1.000)", # 软骨主验证 GSE57218 AUC (95%CI, DeLong 实算 2026-09-16)
  n_mr_tests    = 32,                    # MR 检验组合数
  n_nominal     = 6,                     # 名义显著 MR 信号（via7②：图内措辞 loci→signals）
  drugg_hit     = "5/10",                # 可成药基因组命中
  drugg_universe= "5,928",               # 可成药基因组宇宙大小
  out_dir       = "figures/ClinRheum风格",
  out_file      = "Fig1_研究流程图",
  width_mm      = 174,
  height_mm     = 190                    # via2⑤：200 → 190（顶区压缩）
)

ROOT <- locate_root()
OUT  <- file.path(ROOT, CONFIG$out_dir, CONFIG$out_file)

## ---- 2. 画布坐标系说明 -----------------------------------------------------------
## 以 mm 为逻辑坐标：x ∈ [0,174]，y ∈ [0,190]，与导出尺寸 1:1，布局所见即所得。
## 各区几何参数集中于此，便于微调（改这里即可整体挪动分区）。

## ---- 3. 各面板构造函数 -----------------------------------------------------------

## 3.1 分区底色层（范文同款 pastel 底色，无描边）
##     via2⑤：数据层 150-186 / 左右分析层 62-146 / 底部因果层 6-58
build_zones <- function() {
  data.frame(
    xmin = c(4,  4,  90, 4),
    xmax = c(170, 84, 170, 170),
    ymin = c(150, 62, 62,  6),
    ymax = c(186, 146, 146, 58),
    fill = c(CR_COL$zone_data, CR_COL$zone_analysis,
             CR_COL$zone_method, CR_COL$zone_result),
    stringsAsFactors = FALSE)
}

## 3.2 白色内容框（阈值写进框内；血液层框虚线边=探索性）
##     每框：title=粗体框题（支持 "\n" 两行）；body=多行小字（"• " 项目符，范文同款）；
##     note=框内红色强调行（探索性标注）；dashed=TRUE 时虚线描边
##     via2①：正文/强调行起始 y 按标题行数 n_title 自适应，杜绝两行标题压正文
build_boxes <- function() {
  boxes <- list(
    ## —— 顶层数据框（4 个；via2①：全部两行标题，body 4 行，框高 21mm）——
    list(x = c(8, 44),    y = c(152, 173),
         title = sprintf("GEO bulk transcriptomes\n(n = %d)", CONFIG$n_bulk),
         body  = "OA: cartilage / synovium /\nbone / blood\nMuscle: sarcopenia anchor +\nfrailty / aging proxy cohorts"),
    list(x = c(48, 84),   y = c(152, 173),
         title = sprintf("sc/snRNA-seq\n(n = %d)", CONFIG$n_sc),
         body  = "GSE152805 / GSE169454\n(OA synovium)\nGSE167186 (muscle)\n\u2014 aging proxy \u2014"),
    list(x = c(90, 126),  y = c(152, 173),
         title = sprintf("Outcome GWAS\n(n = %d)", CONFIG$n_gwas),
         body  = "OA- & sarcopenia-related\noutcomes (OpenGWAS)\ndisease GWAS summary\nstatistics"),  # via7③
    list(x = c(130, 168), y = c(152, 173),
         title = "eQTL & druggable\nresources",
         body  = "eQTLGen blood cis-eQTL\nGTEx v8 muscle eQTL\ndruggable genome\n(n = 5,928)"),
    ## —— 左粉区：发现-机制-标志物路径（4 框纵链）——
    list(x = c(8, 80), y = c(124, 139),
         title = "Dual-track shared-gene construction",
         body  = sprintf(paste0(
           "limma/DESeq2 per cohort + RRA integration\n",
           "\u2022 adj.P < 0.05 \u2192 %d shared candidates\n",
           "\u2022 mirror %d/%d (%.1f%%), concordant %d; P = %s"),  # via3⑧：4→3 行防压底边
           CONFIG$n_shared, CONFIG$n_mirror, CONFIG$n_shared,
           CONFIG$mirror_pct, CONFIG$n_conc, CONFIG$p_mirror)),
    list(x = c(8, 80), y = c(105, 120),
         title = "Track-level functional enrichment",
         body  = sprintf(paste0(
           "GO / KEGG / Reactome (ORA)\n",
           "\u2022 FDR < 0.05 \u2192 %d track-level terms\n",
           "\u2022 OXPHOS / ubiquitin / complement axes"), CONFIG$n_terms)),
    list(x = c(8, 80), y = c(86, 101),
         title = "Machine-learning biomarker selection",
         body  = sprintf(paste0(
           "LASSO / SVM-RFE / RF / XGBoost\n",
           "\u2022 LASSO\u2229SVM-RFE primary + RF/XGBoost\n  top-10 voting \u2192 %d biomarkers; nomogram"),  # via8①: 折回 3 行防压底边
           CONFIG$n_biomarkers)),
    list(x = c(8, 80), y = c(64.5, 83.5), dashed = TRUE,   # via8②: 67→64.5（note 压底边修复）
         title = "Three-tier external validation",
         body  = sprintf(paste0(
           "\u2022 cartilage GSE57218 AUC %s\n",
           "\u2022 calibration + DCA\n",
           "\u2022 muscle-proxy tier: exploratory"), CONFIG$auc_ext),   # via7⑤
         note  = "blood layer: hypothesis-generating"),
    ## —— 右蓝区：免疫浸润 / 单细胞 / 调控通讯路径 ——
    list(x = c(94, 166), y = c(122, 139),
         title = "Immune infiltration",
         body  = paste0(
           "CIBERSORT / ssGSEA / xCell\n",
           "\u2022 16 replicated signals, all in the\n  OA tissue layer\n",
           "\u2022 muscle / blood / aging: none")),
    list(x = c(94, 166), y = c(101, 118),
         title = "Single-cell localization",
         body  = paste0(
           "3 datasets / Seurat / Harmony / UCell\n",
           "\u2022 GSE167186 = aging proxy\n",
           "\u2022 biomarkers \u2192 stromal / progenitor\n  compartments")),
    list(x = c(94, 166), y = c(69, 97),
         title = "Upstream regulation & communication",
         body  = paste0(
           "TF-miRNA network (75 edges)\n",
           "\u2022 TRRUST / ChEA_2022 (TF)\n",
           "\u2022 miRTarBase (miRNA), via Enrichr\n",   # via7①: 对齐 Methods（原 RegNetwork+miRWalk 有误）
           "\u2022 CellChat cell-cell communication\n",
           "\u2022 marker-related L-R interactions")),
    ## —— 底青区：因果与转化收束（3 框横排）——
    list(x = c(8, 60),   y = c(16, 52),
         title = "Two-sample MR",
         body  = sprintf(paste0(
           "\u2022 %d exposure-outcome tests\n",
           "\u2022 IVW primary; MR-Egger, weighted\n  median, MR-PRESSO sensitivity\n",
           "\u2022 instrument screen: P < 5e-8,\n  LD r2 < 0.001, F-statistic > 10\n",   # via9: 删 MAF>1%（未执行），补 P<5e-8
           "\u2022 %d nominally significant signals"),   # via7②
           CONFIG$n_mr_tests, CONFIG$n_nominal)),
    list(x = c(64, 116), y = c(16, 52),
         title = "Bayesian colocalization",
         body  = paste0(
           "\u2022 PPH4 > 0.8: high support\n",
           "\u2022 0.5 < PPH4 < 0.8: medium support\n",
           "\u2022 coloc.abf + SuSiE fine-mapping\n  sensitivity\n",
           "\u2022 eQTLGen blood + GTEx v8 muscle")),
    list(x = c(120, 166), y = c(16, 52),
         title = "Druggable genome mapping",
         body  = sprintf(paste0(
           "\u2022 %s biomarkers within the\n  druggable genome (n = %s)\n",
           "\u2022 DSigDB candidates: exploratory\n  (none passed FDR)\n",
           "\u2022 drug repurposing hypotheses"),
           CONFIG$drugg_hit, CONFIG$drugg_universe))
  )
  rects <- do.call(rbind, lapply(boxes, function(b) data.frame(
    xmin = b$x[1], xmax = b$x[2], ymin = b$y[1], ymax = b$y[2],
    linetype = ifelse(isTRUE(b$dashed), "dashed", "solid"),
    stringsAsFactors = FALSE)))
  titles <- do.call(rbind, lapply(boxes, function(b) data.frame(
    x = b$x[1] + 1.6, y = b$y[2] - 2.0, label = b$title,
    stringsAsFactors = FALSE)))
  bodies <- do.call(rbind, lapply(boxes, function(b) {
    n_title <- length(strsplit(b$title, "\n", fixed = TRUE)[[1]])   # via2①
    data.frame(x = b$x[1] + 1.6, y = b$y[2] - 2.0 - n_title * 2.6 - 1.0,
               label = b$body, stringsAsFactors = FALSE)
  }))
  notes <- do.call(rbind, lapply(boxes, function(b) {
    if (is.null(b$note)) return(NULL)
    n_title <- length(strsplit(b$title, "\n", fixed = TRUE)[[1]])   # via2①
    n_lines <- length(strsplit(b$body, "\n", fixed = TRUE)[[1]])
    data.frame(x = b$x[1] + 1.6,
               y = b$y[2] - 2.0 - n_title * 2.6 - 1.0 - n_lines * 2.55,
               label = paste0("\u2022 ", b$note), stringsAsFactors = FALSE)
  }))
  list(rects = rects, titles = titles, bodies = bodies, notes = notes)
}

## 3.3 分区标题与数据层分组括号（范文 Exposure/Outcome 括号同款）
##     via2②：左/右/底区标题居中（hjust=0.5），入区箭头移至两侧干净通道；
##     标题字号随行宽分级（3.3 / 3.1 / 3.0），用 scale_size_identity 逐行控制
build_zone_labels <- function() {
  zone_titles <- data.frame(
    x     = c(87, 44, 130),
    y     = c(182.6, 142.8, 142.8),
    hjust = c(0.5, 0.5, 0.5),
    size  = c(3.3, 3.1, 3.1),
    label = c("Public data resources (21 datasets)",
              "Discovery \u2192 mechanism \u2192 biomarker",
              "Immune / single-cell / regulation"),
    stringsAsFactors = FALSE)
  ## via3⑨：底区标题独立成行，改用 geom_label 区底色底板——
  ## 两条红色收束弧必经中央标题带（几何上无法绕行），加底板后弧从标题"底下"穿过，
  ## 文字完整可读、弧线视觉连续（范文红箭头穿框底同款处理）
  bottom_title <- data.frame(
    x = 8, y = 55.8,
    label = "Causal inference &\ntranslational convergence",
    stringsAsFactors = FALSE)
  sub_titles <- data.frame(
    x = c(46, 129), y = 178.4,
    label = c("Transcriptome data", "Genetic & druggable resources"),
    stringsAsFactors = FALSE)
  brackets <- data.frame(
    x  = c(9, 91), xend = c(83, 167), y = 176, yend = 176)
  ticks <- data.frame(
    x = c(9, 83, 91, 167), xend = c(9, 83, 91, 167),
    y = 176, yend = 174.8)
  list(zone_titles = zone_titles, bottom_title = bottom_title,
       sub_titles = sub_titles, brackets = brackets, ticks = ticks)
}

## 3.4 箭头层：黑色直箭头=主分析流；红色弯箭头=收束至因果层（范文 Stage2 红箭头同款）
##     via2②：数据层入区箭头 x=46/128 → x=10/96（避开居中区标题）
##     via2③：红弧缩短+左右对称（B4 右下→C2 顶左 / R3 左下→C2 顶右），
##            全程不与 x=44 / x=130 黑色直箭交叉、不穿越其他框
##     via2④：直箭 linewidth 0.5→0.85，箭头 1.6→2.2mm；红弧 linewidth 1.1
build_arrows <- function() {
  straight <- data.frame(
    x    = c(10, 96,               # via2②：数据层 → 左/右两区（两侧干净通道）
             44, 44, 44,           # 左区纵链 B1→B2→B3→B4
             130, 130,             # 右区纵链 R1→R2→R3
             56, 130,              # B4→C1、R3→C3（下行入因果层；via5：44→56 避左置标题）
             60.2, 116.2),         # C1→C2、C2→C3（横向推进）
    y    = c(150, 150,
             123.5, 104.5, 85.5,
             121.5, 100.5,
             64, 68.5,             # via8②: B4→C1 起点随框底下移 66.5→64
             34, 34),
    xend = c(10, 96,
             44, 44, 44,
             130, 130,
             56, 130,
             63.8, 119.8),
    yend = c(139.8, 139.8,
             120.5, 101.5, 84,
             118.5, 97.5,
             52.5, 52.5,
             34, 34),
    stringsAsFactors = FALSE)
  ## 红色弯箭头：双路径收束至共定位框 C2（curvature 为 geom_curve 参数而非美学，
  ## 须分两次调用、各给固定值）；via2③：短弧、对称、互不交叉
  curved1 <- data.frame(x = 64,  y = 64.5, xend = 80, yend = 52.6)  # B4 底右 → C2 顶左（via8②: 67→64.5）
  curved2 <- data.frame(x = 126, y = 69, xend = 98, yend = 52.6)  # R3 底左 → C2 顶右
  list(straight = straight, curved1 = curved1, curved2 = curved2)
}

## 3.5 主面板拼装（单面板，无字母标签 —— 范文 Fig.1 同款）
panel_main <- function() {
  z  <- build_zones()
  b  <- build_boxes()
  lb <- build_zone_labels()
  ar <- build_arrows()
  arw <- grid::arrow(length = grid::unit(2.2, "mm"), type = "closed")  # via2④

  ggplot() +
    ## 分区底色
    geom_rect(data = z,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
              color = NA) +
    scale_fill_identity() +
    ## 白色内容框
    geom_rect(data = b$rects,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                  linetype = linetype),
              fill = "white", color = CR_COL$box_border, linewidth = 0.35) +
    scale_linetype_identity() +
    ## 分组括号与端点短竖线（数据层）
    geom_segment(data = lb$brackets, aes(x = x, y = y, xend = xend, yend = yend),
                 color = CR_COL$box_border, linewidth = 0.4) +
    geom_segment(data = lb$ticks, aes(x = x, y = y, xend = xend, yend = yend),
                 color = CR_COL$box_border, linewidth = 0.4) +
    ## 箭头（via2④：加粗）
    geom_segment(data = ar$straight, aes(x = x, y = y, xend = xend, yend = yend),
                 arrow = arw, color = CR_COL$arrow, linewidth = 0.85) +
    geom_curve(data = ar$curved1, aes(x = x, y = y, xend = xend, yend = yend),
               curvature = 0.25, arrow = arw, color = CR_COL$case, linewidth = 1.1) +
    geom_curve(data = ar$curved2, aes(x = x, y = y, xend = xend, yend = yend),
               curvature = -0.25, arrow = arw, color = CR_COL$case, linewidth = 1.1) +
    ## 文字层
    geom_text(data = lb$zone_titles,
              aes(x = x, y = y, label = label, hjust = hjust, size = size),
              vjust = 0.5, fontface = "bold", color = "black") +
    scale_size_identity() +
    geom_label(data = lb$bottom_title, aes(x = x, y = y, label = label),  # via3⑨/via5
               size = 3.0, fontface = "bold", color = "black", hjust = 0,
               lineheight = 1.05,
               fill = CR_COL$zone_result, label.size = 0,
               label.padding = grid::unit(0.12, "lines")) +
    geom_text(data = lb$sub_titles, aes(x = x, y = y, label = label),
              size = 2.6, fontface = "bold", color = "black") +
    geom_text(data = b$titles, aes(x = x, y = y, label = label),
              hjust = 0, vjust = 1, size = 2.4, lineheight = 1.05,
              fontface = "bold", color = "black") +
    geom_text(data = b$bodies, aes(x = x, y = y, label = label),
              hjust = 0, vjust = 1, size = 2.0, lineheight = 1.12, color = "#333333") +
    geom_text(data = b$notes, aes(x = x, y = y, label = label),
              hjust = 0, vjust = 1, size = 2.0, fontface = "bold",
              color = CR_COL$case) +
    coord_cartesian(xlim = c(0, 174), ylim = c(0, 190), expand = FALSE) +
    theme_void(base_family = "sans") +
    theme(plot.margin = margin(2, 2, 2, 2))
}

## ---- 4. 主流程 -------------------------------------------------------------------
fig1 <- panel_main()
save_cr(fig1, OUT, width_mm = CONFIG$width_mm, height_mm = CONFIG$height_mm)
message("[17_fig1 via9] 完成：", OUT, ".pdf / .tiff（174 × 190 mm）")
