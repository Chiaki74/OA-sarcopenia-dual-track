# =============================================================================
# 18_fig2_双轨镜像_富集_ClinRheum风格_20260915_via23.R（via19 派生副本）
# via23（2026-10-02，via22 实跑修正）：via22 把 BCL6 左推 2.0 后撞上 KLF9
#   的自动位（两标签同槽互叠）。改向：BCL6 推向上方中空区——两条 ±1 虚线
#   之间、灰色 NotSig 云上方为验证过的空白带（nudge_x +1.6 / nudge_y +0.6），
#   蓝色标签+引线指向性不变；其余 9 个标签仍走 via20 自动布局。
# via22（2026-10-02，via21 实跑修正）：via21 的全员统一外推导致 SLC19A2
#   与 BCL6 标签互叠。改为：撤掉统一 nudge 与 force 改动（恢复 via20 对
#   其余 9 个标签的良好布局），仅对 BCL6 单基因定点外推（x -2.0 / y +0.8，
#   左上稀疏区），白描边 halo 保留。
# via21（2026-10-02，B 火山图 BCL6 标签可读性终审修复）：via19 已记
#   「BCL6 压点」未根治——该标签落在阈值虚线与密集点云交界处被埋没。
#   修复：按方向外推 nudge_x ±0.9 + force 2/force_pull 0.6 + 白描边 halo
#   （bg.color，需 ggrepel >= 0.9.2；若老版本报错删去 bg.color/bg.r 即可）。
# via20（2026-09-23）：主题引用改为自动探测最高 via —— 配合 16_style_theme via3（TIFF 600dpi，
#   刊规组合彩图要求），本脚本内容逻辑零改动，重跑即得 600dpi 面板/整版。
# via19 修订（2026-09-22，via18 整版预览终审 4 处）：
#   ① B 火山图 ggrepel：box.padding 0.35→0.5 + point.padding 0.3 + seed 锁定
#     （SLC19A2 引线穿字 / BCL6 压点 / RHOB 蹭阈值虚线）
#   ② C 热图 strip 去 Mus 缩写 → 与 E/F 统一 "Muscle-down/up"
#   ③ E 副标题删生产备注 "Color scale unified with panel F"（统一色阶写进图注）
# 功能：Fig 2 双轨共享基因集构建与镜像失调格局 + 分轨富集（6 面板，2 列矩阵）
#   (a) 批次校正前后 PCA 两联（分组色 ctrl/case + stat_ellipse 置信椭圆）
#   (b) 代表队列火山图（volc_* 三轨着色 Up/Down/NotSig，阈值虚线，top 基因标注）
#   (c) 双轨热图（行=99 基因按轨分面，列=队列 logFC，heat_low→mid→high；
#       大字 inset：镜像占比 72.7% (72/99) 与单侧二项 P = 3.5×10⁻⁶ —— R1「统计量入图」）
#   (d) sarcopenia-only 敏感性条形（镜像占比两柱，计数直接标柱上）
#   (e) 分轨富集气泡图（x=轨，y=term，size=Count，color=−log10(p.adj)，FoxO 行标注）
#   (f) compareCluster dotplot（x=GeneRatio，database × track 分面，color=p.adjust）
# 风格锚点：Clin Rheumatol 2026;45 RA 锰代谢×ML 范文 Fig.1（img_p05_x88）——
#   多面板矩阵、白底无网格、左上大写粗体字母、统计量入图
# 输入（均在 CONFIG 块声明，pick_col() 防御性探测列名；相对 locate_root() 根目录）：
#   f_pca_scores  批次校正前后 PCA 得分长表。期望列：sample（样本名）/ PC1 / PC2 /
#                 group（分组，含 ctrl|control|normal 者判为 Control，其余 Case）/
#                 stage（含 before|pre|raw 者判为校正前，其余校正后）[cohort 可选]
#   【a 兜底】expr_before/expr_after 校正前后表达矩阵 rds（行=基因，列=样本）+
#             f_group 分组表 csv（列=sample, group）；缺 PCA 得分表时现场计算
#   f_volcano     代表队列差异表（limma topTable）。期望列：gene / logFC /
#                 P.Value / adj.P.Val。文件缺失时自动取 de_dir 下首个 layer1_DE_*.csv
#   f_heat_long   【c 首选】热图长表。期望列：gene / cohort / logFC；
#                 缺失时由 de_dir 的 layer1_DE_*.csv（逐队列 topTable）现算
#   f_tracks      轨成员表。期望列：gene / track（concordant /
#                 mirror_OAup_muscleDown / mirror_OAdown_muscleUp）；
#                 缺失时由 f_core_con + f_core_mir（列=gene[, direction_OA]）推导
#   f_shared      主分析共享基因表。期望列：gene / direction_OA / direction_muscle
#                 [concordant 逻辑列可选]（layer2_shared_DE_genes.csv）
#   f_sens        sarcopenia-only 敏感性表，列结构同 f_shared
#   f_enrich      分轨富集显著条目表。期望列：track / Description / p.adjust /
#                 GeneRatio（形如 "12/340"）/ Count / database（GO-BP|Reactome|KEGG）
# 输出：figures/ClinRheum风格/Fig2_panel{A..F}_*.pdf/.tiff（逐面板小图，29 号脚本拼接）
# via11（2026-09-16，审图意见 13 条落地版）：
#   ① B 火山 y 轴改 -log10(adj.P)（与着色口径一致，虚线仍 1.3）+ 零点兜底；
#   ② D 重叠标注改 15/16 (93.8%) 并移至 Sarcopenia-only 柱，柱色改中性灰(⑫)；
#   ③ C 热图加组织类型列分面（OA joint tissue vs Skeletal muscle，按 GSE 映射）；
#      ※审图意见正文把 GSE114007(实为OA软骨)/GSE25941(实为骨骼肌)归属写反，按 GEO 事实映射；
#   ④ A 之 Before 按数据集着色（需 27 via2 的 cohort 列，缺失则回退并提示）；
#   ⑤⑥ E 保留但色阶与 F 统一为 p.adjust（low=case/high=ctrl），降强调仅作概览；
#   ⑩ B 点色加深（case/ctrl 实色）、标签放大防重叠、x 轴补 ±2 刻度；
#   ⑪ C 色带改 blue-white-red（#2166AC/#FFFFFF/#B2182B）+ strip 标注基因数并压窄；
#   ⑦⑧⑨⑬ 及椭圆置信水平/离群点说明为图注与正文层面——见决策备忘录 §198。
# via12（2026-09-16）：审图意见⑬试案（用户定规「先尝试压缩画布」）——
#   画布压缩为 A 174×36 / B–E 85×58 / F 174×62 mm，配 29 via6（字母条 48px/
#   字号 36/行缝 8px），整版 ≈232mm ≤ 234mm 版心；面板内容与 via11 完全一致。
# via13（2026-09-16，压缩后实图 QA 修复，画布尺寸不变仍 ≈232mm）：
#   B top 标注 8→5/侧、字号 2.5→2.3（58mm 画布内 16 标签互叠）；
#   D 柱内白字改三行短行 + 柱宽 0.6→0.68（原长行白字两端出柱身落在白底隐形）；
#   E top_n_e 4→3、折行 38→44 字符、y 字号 rel(0.72)→rel(0.65)（12 项互叠）；
#   F top_n_f 6→4（Reactome 行 15+ 术语在 62mm 内互叠）。
# via14（2026-09-16，A 面板 cohort 匹配健壮化）：
#   实图 QA 发现 Before 的 Dataset 图例缺 GSE114007——scores 表样本名与
#   sample_group.csv 的 GSM 号精确匹配失败所致。改双路匹配：先抽 GSM 号互配
#   （容忍前缀/后缀），失败再用样本名自带 GSE 号（队列均值点情形）；
#   并打印 Before/After 分阶段匹配率 + 未匹配样本示例，杜绝静默丢队列。
# via15（2026-09-16，F 标题与列 strip 遮挡微调——仅 panel_f 主题，不涉他图）：
#   标题 size rel(0.9) + margin(b=5pt) 撑开与 strip 行的垂直间隙；
#   plot.title.position="plot" 从整幅图左缘起笔（与 E 面板同款），
#   标题主体落在 y 轴标签区上方，避开中部 strip 密集区。
# via16（2026-09-16，标题排版统一，用户定规）：
#   ① 全部面板 plot.title.position="plot"——加粗标题与图号左缘平齐、
#     同排面板标题同一水平线（A/B/C/D 补齐，E/F 已是）；
#   ② 标题-副标题间距统一：plot.title margin(b=4pt)、plot.subtitle margin(b=5pt)，
#     解决 B/C 标题与下一行相碰；间距由面板区收缩吸收，无裁切风险。
# via18（2026-09-16，A 图例压缩 + rds 路径举证，用户实图反馈 Control 键被裁）：
#   A 为 36mm 矮行，合并图例（Dataset 4-5 键 + Group 2 键 + 2 标题）默认键高
#   超出画布，Control 键被底缘裁掉——压缩 legend 键高/行距/字号（仅作用于
#   panel_a，其余面板不受影响）；pca_from_rds 增加控制台举证行（样本数 +
#   cohort 分布），用于坐实 expr_before_batch.rds 缺 GSE114007 的数据层问题。
# via17（2026-09-16，图号-标题同行化，用户实图反馈「字母与标题不平行」）：
#   统一全部面板 plot.margin 左=14pt、上=6pt（F 顶 8pt 保留防气泡顶裁；
#   C 斜排标签左突 4.6mm < 14pt=4.9mm 仍可容纳）——标题全部从同一 x 起笔；
#   配合 29 via7：废弃 gutter 留白条，字母直接盖在每行左上角标题同行留白区。
# 作者：R 绘图工程师（subagent）  日期：2026-09-15  种子：set.seed(20260915)
# =============================================================================

## ---- 0. 定位并加载共享主题 ---------------------------------------------------------
this_file <- tryCatch(
  normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])),
  error = function(e) NA_character_)
if (is.na(this_file) || !file.exists(this_file))
  this_file <- file.path(getwd(), "18_fig2_双轨镜像_富集_ClinRheum风格_20260915_via4.R")
SCRIPT_DIR <- dirname(this_file)
theme_cand <- c(list.files(SCRIPT_DIR, pattern = "^16_style_theme.*\\.R$", full.names = TRUE),
                list.files(getwd(), pattern = "^16_style_theme.*\\.R$", full.names = TRUE))
if (length(theme_cand) == 0)
  stop("找不到共享主题 16_style_theme 脚本；请将本脚本与主题文件放在同一目录后重跑。")
theme_path <- sort(theme_cand, decreasing = TRUE)[1]   # 文件名最大 = 最高 via（via3 起 TIFF=600dpi）
source(theme_path)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
})
set.seed(20260915)

## ---- 1. CONFIG：输入文件路径与全局参数 ---------------------------------------------
CONFIG <- list(
  ## (a) PCA
  f_pca_scores = "results/moduleA/fig2a_pca_scores_long.csv",
  expr_before  = "results/moduleA/expr_before_batch.rds",
  expr_after   = "results/moduleA/expr_after_batch.rds",
  f_group      = "data/sample_group.csv",
  ## (b) 火山（代表队列；缺失时自动取 de_dir 首个 layer1_DE_*.csv）
  f_volcano    = "results/moduleB/layer1_DE_GSE55235.csv",
  de_dir       = "results/moduleB",
  ## (c) 热图
  f_heat_long  = "results/moduleB/fig2c_heatmap_logFC_long.csv",
  f_tracks     = "results/moduleC/moduleC_gene_tracks.csv",
  f_core_con   = "results/moduleB/CORE_concordant_genes.csv",
  f_core_mir   = "results/moduleB/CORE_mirror_genes.csv",
  ## (d) 敏感性
  f_shared     = "results/moduleB/layer2_shared_DE_genes.csv",
  f_sens       = "results/moduleB/layer2_sensitivity_sarcopenia_only.csv",
  ## (e)(f) 富集
  f_enrich     = "results/moduleC/enrich_ALL_tracks_significant.csv",
  ## 参数与冻结口径（d 面板兜底 / c 面板 inset 缺数据时回退）
  fc_cut       = 1,          # 火山 |logFC| 阈值（仅图上标注，不改分析口径）
  p_cut        = 0.05,       # adj.P 阈值
  top_n_e      = 3,          # (e) 每轨展示条目数（via13：4→3，58mm 画布 12 项必然叠字）
  top_n_f      = 4,          # (f) 每 database×track 展示条目数（via13：6→4，62mm 画布防叠字）
  heat_clip    = 2,          # 热图 logFC 截断
  frozen_n_shared = 99, frozen_n_mirror = 72, frozen_pct = 72.7,  # 冻结口径
  frozen_sens_total = 21, frozen_sens_mirror = 16, frozen_sens_pct = 76.2,
  frozen_overlap = "15/16 (93.8%)",   # via11②：分母=敏感性镜像轨 16，非主分析 72
  ## 输出
  out_dir      = "figures/ClinRheum风格",
  out_file     = "Fig2_双轨镜像_分轨富集",
  width_mm     = 174,
  height_mm    = 215
)

ROOT <- locate_root()
rp <- function(...) file.path(ROOT, ...)   # 根目录相对路径快捷函数

## ---- 2. 小工具 ---------------------------------------------------------------------

## 文件存在性检查（缺失→带候选路径的报错，便于用户一次性补齐）
need_file <- function(path, label) {
  if (!file.exists(path))
    stop(sprintf("找不到%s：\n  %s\n请在 CONFIG 块修正路径后重跑。", label, path),
         call. = FALSE)
  path
}

## 轨标签归一化：识别 moduleB/C 的三轨命名，未识别时保留原值并提示
norm_track <- function(x) {
  x0 <- as.character(x)
  out <- x0
  out[grepl("concord", x0, ignore.case = TRUE)] <- "Concordant"
  out[grepl("mirror.*OAup.*muscleDown|OA.?up.*mus.?down", x0, ignore.case = TRUE)] <-
    "Mirror OA-up / Muscle-down"
  out[grepl("mirror.*OAdown.*muscleUp|OA.?down.*mus.?up", x0, ignore.case = TRUE)] <-
    "Mirror OA-down / Muscle-up"
  unk <- setdiff(unique(out), c("Concordant", "Mirror OA-up / Muscle-down",
                                "Mirror OA-down / Muscle-up"))
  if (length(unk) > 0)
    message("[norm_track] 未识别的轨标签（原样保留）：", paste(unk, collapse = ", "))
  factor(out, levels = c("Concordant", "Mirror OA-up / Muscle-down",
                         "Mirror OA-down / Muscle-up", unk))
}

## 分组标签归一化（PCA/分组表通用）：含 ctrl/control/normal/healthy → Control
norm_group <- function(x) {
  x0 <- as.character(x)
  ifelse(grepl("ctrl|control|normal|healthy", x0, ignore.case = TRUE),
         "Control", "Case")
}

## 从共享/敏感性表统计镜像数（concordant 列或 direction 两列均可）
count_mirror <- function(d, label = "") {
  con_col <- pick_col(d, c("concordant", "is_concordant", "same_direction"),
                      required = FALSE, label = label)
  if (!is.null(con_col)) {
    con <- d[[con_col]]
    con <- if (is.logical(con)) con else grepl("^T|TRUE|yes|1", as.character(con),
                                               ignore.case = TRUE)
    return(list(total = nrow(d), mirror = sum(!con, na.rm = TRUE),
                gene = d[[pick_col(d, c("gene", "Gene", "symbol"), label = label)]],
                is_mirror = !con))
  }
  doa <- d[[pick_col(d, c("direction_OA", "OA_direction", "dir_OA"), label = label)]]
  dmu <- d[[pick_col(d, c("direction_muscle", "direction_mus", "muscle_direction",
                          "dir_muscle"), label = label)]]
  is_mir <- as.character(doa) != as.character(dmu)
  list(total = nrow(d), mirror = sum(is_mir, na.rm = TRUE),
       gene = d[[pick_col(d, c("gene", "Gene", "symbol"), label = label)]],
       is_mirror = is_mir)
}

## GeneRatio "12/340" → 数值
parse_gene_ratio <- function(x) {
  x <- as.character(x)
  ifelse(grepl("/", x),
         sapply(strsplit(x, "/", fixed = TRUE),
                function(z) as.numeric(z[1]) / as.numeric(z[2])),
         suppressWarnings(as.numeric(x)))
}

## ---- 3. 面板 (a)：批次校正前后 PCA 两联 --------------------------------------------
## 返回 patchwork（左=Before，右=After）；作为整体在外层获得一个字母标签
panel_a <- function() {
  f_scores <- rp(CONFIG$f_pca_scores)
  mk_pca_plot <- function(dd, ttl, color_by = "group") {
    xlab <- if (is.na(dd$ve1[1])) "PC1" else sprintf("PC1 (%.1f%%)", dd$ve1[1])
    ylab <- if (is.na(dd$ve2[1])) "PC2" else sprintf("PC2 (%.1f%%)", dd$ve2[1])
    ## via11（审图意见④）：Before 按数据集着色以展示批次效应；
    ## 无 cohort 信息时回退 Case/Control 并提示先跑 27 via2。
    use_cohort <- identical(color_by, "cohort") && "cohort" %in% names(dd) &&
                  any(!is.na(dd$cohort)) &&
                  length(unique(stats::na.omit(dd$cohort))) > 1
    if (use_cohort) {
      p <- ggplot(dd, aes(PC1, PC2, color = cohort)) +
        geom_point(size = 1.6, alpha = 0.85) +
        stat_ellipse(aes(group = group), level = 0.95, linewidth = 0.4,
                     color = "grey40", show.legend = FALSE) +
        scale_color_brewer(palette = "Set2", name = "Dataset")
    } else {
      if (identical(color_by, "cohort"))
        message("[panel_a] Before 阶段无 cohort 列，回退 Case/Control 着色",
                "（先运行 27_出图数据准备 via2 生成含 cohort 的 sample_group.csv）。")
      p <- ggplot(dd, aes(PC1, PC2, color = group)) +
        geom_point(size = 1.6, alpha = 0.85) +
        stat_ellipse(aes(group = group), level = 0.95, linewidth = 0.4,
                     show.legend = FALSE) +
        scale_color_manual(values = c(Control = CR_COL$ctrl, Case = CR_COL$case),
                           name = "Group")   # via14：与 Dataset 图例并置时加标题
    }
    p + labs(title = ttl, x = xlab, y = ylab) +
      theme_cr(base_size = 8) +   # via6：全图统一 8pt 基准，同层级文字同号
      theme(plot.title.position = "plot",          # via16：标题与图号左缘平齐
            plot.title = element_text(margin = margin(b = 4)),
            plot.margin = margin(6, 4, 6, 14),     # via17：统一左/上边距
            legend.key.size = unit(8, "pt"),       # via18：压缩图例总高——36mm
            legend.spacing.y = unit(0.5, "pt"),    #   矮行内默认键高放不下
            legend.margin = margin(0, 0, 0, 0),    #   Dataset+Group 全部键，
            legend.title = element_text(size = rel(0.85), margin = margin(b = 2)),
            legend.text = element_text(size = rel(0.85)))  # Control 曾被底缘裁掉
  }

  if (file.exists(f_scores)) {
    ## 首选：直接读 PCA 得分长表
    d <- read.csv(f_scores, check.names = FALSE)
    d$PC1   <- as.numeric(d[[pick_col(d, c("PC1", "pc1"), label = "PCA-PC1")]])
    d$PC2   <- as.numeric(d[[pick_col(d, c("PC2", "pc2"), label = "PCA-PC2")]])
    d$group <- norm_group(d[[pick_col(d, c("group", "Group", "condition", "phenotype"),
                                      label = "PCA-group")]])
    stage_raw <- as.character(d[[pick_col(d, c("stage", "batch", "correction", "phase"),
                                          label = "PCA-stage")]])
    d$stage <- ifelse(grepl("before|pre|raw|uncorr", stage_raw, ignore.case = TRUE),
                      "Before batch correction", "After batch correction")
    ## via11（审图意见④）：cohort 列缺失时从 sample_group.csv 按样本名补齐
    coh_col <- pick_col(d, c("cohort", "gse", "dataset", "series"),
                        required = FALSE, label = "PCA-cohort")
    if (!is.null(coh_col)) {
      d$cohort <- as.character(d[[coh_col]])
    } else {
      s_col <- pick_col(d, c("sample", "Sample", "gsm", "GSM", "id"),
                        required = FALSE, label = "PCA-sample")
      f_grp2 <- rp(CONFIG$f_group)
      if (!is.null(s_col) && file.exists(f_grp2)) {
        g2 <- read.csv(f_grp2, check.names = FALSE)
        g2s <- pick_col(g2, c("sample", "Sample", "gsm", "GSM", "id"),
                        required = FALSE, label = "group-sample")
        g2c <- pick_col(g2, c("cohort", "gse", "dataset", "series"),
                        required = FALSE, label = "group-cohort")
        if (!is.null(g2s) && !is.null(g2c)) {
          ## via14：双路匹配——先抽 GSM 号互配（容忍前缀/后缀），
          ## 失败再用样本名自带 GSE 号（队列均值点情形）
          extract_pat <- function(x, pat) {
            x <- as.character(x)
            hit <- grepl(pat, x)
            out <- rep(NA_character_, length(x))
            out[hit] <- regmatches(x[hit], regexpr(pat, x[hit]))
            out
          }
          gsm_d <- extract_pat(d[[s_col]], "GSM[0-9]+")
          gsm_g <- extract_pat(g2[[g2s]], "GSM[0-9]+")
          coh <- as.character(g2[[g2c]])[match(gsm_d, gsm_g)]
          gse_self <- extract_pat(d[[s_col]], "GSE[0-9]+")
          fix <- is.na(coh) & !is.na(gse_self)
          coh[fix] <- gse_self[fix]
          d$cohort <- coh
          for (st in unique(d$stage)) {
            ix <- d$stage == st
            message(sprintf("[panel_a] cohort 匹配率 %s: %d/%d",
                            st, sum(!is.na(d$cohort[ix])), sum(ix)))
          }
          if (any(is.na(d$cohort)))
            message("[panel_a] 未匹配样本示例：",
                    paste(head(as.character(d[[s_col]])[is.na(d$cohort)], 5),
                          collapse = ", "))
        }
      }
    }
    ve_col <- pick_col(d, c("ve1", "var1", "pct1"), required = FALSE, label = "PCA-ve1")
    if (!is.null(ve_col)) {
      ve2_col <- pick_col(d, c("ve2", "var2", "pct2"), required = FALSE,
                          label = "PCA-ve2")
      d$ve1 <- as.numeric(d[[ve_col]])
      d$ve2 <- if (!is.null(ve2_col)) as.numeric(d[[ve2_col]]) else NA_real_
    } else {
      d$ve1 <- NA_real_; d$ve2 <- NA_real_
    }
    lst <- split(d, d$stage)
    if (is.null(lst[["After batch correction"]]))
      stop("[panel_a] stage 列无法识别校正后阶段（期望含 before/pre/raw 之外的取值）；",
           "实际取值：", paste(unique(stage_raw), collapse = ", "))
    p_after <- mk_pca_plot(lst[["After batch correction"]], "After batch correction")
    if (is.null(lst[["Before batch correction"]])) {
      message("[panel_a] 得分表无 before 阶段，仅绘制 After。")
      return(p_after)
    }
    p_before <- mk_pca_plot(lst[["Before batch correction"]], "Before batch correction",
                            color_by = "cohort")   # via11④
    return((p_before | p_after) + plot_layout(guides = "collect") &
             theme(legend.position = "right",
                   legend.key.size = unit(8, "pt"),           # via18：收集层同步压缩
                   legend.spacing.y = unit(0.5, "pt"),
                   legend.margin = margin(0, 0, 0, 0),
                   legend.title = element_text(size = rel(0.85), margin = margin(b = 2)),
                   legend.text = element_text(size = rel(0.85))))
  }

  ## 兜底：由表达矩阵 rds + 分组表现算（top 5000 高变基因，prcomp）
  message("[panel_a] 无 PCA 得分表，改由表达矩阵现算：", CONFIG$expr_after)
  f_grp <- need_file(rp(CONFIG$f_group), "样本分组表 f_group")
  grp <- read.csv(f_grp, check.names = FALSE)
  g_sample <- pick_col(grp, c("sample", "Sample", "gsm", "GSM", "id"), label = "group-sample")
  g_group  <- pick_col(grp, c("group", "Group", "condition", "phenotype"), label = "group-group")
  g_cohort <- pick_col(grp, c("cohort", "gse", "dataset", "series"),
                       required = FALSE, label = "group-cohort")   # via11④
  grp <- data.frame(sample = as.character(grp[[g_sample]]),
                    group  = norm_group(grp[[g_group]]),
                    cohort = if (!is.null(g_cohort)) as.character(grp[[g_cohort]])
                             else NA_character_,
                    stringsAsFactors = FALSE)

  pca_from_rds <- function(f_rds, label) {
    if (!file.exists(f_rds)) return(NULL)
    m <- as.matrix(readRDS(f_rds))
    com <- intersect(colnames(m), grp$sample)
    if (length(com) < 4)
      stop("[panel_a] ", label, " 表达矩阵列名与分组表 sample 交集 <4，请核对样本命名。")
    m <- m[, com, drop = FALSE]
    v <- apply(m, 1, var, na.rm = TRUE)
    m2 <- m[order(v, decreasing = TRUE)[seq_len(min(5000, nrow(m)))], , drop = FALSE]
    pc <- prcomp(t(m2), center = TRUE, scale. = FALSE)
    ve <- 100 * (pc$sdev^2 / sum(pc$sdev^2))
    coh_tab <- table(grp$cohort[match(com, grp$sample)])   # via18：举证
    message("[panel_a] ", label, " 矩阵：", length(com), " 样本；cohort 分布：",
            paste(sprintf("%s=%d", names(coh_tab), as.integer(coh_tab)),
                  collapse = ", "))
    data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2],
               group = grp$group[match(com, grp$sample)],
               cohort = grp$cohort[match(com, grp$sample)],   # via11④
               ve1 = ve[1], ve2 = ve[2])
  }
  d_after  <- pca_from_rds(rp(CONFIG$expr_after),  "after")
  d_before <- pca_from_rds(rp(CONFIG$expr_before), "before")
  if (is.null(d_after)) stop("[panel_a] 既无 PCA 得分表也无 expr_after rds，无法绘图。")
  p_after <- mk_pca_plot(d_after, "After batch correction")
  if (is.null(d_before)) {
    message("[panel_a] 无 expr_before rds，仅绘制 After。")
    return(p_after)
  }
  (mk_pca_plot(d_before, "Before batch correction", color_by = "cohort") | p_after) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right",
          legend.key.size = unit(8, "pt"),              # via18：收集层同步压缩
          legend.spacing.y = unit(0.5, "pt"),
          legend.margin = margin(0, 0, 0, 0),
          legend.title = element_text(size = rel(0.85), margin = margin(b = 2)),
          legend.text = element_text(size = rel(0.85)))
}

## ---- 4. 面板 (b)：代表队列火山图 ---------------------------------------------------
panel_b <- function() {
  f_vol <- rp(CONFIG$f_volcano)
  if (!file.exists(f_vol)) {
    alts <- list.files(rp(CONFIG$de_dir), pattern = "^layer1_DE_.*\\.csv$",
                       full.names = TRUE)
    if (length(alts) == 0)
      stop("[panel_b] 找不到代表队列差异表，且 ", CONFIG$de_dir,
           " 下无 layer1_DE_*.csv；请在 CONFIG 修正 f_volcano。")
    f_vol <- alts[1]
    message("[panel_b] f_volcano 不存在，自动改用：", basename(f_vol))
  }
  d <- read.csv(f_vol, check.names = FALSE)
  cohort <- sub("^layer1_DE_(.*)\\.csv$", "\\1", basename(f_vol))
  d$gene  <- as.character(d[[pick_col(d, c("gene", "Gene", "symbol"), label = "volcano-gene")]])
  d$logFC <- as.numeric(d[[pick_col(d, c("logFC", "log2FC"), label = "volcano-logFC")]])
  d$P     <- as.numeric(d[[pick_col(d, c("P.Value", "pvalue", "P"), label = "volcano-P")]])
  d$adjP  <- as.numeric(d[[pick_col(d, c("adj.P.Val", "FDR", "padj", "adjP"),
                                      label = "volcano-adjP")]])
  d <- d[is.finite(d$logFC) & is.finite(d$P) & is.finite(d$adjP), ]
  ## via11（审图意见①）：y 轴与着色统一为 adj.P 口径；adj.P=0 取最小正值 1/10 兜底
  pos_p <- d$adjP[d$adjP > 0]
  d$adjP <- pmax(d$adjP, if (length(pos_p) > 0) min(pos_p) / 10 else 1e-300)
  d$neglog_adj <- -log10(d$adjP)
  d$sig <- ifelse(d$adjP < CONFIG$p_cut & d$logFC >  CONFIG$fc_cut, "Up",
           ifelse(d$adjP < CONFIG$p_cut & d$logFC < -CONFIG$fc_cut, "Down", "NotSig"))
  d$sig <- factor(d$sig, levels = c("Up", "Down", "NotSig"))

  ## 标注 top 基因（上/下调各取 adj.P 最小 5 个；via13：8→5，58mm 画布防互叠）
  lab_up   <- head(d[d$sig == "Up",   ][order(d$adjP[d$sig == "Up"]),   ], 5)
  lab_down <- head(d[d$sig == "Down", ][order(d$adjP[d$sig == "Down"]), ], 5)
  lab <- rbind(lab_up, lab_down)
  n_up <- sum(d$sig == "Up"); n_down <- sum(d$sig == "Down")

  ggplot(d, aes(logFC, neglog_adj, color = sig)) +
    geom_point(size = 0.9, alpha = 0.75) +
    geom_vline(xintercept = c(-CONFIG$fc_cut, CONFIG$fc_cut),
               linetype = "dashed", color = "grey40", linewidth = 0.35) +
    geom_hline(yintercept = -log10(CONFIG$p_cut),   # adj.P = 0.05 → 1.30
               linetype = "dashed", color = "grey40", linewidth = 0.35) +
    geom_text_repel(data = lab, aes(label = gene), size = 2.3, max.overlaps = Inf,
                    box.padding = 0.5, point.padding = 0.3,  # via19①：标签-点间距防引线穿字
                    nudge_x = ifelse(lab$gene == "BCL6",  1.6, 0),  # via23：BCL6 右推至两虚线间空白带
                    nudge_y = ifelse(lab$gene == "BCL6",  0.6, 0),  # via23：微抬，避开点云上缘
                    bg.color = "white", bg.r = 0.12,         # via21：白描边防点云遮字（保留）
                    seed = 20260915,                         # via19①：锁定 ggrepel 随机起点
                    min.segment.length = 0, segment.size = 0.25,
                    show.legend = FALSE) +
    scale_color_manual(values = c(Up = CR_COL$case, Down = CR_COL$ctrl,
                                  NotSig = CR_COL$volc_ns),
                       name = "Expression") +
    scale_x_continuous(breaks = scales::breaks_width(2)) +   # via11⑩：±2 刻度
    labs(title = sprintf("Volcano plot: %s", cohort),
         subtitle = sprintf("adj.P < %g & |log2FC| > %g  (Up %d / Down %d)",
                            CONFIG$p_cut, CONFIG$fc_cut, n_up, n_down),
         x = "log2 Fold Change", y = "-log10 adj.P") +
    theme_cr(base_size = 8) +   # via6：全图统一 8pt 基准，同层级文字同号
    theme(plot.title.position = "plot",            # via16：标题与图号左缘平齐
          plot.title = element_text(margin = margin(b = 4)),
          plot.subtitle = element_text(margin = margin(b = 5)),
          plot.margin = margin(6, 6, 6, 14))       # via17：统一左/上边距
}

## ---- 5. 面板 (c)：双轨热图 + 镜像占比大字 inset ------------------------------------
panel_c <- function() {
  ## 5.1 轨成员表（首选 moduleC_gene_tracks.csv；兜底 CORE_* 两表推导）
  f_tr <- rp(CONFIG$f_tracks)
  if (file.exists(f_tr)) {
    tr <- read.csv(f_tr, check.names = FALSE)
    track_df <- data.frame(
      gene  = as.character(tr[[pick_col(tr, c("gene", "Gene", "symbol"), label = "track-gene")]]),
      track = norm_track(tr[[pick_col(tr, c("track", "Track", "group", "Group"),
                                        label = "track-track")]]),
      stringsAsFactors = FALSE)
  } else {
    message("[panel_c] 无 moduleC_gene_tracks.csv，改由 CORE_concordant/mirror 推导三轨。")
    fc <- need_file(rp(CONFIG$f_core_con), "CORE_concordant_genes.csv")
    fm <- need_file(rp(CONFIG$f_core_mir), "CORE_mirror_genes.csv")
    dc <- read.csv(fc, check.names = FALSE); dm <- read.csv(fm, check.names = FALSE)
    gc <- pick_col(dc, c("gene", "Gene", "symbol"), label = "core_con-gene")
    gm <- pick_col(dm, c("gene", "Gene", "symbol"), label = "core_mir-gene")
    dm_dir <- pick_col(dm, c("direction_OA", "OA_direction", "dir_OA"),
                       required = FALSE, label = "core_mir-direction_OA")
    trk_m <- if (!is.null(dm_dir))
      ifelse(grepl("^up", as.character(dm[[dm_dir]]), ignore.case = TRUE),
             "mirror_OAup_muscleDown", "mirror_OAdown_muscleUp") else "mirror"
    track_df <- data.frame(
      gene  = c(as.character(dc[[gc]]), as.character(dm[[gm]])),
      track = norm_track(c(rep("concordant", nrow(dc)), trk_m)),
      stringsAsFactors = FALSE)
  }
  track_df <- track_df[!duplicated(track_df$gene), ]

  ## 5.2 热图长表（首选现成长表；兜底由 de_dir 逐队列 topTable 现算）
  f_heat <- rp(CONFIG$f_heat_long)
  if (file.exists(f_heat)) {
    h <- read.csv(f_heat, check.names = FALSE)
    heat_long <- data.frame(
      gene   = as.character(h[[pick_col(h, c("gene", "Gene", "symbol"), label = "heat-gene")]]),
      cohort = as.character(h[[pick_col(h, c("cohort", "gse", "dataset", "series"),
                                          label = "heat-cohort")]]),
      logFC  = as.numeric(h[[pick_col(h, c("logFC", "log2FC"), label = "heat-logFC")]]),
      stringsAsFactors = FALSE)
  } else {
    message("[panel_c] 无热图长表，由 ", CONFIG$de_dir, " 的 layer1_DE_*.csv 现算。")
    de_files <- list.files(rp(CONFIG$de_dir), pattern = "^layer1_DE_.*\\.csv$",
                           full.names = TRUE)
    if (length(de_files) == 0)
      stop("[panel_c] 无热图长表且 de_dir 下无 layer1_DE_*.csv，无法构建热图。")
    heat_long <- do.call(rbind, lapply(de_files, function(f) {
      d <- read.csv(f, check.names = FALSE)
      data.frame(
        gene   = as.character(d[[pick_col(d, c("gene", "Gene", "symbol"), label = "DE-gene")]]),
        cohort = sub("^layer1_DE_(.*)\\.csv$", "\\1", basename(f)),
        logFC  = as.numeric(d[[pick_col(d, c("logFC", "log2FC"), label = "DE-logFC")]]),
        stringsAsFactors = FALSE)
    }))
  }
  heat_long <- heat_long[heat_long$gene %in% track_df$gene, ]
  if (nrow(heat_long) == 0)
    stop("[panel_c] 热图 logFC 数据与轨成员表基因无交集——请核对两侧基因命名",
         "（基因列示例：", paste(head(unique(track_df$gene), 3), collapse = ", "), "）")
  heat_long <- merge(heat_long, track_df, by = "gene")
  heat_long$logFC[heat_long$logFC >  CONFIG$heat_clip] <-  CONFIG$heat_clip
  heat_long$logFC[heat_long$logFC < -CONFIG$heat_clip] <- -CONFIG$heat_clip

  ## via11（审图意见③）：组织类型列分面。按 GSE 号映射组织来源——
  ## 注意：审图意见正文将 GSE114007（实为 OA 软骨）与 GSE25941（实为骨骼肌）写反，
  ## 此处按 GEO 事实映射；热图出现未映射队列时立即报错提示补充，绝不静默错标。
  TISSUE_MAP <- c(GSE55235 = "OA joint tissue", GSE55457 = "OA joint tissue",
                  GSE12021 = "OA joint tissue", GSE51588 = "OA joint tissue",
                  GSE114007 = "OA joint tissue",
                  GSE25941 = "Skeletal muscle", GSE111016 = "Skeletal muscle",
                  GSE111010 = "Skeletal muscle", GSE111006 = "Skeletal muscle",
                  GSE1428 = "Skeletal muscle")
  cc_obs <- unique(as.character(heat_long$cohort))
  unmapped <- setdiff(cc_obs, names(TISSUE_MAP))
  if (length(unmapped) > 0)
    stop("[panel_c] 热图队列缺少组织映射：", paste(unmapped, collapse = ", "),
         "——请按 GEO 组织来源补充 panel_c 中的 TISSUE_MAP。")
  message("[panel_c] 队列→组织：",
          paste(sprintf("%s=%s", cc_obs, TISSUE_MAP[cc_obs]), collapse = ", "))
  heat_long$tissue <- factor(TISSUE_MAP[as.character(heat_long$cohort)],
                             levels = c("OA joint tissue", "Skeletal muscle"))

  ## 5.3 行序：轨 → 轨内平均 logFC；列序按数据出现顺序
  tr_lv <- levels(track_df$track)
  tr_lv <- tr_lv[tr_lv %in% unique(as.character(heat_long$track))]
  heat_long$track <- factor(as.character(heat_long$track), levels = tr_lv)
  ord <- do.call(c, lapply(tr_lv, function(tk) {
    g <- unique(heat_long$gene[heat_long$track == tk])
    m <- tapply(heat_long$logFC[heat_long$gene %in% g],
                heat_long$gene[heat_long$gene %in% g], mean, na.rm = TRUE)
    names(sort(m[as.character(g)], decreasing = TRUE))
  }))
  heat_long$gene   <- factor(heat_long$gene, levels = rev(ord))
  heat_long$cohort <- factor(heat_long$cohort, levels = unique(heat_long$cohort))

  ## 5.4 镜像占比统计（数据缺失时回退冻结口径）
  ## via5：统计量移入副标题——原「右侧 2 空列 + inset geom_text」在 85mm 宽内
  ## 长文本必然溢出面板被裁（via4 实图确认），废弃空列方案。
  n_tot <- length(unique(track_df$gene))
  n_mir <- length(unique(track_df$gene[track_df$track != "Concordant"]))
  if (n_tot == 0) { n_tot <- CONFIG$frozen_n_shared; n_mir <- CONFIG$frozen_n_mirror }
  pct <- 100 * n_mir / n_tot
  sub_lab <- sprintf("Mirror proportion = %.1f%% (%d/%d) | One-sided binomial P = 3.5e-6",
                     pct, n_mir, n_tot)

  ## via5：strip 轨名去 "Mirror " 前缀（标题已点明 dual-track/mirror），
  ## 水平 strip 在 85mm 内才放得下，灰块不再挤压热图主体。
  ## via11（审图意见⑪）：strip 两行化并追加各 block 实际入图基因数
  ## （Muscle 缩写为 Mus 控制灰块宽度）。
  tr_n <- vapply(split(as.character(heat_long$gene), as.character(heat_long$track)),
                 function(g) length(unique(g)), integer(1))
  names(tr_n) <- sub("^Mirror ", "", names(tr_n))
  levels(heat_long$track) <- sub("^Mirror ", "", levels(heat_long$track))
  lv_now <- levels(heat_long$track)
  lv_short <- lv_now   # via19②：strip 不再缩写 Mus-，与 E/F 统一 "Muscle-down/up"
                       # （两行 strip 已含 (n=)，宽度代价可接受）
  hit <- lv_now %in% names(tr_n)
  lv_short[hit] <- sprintf("%s\n(n=%d)", lv_short[hit], tr_n[lv_now[hit]])
  levels(heat_long$track) <- lv_short

  ggplot(heat_long, aes(cohort, gene)) +
    geom_tile(aes(fill = logFC), color = NA) +
    facet_grid(track ~ tissue, scales = "free", space = "free") +   # via11③：组织列分面
    scale_fill_gradient2(low = "#2166AC", mid = "#FFFFFF",   # via11⑪：色盲友好
                         high = "#B2182B", midpoint = 0,     # blue-white-red
                         limits = c(-CONFIG$heat_clip, CONFIG$heat_clip),
                         oob = scales::squish, name = "logFC") +
    labs(title = sprintf("Dual-track candidates (%d genes)", n_tot),
         subtitle = sub_lab, x = NULL, y = NULL) +
    theme_cr(base_size = 8) +
    theme(plot.title.position = "plot",   # via16：标题与图号左缘平齐
          plot.title = element_text(margin = margin(b = 4)),
          axis.text.y = element_blank(),   # 99 行基因名在双栏矩阵中不可读，隐藏；
          axis.ticks.y = element_blank(),  # 需基因名时改回 element_text(size = 3)
          ## via8：按用户指定改回 45° 斜排（参考范文样式：斜行且完整）。
          ## 关键在左侧距——45°+hjust=1 时最左标签向左下探出
          ## ≈ cos45°×标签宽 ≈ 4.6mm（"GSE55235"），左边距 4→16pt(≈5.6mm)
          ## 即可完整容纳，不会再被画布左缘裁切（via7 的 90° 方案废止）。
          axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.72)),
          strip.text.y = element_text(size = rel(0.62), angle = 0, lineheight = 0.95,
                                      margin = margin(1, 1, 1, 1)),   # via11⑪：压窄灰块
          strip.text.x = element_text(size = rel(0.7), face = "bold",
                                      margin = margin(1.5, 1, 1.5, 1)),   # via11③：组织条
          strip.background = element_rect(fill = "grey92", color = NA),
          plot.subtitle = element_text(size = rel(0.72), margin = margin(b = 5)),   # via16
          plot.margin = margin(6, 8, 6, 14),   # via17：统一边距（斜排标签 4.6mm<4.9mm 仍容纳）
          legend.position = "right")
}

## ---- 6. 面板 (d)：sarcopenia-only 敏感性条形 ----------------------------------------
panel_d <- function() {
  ## 主分析镜像统计（文件缺失/列不识→冻结口径兜底并 message 提示）
  stat_main <- tryCatch({
    d <- read.csv(need_file(rp(CONFIG$f_shared), "主分析共享基因表 f_shared"),
                  check.names = FALSE)
    count_mirror(d, "shared")
  }, error = function(e) {
    message("[panel_d] 主分析表读取失败，回退冻结口径 72/99：", conditionMessage(e))
    list(total = CONFIG$frozen_n_shared, mirror = CONFIG$frozen_n_mirror,
         gene = NULL, is_mirror = NULL)
  })
  stat_sens <- tryCatch({
    d <- read.csv(need_file(rp(CONFIG$f_sens), "敏感性表 f_sens"), check.names = FALSE)
    count_mirror(d, "sens")
  }, error = function(e) {
    message("[panel_d] 敏感性表读取失败，回退冻结口径 16/21：", conditionMessage(e))
    list(total = CONFIG$frozen_sens_total, mirror = CONFIG$frozen_sens_mirror,
         gene = NULL, is_mirror = NULL)
  })

  df <- data.frame(
    analysis = factor(c("Primary\n(all muscle cohorts)", "Sarcopenia-only"),
                      levels = c("Primary\n(all muscle cohorts)", "Sarcopenia-only")),
    total  = c(stat_main$total,  stat_sens$total),
    mirror = c(stat_main$mirror, stat_sens$mirror),
    stringsAsFactors = FALSE)
  df$pct <- 100 * df$mirror / df$total
  df$lab <- sprintf("%d/%d (%.1f%%)", df$mirror, df$total, df$pct)

  ## via11（审图意见②）：重合度分母应为敏感性镜像轨（16）而非主分析（72）；
  ## 标注移至 Sarcopenia-only 柱（x=2）；柱色改中性灰（审图意见⑫，
  ## 避免与 Case/Up 红色语义冲突）。
  if (!is.null(stat_main$gene) && !is.null(stat_sens$gene)) {
    ov <- length(intersect(stat_sens$gene[stat_sens$is_mirror],
                           stat_main$gene[stat_main$is_mirror]))
    overlap_lab <- sprintf("Replicated in\nprimary mirror:\n%d/%d (%.1f%%)",
                           ov, stat_sens$mirror, 100 * ov / stat_sens$mirror)
  } else {
    overlap_lab <- sprintf("Replicated in\nprimary mirror:\n%s",
                           CONFIG$frozen_overlap)
  }

  ggplot(df, aes(analysis, pct)) +
    geom_col(fill = "grey65", width = 0.68) +   # via13：加宽以容纳三行白字
    geom_text(aes(label = lab), vjust = -0.4, size = 2.8, fontface = "bold") +
    annotate("text", x = 2, y = df$pct[2] * 0.5, label = overlap_lab,
             size = 2.2, color = "white", fontface = "bold", lineheight = 0.95) +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.08))) +
    labs(title = "Sensitivity: sarcopenia-only reconstruction",
         x = NULL, y = "Mirror-dysregulated candidates (%)") +
    theme_cr(base_size = 8) +   # via6：全图统一 8pt 基准，同层级文字同号
    theme(plot.title.position = "plot",            # via16：标题与图号左缘平齐
          plot.title = element_text(margin = margin(b = 4)),
          plot.margin = margin(6, 6, 6, 14))       # via17：统一左/上边距
}

## ---- 7. 富集表读取（e/f 共用） ------------------------------------------------------
load_enrich <- function() {
  f_en <- need_file(rp(CONFIG$f_enrich), "分轨富集显著条目表 f_enrich")
  d <- read.csv(f_en, check.names = FALSE)
  out <- data.frame(
    track = norm_track(d[[pick_col(d, c("track", "Track", "group", "Group", "Cluster"),
                                      label = "enrich-track")]]),
    term  = as.character(d[[pick_col(d, c("Description", "description", "Term", "term", "ID"),
                                        label = "enrich-term")]]),
    padj  = as.numeric(d[[pick_col(d, c("p.adjust", "p_adj", "FDR", "qvalue", "padj"),
                                      label = "enrich-padj")]]),
    stringsAsFactors = FALSE)
  gr_col  <- pick_col(d, c("GeneRatio", "gene_ratio", "RichFactor"),
                      required = FALSE, label = "enrich-GeneRatio")
  cnt_col <- pick_col(d, c("Count", "count", "n"), required = FALSE, label = "enrich-Count")
  db_col  <- pick_col(d, c("database", "Database", "source", "Source", "ONTOLOGY"),
                      required = FALSE, label = "enrich-db")
  out$gene_ratio <- if (!is.null(gr_col)) parse_gene_ratio(d[[gr_col]]) else NA_real_
  out$count      <- if (!is.null(cnt_col)) as.numeric(d[[cnt_col]]) else NA_real_
  out$database   <- if (!is.null(db_col)) as.character(d[[db_col]]) else "all"
  out <- out[is.finite(out$padj), ]
  out$neglog <- -log10(out$padj)
  out
}

## ---- 8. 面板 (e)：分轨富集气泡图（x=轨，y=term，size=Count，color=−log10(p.adj)）----
panel_e <- function(en) {
  df <- en[order(en$padj), ]
  df <- do.call(rbind, lapply(split(df, df$track), function(x) head(x, CONFIG$top_n_e)))
  ## y 轴条目顺序：轨 → padj 降序（最显著条目置于底部，dotplot 惯例）
  df <- df[order(as.numeric(df$track), df$padj, decreasing = c(FALSE, TRUE), method = "radix"), ]
  df$term <- factor(df$term, levels = unique(df$term))
  df <- df[!is.na(df$term), ]

  ## via5：x 轴轨名去 "Mirror " 前缀（与 panel_c 统一口径）；
  ## 长通路名按 38 字符自动折行——85mm 宽内 20+ 条全长 Reactome 名必然溢出
  df$track_disp <- sub("^Mirror ", "", as.character(df$track))
  df$track_disp <- factor(df$track_disp, levels = unique(df$track_disp))

  ## via11（审图意见⑤⑥，用户定规「保留 E，仅做区分」）：
  ## 色阶与 F 完全统一为 p.adjust（low=case/high=ctrl）；
  ## 降透明度与点径作视觉降级——E 仅保留为轨级概览。
  p <- ggplot(df, aes(track_disp, term)) +
    geom_point(aes(size = count, color = padj), alpha = 0.78) +
    scale_color_gradient(low = CR_COL$case, high = CR_COL$ctrl, name = "p.adjust") +
    scale_size_continuous(range = c(1.2, 4), name = "Count") +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, 44)) +   # via13：少折一行
    labs(title = "Track-level enrichment (ORA, FDR < 0.05)",
         x = NULL, y = NULL) +   # via19③：删生产备注副标题（色阶统一口径移入图注）
    theme_cr(base_size = 8) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = rel(0.72)),
          axis.text.y = element_text(size = rel(0.65), lineheight = 1.05),   # via13
          plot.title = element_text(size = rel(0.95), margin = margin(b = 4)),   # via16
          plot.subtitle = element_text(margin = margin(b = 5)),                  # via16
          plot.title.position = "plot",   # via10：标题从整幅图左缘起笔，
                                          # 不再被 y 轴长标签挤出版心右缘
          legend.key.height = grid::unit(3, "mm"),
          plot.margin = margin(6, 8, 6, 14))   # via17：统一边距

  ## FoxO 行标注（机制检查点；真实存在才标）
  ## via5：nudge_x 由 0.55 改 -0.3——原值把标注推出面板右缘被裁（via4 实图确认）
  ix <- grepl("FoxO", df$term, ignore.case = TRUE)
  if (any(ix)) {
    fx <- df[which(ix)[1], ]
    fx$lab <- sprintf("FoxO signaling\np.adj = %.3f", fx$padj)
    p <- p + geom_text_repel(data = fx, aes(label = lab), size = 2.0,
                             fontface = "bold", color = CR_COL$case,
                             nudge_x = -0.3, direction = "y",
                             min.segment.length = 0, segment.size = 0.25,
                             show.legend = FALSE, inherit.aes = TRUE)
  }
  p
}

## ---- 9. 面板 (f)：compareCluster dotplot（database × track 分面）-------------------
panel_f <- function(en) {
  df <- en[order(en$padj), ]
  df <- do.call(rbind, lapply(split(df, list(df$database, df$track), drop = TRUE),
                              function(x) head(x, CONFIG$top_n_f)))
  df <- df[order(df$database, as.numeric(df$track), df$padj,
                 decreasing = c(FALSE, FALSE, TRUE), method = "radix"), ]
  df$term <- factor(df$term, levels = unique(df$term))
  df <- df[!is.na(df$term), ]
  if (all(is.na(df$gene_ratio))) {
    message("[panel_f] 无 GeneRatio 列，x 轴改用 -log10(p.adj)。")
    df$gene_ratio <- df$neglog
    x_lab <- expression(-log[10](p.adj))
  } else {
    x_lab <- "GeneRatio"
  }
  ## via6：分面行序固定 GO-BP → Reactome → KEGG，与标题声明顺序一致
  db_lv <- unique(as.character(df$database))
  db_rank <- ifelse(grepl("^GO", db_lv, ignore.case = TRUE), 1,
             ifelse(grepl("reactome", db_lv, ignore.case = TRUE), 2,
             ifelse(grepl("kegg", db_lv, ignore.case = TRUE), 3, 4)))
  df$database <- factor(as.character(df$database), levels = db_lv[order(db_rank, db_lv)])

  ggplot(df, aes(gene_ratio, term)) +
    geom_point(aes(size = count, color = padj), alpha = 0.9) +
    facet_grid(database ~ track, scales = "free", space = "free") +
    ## via6 三连修：① y 扩 1.2 行距 + clip="off" + 上边距 8pt——Count=7 大气泡
    ## 半径超出 0.9 行距仍被 strip 顶裁（用户实图反馈）；② 右边距 18pt——
    ## "Count" 图例标题贴边被裁；③ 行序见上。
    scale_y_discrete(expand = expansion(add = 1.2)) +
    coord_cartesian(clip = "off") +
    scale_color_gradient(low = CR_COL$case, high = CR_COL$ctrl, name = "p.adjust") +
    ## via10：右侧图例盒总高超过画布可用高度，上下对称溢出——"Count" 图例
    ## 标题顶裁（用户实图反馈）。压缩：Count 仅 2/4/6 三档、色条限高 16mm。
    scale_size_continuous(range = c(1.2, 4.5), breaks = c(2, 4, 6), name = "Count") +
    guides(size = guide_legend(order = 1),
           color = guide_colorbar(barheight = grid::unit(16, "mm"),
                                  barwidth = grid::unit(2.5, "mm"), order = 2)) +
    labs(title = "compareCluster dotplot (GO-BP / Reactome / KEGG)",
         x = x_lab, y = NULL) +
    theme_cr(base_size = 8) +
    theme(strip.text = element_text(size = rel(0.75)),
          axis.text.y = element_text(size = rel(0.7)),
          axis.text.x = element_text(size = rel(0.7)),
          ## via15：标题与列 strip 文字相碰（用户实图反馈）——缩小一号 +
          ## 下边距 5pt 撑开垂直间隙；position="plot" 从整幅图左缘起笔，
          ## 避开 strip 密集区（与 panel_e via10 同款处理）。
          plot.title = element_text(size = rel(0.9), margin = margin(b = 5)),
          plot.title.position = "plot",
          ## via9：下边距 4→14pt——p.adjust 色条直抵画布下缘易被误认为裁切
          ## 不全（用户实图反馈），色条底部留白悬空。
          plot.margin = margin(8, 18, 14, 14))   # via17：统一左边距（顶 8pt 保留防气泡顶裁）
}

## ---- 10. 主流程（via5 化整为零版）：6 面板各自独立出小图，拼接后续专门处理 -------
## 背景：整机 patchwork 在 RStudio 反复触发 "Plots window may be too small"（设备层
## 问题，非代码问题）。via4 改为逐面板 save_cr，每个面板独立 tryCatch——单个面板
## 失败不影响其它面板，末尾打印成败清单。拼接脚本随后另出。
.pdf_null <- grDevices::pdf(NULL)   # 空设备兜底：屏幕绘制一律丢弃

out_dir <- rp(CONFIG$out_dir)
panel_results <- character(0)

## make_p 惰性求值：面板构建也在 tryCatch 内，构建/渲染任一环节报错都被隔离
save_panel <- function(make_p, tag, w, h) {
  fn <- file.path(out_dir, paste0("Fig2_panel", tag))
  tryCatch({
    p <- make_p()
    ## via6：字母不再 R 端盖章——tag 落在标题/副标题上造成重叠（用户实图反馈）。
    ## 按用户定规：小图保持无字母，拼接完成后由 29 号脚本在每行上方
    ## 专属留白条内、对齐各面板左缘加序号字母。
    save_cr(p, fn, width_mm = w, height_mm = h)
    panel_results <<- c(panel_results, paste0("✅ panel ", tag))
  }, error = function(e) {
    panel_results <<- c(panel_results,
                        paste0("❌ panel ", tag, "：", conditionMessage(e)))
  })
}

## via12 画布（审图意见⑬试案）：整版高度压缩至 ≈232mm——
## A 44→36，中段四面版 72→58，F 80→62；配合 29 via6 字母条 48px + 行缝 8px。
## 面板内部字号不变（8pt），若出现拥挤/裁切再回退或改拆 F 方案。
save_panel(panel_a, "A_批次校正前后PCA", 174, 36)
save_panel(panel_b, "B_代表队列火山图",  85, 58)
save_panel(panel_c, "C_双轨镜像热图",    85, 58)
save_panel(panel_d, "D_镜像占比条形图",  85, 58)
en <- load_enrich()
save_panel(function() panel_e(en), "E_分轨富集气泡图", 85, 58)
save_panel(function() panel_f(en), "F_富集compareCluster点图", 174, 62)

invisible(grDevices::dev.off())   # 关闭空设备
cat("\n========== Fig 2 逐面板出图清单 ==========\n")
cat(paste(panel_results, collapse = "\n"), "\n")
message("[18_fig2] 面板小图输出目录：", out_dir,
        "（拼接脚本随后另出，先把 6 张小图拿到手）")
