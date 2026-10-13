# =============================================================================
# 20_fig4_外部验证_ClinRheum风格_20260915_via9.R
# via9（2026-10-02，B 面板误差棒视觉强化，纯样式零统计改动）：
#   用户反馈部分误差棒看似显示不全——核查确认：观测率=100% 的分箱其精确
#   CI 上界=1.0 与点估计重合，上半段数学上为零长（非渲染缺陷）；本版加粗
#   误差棒（linewidth 0.4→0.9）、加宽帽线（width 0.02→0.045）、加深颜色
#   （grey40→grey20），让已有的下半段在印刷尺寸下清晰可辨。
# via8（2026-09-30，B 面板两处用户反馈修复）：
#   ① CI 算法 Wald → 精确二项（Clopper-Pearson, binom.test）：高患病率队列顶部
#     分箱观测 p=1.0 时 Wald SE=0 致误差棒退化为零长不可见；精确 CI 恒非退化
#     （如 8/8 → [0.63,1.00]），小样本校准图标准做法；
#   ② 边缘呼吸空间：coord_cartesian(0–1 硬边界) 改 scale limits + expand 4%，
#     落在 x=0 或 (1,1) 边界的分箱点不再被面板边缘裁切一半。
# via7（2026-09-23）：主题引用改为自动探测最高 via —— 配合 16_style_theme via3（TIFF 600dpi，
#   刊规组合彩图要求），本脚本内容逻辑零改动，重跑即得 600dpi 面板/整版。
# via6 修订（2026-09-22，用户三条全局规范 + 范文 Fig3 配色核对）：
#   ① 规格统一：A ROC 弃 base plot 改 ggplot geom_step——四面板同 ggplot+
#     theme_cr(base_size=8)，直角坐标系规格/字体层级一致；B y 刻度改 0.2 步长
#     一位小数，与 A/C 轴标签同宽（绘图区宽度对齐）
#   ② 配色对齐范文 Fig3：DCA 全实线 model=红/All=grey70 浅灰/None=grey30 深灰
#     （范文 None 为深灰实线，非虚线）；ROC 主曲线红+细灰对角线（范文同款）
#   ③ 字体层级：轴题/轴字=theme_cr 8pt 层级；图例 rel(0.85)；图内标注 2.4–2.6mm
# via5 修订（2026-09-22，B 面板终审，用户拍板）：
#   ② B 校准：删 loess 平滑线（中段 >1 被裁致断裂 + 0.8 处伪下凹，
#     loess 用于二分类小样本固有毛病）——只留分箱点+误差棒+对角线+inset
#   ③ 新增控制台实算校验：从校准数据重拟 logistic 再校准
#     （intercept: glm(y~1, offset=qlogis(p))；slope: glm(y~qlogis(p))），
#     与冻结常数 CAL_INTERCEPT=5.45 / CAL_SLOPE=0.268 比对打印
# via4 修订（2026-09-22，via3 实跑确认 A 面板翻车根因）：
#   ① A ROC：pROC legacy.axes 本质是"反向 usr + 标签翻转"，axes=FALSE 后
#     手工轴在反向 usr 下反向、legend 关键词定位持续漂移 → 弃 pROC 排版，
#     coords() 取 (specificity, sensitivity) 后 base plot 正向 usr 手工绘制，
#     轴/对角线/图例全受控（type="s" 阶梯线=ROC 惯例）
# via3 修订（2026-09-22，via2 实跑预览 2 处）：
#   ① A ROC：pROC 自动刻度在 wrap_elements 下退化（科学计数法 -2e-01/8e-01）
#     且 legend(bottomright) 定位漂移 → axes=FALSE + 手工 axis/box +
#     固定坐标手写图例（色块短线 + text），完全脱离 pROC 自动排版
#   ② D：hypothesis-generating 角标与末行柱标签同行歧义 → 改红色副标题
# via2 修订（2026-09-22，v1 实跑预览 4 项整改，用户已拍板）：
#   ① A ROC：legacy.axes 实际画的是 FPR 轴 → xlab 改 "1 − Specificity"
#   ② B 校准：分箱点+二项 95%CI 误差棒+细平滑线（se=FALSE 去粉色大 CI 带）；
#     inset（Slope/Intercept/MAE）移至右下角空白区，字号 3→2.6
#   ③ C DCA：弃 dcurves 原生灰网格 → as_tibble 取数后 ggplot+theme_cr 重绘；
#     Model=红实线 / Treat All=灰实线 / Treat None=灰虚线；x 限 0–80%
#     （阈值>0.8 小样本尖刺伪影区不展示，DCA 惯例临床相关区间）
#   ④ D 稳定性：36 基因→Top 20（用户定）；末端 % 标签统一右偏+缩至 2.4；
#     「hypothesis-generating」红角标从右上（压 100% 柱群）移至右下空白区
# -----------------------------------------------------------------------------
# 功能：主文 Fig 4「三层架构外部验证」出图（4 面板 A–D）
#   (a) 软骨主验证 ROC：pROC 原生风（legacy.axes 与范文一致），
#       GSE57218 主曲线（红）+ GSE168505 探索性降权（灰虚线），
#       AUC + 95%CI（DeLong）写入图内图例，如 "GSE57218: 0.818 (0.614–1.000)"
#   (b) 校准曲线：分箱观测点 + loess 平滑 + 对角虚线；图内 inset 标注
#       bootstrap B=200、slope=0.268（外部实算）、MAE（范文校准图内嵌参数同款）
#   (c) DCA：dcurves 原生风，treat-all / treat-none 灰线
#   (d) 血液子签名折内特征选择稳定性条形图，右上角红字
#       「hypothesis-generating」角标
# 风格锚点：Clin Rheumatol 2026;45 RA 锰代谢×ML 范文（DOI 10.1007/s10067-026-08081-3）
#   Fig3（ROC 图内图例 / 校准曲线内嵌 bootstrap 参数 / DCA 灰线）
#   锚点图版：稿件工作区/参考文章原文/RA_manganese_ML_figs/img_p07_x105.png
# 内容依据：论文初稿_中文主稿_v6（Fig 4 图注，方案 A）
# 输入（路径见 CONFIG 块；候选位置自动探测，找不到则报错并列出已尝试路径）：
#   1) 外部验证预测值表（用户本地，主路长表 extval_predictions.csv）：
#        cohort（含 GSE57218 / GSE168505 字样）, prob（模型预测概率）, label（OA/Control 或 0/1）
#      备路：pred_GSE57218.csv + pred_GSE168505.csv 两个单队列文件
#        （各含 prob + label 两列；GSE168505 可缺省 → 灰线自动跳过并提示）
#   2) 校准/DCA 数据：默认复用 GSE57218 预测值；亦可用 calib_data.csv /
#      dca_data.csv 单独提供（列同 prob + label）
#   3) 血液子签名稳定性表 stability_blood.csv（用户本地）：
#        gene + freq（0–1 或 0–100 百分比均可）；
#        或 gene + n_selected + n_folds（脚本换算 freq = n_selected / n_folds）
# 输出：figures/ClinRheum风格/Fig4_三层外部验证.pdf / .tiff（174×140 mm）
# 作者：稿件出图流水线（脚本 20）  日期：2026-09-15  种子：20260915
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(pROC)      # ROC 原生风 + DeLong 95%CI
  library(dcurves)   # DCA 原生风（treat-all/none 灰线）
  library(scales)    # alpha 透明度
})

set.seed(20260915)

## ---------------------------------------------------------------------------
## 0. 定位脚本目录 → source 共享主题
## ---------------------------------------------------------------------------
this_file <- tryCatch(
  normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])),
  error = function(e) NA_character_
)
script_dir <- if (length(this_file) == 0 || is.na(this_file) || !file.exists(this_file)) {
  getwd()
} else {
  dirname(this_file)
}
theme_cand <- c(list.files(script_dir, pattern = "^16_style_theme.*\\.R$", full.names = TRUE),
                list.files(getwd(), pattern = "^16_style_theme.*\\.R$", full.names = TRUE))
if (length(theme_cand) == 0)
  stop("找不到共享主题 16_style_theme 脚本；请将本脚本与主题文件放在同一目录后重跑。")
theme_path <- sort(theme_cand, decreasing = TRUE)[1]   # 文件名最大 = 最高 via（via3 起 TIFF=600dpi）
source(theme_path)

ROOT <- locate_root()
OUT_DIR  <- file.path(ROOT, "figures", "ClinRheum风格")
OUT_BASE <- file.path(OUT_DIR, "Fig4_三层外部验证")

## ---------------------------------------------------------------------------
## 1. CONFIG —— 输入路径（候选制）+ 论文既定常数（必须出现在图内）
## ---------------------------------------------------------------------------

find_input <- function(cands, label) {
  tried <- character(0)
  for (cnd in cands) {
    for (p in c(cnd, file.path(script_dir, cnd), file.path(ROOT, cnd))) {
      tried <- c(tried, p)
      if (file.exists(p)) {
        message("[CONFIG] ", label, " -> ", normalizePath(p))
        return(normalizePath(p))
      }
    }
  }
  stop("[CONFIG] 找不到输入文件【", label, "】。\n",
       "请将文件放到下列任一位置（或修改脚本 CONFIG 块的候选路径）后重跑：\n  ",
       paste(unique(tried), collapse = "\n  "), call. = FALSE)
}

find_optional <- function(cands) {
  for (cnd in cands) {
    for (p in c(cnd, file.path(script_dir, cnd), file.path(ROOT, cnd))) {
      if (file.exists(p)) return(normalizePath(p))
    }
  }
  NA_character_
}

CFG <- list(
  ## 主路：长表（cohort + prob + label）
  pred_long = find_optional(c(
    "results/moduleE/extval_predictions.csv",
    "results/moduleI/extval_predictions.csv",
    "moduleI_data/extval_predictions.csv",
    "extval_predictions.csv"
  )),
  ## 备路：单队列文件（prob + label）
  pred_57218  = c("复核输出_20260916/PartB0_GSE57218预测概率.csv",  # 2026-09-16 复核实算（sample/prob/label）
                  "results/moduleE/pred_GSE57218.csv",  "moduleI_data/pred_GSE57218.csv",  "pred_GSE57218.csv"),
  pred_168505 = c("results/moduleE/pred_GSE168505.csv", "moduleI_data/pred_GSE168505.csv", "pred_GSE168505.csv"),
  ## 校准/DCA 独立数据（可缺省 → 默认复用 GSE57218 预测值）
  calib_csv = c("results/moduleE/calib_data.csv", "moduleI_data/calib_data.csv", "calib_data.csv"),
  dca_csv   = c("results/moduleE/dca_data.csv",   "moduleI_data/dca_data.csv",   "dca_data.csv"),
  ## 血液子签名稳定性表（必需）
  stab_csv = find_input(c(
    "results/moduleD_v2/blood_exploratory/EXPLORATORY_blood_fold_selection_stability.csv",  # 05c 实际产出（gene + folds_selected + stability）
    "results/moduleE/stability_blood.csv",
    "results/moduleD/stability_blood.csv",
    "moduleI_data/stability_blood.csv",
    "stability_blood.csv"
  ), "血液子签名稳定性表（期望列：gene + freq，或 gene + n_selected + n_folds）")
)

## 论文既定常数（写入图内 inset / 用于控制台一致性校验）
BOOT_B   <- 200L        # bootstrap 重复次数（内部 validate 口径）
CAL_SLOPE <- 0.268      # 外部校准斜率（GSE57218 实算，2026-09-16 复核；图内 inset 必须出现）
CAL_INTERCEPT <- 5.45   # 外部校准截距（同次复核）
EXPECT <- list(
  auc_57218   = 0.818,                    # GSE57218 主验证 AUC
  ci_57218    = c(0.614, 1.000),          # 95%CI（DeLong 实算 2026-09-16，重建模型预测）
  auc_168505  = 0.750,                    # GSE168505 探索性 AUC（n=7，降权灰线）
  stab_top5   = c("CALML4", "DCLRE1C", "GADD45A", "KLF2", "NFIC")  # 折内稳定前 5 基因
)

## ---------------------------------------------------------------------------
## 2. 小工具
## ---------------------------------------------------------------------------

## 分组标签归一化为 0/1（1=OA 病例）；无法识别即明确报错
norm_label <- function(x, label = "label") {
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x) && all(na.omit(unique(x)) %in% c(0, 1))) return(as.integer(x))
  s <- tolower(trimws(as.character(x)))
  is_case <- s %in% c("oa", "case", "cases", "patient", "disease", "1", "oa_case")
  is_ctrl <- s %in% c("control", "ctrl", "healthy", "normal", "hc", "0", "non-oa")
  if (all(is_case | is_ctrl)) return(as.integer(is_case))
  stop("分组标签无法归一化为 0/1 [", label, "]：唯一值 = ",
       paste(unique(as.character(x)), collapse = ", "),
       "（期望 OA/Control 或 0/1）", call. = FALSE)
}

## 读取单队列 prob+label 文件
read_pred_file <- function(path, cohort_name) {
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  pc <- pick_col(df, c("prob", "pred", "prediction", "score", "risk", "fitted", "probability"),
                 label = paste0(cohort_name, " 预测概率列"))
  lc <- pick_col(df, c("label", "group", "status", "y", "outcome", "class"),
                 label = paste0(cohort_name, " 标签列"))
  data.frame(cohort = cohort_name,
             prob   = as.numeric(df[[pc]]),
             label  = norm_label(df[[lc]], label = cohort_name))
}

## 概率取值域校验（应为 0–1）
check_prob <- function(df, label) {
  if (any(df$prob < 0 | df$prob > 1, na.rm = TRUE)) {
    warning("[校验] ", label, " 的预测概率超出 [0,1]，请确认传入的是概率而非线性预测值")
  }
  invisible(df)
}

## ---------------------------------------------------------------------------
## 3. 数据读取（防御性）
## ---------------------------------------------------------------------------

## 3.1 外部验证预测值：长表主路 / 单队列文件备路；GSE168505 可缺省
load_predictions <- function() {
  if (!is.na(CFG$pred_long)) {
    df <- read.csv(CFG$pred_long, check.names = FALSE, stringsAsFactors = FALSE)
    co <- pick_col(df, c("cohort", "dataset", "gse", "series", "study"),
                   label = "队列举")
    pc <- pick_col(df, c("prob", "pred", "prediction", "score", "risk", "fitted", "probability"),
                   label = "预测概率列")
    lc <- pick_col(df, c("label", "group", "status", "y", "outcome", "class"),
                   label = "标签列")
    out <- data.frame(cohort = toupper(as.character(df[[co]])),
                      prob   = as.numeric(df[[pc]]),
                      label  = norm_label(df[[lc]]))
    ## 队列名宽松匹配（允许 "GSE57218" / "57218" / "gse57218" 等写法）
    out$cohort <- ifelse(grepl("57218", out$cohort), "GSE57218",
                  ifelse(grepl("168505", out$cohort), "GSE168505", out$cohort))
  } else {
    f1 <- find_optional(CFG$pred_57218)
    if (is.na(f1)) {
      stop("[CONFIG] 外部验证预测值缺失：既无 extval_predictions.csv 长表，",
           "也无 pred_GSE57218.csv。\n期望列：prob（预测概率）+ label（OA/Control 或 0/1）；",
           "长表另需 cohort 列（含 GSE57218 / GSE168505 字样）。", call. = FALSE)
    }
    out <- read_pred_file(f1, "GSE57218")
    f2 <- find_optional(CFG$pred_168505)
    if (!is.na(f2)) out <- rbind(out, read_pred_file(f2, "GSE168505"))
  }
  if (!any(out$cohort == "GSE57218")) {
    stop("[CONFIG] 预测值表中未找到 GSE57218 主验证队列记录。", call. = FALSE)
  }
  if (!any(out$cohort == "GSE168505")) {
    message("[预测值] 未提供 GSE168505 → 面板 A 仅画 GSE57218 主曲线（无探索性灰线）")
  }
  check_prob(out, "外部验证预测值")
  message("[预测值] ", paste(sprintf("%s n=%d", unique(out$cohort),
          as.integer(table(out$cohort)[unique(out$cohort)])), collapse = "；"))
  out
}

## 3.2 校准/DCA 数据：独立文件优先，否则复用 GSE57218 预测值
load_calib_dca <- function(pred) {
  d57218 <- pred[pred$cohort == "GSE57218", ]
  f_cal <- find_optional(CFG$calib_csv)
  calib <- if (!is.na(f_cal)) read_pred_file(f_cal, "calibration") else d57218
  f_dca <- find_optional(CFG$dca_csv)
  dca   <- if (!is.na(f_dca)) read_pred_file(f_dca, "DCA") else d57218
  check_prob(calib, "校准数据"); check_prob(dca, "DCA 数据")
  list(calib = calib, dca = dca)
}

## 3.3 血液子签名稳定性表：freq（0–1 或百分比）或 n_selected/n_folds 换算
load_stability <- function(path) {
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  gc <- pick_col(df, c("gene", "symbol", "feature"), label = "稳定性表 gene 列")
  fc <- pick_col(df, c("freq", "frequency", "stability", "selection_freq",
                       "prop", "proportion", "rate"), required = FALSE,
                 label = "稳定性频率列")
  if (!is.null(fc)) {
    freq <- as.numeric(df[[fc]])
    if (max(freq, na.rm = TRUE) > 1.0001) {          # 兼容 0–100 百分比写法
      message("[稳定性] 频率列最大值 >1，按 0–100 百分比处理（除以 100）")
      freq <- freq / 100
    }
  } else {
    nc <- pick_col(df, c("n_selected", "times_selected", "count", "n", "times"),
                   label = "入选次数列")
    kc <- pick_col(df, c("n_folds", "folds", "n_fold", "n_cv"), label = "折数列")
    freq <- as.numeric(df[[nc]]) / as.numeric(df[[kc]])
  }
  out <- data.frame(gene = as.character(df[[gc]]), freq = freq)
  out <- out[stats::complete.cases(out), ]
  if (any(out$freq < 0 | out$freq > 1)) {
    stop("[稳定性] 频率换算后超出 [0,1]，请复核 n_selected / n_folds 列", call. = FALSE)
  }
  out <- out[order(-out$freq), ]
  ## 与论文「折内稳定前 5 基因」交叉校验（仅提示，不阻断）
  top5 <- head(out$gene, 5)
  if (!identical(sort(top5), sort(EXPECT$stab_top5))) {
    warning("[校验] 稳定性前 5 基因为 ", paste(top5, collapse = "/"),
            "，与论文值 ", paste(EXPECT$stab_top5, collapse = "/"), " 不一致，请复核输入")
  }
  message("[稳定性] ", nrow(out), " 基因；前 5：", paste(top5, collapse = "/"))
  out
}

## ---------------------------------------------------------------------------
## 4. 面板函数（图内文字全英文；统计量直接入图——范文同款）
## ---------------------------------------------------------------------------

## (a) 软骨主验证 ROC：pROC 原生风（legacy.axes，与范文 Fig3A 一致）；
##     AUC + DeLong 95%CI 写入图内图例；GSE168505 探索性降权 = 灰虚线
panel_a <- function(pred) {
  d1 <- pred[pred$cohort == "GSE57218", ]
  roc1 <- pROC::roc(factor(d1$label, levels = c(0, 1), labels = c("Control", "OA")),
                    d1$prob, quiet = TRUE)
  ci1 <- pROC::ci.auc(roc1, method = "delong")   # c(5%, AUC, 95%)
  txt1 <- sprintf("GSE57218: %.3f (%.3f\u2013%.3f)", ci1[2], ci1[1], ci1[3])
  ## 与论文既定数字一致性校验（图内数字始终以数据实算为准）
  if (abs(ci1[2] - EXPECT$auc_57218) > 0.01 ||
      abs(ci1[1] - EXPECT$ci_57218[1]) > 0.02 || abs(ci1[3] - EXPECT$ci_57218[2]) > 0.02) {
    warning("[校验] GSE57218 实算 ", txt1, " 与论文值 0.818 (0.614–1.000) 不一致")
  }
  has2 <- any(pred$cohort == "GSE168505")
  roc2 <- NULL; txt2 <- NULL
  if (has2) {
    d2 <- pred[pred$cohort == "GSE168505", ]
    roc2 <- pROC::roc(factor(d2$label, levels = c(0, 1), labels = c("Control", "OA")),
                      d2$prob, quiet = TRUE)
    txt2 <- sprintf("GSE168505: %.3f (exploratory)", as.numeric(pROC::auc(roc2)))
    if (abs(as.numeric(pROC::auc(roc2)) - EXPECT$auc_168505) > 0.02) {
      warning("[校验] GSE168505 实算 AUC 与论文值 ", EXPECT$auc_168505, " 不一致")
    }
  }
  ## via4①：coords() 取坐标 → 正向 usr 手工绘制（弃 pROC 自动排版）
  co1 <- pROC::coords(roc1, "all", ret = c("specificity", "sensitivity"),
                      transpose = FALSE)
  fpr1 <- 1 - co1$specificity; sen1 <- co1$sensitivity
  o1 <- order(fpr1, sen1); fpr1 <- fpr1[o1]; sen1 <- sen1[o1]
  if (has2) {
    co2 <- pROC::coords(roc2, "all", ret = c("specificity", "sensitivity"),
                        transpose = FALSE)
    fpr2 <- 1 - co2$specificity; sen2 <- co2$sensitivity
    o2 <- order(fpr2, sen2); fpr2 <- fpr2[o2]; sen2 <- sen2[o2]
  }
  ## via6①：ggplot geom_step（hv 阶梯 = ROC 惯例）——与 B/C/D 同 theme_cr 规格
  roc_df <- data.frame(fpr = fpr1, sens = sen1, series = txt1)
  if (has2) roc_df <- rbind(roc_df, data.frame(fpr = fpr2, sens = sen2, series = txt2))
  lev <- if (has2) c(txt1, txt2) else txt1
  roc_df$series <- factor(roc_df$series, levels = lev)
  ggplot(roc_df, aes(fpr, sens, color = series, linetype = series)) +
    geom_abline(slope = 1, intercept = 0, color = "grey85", linewidth = 0.5) +  # 范文细灰对角线
    geom_step(linewidth = 0.7) +
    scale_color_manual(values = stats::setNames(c(CR_COL$case, "grey60")[seq_along(lev)], lev),
                       name = NULL) +
    scale_linetype_manual(values = stats::setNames(c("solid", "dashed")[seq_along(lev)], lev),
                          name = NULL) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "1 \u2212 Specificity", y = "Sensitivity") +
    theme_cr(base_size = 8) +
    theme(legend.position = c(0.98, 0.03), legend.justification = c(1, 0),
          legend.text = element_text(size = rel(0.85)),
          legend.key.width = grid::unit(5, "mm"))
}

## (b) 校准曲线：分箱观测点 + loess 平滑 + 对角虚线；
##     图内 inset：bootstrap B=200、slope=0.268（外部实算）、MAE（范文内嵌参数同款）
panel_b <- function(calib) {
  df <- data.frame(prob = calib$prob, y = calib$label)
  ## 分箱数防御：n≈40 → 5 箱；样本量过小自动降档，保证每箱非空
  g <- max(3L, min(5L, floor(nrow(df) / 8L)))
  br <- unique(stats::quantile(df$prob, probs = seq(0, 1, length.out = g + 1L),
                               na.rm = TRUE))
  bin <- cut(df$prob, breaks = br, include.lowest = TRUE)
  agg <- data.frame(pred = as.numeric(tapply(df$prob, bin, mean)),
                    obs  = as.numeric(tapply(df$y,    bin, mean)),
                    n    = as.numeric(tapply(df$y,    bin, length)))
  ## via2②：分箱点 + 95%CI 误差棒（范文主流校准图样式，替代大 CI 带）
  ## via9：obs=1.0 分箱的 CI 上界=1.0 与点重合 → 上半段零长为数学必然，
  ##         图注已加说明句；此处仅强化已有棒段的视觉可读性
  ## via8①：Wald → 精确二项（Clopper-Pearson）——obs=0 或 1 时 Wald SE=0
  ##         会退化成零长不可见误差棒；精确法恒给出非退化区间
  ci <- t(mapply(function(k, n) stats::binom.test(k, n)$conf.int,
                 round(agg$obs * agg$n), agg$n))
  agg$lo <- ci[, 1]
  agg$hi <- ci[, 2]
  mae <- mean(abs(agg$obs - agg$pred), na.rm = TRUE)
  inset <- sprintf("Slope = %.3f\nIntercept = %.2f\nMAE = %.3f",
                   CAL_SLOPE, CAL_INTERCEPT, mae)  # 外部校准实算参数（2026-09-16 复核）
  ggplot(df, aes(prob, y)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey60") +
    ## via5②：不再画 loess 平滑线（中段冲出 y=1 被裁致断裂 + 小样本伪下凹）
    geom_errorbar(data = agg, aes(x = pred, ymin = lo, ymax = hi),
                  width = 0.045, color = "grey20", linewidth = 0.9,
                  inherit.aes = FALSE) +                     # via9：加粗加深+宽帽线
    geom_point(data = agg, aes(pred, obs), size = 2, color = "black") +
    annotate("text", x = 0.98, y = 0.02, label = inset,          # via2②：inset 移至右下角空白区
             hjust = 1, vjust = 0, size = 2.6, lineheight = 1.25) +
    ## via8②：limits 仍在 0–1（刻度不变），expand 4% 给边界点留呼吸空间，
    ##         边界分箱点（x=0 或角点 (1,1)）不再被面板边缘裁半
    scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0.04)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1),
                       expand = expansion(mult = 0.04)) +   # via6③：y 刻度 0.2 步长一位小数，与 A/C 轴标签同宽
    labs(x = "Predicted Probability", y = "Observed Probability") +
    theme_cr(base_size = 8)                          # via6③：全图统一 8pt 基准
}

## (c) DCA：dcurves 计算净获益（计算引擎不变），as_tibble 取数后
##     ggplot + theme_cr 重绘（via2③：弃原生灰网格，与全图主题统一）；
##     x 轴限 0–80% 临床相关阈值区间（>0.8 段为 n≈40 小样本尖刺伪影）
panel_c <- function(dca_df) {
  dd <- data.frame(outcome = dca_df$label, Model = dca_df$prob)
  d <- dcurves::dca(outcome ~ Model, data = dd)
  tb <- tryCatch(tibble::as_tibble(d), error = function(e) d$dca)  # 对象结构防御
  idc <- if ("label" %in% names(tb)) "label" else "variable"       # 曲线标识列名防御
  tb <- tb[!is.na(tb$net_benefit), ]
  nm <- as.character(tb[[idc]])
  nm <- ifelse(grepl("model", nm, ignore.case = TRUE), "Model",
        ifelse(grepl("all",   nm, ignore.case = TRUE), "Treat All",
        ifelse(grepl("none",  nm, ignore.case = TRUE), "Treat None", nm)))
  tb$series <- factor(nm, levels = c("Model", "Treat All", "Treat None"))
  ggplot(tb, aes(threshold, net_benefit, color = series)) +
    geom_line(linewidth = 0.6) +
    scale_color_manual(values = c(Model = CR_COL$case, "Treat All" = "grey70",
                                  "Treat None" = "grey30"), name = NULL) +   # via6②：范文 DCA 配色（红/浅灰/深灰全实线）
    coord_cartesian(xlim = c(0, 0.8), ylim = c(-0.1, NA)) +        # via2③
    scale_x_continuous(labels = function(z) sprintf("%.0f%%", z * 100)) +
    labs(x = "Threshold Probability", y = "Net Benefit") +
    theme_cr(base_size = 8) +                        # via6③
    theme(legend.position = c(0.98, 0.72), legend.justification = c(1, 0.5),
          legend.text = element_text(size = rel(0.85)))
}

## (d) 血液子签名折内特征选择稳定性条形（via2④ 重画版）：
##     Top 20（用户定 2026-09-22，36 基因在半高面板必然互叠）；
##     基因名走 y 轴外侧、末端 % 标签统一右偏 2.4 号；
##     「hypothesis-generating」红角标移右下空白区（低频区右侧无柱）
panel_d <- function(stab) {
  df <- head(stab, 20)                                # via2④：Top 20（stab 已按 freq 降序）
  df$gene <- factor(df$gene, levels = rev(df$gene))   # 频次最高者置顶
  ggplot(df, aes(freq, gene)) +
    geom_col(fill = CR_COL$case, width = 0.7) +
    geom_text(aes(label = sprintf("%.0f%%", freq * 100)),
              hjust = -0.15, size = 2.4) +            # via2④：末端标签统一右偏、缩字号
    scale_x_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.25),
                       labels = function(z) sprintf("%.0f%%", z * 100),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Selection frequency across CV folds", y = NULL,
         subtitle = "hypothesis-generating") +        # via3②：红角标改副标题，脱离柱区
    theme_cr(base_size = 8) +                        # via6③
    theme(plot.subtitle = element_text(color = CR_COL$case, face = "bold",
                                       hjust = 1))    # via3②：红色右对齐
}

## ---------------------------------------------------------------------------
## 4.5 via5③：logistic 再校准实算校验（与冻结常数比对，仅打印不阻断）
##     intercept = glm(y ~ 1, offset = qlogis(p))；slope = glm(y ~ qlogis(p))
##     p=0/1 边界值截断至 (1e-4, 1-1e-4) 防 qlogis  Inf
check_calib_constants <- function(calib) {
  p <- as.numeric(calib$prob); y <- as.numeric(calib$label)
  n_clip <- sum(p <= 0 | p >= 1, na.rm = TRUE)
  p <- pmin(pmax(p, 1e-4), 1 - 1e-4)
  lp <- qlogis(p)
  fit_i <- glm(y ~ 1,  offset = lp, family = binomial)
  fit_s <- glm(y ~ lp, family = binomial)
  est_int <- unname(coef(fit_i)[1])
  est_slp <- unname(coef(fit_s)[2])
  message(sprintf("[校验] logistic 再校准实算：intercept = %.3f（冻结 %.2f）；slope = %.3f（冻结 %.3f）%s",
                  est_int, CAL_INTERCEPT, est_slp, CAL_SLOPE,
                  if (n_clip > 0) sprintf("；%d 个 p=0/1 已截断", n_clip) else ""))
  if (abs(est_int - CAL_INTERCEPT) > 0.5 || abs(est_slp - CAL_SLOPE) > 0.05) {
    warning("[校验] 实算再校准参数与冻结常数差异较大，请复核 CAL_SLOPE / CAL_INTERCEPT 口径",
            call. = FALSE)
  }
  invisible(list(intercept = est_int, slope = est_slp))
}

## ---------------------------------------------------------------------------
## 5. 主流程
## ---------------------------------------------------------------------------
main <- function() {
  message("========== Fig 4 出图开始 ==========")
  pred <- load_predictions()
  cd   <- load_calib_dca(pred)
  stab <- load_stability(CFG$stab_csv)
  check_calib_constants(cd$calib)   # via5③：控制台实算比对冻结常数

  fig <- (panel_a(pred) | panel_b(cd$calib)) /
         (panel_c(cd$dca) | panel_d(stab)) +
    plot_annotation(tag_levels = "A") & CR_TAG      # 左上角大写粗体字母标签
  save_cr(fig, OUT_BASE, width_mm = 174, height_mm = 140)
  message("========== Fig 4 出图完成 ==========")
  invisible(fig)
}

main()
