# =============================================================================
# 19_fig3_机器学习_诊断模型_ClinRheum风格_20260915_via12.R
# via12（2026-09-23）：主题引用改为自动探测最高 via —— 配合 16_style_theme via3（TIFF 600dpi，
#   刊规组合彩图要求），本脚本内容逻辑零改动，重跑即得 600dpi 面板/整版。
# via11（2026-09-21，字体统一）：
#   ⑫ B Venn：fontfamily/cat.fontfamily="sans"——VennDiagram 默认衬线体
#      （Times），与全图无衬线（Helvetica）及范文 Venn 不一致；纯字体统一，
#      布局参数不动。
# via10（2026-09-21，方案 A 第三轮：范文像素级对照后微调）：
#   ① lmgp 0.5→0.15——刻度数字由线下 4.0mm 上提至 ~2.6mm（范文 2.9mm），
#      数字归属自己的尺而非下一条；TP 数字与 LP 尺线间隙同时拉开；
#   ② xfrac 0.19→0.16 + cex.var 0.60→0.55——标签区 33→27mm（范文 25.6mm），
#      Linear Predictor 6.6pt 实测可放下；cex.axis 0.50、label.every=2 不变。
# via9（2026-09-21，方案 A 第二轮：按范文锚点图像素比例定字号）：
#   ⑩ F 列线图：实测 img_p07_x105——范文刻度数字仅约 4.8pt、变量名约 5.8pt，
#      via8 字号相对轴间距（6.3mm）过大才是压线/互叠主因
#      → cex.axis 0.50/cex.var 0.60（对齐范文比例）；
#      label.every=2 稀疏短轴刻度（治 GDE1/STEAP1/CBR3/RNF14/BNIP3 数字互叠）；
#      xfrac 0.17→0.19（最长标签 Linear Predictor 与首位刻度数字留净距）；
#      lmgp 回 0.5。
# via8（2026-09-21，列线图方案抉择=A：rms 原生刻度尺精修，范文 Fig3B 像素级同款）：
#   ⑨ F 列线图排版精修（rms 官方参数对症）：
#      xfrac 0.35→0.17——标签区约 61→30mm，基因名贴近各自轴线（意见②）；
#      面板 76→88mm + lmgp 0.4→0.55——轴间距 +16%，刻度数字与轴线分离（意见①）；
#      fun.at 加 0.01——LP 轴跨 -25~15，概率标签沿轴散开不再挤叠；
#      附：GADD45A 原无数字 = rms force.label=FALSE 碰撞自动省略，空间足即恢复；
#      A/B/C 70→64mm、D/E 70→66mm——整版约 232mm ≤ 234mm 预算再平衡。
# via7（2026-09-21，二轮目检 3 条：Venn 集合名压线 + 列线图版式重构）：
#   ⑥ B Venn：cat.dist c(0.28,0.28,0.14,0.14) + margin 0.17 + cat.cex 0.85——
#      LASSO/SVM_RFE 集合名原压在椭圆线上（文字-图形互叠）；
#   ⑦ 列线图 58mm 窄格放不下 14 条轴（刻度朝向看似不一、刻度文字互叠串行）→
#      版式重构：列线图升整宽底排面板 F（174×76mm），原 F 稳健性升 E；
#   ⑧ 版式配套：D/E 升 85×70mm、A/B/C 58×70mm，整版 174×约 229mm ≤ 234mm。
# via6（2026-09-16，整版目检 5 条审图意见重置，沿用 Fig2 定版经验）：
#   ① B Venn：margin 0.06→0.14、cex/cat.cex 1.1→0.9——集合名 LASSO/SVM_RF
#      原被画布边缘裁切（同 Fig2 A 图例裁切病）；
#   ② C 柱端数值：expand mult 0.14→0.22——GADD45A 数值原被右缘裁切；
#   ③ E 列线图：par(mai) 收紧 + cex.var 0.6/cex.axis 0.5——原全部轴标签
#      在 58×75mm 画布上互相叠印不可读；
#   ④ F 附注：图内 annotate（与图例互叠、被右缘裁切）→ labs(caption)
#      移至面板外下方（同 Fig2 标题/图例分区经验）；
#   ⑤ theme_cr() 裸调 9pt → theme_cr(base_size = 8) 全图统一（Fig2 铁律）。
# -----------------------------------------------------------------------------
# 功能：主文 Fig 3「机器学习标志物筛选与诊断模型」出图（6 面板 A–F）
#   (a) 四算法特征选择投票条形图（x=gene，y=入选算法数，单色 CR_COL$case）
#   (b) 四算法交集 Venn 图（VennDiagram 原生，粉/绿/蓝/紫 pastel 四圆）
#   (c) mean|SHAP| 全局重要性条形图（降序、单色、数值入图）
#   (d) SHAP 蜜蜂图（shapviz 原生，特征值 低→高 = 蓝→红）
#   (e) 10 标志物 logistic 列线图（rms::nomogram + base plot 原生输出）
#   (f) 七分类器稳健性（内部 CV AUC 与 GSE57218 外部 AUC 双系列点图，数值入图）
# 风格锚点：Clin Rheumatol 2026;45 RA 锰代谢×ML 范文（DOI 10.1007/s10067-026-08081-3）
#   Fig2（LASSO/SVM/RF/Venn/统计量入图）+ Fig3（列线图原生输出）
#   锚点图版：稿件工作区/参考文章原文/RA_manganese_ML_figs/img_p06_x93.png
#             稿件工作区/参考文章原文/RA_manganese_ML_figs/img_p07_x105.png
# 内容依据：论文初稿_中文主稿_v6（Fig 3 图注，方案 A）
# 输入（路径见下方 CONFIG 块；候选位置自动探测，全部找不到则报错并列出已尝试路径）：
#   1) moduleI_classifier_robustness.csv  【本仓库已存在】
#      期望列：classifier, cv_AUC, GSE57218_AUC, seed
#      （8 行；xgbTree 行 cv_AUC/GSE57218_AUC 为 NA，脚本自动剔除）
#   2) moduleI_single_gene_AUC.csv        【本仓库已存在】
#      期望列：gene, train_AUC, case_direction, GSE57218_AUC, GSE57218_CI
#      （10 行；CI 形如 "0.153-1.000"，脚本拆成 ci_low/ci_up 用于一致性校验）
#   3) 特征投票表（用户本地）：列 gene + 每个算法一列 0/1（列名即算法名，
#      如 LASSO / SVM_RFE / RF / XGBoost；3–5 个算法列均可，4 列为设计情形）
#   4) SHAP 输入（用户本地，二选一）：
#      - shap_object.rds：shapviz 对象 / kernelshap 对象 / 含 S+X 槽位的 list
#        （防御性读取：list 时探测 shapviz/kernelshap 元素及 S/shap_values/phi
#         与 X/feature_values 槽位；均无则明确报错）
#      - shap_long.csv：长表，列 gene + shap（+ 可选 feature_value、sample）；
#        有 feature_value 时蜜蜂图按特征值着色，否则退化为单色抖动点并提示
#   5) 列线图输入（用户本地，主+兜底两路）：
#      - 主路 nomogram_fit.rds：lrm/rms 拟合对象（或 list 包裹，自动探测）；
#        若对象缺 datadist 导致 nomogram() 失败，自动转入兜底路
#      - 兜底 train_expr_z.csv（行=样本，列=10 个标志物基因，z 分数表达，
#        可选首列 sample）+ train_group.csv（列 sample[可选] + label/group）
#        → 以 rms::lrm 重新拟合后建 datadist 再画列线图
# 输出：figures/ClinRheum风格/Fig3_机器学习与诊断模型.pdf / .tiff（174×150 mm）
# 作者：稿件出图流水线（脚本 19）  日期：2026-09-15  种子：20260915
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(VennDiagram)   # 四算法 Venn（原生 pastel 圆）
  library(shapviz)      # SHAP 蜜蜂图（sv_importance 原生）
  library(rms)          # 列线图（nomogram 原生输出）
})

set.seed(20260915)

## ---------------------------------------------------------------------------
## 0. 定位脚本目录 → source 共享主题（CR_COL/theme_cr/CR_TAG/cr_stars/save_cr/
##    pick_col/locate_root）；脚本与主题同目录，Rscript --file= 探测，回退 getwd
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

ROOT <- locate_root()                       # 向上找含 data/results/01_raw/GEO 的根
OUT_DIR  <- file.path(ROOT, "figures", "ClinRheum风格")
OUT_BASE <- file.path(OUT_DIR, "Fig3_机器学习与诊断模型")

## ---------------------------------------------------------------------------
## 1. CONFIG —— 集中声明全部输入路径（候选制：逐一探测，找不到则一次性报错）
## ---------------------------------------------------------------------------

## 输入文件定位器：按候选相对路径在「当前目录 / 脚本目录 / ROOT」三处探测
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

## 可选文件探测（不报错，返回 NA_character_）
find_optional <- function(cands) {
  for (cnd in cands) {
    for (p in c(cnd, file.path(script_dir, cnd), file.path(ROOT, cnd))) {
      if (file.exists(p)) return(normalizePath(p))
    }
  }
  NA_character_
}

CFG <- list(
  ## —— 已知存在于本仓库的两个 moduleI 中间结果 ——
  robust_csv = find_input(c(
    "results/moduleI/moduleI_classifier_robustness.csv",
    "moduleI_data/moduleI_classifier_robustness.csv",
    "备份/20260913_模块I完成_v5稿件整合/moduleI_classifier_robustness.csv"
  ), "七分类器稳健性表 moduleI_classifier_robustness.csv"),

  gene_auc_csv = find_input(c(
    "results/moduleI/moduleI_single_gene_AUC.csv",
    "moduleI_data/moduleI_single_gene_AUC.csv",
    "备份/20260913_模块I完成_v5稿件整合/moduleI_single_gene_AUC.csv"
  ), "单基因 AUC 表 moduleI_single_gene_AUC.csv"),

  ## —— 用户本地文件（缺失即在 CONFIG 阶段给出带期望列说明的报错）——
  vote_csv = find_input(c(
    "results/moduleD_v2/feature_votes_by_algorithm.csv",
    "results/moduleD/feature_votes.csv",
    "results/moduleD/moduleD_feature_votes.csv",
    "moduleI_data/feature_votes.csv",
    "feature_votes.csv"
  ), "四算法特征投票表（期望列：gene + 每算法一列 0/1，列名=算法名）"),

  ## SHAP / 列线图 rds 为「主路」，缺失不报错——load_shap()/load_nomogram()
  ## 内部会自动转入备路（shap_long.csv / train_expr_z+train_group 重建）
  shap_rds = find_optional(c(
    "results/moduleD/shap_object.rds",
    "moduleI_data/shap_object.rds",
    "shap_object.rds"
  )),

  nomo_rds = find_optional(c(
    "results/moduleI/nomogram_fit.rds",
    "moduleI_data/nomogram_fit.rds",
    "nomogram_fit.rds"
  )),

  ## 兜底重建列线图用（仅当 nomo_rds 缺失或缺 datadist 时使用；可缺省）
  train_expr  = c("results/moduleI/train_expr_z.csv",  "results/moduleD/train_expr_z.csv",  "moduleI_data/train_expr_z.csv",  "train_expr_z.csv"),
  train_group = c("results/moduleI/train_group.csv",   "results/moduleD/train_group.csv",   "moduleI_data/train_group.csv",   "train_group.csv")
)

## 已知关键数字（用于控制台一致性校验；图内数字一律由数据实时计算）
EXPECT <- list(
  n_markers   = 10L,                        # 10 个四算法交集标志物
  glm_ext_auc = 0.818,                      # 主模型 logistic 的 GSE57218 AUC
  cv_note     = "Internal CV is for reference only\n(feature selection performed prior to CV)"
)

## ---------------------------------------------------------------------------
## 2. 小工具（防御性）
## ---------------------------------------------------------------------------

## 0/1 归一化：logical/数值/字符（yes/true/selected…）统一成 0/1
to01 <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) return(as.integer(!is.na(x) & x != 0))
  s <- tolower(trimws(as.character(x)))
  as.integer(s %in% c("1", "yes", "y", "true", "selected", "sel", "+"))
}

## 分组标签归一化为 0/1（1=病例 OA，0=对照）；无法识别即明确报错
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

## ---------------------------------------------------------------------------
## 3. 数据读取（全部防御性：列名用 pick_col 探测，结构不符给中文报错）
## ---------------------------------------------------------------------------

## 3.1 特征投票表：gene + 每算法一列 0/1
load_votes <- function(path) {
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  gcol <- pick_col(df, c("gene", "symbol", "feature", "Gene"), label = "投票表 gene 列")
  names(df)[names(df) == gcol] <- "gene"
  ## 算法列 = 除 gene 外、内容可解释为 0/1 的列（3–5 列，设计为 4 列）
  cand <- setdiff(names(df), "gene")
  alg_cols <- cand[vapply(cand, function(cl) {
    v <- df[[cl]]
    if (is.logical(v)) return(TRUE)
    if (is.numeric(v)) return(all(na.omit(unique(v)) %in% c(0, 1)))
    all(na.omit(tolower(trimws(as.character(v)))) %in%
          c("0", "1", "yes", "no", "y", "n", "true", "false", "selected", ""))
  }, logical(1))]
  if (length(alg_cols) < 3L || length(alg_cols) > 5L) {
    stop("投票表算法列数异常（探测到 ", length(alg_cols), " 列：",
         paste(alg_cols, collapse = ", "), "）。期望 3–5 个 0/1 算法列（设计为 4：",
         "LASSO / SVM_RFE / RF / XGBoost）。实际列：",
         paste(names(df), collapse = ", "), call. = FALSE)
  }
  for (cl in alg_cols) df[[cl]] <- to01(df[[cl]])
  df$n_votes <- rowSums(df[, alg_cols, drop = FALSE], na.rm = TRUE)
  attr(df, "alg_cols") <- alg_cols
  message("[投票表] ", nrow(df), " 基因 × ", length(alg_cols), " 算法（",
          paste(alg_cols, collapse = "/"), "）；",
          sum(df$n_votes == length(alg_cols)), " 个基因为全票（交集）")
  df
}

## 3.2 单基因 AUC 表：拆 CI、校验、输出 10 标志物 canonical 顺序
load_gene_auc <- function(path) {
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  gcol  <- pick_col(df, c("gene", "symbol", "feature"), label = "单基因表 gene 列")
  names(df)[names(df) == gcol] <- "gene"
  auc_col <- pick_col(df, c("GSE57218_AUC", "gse57218_auc", "ext_auc"), label = "外部 AUC 列")
  ci_col  <- pick_col(df, c("GSE57218_CI", "gse57218_ci", "CI", "ci"), label = "CI 列")
  ## CI 形如 "0.153-1.000" → 拆成 ci_low / ci_up 两列
  parts <- strsplit(as.character(df[[ci_col]]), "-", fixed = TRUE)
  bad <- lengths(parts) != 2L
  if (any(bad)) stop("CI 列拆分失败（期望 'lower-upper' 格式）：",
                     paste(df[[ci_col]][bad], collapse = ", "), call. = FALSE)
  df$ci_low <- as.numeric(vapply(parts, `[`, character(1), 1L))
  df$ci_up  <- as.numeric(vapply(parts, `[`, character(1), 2L))
  ## 一致性校验：点估计应落在 CI 内（容忍 1e-3 浮点误差）
  outside <- df[[auc_col]] < df$ci_low - 1e-3 | df[[auc_col]] > df$ci_up + 1e-3
  if (any(outside, na.rm = TRUE)) {
    warning("[校验] 以下基因 GSE57218_AUC 落在 CI 之外，请复核：",
            paste(df$gene[outside], collapse = ", "))
  }
  message("[单基因表] ", nrow(df), " 基因；CI 已拆分为 ci_low/ci_up")
  df
}

## 3.3 七分类器稳健性表：剔除 NA 行（xgbTree），支持每分类器多种子聚合
load_robust <- function(path) {
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  ccol <- pick_col(df, c("classifier", "model", "method", "algo"), label = "分类器列")
  vcol <- pick_col(df, c("cv_AUC", "cv_auc", "internal_cv", "cv"), label = "CV AUC 列")
  ecol <- pick_col(df, c("GSE57218_AUC", "gse57218_auc", "ext_auc", "external_auc"),
                   label = "GSE57218 AUC 列")
  out <- data.frame(
    classifier = as.character(df[[ccol]]),
    cv_AUC     = suppressWarnings(as.numeric(df[[vcol]])),
    ext_AUC    = suppressWarnings(as.numeric(df[[ecol]])),
    stringsAsFactors = FALSE
  )
  n0 <- nrow(out)
  out <- out[stats::complete.cases(out), ]        # 剔除 xgbTree 等 NA 行
  if (nrow(out) < n0) message("[稳健性] 剔除 ", n0 - nrow(out), " 行含 NA 的记录（如 xgbTree）")
  ## 多种子时聚合为 mean±sd（单种子时 sd=NA → 面板不画误差线）
  agg <- do.call(rbind, lapply(split(out, out$classifier), function(z) {
    data.frame(
      classifier = z$classifier[1],
      cv_mean  = mean(z$cv_AUC),  cv_sd  = if (nrow(z) > 1) stats::sd(z$cv_AUC)  else NA_real_,
      ext_mean = mean(z$ext_AUC), ext_sd = if (nrow(z) > 1) stats::sd(z$ext_AUC) else NA_real_,
      n_seed   = nrow(z), stringsAsFactors = FALSE
    )
  }))
  rownames(agg) <- NULL
  glm_ext <- agg$ext_mean[agg$classifier %in% c("glm", "logistic", "lr")]
  if (length(glm_ext) && abs(glm_ext[1] - EXPECT$glm_ext_auc) > 0.01) {
    warning("[校验] logistic 外部 AUC=", round(glm_ext[1], 3),
            "，与论文值 ", EXPECT$glm_ext_auc, " 不一致，请复核输入表")
  }
  message("[稳健性] ", nrow(agg), " 个分类器；外部 AUC 范围 ",
          sprintf("%.3f–%.3f", min(agg$ext_mean), max(agg$ext_mean)))
  agg
}

## 3.4 SHAP：rds 对象（shapviz/kernelshap/list 探测）或长表重建
as_shapviz <- function(obj, src = "") {
  ## (1) 已是 shapviz 对象
  if (inherits(obj, "shapviz")) return(obj)
  ## (2) kernelshap 对象 → shapviz 原生转换
  if (inherits(obj, "kernelshap")) return(shapviz::shapviz(obj))
  ## (3) list：先找包裹的 shapviz/kernelshap 元素，再探测 S/X 槽位
  if (is.list(obj)) {
    for (el in obj) if (inherits(el, "shapviz"))   return(el)
    for (el in obj) if (inherits(el, "kernelshap")) return(shapviz::shapviz(el))
    nm <- names(obj); if (is.null(nm)) nm <- ""
    s_hit <- grep("^(S|shap_?values?|phi)$", nm, ignore.case = TRUE, value = TRUE)[1]
    x_hit <- grep("^(X|feature_?values?)$", nm, ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(s_hit)) {
      S <- as.matrix(obj[[s_hit]])
      X <- if (!is.na(x_hit)) as.data.frame(obj[[x_hit]]) else NULL
      return(shapviz::shapviz(S, X = X))
    }
  }
  stop("SHAP rds 结构无法识别：", src, "\n",
       "已尝试：shapviz 对象 / kernelshap 对象 / list 内嵌 shapviz|kernelshap / ",
       "list 槽位 S|shap_values|phi (+X|feature_values)。\n",
       "请改用 shapviz::shapviz() 构造后保存，或提供 shap_long.csv 长表。", call. = FALSE)
}

load_shap <- function() {
  ## 主路：rds 对象
  if (!is.na(CFG$shap_rds) && file.exists(CFG$shap_rds)) {
    shp <- as_shapviz(readRDS(CFG$shap_rds), CFG$shap_rds)
    message("[SHAP] rds 读取成功：", nrow(shp$S), " 样本 × ", ncol(shp$S), " 特征",
            if (is.null(shp$X)) "（无特征值矩阵 X，蜜蜂图将退化为单色）" else "")
    return(shp)
  }
  ## 备路：长表（gene + shap [+ feature_value]）→ 宽矩阵 → shapviz
  csv <- find_optional(c("results/moduleD_v2/shap_long.csv",
                         "results/moduleD/shap_long.csv", "moduleI_data/shap_long.csv",
                         "shap_long.csv"))
  if (!is.na(csv)) {
    lg <- read.csv(csv, check.names = FALSE, stringsAsFactors = FALSE)
    gc <- pick_col(lg, c("gene", "feature", "variable"), label = "SHAP 长表 gene 列")
    sc <- pick_col(lg, c("shap", "shap_value", "phi", "value"), label = "SHAP 长表 shap 列")
    vc <- pick_col(lg, c("feature_value", "expr", "expression", "x"),
                   required = FALSE, label = "SHAP 长表 feature_value 列")
    spl <- split(lg, lg[[gc]])
    len <- vapply(spl, nrow, integer(1))   # 注意：lengths() 对 data.frame 返回列数
    if (length(unique(len)) != 1L) {
      stop("SHAP 长表各基因行数不一致（", paste(unique(len), collapse = "/"),
           "），无法还原为 样本×基因 矩阵；请提供 shapviz 对象 rds。", call. = FALSE)
    }
    S <- vapply(spl, function(z) as.numeric(z[[sc]]), numeric(len[1]))
    X <- if (!is.null(vc)) as.data.frame(vapply(spl, function(z) as.numeric(z[[vc]]),
                                                numeric(len[1]))) else NULL
    message("[SHAP] 长表重建 shapviz：", nrow(S), " 样本 × ", ncol(S), " 特征")
    return(shapviz::shapviz(as.matrix(S), X = X))
  }
  stop("[CONFIG] SHAP 输入缺失：既未找到 shap_object.rds，也未找到 shap_long.csv。\n",
       "期望：shapviz/kernelshap rds（或含 S+X 槽位的 list），",
       "或长表（列 gene, shap, feature_value[可选]）。\n",
       "候选路径：results/moduleD/、moduleI_data/ 或脚本同目录。", call. = FALSE)
}

## 3.5 列线图拟合对象：rds 主路（含 datadist 失败自动降级）+ lrm 重建兜底
find_rms_fit <- function(obj) {
  ## 两轮探测：优先 rms/lrm，其次 glm（glm 亦可喂给 rms::nomogram 的兼容路径）
  if (inherits(obj, c("rms", "lrm", "ols"))) return(obj)
  if (is.list(obj)) for (el in obj) if (inherits(el, c("rms", "lrm", "ols"))) return(el)
  if (inherits(obj, c("glm", "lm"))) return(obj)
  if (is.list(obj)) for (el in obj) if (inherits(el, c("glm", "lm"))) return(el)
  NULL
}

load_nomogram <- function(genes10) {
  nomo_fun <- function(fit) rms::nomogram(
    fit, fun = stats::plogis, fun.at = c(0.01, 0.1, 0.5, 0.9, 0.99),
    funlabel = "Prob of OA"            # 范文同款标签；via8⑨：加 0.01 沿 LP 轴散开
  )
  ## —— 主路：rds 拟合对象 ——
  if (!is.na(CFG$nomo_rds) && file.exists(CFG$nomo_rds)) {
    fit <- find_rms_fit(readRDS(CFG$nomo_rds))
    if (!is.null(fit)) {
      nom <- tryCatch(nomo_fun(fit), error = function(e) {
        message("[列线图] rds 拟合对象直接建模失败（多为缺 datadist）：",
                conditionMessage(e), "\n         → 转入表达矩阵重建兜底")
        NULL
      })
      if (!is.null(nom)) {
        message("[列线图] 使用 rds 拟合对象：", class(fit)[1])
        return(nom)
      }
    } else {
      message("[列线图] rds 内未找到 lrm/rms/glm 拟合对象 → 转入表达矩阵重建兜底")
    }
  } else {
    message("[列线图] 未找到 nomogram_fit.rds → 转入表达矩阵重建兜底")
  }
  ## —— 兜底：训练表达矩阵（z 分数）+ 分组 → rms::lrm 重建 ——
  f_expr  <- find_optional(CFG$train_expr)
  f_group <- find_optional(CFG$train_group)
  if (is.na(f_expr) || is.na(f_group)) {
    stop("[CONFIG] 无法构建列线图：nomogram_fit.rds 不可用，且兜底文件缺失。\n",
         "请提供以下任一组合：\n",
         "  (i) nomogram_fit.rds：rms::lrm 拟合对象（建模时保留 datadist）；或\n",
         "  (ii) train_expr_z.csv（行=样本，列=10 标志物基因，z 分数；可选首列 sample）",
         " + train_group.csv（列 sample[可选], label[OA/Control 或 0/1]）。", call. = FALSE)
  }
  expr <- read.csv(f_expr, check.names = FALSE, stringsAsFactors = FALSE)
  grp  <- read.csv(f_group, check.names = FALSE, stringsAsFactors = FALSE)
  lcol <- pick_col(grp, c("label", "group", "status", "y", "outcome"), label = "分组列")
  y <- norm_label(grp[[lcol]], label = "train_group")
  genes <- intersect(genes10, names(expr))
  if (length(genes) < 3L) {
    stop("train_expr_z.csv 与 10 标志物交集 <3 列（命中：",
         paste(genes, collapse = ", "), "）。实际列：",
         paste(names(expr), collapse = ", "), call. = FALSE)
  }
  if (length(genes) < length(genes10)) {
    warning("[列线图] 表达矩阵仅命中 ", length(genes), "/", length(genes10),
            " 个标志物，列线图将用命中基因重建")
  }
  df <- data.frame(y = y, expr[, genes, drop = FALSE], check.names = FALSE)
  ## via5 关键修复：datadist 必须先于 lrm 创建、且必须显式入全局环境。
  ## 机理（已核 rms 源码 rms.s::Design）：lrm() 拟合时刻内部即按
  ## options(datadist) 的名字、沿词法作用域链（rms 命名空间 → .GlobalEnv）
  ## 查找 datadist——函数局部变量（含 <<- 前的局部态）永远不可见；且
  ## options(datadist) 全会话持久，上一轮残留 "dd" 而全局无 dd 时，
  ## lrm() 会在拟合时刻直接报 "dataset dd not found for options(datadist=)"，
  ## 拟合之后再赋值根本来不及（via4 失败原因）。先建后拟合，
  ## limits 随 fit$Design 一并嵌入，nomogram() 的 Getlim 亦可从全局命中。
  assign("dd", rms::datadist(df), envir = .GlobalEnv)
  options(datadist = "dd")
  fit <- rms::lrm(y ~ ., data = df, x = TRUE, y = TRUE)
  message("[列线图] lrm 重建完成：", nrow(df), " 样本 × ", length(genes), " 基因")
  nomo_fun(fit)
}

## ---------------------------------------------------------------------------
## 4. 面板函数（每面板一个函数；图内文字全英文；统计量直接入图——范文同款）
## ---------------------------------------------------------------------------

## (a) 四算法特征选择投票条形图：x=gene（封板 10 标志物，按票数降序），y=入选算法数
##     口径与 05b 封板规则一致（LASSO∩SVM-RFE 交集为主、RF/XGBoost 前10 投票为辅），
##     故展示对象为 genes10（单基因 AUC 表 canonical 顺序）而非全表票数 top10——
##     后者会把非交集的 3 票基因（如 ZMPSTE24）排进前 10、挤出 2 票的封板基因 CBR3。
panel_a <- function(vote, genes10) {
  alg_cols <- attr(vote, "alg_cols")
  top <- vote[vote$gene %in% genes10, , drop = FALSE]
  if (nrow(top) != length(genes10))
    stop("[投票表] 封板标志物在投票表中仅命中 ", nrow(top), "/", length(genes10),
         "；缺失：", paste(setdiff(genes10, top$gene), collapse = ", "),
         "。请核对投票表与单基因 AUC 表的基因集。")
  if (any(top$n_votes < 2L))
    warning("[校验] 以下封板标志物票数 <2，请复核：",
            paste(top$gene[top$n_votes < 2L], collapse = ", "))
  top <- top[order(-top$n_votes, top$gene), ]
  top$gene <- factor(top$gene, levels = top$gene)
  ggplot(top, aes(gene, n_votes)) +
    geom_col(fill = CR_COL$case, width = 0.72) +
    geom_text(aes(label = n_votes), vjust = -0.4, size = 2.8) +   # 票数入图
    scale_y_continuous(limits = c(0, length(alg_cols) + 0.4),
                       breaks = seq_len(length(alg_cols))) +
    labs(x = NULL, y = "Number of algorithms") +
    theme_cr(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))
}

## (b) 四算法交集 Venn：VennDiagram 原生 pastel 四圆（粉/蓝/绿/紫，范文 Fig2E 同款）
panel_b <- function(vote) {
  alg_cols <- attr(vote, "alg_cols")
  sets <- lapply(alg_cols, function(cl) vote$gene[vote[[cl]] == 1L])
  names(sets) <- alg_cols
  pastel <- c("#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6")  # pastel 五色
  invisible(capture.output(                      # 屏蔽 VennDiagram 控制台日志
    v <- VennDiagram::venn.diagram(
      sets, filename = NULL, disable.logging = TRUE,
      fill = pastel[seq_along(sets)], alpha = 0.55,
      cex = 0.9, cat.cex = 0.85, cat.fontface = "bold",  # via7⑥：集合名再缩
      fontfamily = "sans", cat.fontfamily = "sans",      # via11⑫：统一无衬线
      cat.dist = c(0.28, 0.28, 0.14, 0.14),              # via7⑥：推离椭圆线
      margin = 0.17                                      # 防集合名边缘裁切
    )
  ))
  patchwork::wrap_elements(panel = grid::grobTree(v))
}

## (c) mean|SHAP| 全局重要性条形：降序、单色、数值入图
panel_c <- function(shp) {
  imp <- sort(colMeans(abs(as.matrix(shp$S))), decreasing = TRUE)
  df <- data.frame(gene = factor(names(imp), levels = rev(names(imp))),
                   val  = as.numeric(imp))
  ggplot(df, aes(val, gene)) +
    geom_col(fill = CR_COL$case, width = 0.72) +
    geom_text(aes(label = sprintf("%.2f", val)), hjust = -0.12, size = 2.8) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +   # via6②：柱端数值防裁
    labs(x = "mean(|SHAP value|)", y = NULL) +
    theme_cr(base_size = 8)
}

## (d) SHAP 蜜蜂图：shapviz 原生 sv_importance(kind="beeswarm")，
##     特征值低→高 = 蓝→红渐变；无特征值矩阵时退化为单色抖动点并提示
panel_d <- function(shp) {
  if (!is.null(shp$X)) {
    p <- shapviz::sv_importance(shp, kind = "beeswarm", max_display = 10L) +
      scale_color_gradient(low = "#2166AC", high = "#B2182B",   # SHAP 惯例 蓝→红
                           name = "Feature value") +
      theme_cr(base_size = 8) +
      theme(legend.key.height = grid::unit(3, "mm"))
    return(p)
  }
  message("[面板 d] shapviz 对象无特征值矩阵 X → 单色抖动点兜底")
  S <- as.data.frame(as.matrix(shp$S))
  lg <- utils::stack(S)
  lg$feature <- factor(lg$ind, levels = names(sort(colMeans(abs(as.matrix(shp$S))),
                                                  decreasing = TRUE)))
  ggplot(lg, aes(values, feature)) +
    geom_point(position = position_jitter(height = 0.18, seed = 20260915),
               size = 0.9, alpha = 0.6, color = CR_COL$case) +
    labs(x = "SHAP value", y = NULL) +
    theme_cr(base_size = 8)
}

## (f) 列线图：rms::nomogram + base plot 原生输出（范文 Fig3B 像素级同款）
##     via7⑦：升整宽底排 174×76mm——14 条轴在 58mm 窄格必然互叠
##     via10⑪：范文像素级对照后微调——lmgp 0.15（数字上提、归属自己的尺）；
##     xfrac 0.16 + cex.var 0.55（标签区收窄贴近范文）；cex.axis 0.50 与
##     label.every=2 沿用 via9（锚点图实测：刻度数字≈4.8pt、变量名≈5.8pt）
panel_e <- function(nom) {
  patchwork::wrap_elements(full = ~ {
    op <- par(mai = c(0.15, 0.15, 0.15, 0.15))
    plot(nom, cex.var = 0.55, cex.axis = 0.50,
         xfrac = 0.16, lmgp = 0.15, label.every = 2)   # via10⑪
    par(op)
  })
}

## (f) 七分类器稳健性：内部 CV 与 GSE57218 外部 AUC 双系列点图，数值入图；
##     多种子时加 mean±sd 误差线；图内附注内部 CV 选择泄漏声明
panel_f <- function(agg) {
  agg <- agg[order(-agg$ext_mean), ]
  agg$classifier <- factor(agg$classifier, levels = agg$classifier)
  long <- data.frame(
    classifier = rep(agg$classifier, 2),
    series = factor(rep(c("Internal CV", "GSE57218 (external)"), each = nrow(agg)),
                    levels = c("Internal CV", "GSE57218 (external)")),
    auc = c(agg$cv_mean, agg$ext_mean),
    sd  = c(agg$cv_sd,  agg$ext_sd)
  )
  pd <- position_dodge(width = 0.32)
  p <- ggplot(long, aes(classifier, auc, color = series, group = series)) +
    geom_point(position = pd, size = 2) +
    geom_text(aes(label = sprintf("%.3f", auc)), position = pd,
              vjust = -0.9, size = 2.3, show.legend = FALSE) +
    scale_color_manual(values = c("Internal CV" = CR_COL$ctrl,
                                  "GSE57218 (external)" = CR_COL$case),
                       name = NULL) +
    scale_y_continuous(limits = c(min(long$auc, na.rm = TRUE) - 0.10, 1.03),
                       breaks = seq(0.7, 1, 0.1)) +
    labs(x = NULL, y = "AUC", caption = EXPECT$cv_note) +   # via6④：附注移出面板
    theme_cr(base_size = 8) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
          plot.caption = element_text(size = rel(0.72), color = "grey30",
                                      hjust = 0, lineheight = 1.15,
                                      margin = margin(t = 3)))
  ## 多种子（sd 非全 NA）时补误差线
  if (any(!is.na(long$sd))) {
    p <- p + geom_errorbar(aes(ymin = auc - sd, ymax = auc + sd),
                           position = pd, width = 0.18, linewidth = 0.4,
                           na.rm = TRUE)
  }
  p
}

## ---------------------------------------------------------------------------
## 5. 主流程
## ---------------------------------------------------------------------------
main <- function() {
  message("========== Fig 3 出图开始 ==========")
  vote     <- load_votes(CFG$vote_csv)
  gene_auc <- load_gene_auc(CFG$gene_auc_csv)
  robust   <- load_robust(CFG$robust_csv)
  shp      <- load_shap()
  genes10  <- gene_auc$gene                       # 10 标志物 canonical 顺序
  if (length(genes10) != EXPECT$n_markers) {
    warning("[校验] 单基因表基因数=", length(genes10),
            "，与预期 ", EXPECT$n_markers, " 不一致")
  }
  ## 交叉校验：SHAP 特征与 10 标志物应对齐
  miss <- setdiff(genes10, colnames(shp$S))
  if (length(miss)) warning("[校验] SHAP 矩阵缺少标志物：", paste(miss, collapse = ", "))

  ## via5 化整为零：6 面板各自独立出小图（同 18_fig2 via4 策略），拼接后续统一处理
  .pdf_null <- grDevices::pdf(NULL)
  panel_results <- character(0)
  save_panel <- function(make_p, tag, w, h) {
    fn <- file.path(OUT_DIR, paste0("Fig3_panel", tag))
    tryCatch({
      p <- make_p()
      save_cr(p, fn, width_mm = w, height_mm = h)
      panel_results <<- c(panel_results, paste0("✅ panel ", tag))
    }, error = function(e) {
      panel_results <<- c(panel_results,
                          paste0("❌ panel ", tag, "：", conditionMessage(e)))
    })
  }
  ## via7⑦⑧：三排版式（ABC / DE / F 整宽）——列线图升 F 整宽 174x76mm；
  ## 稳健性升 E（85x70mm）；D 同步加宽 85x70mm；A/B/C 58x70mm
  save_panel(function() panel_a(vote, genes10), "A_四算法投票", 58, 64)
  save_panel(function() panel_b(vote),          "B_四算法Venn", 58, 64)
  save_panel(function() panel_c(shp),           "C_SHAP全局重要性", 58, 64)
  save_panel(function() panel_d(shp),           "D_SHAP蜜蜂图", 85, 66)
  save_panel(function() panel_f(robust),        "E_七分类器稳健性", 85, 66)
  save_panel(function() panel_e(load_nomogram(genes10)), "F_列线图", 174, 88)
  invisible(grDevices::dev.off())
  cat("\n========== Fig 3 逐面板出图清单 ==========\n")
  cat(paste(panel_results, collapse = "\n"), "\n")
  message("========== Fig 3 出图完成 ==========")
  invisible(panel_results)
}

main()
