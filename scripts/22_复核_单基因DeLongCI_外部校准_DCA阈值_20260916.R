## =============================================================================
## 22_复核_单基因DeLongCI_外部校准_DCA阈值_20260916.R  （v2 自给自足版）
## 功能：v7 终稿「待确认事项」三项复核
##   Part A ：10 标志物单基因外部 AUC 的 95%CI 用 pROC DeLong 法重算
##   Part B0：预测值获取链（三级回退）：
##             ① 现成预测值文件 → ② 系数表 → ③ 训练集矩阵现场重建 glm 冻结系数
##   Part B ：GSE57218 外部校准斜率/截距（替换正文中误植的内部值 0.394）
##   Part C ：Fig 4c DCA 正净获益阈值概率区间提取
## 运行：把本文件放在项目根目录（含 processed/ results/ data/ 的那层），然后：
##       source("22_复核_单基因DeLongCI_外部校准_DCA阈值_20260916.R")
## 包依赖：pROC（必需）、dcurves（Part C 用，缺则自动跳过）
## =============================================================================

set.seed(20260916)

## ---- 0. 自动定位项目根目录并 setwd（脚本须放在根目录下）---------------------
.script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- grep("--file=", a, value = TRUE)
  if (length(f) >= 1) return(dirname(normalizePath(sub("--file=", "", f[1]))))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable() &&
      nzchar(rstudioapi::getActiveDocumentContext()$path))
    return(dirname(normalizePath(rstudioapi::getActiveDocumentContext()$path)))
  getwd()
}
.d <- .script_dir()
for (i in 1:8) {
  if (any(dir.exists(file.path(.d, c("processed", "results", "data"))))) break
  if (dirname(.d) == .d) break
  .d <- dirname(.d)
}
setwd(.d)
message("[setwd] 项目根目录：", .d)

## ---- CONFIG（已按 2026-09-16 本地文件清单填好；训练集路径若不同请改）-------
CFG <- list(
  ## Part A 输入（已确认存在）
  expr_file  = "processed/GSE57218/GSE57218_expr_gene.rds",
  group_file = "processed/GSE57218/GSE57218_pheno.csv",
  genes_file = "results/moduleI/moduleI_single_gene_AUC.csv",
  ## Part B0 三级回退输入
  pred_file  = "results/moduleI/GSE57218_predictions.csv",   ## ① 若有现成预测值
  coef_file  = "results/moduleI/moduleI_coefficients.csv",   ## ② 若有系数表(gene,coef[,intercept])
  train_expr = "moduleI_data/train_expr_z.csv",              ## ③ 训练集 z 矩阵(样本×基因)
  train_grp  = "moduleI_data/train_group.csv",               ## ③ 训练集分组(sample,group)
  dca_file   = "",                                            ## 留空则复用 Part B0 得到的预测值
  out_dir    = "复核输出_20260916"
)
dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 通用小工具 --------------------------------------------------------------
read_any <- function(path) {
  if (!file.exists(path)) stop("文件不存在：", path, call. = FALSE)
  if (grepl("\\.rds$", path, ignore.case = TRUE)) return(readRDS(path))
  if (grepl("\\.tsv$|\\.txt$", path, ignore.case = TRUE))
    return(utils::read.delim(path, check.names = FALSE))
  utils::read.csv(path, check.names = FALSE)
}

pick_col <- function(df, candidates, label) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) >= 1) return(hit[1])
  stop("在「", label, "」中找不到候选列：", paste(candidates, collapse = "/"),
       "；实际列为：", paste(colnames(df), collapse = ", "), call. = FALSE)
}

norm_label <- function(x) {  ## 归一化 1=OA 病例, 0=对照（正则模糊匹配）
  y <- tolower(trimws(as.character(x)))
  is_case <- grepl("^(oa|case|1|disease|oa_patient)$|osteoarthritis", y)
  is_ctrl <- grepl("^(control|ctrl|normal|healthy|0|nc|hc|donor)$|control|healthy|normal", y)
  out <- ifelse(is_case & !is_ctrl, 1L, ifelse(is_ctrl & !is_case, 0L, NA))
  out[grepl("preserved", y)] <- 2L   ## RAAK 保留侧软骨：标记为 2，主验证剔出
  if (any(is.na(out)))
    stop("分组标签无法归一化，原始取值：", paste(unique(x), collapse = " | "),
         "\n请把上述取值发回，补充映射规则。", call. = FALSE)
  out
}

as_gene_x_sample <- function(expr, label_txt) {
  ## 防御性整理表达矩阵为 基因×样本；兼容 data.frame/matrix/list(rds)
  if (is.list(expr) && !is.data.frame(expr)) {
    cand <- intersect(c("expr", "exprs", "mat", "data", "expr_gene"), names(expr))
    if (length(cand) >= 1) expr <- expr[[cand[1]]] else
      stop("rds list 中找不到表达矩阵槽位（", label_txt, "）")
  }
  expr <- as.data.frame(expr, check.names = FALSE)
  if (is.null(rownames(expr)) || rownames(expr)[1] %in% as.character(seq_len(nrow(expr)))) {
    rownames(expr) <- make.unique(as.character(expr[[1]]))
    expr <- expr[, -1, drop = FALSE]
  }
  m <- as.matrix(expr)
  storage.mode(m) <- "numeric"
  m
}

## =============================================================================
## Part A：单基因外部 AUC 的 DeLong 95%CI 重算
## =============================================================================
if (!requireNamespace("pROC", quietly = TRUE)) stop("请先 install.packages('pROC')")

## -- A1. 10 个标志物及 OA 病例侧方向 --
gi <- read_any(CFG$genes_file)
g_gene <- pick_col(gi, c("gene", "Gene", "symbol"), "单基因结果表")
dir_col <- intersect(c("direction_OA", "OA_direction", "oa_direction", "oa_dir"),
                     colnames(gi))
if (length(dir_col) >= 1) {
  ten_dir <- stats::setNames(tolower(as.character(gi[[dir_col[1]]])), gi[[g_gene]])
  message("[A1] 方向来源：", dir_col[1])
} else {
  ## 手工兜底（v7 结果节：GADD45A/BNIP3/EGFR 为 OA 下调，其余 7 个 OA 上调）
  ten_dir <- stats::setNames(rep("up", nrow(gi)), gi[[g_gene]])
  ten_dir[c("GADD45A", "BNIP3", "EGFR")] <- "down"
  message("[A1] 未找到方向列 → 手工兜底：GADD45A/BNIP3/EGFR=down，其余=up（请核对）")
}
## 注意：表中 case_direction 全为 up_in_case 是 pROC direction 修正伪影，不可用
biomarkers <- gi[[g_gene]]

## -- A2. 表达矩阵与分组 --
expr <- as_gene_x_sample(read_any(CFG$expr_file), "GSE57218 表达")
grp  <- read_any(CFG$group_file)
g_s  <- pick_col(grp, c("sample", "Sample", "gsm", "GSM", "geo_accession", "id"), "分组表")
g_g  <- pick_col(grp, c("group", "Group", "label", "phenotype", "status",
                        "condition", "disease", "diagnosis", "characteristics",
                        "condition_raw"), "分组表")

common <- intersect(colnames(expr), as.character(grp[[g_s]]))
if (length(common) < 10)
  stop("表达矩阵列名与分组表样本名交集 <10（", length(common), "）。",
       "\n矩阵列名示例：", paste(head(colnames(expr)), collapse = ", "),
       "\n分组样本示例：", paste(head(as.character(grp[[g_s]])), collapse = ", "), call. = FALSE)
grp  <- grp[match(common, as.character(grp[[g_s]])), ]
expr <- expr[, common, drop = FALSE]
label <- norm_label(grp[[g_g]])
if (any(label == 2L)) {  ## RAAK OA_preserved（保留侧）剔出主验证，与正文 7/33 口径一致
  keep <- label != 2L
  message("[A2] OA_preserved ", sum(!keep), " 例剔出主验证（模块 E 口径：仅供体内配对分析）")
  label <- label[keep]
  expr  <- expr[, keep, drop = FALSE]
}
message("[A2] 样本 ", length(label), "（病例 ", sum(label == 1), " / 对照 ", sum(label == 0), "）")
if (!(sum(label == 1) == 33 && sum(label == 0) == 7))
  warning("病例/对照数与正文口径（33/7）不符，请核对分组映射！")

## -- A3. 逐基因 DeLong CI --
resA <- lapply(biomarkers, function(g) {
  if (!g %in% rownames(expr)) { warning("基因 ", g, " 不在表达矩阵，跳过"); return(NULL) }
  d <- ten_dir[g]; if (is.na(d)) d <- "up"
  ## OA 上调基因：病例侧表达更高 → direction="<"；下调基因 → ">"
  roc_o <- pROC::roc(response = label, predictor = expr[g, ], quiet = TRUE,
                     direction = if (d == "down") ">" else "<", levels = c(0, 1))
  ci <- as.numeric(pROC::ci.auc(roc_o, method = "delong"))
  data.frame(gene = g, direction_OA = d, AUC = ci[2], CI_low = ci[1], CI_high = ci[3])
})
resA <- do.call(rbind, resA)
resA[, c("AUC", "CI_low", "CI_high")] <- round(resA[, c("AUC", "CI_low", "CI_high")], 3)
resA$CI_crosses_0.5 <- resA$CI_low < 0.5
write.csv(resA, file.path(CFG$out_dir, "PartA_单基因DeLongCI_重算.csv"), row.names = FALSE)
cat("\n===== Part A：单基因 DeLong 95%CI 重算 =====\n"); print(resA)
cat("→ 与 S11 旧值对照（旧：GADD45A/EGFR CI 0.153-1.000）；",
    "若下限明显 >0.5，正文保守口径与 S11 需联动修订。\n")

## =============================================================================
## Part B0：预测值获取链（① 预测值文件 → ② 系数表 → ③ 训练集重建 glm）
## =============================================================================
pred <- NULL   ## 目标：data.frame(prob, label)

## ① 现成预测值文件
if (nzchar(CFG$pred_file) && file.exists(CFG$pred_file)) {
  pr <- read_any(CFG$pred_file)
  p_l <- pick_col(pr, c("label", "group", "y", "obs"), "预测值表")
  p_p <- pick_col(pr, c("prob", "prob_OA", "pred", "p"), "预测值表")
  pred <- data.frame(prob = pr[[p_p]], label = norm_label(pr[[p_l]]))
  message("[B0] 预测值来源 ①：", CFG$pred_file)
}

## ② 系数表
if (is.null(pred) && nzchar(CFG$coef_file) && file.exists(CFG$coef_file)) {
  cf <- read_any(CFG$coef_file)
  c_g <- pick_col(cf, c("gene", "Gene", "term", "variable"), "系数表")
  c_c <- pick_col(cf, c("coef", "coefficient", "estimate", "beta"), "系数表")
  co  <- stats::setNames(cf[[c_c]], cf[[c_g]])
  int <- if ("(Intercept)" %in% names(co)) co["(Intercept)"] else
         if ("intercept" %in% names(co)) co["intercept"] else
         stop("系数表中找不到截距行 (Intercept)")
  gs  <- setdiff(names(co), c("(Intercept)", "intercept"))
  miss <- setdiff(gs, rownames(expr))
  if (length(miss) > 0) stop("系数表基因不在表达矩阵：", paste(miss, collapse = ", "))
  lp <- int + as.numeric(co[gs] %*% expr[gs, , drop = FALSE])
  pred <- data.frame(prob = plogis(lp), label = label)
  message("[B0] 预测值来源 ②：系数表（", length(gs), " 基因）")
}

## ③ 训练集现场重建 glm（冻结系数后应用于外部队列）
if (is.null(pred) && file.exists(CFG$train_expr) && file.exists(CFG$train_grp)) {
  te <- read_any(CFG$train_expr)
  tg <- read_any(CFG$train_grp)
  t_s <- pick_col(tg, c("sample", "Sample", "gsm", "GSM", "id"), "训练分组表")
  t_g <- pick_col(tg, c("group", "Group", "label", "phenotype", "status",
                        "condition", "condition_raw"), "训练分组表")
  ## 训练矩阵整理为 样本×基因
  te <- as.data.frame(te, check.names = FALSE)
  first_col <- as.character(te[[1]])
  biom_in_col <- any(biomarkers %in% colnames(te))
  if (!biom_in_col && any(biomarkers %in% first_col)) {
    rownames(te) <- make.unique(first_col); te <- te[, -1, drop = FALSE]
    te <- as.data.frame(t(te), check.names = FALSE)          ## 转置为样本×基因
  } else if (biom_in_col) {
    if (!all(as.character(tg[[t_s]]) %in% first_col) && any(tg[[t_s]] %in% first_col))
      rownames(te) <- first_col
  }
  use_g <- intersect(biomarkers, colnames(te))
  if (length(use_g) < length(biomarkers))
    warning("训练矩阵缺基因：", paste(setdiff(biomarkers, use_g), collapse = ", "))
  tr <- data.frame(y = norm_label(tg[[t_g]])[match(rownames(te), as.character(tg[[t_s]]))],
                   te[, use_g, drop = FALSE], check.names = FALSE)
  tr <- tr[!is.na(tr$y), ]
  fit <- stats::glm(y ~ ., data = tr, family = stats::binomial())
  co  <- stats::coef(fit)
  message("[B0] 预测值来源 ③：训练集重建 glm（", length(use_g), " 基因，",
          "训练 n=", nrow(tr), "）")
  ## 外部队列按自身均值/SD 行 z-score（逐队列独立 z 化，与正文口径一致；
  ## 若当时管线是用训练集均值/SD 换算，请改此处）
  ez <- t(scale(t(expr[use_g, , drop = FALSE])))
  lp <- co["(Intercept)"] + as.numeric(co[use_g] %*% t(ez))
  pred <- data.frame(prob = plogis(lp), label = label)
  ## 系数存盘备查（② 路径下次可直接用）
  write.csv(data.frame(gene = names(co), coef = as.numeric(co)),
            file.path(CFG$out_dir, "PartB0_重建模型系数.csv"), row.names = FALSE)
}

if (is.null(pred)) {
  message("\n[B0] 三条路径均未命中，Part B/C 跳过。",
          "\n请提供以下任一：① GSE57218 预测值文件；② 模型系数表；",
          "\n③ 训练集矩阵（修改 CFG$train_expr / CFG$train_grp 指向实际文件）。\n")
}

## =============================================================================
## Part B：外部校准斜率/截距（slope = glm(y~lp)；intercept = glm(y~offset(lp))）
## =============================================================================
if (!is.null(pred)) {
  pv <- pmin(pmax(pred$prob, 1e-6), 1 - 1e-6)
  lp <- qlogis(pv); y <- pred$label
  slope_ext     <- unname(stats::coef(stats::glm(y ~ lp, family = stats::binomial()))[2])
  intercept_ext <- unname(stats::coef(stats::glm(y ~ 1, offset = lp,
                                                 family = stats::binomial()))[1])
  cat("\n===== Part B：外部校准（GSE57218, n=", length(y), "）=====",
      "\n校准斜率 slope     =", round(slope_ext, 3),
      "\n校准截距 intercept =", round(intercept_ext, 3),
      "\n→ 正文「外部校准显示校准斜率 0.394」替换为以上实算值；",
      "\n→ 20_fig4 脚本 CONFIG 的 slope_ext 同步替换。\n")
  write.csv(data.frame(metric = c("slope", "intercept"),
                       value = c(slope_ext, intercept_ext)),
            file.path(CFG$out_dir, "PartB_外部校准斜率截距.csv"), row.names = FALSE)
  write.csv(pred, file.path(CFG$out_dir, "PartB0_GSE57218预测概率.csv"), row.names = FALSE)
}

## =============================================================================
## Part C：DCA 正净获益阈值概率区间提取
## =============================================================================
dca_src <- NULL
if (nzchar(CFG$dca_file) && file.exists(CFG$dca_file)) {
  dd <- read_any(CFG$dca_file)
  d_l <- pick_col(dd, c("label", "group", "y", "obs"), "DCA 数据")
  d_p <- pick_col(dd, c("prob", "prob_OA", "pred", "p"), "DCA 数据")
  dca_src <- data.frame(prob = dd[[d_p]], label = norm_label(dd[[d_l]]))
} else if (!is.null(pred)) {
  dca_src <- pred
}
if (!is.null(dca_src) && requireNamespace("dcurves", quietly = TRUE)) {
  dca_o <- dcurves::dca(label ~ prob, data = dca_src)
  nd    <- as.data.frame(dca_o$dca)
  nm    <- unique(nd$variable)
  v_mod <- nd[nd$variable == "prob", ]
  v_all <- nd[nd$variable %in% c("all", "treat_all", "all_treat"), ]
  v_non <- nd[nd$variable %in% c("none", "treat_none"), ]
  if (nrow(v_all) == 0) v_all <- nd[grepl("all", nd$variable, ignore.case = TRUE), ]
  m <- merge(v_mod[, c("threshold", "net_benefit")],
             v_all[, c("threshold", "net_benefit")], by = "threshold",
             suffixes = c("_model", "_all"))
  m$nb_none <- if (nrow(v_non) > 0)
    stats::approx(v_non$threshold, v_non$net_benefit, xout = m$threshold,
                  rule = 2)$y else 0
  win <- m[m$nb_model > m$nb_all & m$nb_model > m$nb_none, ]
  if (nrow(win) > 0) {
    ## 连续段切分
    br <- c(0, which(diff(win$threshold) > 1.5 * stats::median(diff(m$threshold))), nrow(win))
    cat("\n===== Part C：DCA 正净获益阈值区间 =====\n")
    for (k in seq_len(length(br) - 1)) {
      seg <- win[(br[k] + 1):br[k + 1], ]
      cat(sprintf("  区间 %d：%.3f – %.3f\n", k, min(seg$threshold), max(seg$threshold)))
    }
    cat("→ 正文「在常用阈值概率范围内呈正净获益」替换为实际区间，",
        "并注明计算于外部验证层（GSE57218, n=", nrow(dca_src), "）。\n", sep = "")
  } else {
    cat("\n===== Part C：警告 =====",
        "\n模型在任何阈值下均未同时优于 treat-all 与 treat-none（高患病率队列常见）。",
        "\n正文需据实改写，不可写「呈正净获益」。\n")
  }
  write.csv(m, file.path(CFG$out_dir, "PartC_DCA净获益逐阈值明细.csv"), row.names = FALSE)
} else if (!is.null(dca_src)) {
  message("\n[Part C] 需要 dcurves 包：install.packages('dcurves') 后重跑。")
}

cat("\n全部复核输出见：", normalizePath(CFG$out_dir), "\n")
