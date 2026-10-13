# =============================================================================
# 43_重建_GSE168505预测概率_补齐Fig4A灰线_20261010_via3.R
# 功能：重建 GSE168505（n=7）预测概率，补齐 Fig 4A 探索性灰虚线
# 背景：主图审查发现 Fig 4A 只有 GSE57218 单曲线，图注/正文均引用 GSE168505
#   （AUC = 0.750, n = 7, exploratory）；根因 = pred_GSE168505.csv 从未落盘。
# via3 关键更正：当年 PartB0 并非用锁定系数表计算，而是按 23_重建预测 脚本
#   现场重建 glm（5 训练队列各自 z 分后堆叠 117=65/52 -> glm -> 预测）。
#   因此本脚本【逐行照搬 23 号的 assemble_cohort / resolve_case / CONFIG】，
#   保证管线与当年一字不差；锁定系数表不再使用。
# 自验证（三道硬闸门，任一不过即终止）：
#   1) 训练层构成硬检：117 样本 = 65 case / 52 control（23 号同款封板值）；
#   2) GSE57218 重建预测 vs 归档 PartB0_GSE57218预测概率.csv：相关 > 0.999
#      且实算 AUC = 0.818 ± 0.005（冻结值）；
#   3) GSE168505：n = 7 且实算 AUC = 0.750 ± 0.01（封存值）。
# 输出：results/moduleE/pred_GSE168505.csv（sample, prob, label；1=OA, 0=Control）
# 作者：刘晟宇    日期：2026-10-10
# =============================================================================

message("== 脚本版本：43_重建_GSE168505预测概率_20261010_via3 ==")
if (!requireNamespace("pROC", quietly = TRUE)) stop("请先 install.packages('pROC')")

set.seed(20260913)   ## 与 23 号一致（glm 确定无关种子，仅为纪律）

## ---- 0. 项目根定位（同 23 号）-----------------------------------------------
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
setwd(.d); message("[43] 项目根目录：", .d)

## ---- CONFIG（训练层部分与 23 号完全同口径，勿改）-----------------------------
CONFIG <- list(
  IN_DIR      = "processed",
  TRAIN_SETS  = c("GSE114007", "GSE55235", "GSE55457", "GSE12021", "GSE51588"),
  TRAIN_LABEL = "OA",
  VAL_SET     = "GSE57218",
  VAL_CASE_CANDIDATES = c("OA_affected", "OA"),
  NEW_SET     = "GSE168505",
  NEW_CASE_CANDIDATES = c("OA", "OA_affected", "case"),   ## 按命中顺序取第一个
  BIOMARKERS  = c("NELL1", "STEAP1", "GADD45A", "FBLN5", "HEMK1",
                  "GDE1", "RNF14", "BNIP3", "EGFR", "CBR3"),
  AUC_FROZEN_57218  = 0.818, AUC_TOL_57218  = 0.005,
  AUC_FROZEN_168505 = 0.750, AUC_TOL_168505 = 0.01,
  ARC_FILE    = "复核输出_20260916/PartB0_GSE57218预测概率.csv",
  OUT_FILE    = "results/moduleE/pred_GSE168505.csv"     ## 脚本 20 via9 第一候选路径
)

## ---- 装配函数（逐行抄自 23 号，纪律：勿改）-----------------------------------
assemble_cohort <- function(gse, case_label, ctrl_label = "control") {
  expr  <- readRDS(file.path(CONFIG$IN_DIR, gse, paste0(gse, "_expr_gene.rds")))
  pheno <- read.csv(file.path(CONFIG$IN_DIR, gse, paste0(gse, "_pheno.csv")),
                    stringsAsFactors = FALSE)
  pheno <- pheno[match(colnames(expr), pheno$sample), ]
  if (anyNA(pheno$sample) || anyNA(pheno$condition))
    stop(gse, ": pheno 缺少 expr 中的样本或 condition 缺失——请回查模块 A 输出")
  keep <- pheno$condition %in% c(ctrl_label, case_label)
  if (sum(keep) == 0)
    stop(gse, ": 无 case/control 样本——case 标签 '", case_label,
         "' 不在 condition 取值（", paste(unique(pheno$condition), collapse = ","), "）中")
  expr <- expr[, keep, drop = FALSE]; pheno <- pheno[keep, ]
  if ("donor_id" %in% names(pheno) && anyDuplicated(pheno$donor_id)) {
    if (anyNA(pheno$donor_id)) stop(gse, ": donor_id 含 NA，请人工核对")
    keep_idx <- !duplicated(pheno$donor_id)
    expr <- expr[, keep_idx, drop = FALSE]; pheno <- pheno[keep_idx, ]
    message("  ", gse, ": donor_id 去重 -> ", ncol(expr), " 样本")
  }
  x <- scale(t(expr))
  if (anyNA(x)) {
    bad <- colnames(x)[apply(x, 2, anyNA)]
    message("  ", gse, ": 剔除 ", length(bad), " 个零方差基因（z-score NaN）")
    x <- x[, !colnames(x) %in% bad, drop = FALSE]
  }
  y <- factor(ifelse(pheno$condition == case_label, "case", "control"),
              levels = c("control", "case"))
  if (sum(y == "case") == 0 || sum(y == "control") == 0)
    stop(gse, ": 病例/对照有一侧为 0")
  list(x = x, y = y)
}

resolve_case <- function(gse, candidates) {
  pheno <- read.csv(file.path(CONFIG$IN_DIR, gse, paste0(gse, "_pheno.csv")),
                    stringsAsFactors = FALSE)
  avail <- unique(pheno$condition)
  hit <- intersect(candidates, avail)
  if (length(hit) == 0)
    stop(gse, ": case 候选均不在 condition 实测取值（",
         paste(avail, collapse = ","), "）中——请把取值发回补映射")
  message("  ", gse, ": case 标签 = '", hit[1], "'（实测取值：",
          paste(avail, collapse = ","), "）")
  hit[1]
}

## ---- 闸门 1：训练层装配 + 构成硬检（117 = 65/52）------------------------------
tr <- lapply(CONFIG$TRAIN_SETS, assemble_cohort, case_label = CONFIG$TRAIN_LABEL)
feat <- Reduce(intersect, c(list(CONFIG$BIOMARKERS), lapply(tr, function(cc) colnames(cc$x))))
if (length(feat) != 10) stop("训练层可检标志物 != 10（实测 ", length(feat), "）")
x_tr <- do.call(rbind, lapply(tr, function(cc) cc$x[, feat, drop = FALSE]))
y_tr <- factor(unlist(lapply(tr, function(cc) as.character(cc$y))),
               levels = c("control", "case"))
message("[43] 训练层：", nrow(x_tr), " 样本（case=", sum(y_tr == "case"),
        " / control=", sum(y_tr == "control"), "）")
if (nrow(x_tr) != 117 || sum(y_tr == "case") != 65)
  stop("[闸门1] 训练层样本量偏离封板（117=65/52）——实得 ", nrow(x_tr), "=",
       sum(y_tr == "case"), "/", sum(y_tr == "control"), "，请回查")

## ---- 重建 glm（同 23 号）------------------------------------------------------
fit <- glm(y ~ ., data = data.frame(y = y_tr, x_tr[, feat, drop = FALSE]),
           family = binomial())
message("[43] glm 收敛，系数 ", length(coef(fit)), " 个（含截距）")

## ---- 闸门 2：GSE57218 重建预测，对冻结值 + 归档文件双校验 ----------------------
val_case <- resolve_case(CONFIG$VAL_SET, CONFIG$VAL_CASE_CANDIDATES)
va <- assemble_cohort(CONFIG$VAL_SET, case_label = val_case)
if (!all(feat %in% colnames(va$x))) stop("验证层缺标志物")
prob57 <- as.numeric(predict(fit, newdata = data.frame(va$x[, feat, drop = FALSE]),
                             type = "response"))
names(prob57) <- rownames(va$x)
y57 <- ifelse(va$y == "case", 1L, 0L)
auc57 <- as.numeric(pROC::auc(pROC::roc(y57, prob57, quiet = TRUE,
                                        levels = c(0, 1), direction = "<")))
message("[闸门2] GSE57218 重建 AUC = ", round(auc57, 3),
        "（冻结值 ", CONFIG$AUC_FROZEN_57218, "）")
if (abs(auc57 - CONFIG$AUC_FROZEN_57218) > CONFIG$AUC_TOL_57218)
  stop("[闸门2] 重建 AUC 与冻结值偏差 >", CONFIG$AUC_TOL_57218, "——管线口径不一致，终止")

if (file.exists(CONFIG$ARC_FILE)) {
  arc <- read.csv(CONFIG$ARC_FILE, check.names = FALSE, stringsAsFactors = FALSE)
  common <- intersect(arc$sample, names(prob57))
  r <- suppressWarnings(cor(arc$prob[match(common, arc$sample)], prob57[common]))
  message("[闸门2] 与归档 PartB0 逐样本相关 = ", round(r, 4),
          "（", length(common), "/", nrow(arc), " 样本匹配）")
  if (length(common) < nrow(arc) * 0.9 || !is.finite(r) || r <= 0.999)
    stop("[闸门2] 与归档 PartB0 不一致（相关 ", round(r, 4), "）——管线口径不一致，终止")
} else {
  warning("[闸门2] 未找到归档 ", CONFIG$ARC_FILE, "，仅做了 AUC 冻结值校验")
}
message("[43] 管线验证通过：与当年 23 号口径一致")

## ---- GSE168505 装配 + 预测 ---------------------------------------------------
new_case <- resolve_case(CONFIG$NEW_SET, CONFIG$NEW_CASE_CANDIDATES)
nw <- assemble_cohort(CONFIG$NEW_SET, case_label = new_case)
if (!all(feat %in% colnames(nw$x)))
  stop("GSE168505 缺标志物：", paste(setdiff(feat, colnames(nw$x)), collapse = ", "))
prob16 <- as.numeric(predict(fit, newdata = data.frame(nw$x[, feat, drop = FALSE]),
                             type = "response"))
y16 <- ifelse(nw$y == "case", 1L, 0L)
out <- data.frame(sample = rownames(nw$x), prob = prob16, label = y16)
message("[43] GSE168505 n = ", nrow(out),
        "（case ", sum(y16), " / control ", sum(1 - y16), "）")

## ---- 闸门 3：n 与 AUC 封存值硬验收 --------------------------------------------
if (nrow(out) != 7) stop("[闸门3] GSE168505 n != 7（实得 ", nrow(out), "）——请回查 pheno 口径")
auc16 <- as.numeric(pROC::auc(pROC::roc(y16, prob16, quiet = TRUE,
                                        levels = c(0, 1), direction = "<")))
message("[闸门3] GSE168505 实算 AUC = ", round(auc16, 3),
        "（封存值 ", CONFIG$AUC_FROZEN_168505, "）")
if (abs(auc16 - CONFIG$AUC_FROZEN_168505) > CONFIG$AUC_TOL_168505)
  stop("[闸门3] AUC 与封存值 ", CONFIG$AUC_FROZEN_168505, " 不符（差 ",
       round(abs(auc16 - CONFIG$AUC_FROZEN_168505), 4), "）——终止，请回查")

## ---- 落盘 ---------------------------------------------------------------------
dir.create(dirname(CONFIG$OUT_FILE), showWarnings = FALSE, recursive = TRUE)
write.csv(out, CONFIG$OUT_FILE, row.names = FALSE, quote = TRUE)
message("[43] 已写出：", normalizePath(CONFIG$OUT_FILE))
message("[43] 完成。下一步：重跑 20_fig4 via9 -> Fig 4A 将自动加入 GSE168505 灰虚线")
