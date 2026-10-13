## =============================================================================
## 23_重建预测_外部校准_DCA_20260916.R
## 功能：按模块 I（12_moduleI_nomogram_...）逐行同口径重建训练层 logistic 模型，
##       对 GSE57218 重算预测概率，然后完成：
##   Part B：外部校准斜率/截距（替换正文误植的内部值 0.394）
##   Part C：Fig 4c DCA 正净获益阈值概率区间（替换正文占位话术）
## 自检：重建预测的外部 AUC 须等于冻结值 0.818（容差 0.005），不符即硬停。
## 运行：放在项目根目录（含 processed/ results/ 的那层），
##       source("23_重建预测_外部校准_DCA_20260916.R")
## 包依赖：pROC（必需）、dcurves（Part C，缺则跳过）、rms（可选，不用）
## =============================================================================

set.seed(20260913)   ## 与模块 I SEED_MAIN 一致（glm 本身确定无关种子，仅为纪律）

## ---- 0. 根目录定位 -----------------------------------------------------------
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
setwd(.d); message("[setwd] 项目根目录：", .d)

## ---- CONFIG（与 12_moduleI 完全同口径，勿改）---------------------------------
CONFIG <- list(
  IN_DIR     = "processed",
  TRAIN_SETS = c("GSE114007", "GSE55235", "GSE55457", "GSE12021", "GSE51588"),
  TRAIN_LABEL = "OA",
  VAL_SET    = "GSE57218",
  VAL_CASE_CANDIDATES = c("OA_affected", "OA"),
  BIOMARKERS = c("NELL1", "STEAP1", "GADD45A", "FBLN5", "HEMK1",
                 "GDE1", "RNF14", "BNIP3", "EGFR", "CBR3"),   ## 模块 I 日志封板
  AUC_FROZEN = 0.818, AUC_TOL = 0.005,
  OUT_DIR    = "复核输出_20260916"
)
dir.create(CFG <- CONFIG$OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ---- 装配函数（逐行抄自 12_moduleI，纪律：勿改）------------------------------
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
    message("  ", gse, ": donor_id 去重 → ", ncol(expr), " 样本")
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
         paste(avail, collapse = ","), "）中")
  hit[1]
}

## ---- 训练层装配 + 硬自检（117 = 65/52）---------------------------------------
tr <- lapply(CONFIG$TRAIN_SETS, assemble_cohort, case_label = CONFIG$TRAIN_LABEL)
feat <- Reduce(intersect, c(list(CONFIG$BIOMARKERS), lapply(tr, function(cc) colnames(cc$x))))
if (length(feat) != 10) stop("训练层可检标志物 != 10（实测 ", length(feat), "）")
x_tr <- do.call(rbind, lapply(tr, function(cc) cc$x[, feat, drop = FALSE]))
y_tr <- factor(unlist(lapply(tr, function(cc) as.character(cc$y))),
               levels = c("control", "case"))
message("训练层：", nrow(x_tr), " 样本（case=", sum(y_tr == "case"),
        " / control=", sum(y_tr == "control"), "）")
if (nrow(x_tr) != 117 || sum(y_tr == "case") != 65)
  stop("训练层样本量偏离封板（117=65/52）——实得 ", nrow(x_tr), "=",
       sum(y_tr == "case"), "/", sum(y_tr == "control"), "，请回查")

## ---- 验证层装配 + 硬自检（33/7）----------------------------------------------
val_case <- resolve_case(CONFIG$VAL_SET, CONFIG$VAL_CASE_CANDIDATES)
va <- assemble_cohort(CONFIG$VAL_SET, case_label = val_case)
if (!all(feat %in% colnames(va$x))) stop("验证层缺标志物")
if (sum(va$y == "case") != 33 || sum(va$y == "control") != 7)
  warning("GSE57218 构成偏离封板（33/7）——实得 ",
          sum(va$y == "case"), "/", sum(va$y == "control"))
message("验证层 GSE57218：", nrow(va$x), " 样本（33/7 口径，OA_preserved 未纳入）")

## ---- 重建 glm + 预测 + AUC 自检 ----------------------------------------------
if (!requireNamespace("pROC", quietly = TRUE)) stop("请先 install.packages('pROC')")
fit <- glm(y ~ ., data = data.frame(y = y_tr, x_tr[, feat, drop = FALSE]),
           family = binomial())
co  <- coef(fit)
message("glm 收敛，系数 ", length(co), " 个（含截距）")
write.csv(data.frame(term = names(co), coef = as.numeric(co)),
          file.path(CFG, "PartB0_重建模型系数.csv"), row.names = FALSE)

prob <- as.numeric(predict(fit, newdata = data.frame(va$x[, feat, drop = FALSE]),
                           type = "response"))
y_va <- ifelse(va$y == "case", 1L, 0L)
auc_va <- as.numeric(pROC::auc(pROC::roc(y_va, prob, quiet = TRUE,
                                         levels = c(0, 1), direction = "<")))
message("重建预测的外部 AUC = ", round(auc_va, 3), "（冻结值 ", CONFIG$AUC_FROZEN, "）")
if (abs(auc_va - CONFIG$AUC_FROZEN) > CONFIG$AUC_TOL)
  stop("重建 AUC 与冻结值偏差 >", CONFIG$AUC_TOL,
       "——装配口径与模块 I 不一致，请勿使用后续输出，先回查！")
message("✓ AUC 自检通过：重建口径与模块 I 一致")

pred <- data.frame(sample = rownames(va$x), prob = prob, label = y_va)
write.csv(pred, file.path(CFG, "PartB0_GSE57218预测概率.csv"), row.names = FALSE)

## =============================================================================
## Part B：外部校准斜率/截距
## =============================================================================
lp <- qlogis(pmin(pmax(prob, 1e-6), 1 - 1e-6))
slope_ext     <- unname(coef(glm(y_va ~ lp, family = binomial()))[2])
intercept_ext <- unname(coef(glm(y_va ~ 1, offset = lp, family = binomial()))[1])
cat("\n===== Part B：外部校准（GSE57218, n=", length(y_va), "）=====",
    "\n校准斜率 slope     =", round(slope_ext, 3),
    "\n校准截距 intercept =", round(intercept_ext, 3),
    "\n→ 正文「外部校准显示校准斜率 0.394」应替换为以上实算值；",
    "\n→ 20_fig4 脚本 CONFIG 的 slope_ext 同步替换；",
    "\n（旧值 0.394 实为训练层内部 bootstrap 校正值，见 moduleI_nomogram_validate.txt）\n")
write.csv(data.frame(metric = c("slope", "intercept"),
                     value  = c(slope_ext, intercept_ext)),
          file.path(CFG, "PartB_外部校准斜率截距.csv"), row.names = FALSE)

## =============================================================================
## Part C：DCA 正净获益阈值概率区间
## =============================================================================
if (requireNamespace("dcurves", quietly = TRUE)) {
  dca_o <- dcurves::dca(label ~ prob, data = pred)
  nd    <- as.data.frame(dca_o$dca)
  v_mod <- nd[nd$variable == "prob", ]
  v_all <- nd[grepl("all",  nd$variable, ignore.case = TRUE), ]
  v_non <- nd[grepl("none", nd$variable, ignore.case = TRUE), ]
  m <- merge(v_mod[, c("threshold", "net_benefit")],
             v_all[, c("threshold", "net_benefit")], by = "threshold",
             suffixes = c("_model", "_all"))
  m$nb_none <- if (nrow(v_non) > 0)
    stats::approx(v_non$threshold, v_non$net_benefit, xout = m$threshold, rule = 2)$y else 0
  win <- m[m$nb_model > m$nb_all & m$nb_model > m$nb_none, ]
  if (nrow(win) > 0) {
    step <- stats::median(diff(m$threshold))
    br   <- c(0, which(diff(win$threshold) > 1.5 * step), nrow(win))
    cat("\n===== Part C：DCA 正净获益阈值区间（外部层 GSE57218, n=", nrow(pred), "）=====\n")
    for (k in seq_len(length(br) - 1)) {
      seg <- win[(br[k] + 1):br[k + 1], ]
      cat(sprintf("  区间 %d：%.3f – %.3f\n", k, min(seg$threshold), max(seg$threshold)))
    }
    cat("→ 正文「在常用阈值概率范围内呈正净获益」替换为以上区间，",
        "并注明计算于外部验证层。\n")
  } else {
    cat("\n===== Part C：警告 =====",
        "\n模型在任何阈值下均未同时优于 treat-all 与 treat-none（患病率 82.5% 队列常见）。",
        "\n正文需据实改写，不可写「呈正净获益」。\n")
  }
  write.csv(m, file.path(CFG, "PartC_DCA净获益逐阈值明细.csv"), row.names = FALSE)
} else {
  message("\n[Part C] 需要 dcurves 包：install.packages('dcurves') 后重跑本脚本。")
}

cat("\n全部输出见：", normalizePath(CFG), "\n")
