# ============================================================================
# 03_moduleB_shared_genes.R
# 模块 B｜OA × 肌少症共享基因识别（三层交叉框架，计划书 v1.4a §3.3.2）
# ----------------------------------------------------------------------------
# 三层设计：
#   第 1 层 差异表达：逐数据集 limma（模块 A 输出均为 log 尺度基因矩阵；
#           RNA-seq 已经 VST；若需 DESeq2 原始计数路径见 DESeq2_OPTIONAL）
#           —— 设计公式读取模块 A 写入 pheno 表的 design_flag：
#              paired_by_donor / paired_cartilage → block 因子
#              family_relatedness → duplicateCorrelation
#              POOLED_SAMPLES → 直接排除
#   第 2 层 秩和整合：RobustRankAggreg 分别在 OA 侧与肌肉侧整合；
#           共享判定采用双轨框架（2026-08-22 适配，备忘录 §八-24）：
#           OA 侧 adj<RRA_CUT ∩ 肌肉侧名义 RRA_score<MUSCLE_NOMINAL_CUT，
#           按方向一致性拆为「同向轨」与「镜像轨」两个核心集
#   第 3 层 共表达网络：WGCNA 在各侧锚点数据集上构建（默认 OA=GSE114007，
#           肌肉=GSE111016），作为核心集基因的置信度注释（不再做硬门槛）
# ----------------------------------------------------------------------------
# 前置：已运行 02_moduleA_preprocess.R（processed/<GSE>/ 下有 *_expr_gene.rds
#       与 *_pheno.csv，且 condition 列已人工核对——module A 对不确定样本
#       标 "other"，进入本脚本前必须处理掉）
# 运行：source("03_moduleB_shared_genes.R")；结果输出 results/moduleB/
# ============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(RobustRankAggreg)
  library(WGCNA)
  library(ggplot2)
  library(data.table)
})
allowWGCNAThreads()

# ----------------------------- CONFIG --------------------------------------
CONFIG <- list(
  IN_DIR   = "processed",
  OUT_DIR  = "results/moduleB",
  LOG_FILE = sprintf("moduleB_log_%s.txt", Sys.Date()),
  DE_FDR_CUT  = 0.05,     # 单层 DE 显著性阈值（各数据集内部）
  DE_LFC_CUT  = 0.585,    # |log2FC| ≥ 1.5 倍（rlog/VST 尺度近似；敏感性分析用 0/1.0）
  RRA_CUT     = 0.05,     # RRA 校正后阈值（OA 侧）
  MUSCLE_NOMINAL_CUT = 0.05,  # 肌肉侧名义 RRA 分数阈值（双轨框架；弱信号适配，备忘录 §八-24）
  AGING_PROXY = c("GSE1428", "GSE25941"),  # 衰老代理队列（敏感性分析用）
  WGCNA_DATASETS = c(OA = "GSE114007", muscle = "GSE111016"),  # 网络构建锚点集（各侧训练队列中样本最大/表型最确证者）
  WGCNA_POWER_OA = NA,    # NA = 自动 pickSoftThreshold
  WGCNA_POWER_MUSCLE = NA,
  MODULE_KME_CUT = 0.7,   # 模块成员度阈值
  GS_COR_CUT  = 0.3,      # 基因-表型相关阈值
  EXCLUDE_FROM_DE = c("GSE169077")   # pooled，仅方向参考
)

## 数据集分组（与 01_GEO_download_check.R MANIFEST 的 role 列严格对齐；备忘录 §八-14）
## 仅训练/发现侧进入模块 B；验证侧（OA: GSE169077/GSE168505/GSE57218/GSE206848/GSE48556，
## 肌肉: GSE136344/GSE117525/GSE144304）留给模块 D/E，避免信息泄漏。
## ★ GSE48556 是 OA 外周血队列（GARP），不是肌少症队列——严禁归入肌肉侧（2026-08-22 更正）
OA_SETS     <- c("GSE114007","GSE55235","GSE55457","GSE12021","GSE51588")
MUSCLE_SETS <- c("GSE111016","GSE111010","GSE111006","GSE1428","GSE25941")
DISEASE_LABEL <- c(OA = "OA", muscle = "sarcopenia")
## 逐数据集对比口径覆盖（衰老代理队列：case=aging, ctrl=young）
CONTRAST_OVERRIDE <- list(
  GSE1428  = c(case = "aging", ctrl = "young"),
  GSE25941 = c(case = "aging", ctrl = "young")
)
get_contrast <- function(g, side) {
  if (g %in% names(CONTRAST_OVERRIDE))
    return(c(case = unname(CONTRAST_OVERRIDE[[g]]["case"]), ctrl = unname(CONTRAST_OVERRIDE[[g]]["ctrl"])))
  c(case = unname(DISEASE_LABEL[[side]]), ctrl = "control")
}

# ----------------------------- 工具函数 -------------------------------------
log_msg <- function(...) {
  txt <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = ""))
  cat(txt, "\n"); cat(txt, "\n", file = CONFIG$LOG_FILE, append = TRUE)
}

load_processed <- function(gse) {
  expr_file  <- file.path(CONFIG$IN_DIR, gse, paste0(gse, "_expr_gene.rds"))
  pheno_file <- file.path(CONFIG$IN_DIR, gse, paste0(gse, "_pheno.csv"))
  if (!file.exists(expr_file) || !file.exists(pheno_file))
    stop(gse, ": 缺少模块 A 输出，请先运行 02_moduleA_preprocess.R")
  expr  <- readRDS(expr_file)
  pheno <- read.csv(pheno_file, stringsAsFactors = FALSE)
  list(expr = expr, pheno = pheno)
}

# ----------------------------- 第 1 层：逐数据集 limma ----------------------
# 依据 design_flag 构建设计：
#   - paired_*        : 样本名中含供体/患者编号的需在 pheno 手工补 block_id 列；
#                       若无 block_id 则退回非配对并记 WARNING
#   - family_relatedness: pheno 需有 family_id 列 → duplicateCorrelation
run_limma_one <- function(gse, disease_label, ctrl_label = "control") {
  obj <- load_processed(gse)
  expr <- obj$expr; pheno <- obj$pheno
  pheno <- pheno[match(colnames(expr), pheno$sample), ]
  if (any(is.na(pheno$sample)))
    stop(gse, ": 样本名对齐失败（表达矩阵列名 ↔ pheno$sample）——",
         "回查模块 A 输出的 ", gse, "_pheno.csv")
  if (any(pheno$condition == "other"))
    log_msg(gse, ": ⚠ ", sum(pheno$condition == "other"), " 个样本 condition=other，",
            "将被剔除——请先在模块 A 输出中人工核对补全")
  keep <- pheno$condition %in% c(ctrl_label, disease_label,
                                 if (gse == "GSE57218") c("OA_affected") else NULL)
  # GSE57218 主比较：OA_affected vs control（受累 vs 健康；配对分析作为补充）
  if (gse == "GSE57218") pheno$condition[pheno$condition == "OA_affected"] <- "OA"
  expr <- expr[, keep, drop = FALSE]; pheno <- pheno[keep, ]
  grp <- factor(pheno$condition, levels = c(ctrl_label, disease_label))
  if (nlevels(grp) < 2 || any(table(grp) < 3)) {
    log_msg(gse, ": ⚠ 组样本不足（", paste(names(table(grp)), table(grp), collapse = " vs "),
            "），跳过该数据集"); return(NULL)
  }
  flag <- pheno$design_flag[1]
  corfit <- NULL; block_var <- NULL
  if (grepl("family_relatedness", flag) && "family_id" %in% colnames(pheno)) {
    design <- model.matrix(~ grp)
    block_var <- pheno$family_id
    corfit <- duplicateCorrelation(expr, design, block = block_var)
    log_msg(gse, ": duplicateCorrelation(family) consensus cor = ", round(corfit$consensus, 3))
  } else if (grepl("donor_correlation", flag) && "donor_id" %in% colnames(pheno)) {
    design <- model.matrix(~ grp)
    block_var <- pheno$donor_id
    corfit <- duplicateCorrelation(expr, design, block = block_var)
    log_msg(gse, ": duplicateCorrelation(donor) consensus cor = ", round(corfit$consensus, 3))
  } else if (grepl("paired", flag) && "block_id" %in% colnames(pheno)) {
    design <- model.matrix(~ grp + factor(pheno$block_id))
    log_msg(gse, ": 配对设计（block_id）")
  } else {
    if (grepl("paired|family|donor_correlation", flag))
      log_msg(gse, ": ⚠ WARNING 设计标注为 ", flag, " 但 pheno 缺对应列，退回非配对设计——请人工补列后重跑")
    design <- model.matrix(~ grp)
  }
  fit <- lmFit(expr, design, block = block_var,
               correlation = if (!is.null(corfit)) corfit$consensus else NULL)
  fit <- eBayes(fit)
  coef_name <- paste0("grp", disease_label)
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  log_msg(gse, ": DE 完成（", disease_label, " vs ", ctrl_label, "），FDR<",
          CONFIG$DE_FDR_CUT, " 且 |logFC|≥", CONFIG$DE_LFC_CUT,
          " 的基因数 = ", sum(tt$adj.P.Val < CONFIG$DE_FDR_CUT & abs(tt$logFC) >= CONFIG$DE_LFC_CUT))
  tt[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
}

# ----------------------------- 第 2 层：RRA 整合 ----------------------------
run_rra <- function(de_list, side_name) {
  de_list <- de_list[!vapply(de_list, is.null, logical(1))]
  if (length(de_list) < 2) stop(side_name, " 侧有效数据集 < 2，无法 RRA")
  # 上调/下调分别按 P 值排序（带方向，双侧列表）
  up_lists   <- lapply(de_list, function(tt) tt$gene[order(tt$P.Value, -abs(tt$logFC))][tt$logFC[order(tt$P.Value, -abs(tt$logFC))] > 0])
  down_lists <- lapply(de_list, function(tt) tt$gene[order(tt$P.Value, -abs(tt$logFC))][tt$logFC[order(tt$P.Value, -abs(tt$logFC))] < 0])
  universe <- Reduce(union, lapply(de_list, `[[`, "gene"))
  up   <- aggregateRanks(glist = up_lists,   N = length(universe))
  down <- aggregateRanks(glist = down_lists, N = length(universe))
  colnames(up)   <- c("gene", "RRA_score_up")
  colnames(down) <- c("gene", "RRA_score_down")
  m <- merge(up, down, by = "gene")
  # 方向一致性：取较小分数方向为整合方向；用 BH 校正（aggregateRanks 分数为 p 值量级）
  m$direction <- ifelse(m$RRA_score_up < m$RRA_score_down, "up", "down")
  m$RRA_score <- pmin(m$RRA_score_up, m$RRA_score_down)
  m$RRA_adj   <- p.adjust(m$RRA_score, method = "BH")
  m <- m[order(m$RRA_adj), ]
  log_msg(side_name, " 侧 RRA：", sum(m$RRA_adj < CONFIG$RRA_CUT), " 个基因 adj<", CONFIG$RRA_CUT)
  m
}

# ----------------------------- 第 3 层：WGCNA -------------------------------
run_wgcna <- function(gse, disease_label, power_override) {
  obj <- load_processed(gse)
  expr <- obj$expr; pheno <- obj$pheno
  pheno <- pheno[match(colnames(expr), pheno$sample), ]
  if (any(is.na(pheno$sample)))
    stop(gse, ": 样本名对齐失败——回查模块 A 输出的 ", gse, "_pheno.csv")
  keep <- pheno$condition %in% c("control", disease_label)
  expr <- expr[, keep, drop = FALSE]; pheno <- pheno[keep, ]
  datExpr <- t(expr)
  # 方差过滤提速（保留 top 8000 高变异基因）
  v <- apply(datExpr, 2, var, na.rm = TRUE)
  datExpr <- datExpr[, order(v, decreasing = TRUE)[1:min(8000, ncol(datExpr))]]
  gsg <- goodSamplesGenes(datExpr, verbose = 0)
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  trait <- as.numeric(pheno$condition[match(rownames(datExpr), pheno$sample)] == disease_label)
  if (is.na(power_override)) {
    sft <- pickSoftThreshold(datExpr, powerVector = 1:20, verbose = 0,
                             networkType = "signed", corFnc = "bicor")
    power <- ifelse(is.na(sft$powerEstimate), 6, sft$powerEstimate)
  } else power <- power_override
  log_msg(gse, ": WGCNA soft power = ", power, "（", ncol(datExpr), " 基因 × ",
          nrow(datExpr), " 样本）")
  net <- blockwiseModules(datExpr, power = power, networkType = "signed",
                          TOMType = "signed", minModuleSize = 30,
                          mergeCutHeight = 0.25, numericLabels = TRUE,
                          pamRespectsDendro = FALSE, verbose = 0)
  MEs <- moduleEigengenes(datExpr, colors = labels2colors(net$colors))$eigengenes
  mod_trait_cor <- cor(MEs, trait, use = "p")
  mod_trait_p   <- corPvalueStudent(mod_trait_cor, nrow(datExpr))
  gene_module <- data.frame(gene = colnames(datExpr),
                            module = labels2colors(net$colors),
                            stringsAsFactors = FALSE)
  kME <- signedKME(datExpr, MEs)
  gene_module$kME_own <- vapply(seq_len(nrow(gene_module)), function(i)
    kME[i, paste0("kME", gene_module$module[i])], numeric(1))
  gene_module$GS <- abs(as.numeric(cor(datExpr, trait, use = "p")))
  list(gene_module = gene_module, mod_trait_cor = mod_trait_cor,
       mod_trait_p = mod_trait_p, power = power)
}

# ----------------------------- 主流程 ---------------------------------------
dir.create(CONFIG$OUT_DIR, recursive = TRUE, showWarnings = FALSE)
log_msg("模块 B 启动")

## ---- 第 1 层：逐数据集 DE ----
de_oa <- setNames(vector("list", length(OA_SETS)), OA_SETS)
for (g in setdiff(OA_SETS, CONFIG$EXCLUDE_FROM_DE)) {
  ct <- get_contrast(g, "OA")
  de_oa[[g]] <- tryCatch(run_limma_one(g, ct[["case"]], ct[["ctrl"]]),
                         error = function(e) { log_msg(g, ": ❌ ", conditionMessage(e)); NULL })
}
de_mu <- setNames(vector("list", length(MUSCLE_SETS)), MUSCLE_SETS)
for (g in MUSCLE_SETS) {
  ct <- get_contrast(g, "muscle")
  de_mu[[g]] <- tryCatch(run_limma_one(g, ct[["case"]], ct[["ctrl"]]),
                         error = function(e) { log_msg(g, ": ❌ ", conditionMessage(e)); NULL })
}
saveRDS(list(OA = de_oa, muscle = de_mu),
        file.path(CONFIG$OUT_DIR, "layer1_DE_all_datasets.rds"))
invisible(lapply(names(de_oa), function(g) if (!is.null(de_oa[[g]]))
  write.csv(de_oa[[g]], file.path(CONFIG$OUT_DIR, paste0("layer1_DE_", g, ".csv")), row.names = FALSE)))
invisible(lapply(names(de_mu), function(g) if (!is.null(de_mu[[g]]))
  write.csv(de_mu[[g]], file.path(CONFIG$OUT_DIR, paste0("layer1_DE_", g, ".csv")), row.names = FALSE)))

## ---- 第 2 层：RRA ----
rra_oa <- run_rra(de_oa, "OA")
rra_mu <- run_rra(de_mu, "muscle")
write.csv(rra_oa, file.path(CONFIG$OUT_DIR, "layer2_RRA_OA.csv"), row.names = FALSE)
write.csv(rra_mu, file.path(CONFIG$OUT_DIR, "layer2_RRA_muscle.csv"), row.names = FALSE)

## ---- 第 2 层共享判定：双轨框架（2026-08-22 适配，备忘录 §八-24）----
## 适配依据：肌肉侧基因级信号弱且跨族群异质（GSE111016 全队列仅 1 个 FDR 显著基因；
## PMID 31862890 确立"基因级弱、通路级汇聚"为该领域已知特性）；
## 且 OA 显著基因与肌肉侧呈显著镜像（99 候选中 72 反向=73%，单侧二项 p=3.5e-6，
## GSE111016 单队列独立复现）。故 OA 侧维持 adj<RRA_CUT，肌肉侧改用名义
## RRA_score<MUSCLE_NOMINAL_CUT，按方向一致性拆为同向轨/镜像轨双轨输出。
shared_all <- merge(rra_oa[rra_oa$RRA_adj < CONFIG$RRA_CUT, c("gene", "direction", "RRA_adj")],
                    rra_mu[rra_mu$RRA_score < CONFIG$MUSCLE_NOMINAL_CUT, c("gene", "direction", "RRA_score")],
                    by = "gene", suffixes = c("_OA", "_muscle"))
if (nrow(shared_all) == 0) stop("双轨共享基因为 0——请回查两侧 RRA 表与阈值设置")
shared_all$concordant <- shared_all$direction_OA == shared_all$direction_muscle
core_concordant <- shared_all[shared_all$concordant, ]
core_mirror     <- shared_all[!shared_all$concordant, ]
bt <- binom.test(nrow(core_mirror), nrow(shared_all), p = 0.5, alternative = "greater")
log_msg("第 2 层双轨共享基因：同向轨 ", nrow(core_concordant), " 个 / 镜像轨 ", nrow(core_mirror),
        " 个（镜像占比 ", round(100 * nrow(core_mirror) / nrow(shared_all), 1),
        "%，单侧二项检验 p = ", signif(bt$p.value, 2), "）")
write.csv(shared_all, file.path(CONFIG$OUT_DIR, "layer2_shared_DE_genes.csv"), row.names = FALSE)

## ---- 敏感性分析：肌肉侧剔除衰老代理队列，仅 3 个 sarcopenia 队列重跑 RRA ----
de_mu_sarc <- de_mu[setdiff(names(de_mu), CONFIG$AGING_PROXY)]
rra_mu_sarc <- tryCatch(run_rra(de_mu_sarc, "muscle-sarcopenia-only"),
                        error = function(e) { log_msg("sarcopenia-only RRA 失败：", conditionMessage(e)); NULL })
if (!is.null(rra_mu_sarc)) {
  write.csv(rra_mu_sarc, file.path(CONFIG$OUT_DIR, "layer2_RRA_muscle_sarcopenia_only.csv"), row.names = FALSE)
  sens <- merge(rra_oa[rra_oa$RRA_adj < CONFIG$RRA_CUT, c("gene", "direction")],
                rra_mu_sarc[rra_mu_sarc$RRA_score < CONFIG$MUSCLE_NOMINAL_CUT, c("gene", "direction")],
                by = "gene", suffixes = c("_OA", "_muscle"))
  if (nrow(sens) > 0) {
    sens$concordant <- sens$direction_OA == sens$direction_muscle
    log_msg("敏感性（sarcopenia-only）：候选 ", nrow(sens), " 个，镜像占比 ",
            round(100 * mean(!sens$concordant), 1), "%；与主分析镜像轨重合 ",
            length(intersect(sens$gene[!sens$concordant], core_mirror$gene)), "/", nrow(core_mirror))
    write.csv(sens, file.path(CONFIG$OUT_DIR, "layer2_sensitivity_sarcopenia_only.csv"), row.names = FALSE)
  }
}

## ---- 第 3 层：WGCNA ----
wg <- list()
for (side in names(CONFIG$WGCNA_DATASETS)) {
  gse <- CONFIG$WGCNA_DATASETS[[side]]
  wg[[side]] <- tryCatch(
    run_wgcna(gse, DISEASE_LABEL[[side]],
              if (side == "OA") CONFIG$WGCNA_POWER_OA else CONFIG$WGCNA_POWER_MUSCLE),
    error = function(e) { log_msg(gse, ": ❌ WGCNA ", conditionMessage(e)); NULL })
}
saveRDS(wg, file.path(CONFIG$OUT_DIR, "layer3_WGCNA_networks.rds"))

# 关键模块 = |cor(ME, trait)| 最大且 P<0.05 的模块；模块基因再按 kME/GS 过滤
key_module_genes <- list()
for (side in names(wg)) {
  if (is.null(wg[[side]])) next
  cor_v <- wg[[side]]$mod_trait_cor[, 1]; p_v <- wg[[side]]$mod_trait_p[, 1]
  key_mods <- names(which(abs(cor_v) > 0.3 & p_v < 0.05))
  key_mods <- sub("^ME", "", key_mods)
  gm <- wg[[side]]$gene_module
  sel <- gm$module %in% key_mods & gm$kME_own >= CONFIG$MODULE_KME_CUT & gm$GS >= CONFIG$GS_COR_CUT
  key_module_genes[[side]] <- gm$gene[sel]
  write.csv(gm, file.path(CONFIG$OUT_DIR, paste0("layer3_WGCNA_modules_", side, ".csv")), row.names = FALSE)
  log_msg(side, " 侧关键模块：", paste(key_mods, collapse = ","), " → 高置信模块基因 ",
          length(key_module_genes[[side]]), " 个")
}

## ---- 三层整合：WGCNA 高置信模块成员身份作为核心集注释（双轨；不再做硬门槛）----
wgcna_union <- unique(unlist(key_module_genes))
core_concordant$in_WGCNA_key <- core_concordant$gene %in% wgcna_union
core_mirror$in_WGCNA_key     <- core_mirror$gene %in% wgcna_union
write.csv(core_concordant, file.path(CONFIG$OUT_DIR, "CORE_concordant_genes.csv"), row.names = FALSE)
write.csv(core_mirror,     file.path(CONFIG$OUT_DIR, "CORE_mirror_genes.csv"),     row.names = FALSE)
log_msg("✅ 三层整合完成（WGCNA 注释）：同向轨 ", nrow(core_concordant), " 个（WGCNA 高置信 ",
        sum(core_concordant$in_WGCNA_key), "）→ CORE_concordant_genes.csv；镜像轨 ",
        nrow(core_mirror), " 个（WGCNA 高置信 ", sum(core_mirror$in_WGCNA_key),
        "）→ CORE_mirror_genes.csv")

cat("\n模块 B 完成（双轨框架）。核心产物：\n",
    "  ", CONFIG$OUT_DIR, "/layer1_DE_*.csv          —— 各数据集 DE 表\n",
    "  ", CONFIG$OUT_DIR, "/layer2_shared_DE_genes.csv   —— 双轨共享 DE 全表（含方向一致性）\n",
    "  ", CONFIG$OUT_DIR, "/layer2_sensitivity_sarcopenia_only.csv —— 剔除衰老代理的敏感性分析\n",
    "  ", CONFIG$OUT_DIR, "/layer3_WGCNA_modules_*.csv   —— 模块归属\n",
    "  ", CONFIG$OUT_DIR, "/CORE_concordant_genes.csv    —— 同向轨核心集（下游模块 C-I 输入）\n",
    "  ", CONFIG$OUT_DIR, "/CORE_mirror_genes.csv        —— 镜像轨核心集（下游模块 C-I 输入）\n",
    "检查日志 ", CONFIG$LOG_FILE, " 中的 WARNING（配对/家系列缺失、样本剔除）后再继续。\n", sep = "")
