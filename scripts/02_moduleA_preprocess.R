# ============================================================================
# 02_moduleA_preprocess.R
# 模块 A｜数据获取与预处理（per-dataset：探针→基因映射、标准化、QC）
# 项目：OA × 肌少症共享分子特征整合生信研究（计划书 v1.4a）
# ----------------------------------------------------------------------------
# 设计原则（与计划书 §3.3.1 一致）：
#   1) 逐数据集独立预处理（不跨平台合并原始矩阵；跨数据集整合交给模块 B 的
#      RRA 秩和层 + 逐数据集 WGCNA，从源头规避跨平台批次）；
#   2) 已知批次的唯一例外：GSE114007 双平台计数合并（GPL11154+GPL18573），
#      VST 后用 sva::ComBat 以"平台"为已知批次校正（计划书 3.2.1 附注）；
#   3) 探针→基因：默认用 GPL 官方注释（getGEO 平台表），同一基因多探针
#      取行均值最大者（按均值降序后去重）；Ensembl ID（RNA-seq）
#      经 org.Hs.eg.db 转 symbol；
#   4) 每个数据集输出：基因级表达矩阵 RDS + 样本表 CSV + QC PDF
#      （标准化前后箱线/密度、PCA 按表型着色）；
#   5) 特殊数据集按 CONFIG 中的 flag 处理（pooled / 家族相关 / 配对设计 /
#      双平台），flag 会写入样本表供模块 B 读取。
# ----------------------------------------------------------------------------
# 前置：已运行 01_GEO_download_check.R（生成 data/GEO/<GSE>/<GSE>_eset.rds
#       或 _series_matrix.txt.gz；RNA-seq 原始计数需手动放入对应目录，
#       见 CONFIG$seq_count_file 说明）
# 运行：source("02_moduleA_preprocess.R") ；结果输出至 processed/<GSE>/
# ============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(limma)        # avereps / neqc / plotMDS
  library(DESeq2)      # vst（RNA-seq）
  library(sva)         # ComBat（GSE114007 平台批次）
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  # 注：Affymetrix 的 RMA 需原始 CEL；GEO series matrix 多已 RMA/MAS5，
  #     本脚本对 series matrix 做 log 判定 + quantile 兜底（见函数注释）。
  library(ggplot2)
  library(data.table)
})

# ----------------------------- CONFIG --------------------------------------
CONFIG <- list(
  DATA_DIR  = "data/GEO",
  OUT_DIR   = "processed",
  LOG_FILE  = sprintf("moduleA_log_%s.txt", Sys.Date()),
  # RNA-seq 原始计数文件名约定（需手动从 GEO supplementary 下载后放入）：
  #   data/GEO/<GSE>/<GSE>_raw_counts.txt(.gz)  —— 行=基因(Ensembl 或 symbol)，列=样本
  SEQ_COUNT_FILENAME = "raw_counts",
  MAKE_QC_PDF = TRUE
)

# 数据集清单 + 处理配置（与 01_GEO_download_check.R 的 MANIFEST 对齐；
# platform_type: affy / ilmn(甲基化级 illumina 表达芯片) / seq / agilent
# flag: 设计警示，写入样本表供下游读取）
DATASETS <- data.frame(
  gse = c("GSE114007","GSE55235","GSE55457","GSE12021","GSE51588",
          "GSE169077","GSE168505","GSE57218","GSE206848","GSE48556",
          "GSE1428","GSE25941","GSE111016","GSE111010","GSE111006",
          "GSE136344","GSE117525","GSE144304",
          "GSE152805","GSE169454","GSE167186"),
  tissue_group = c(rep("OA_tissue",9),"OA_blood",
                   rep("muscle",8), rep("singlecell",3)),
  platform_type = c("seq","affy","affy","affy","agilent",        # OA 训练
                    "affy","seq","ilmn","affy","ilmn",           # OA 验证（GSE48556=血液卖点）
                    "affy","affy","seq","seq","seq",             # 肌肉训练
                    "affy","affy","seq",                         # 肌肉验证
                    "seq","seq","seq"),                          # 单细胞（模块 A 只登记，不处理）
  # 注：agilent 走 ilmn 标准化分支（series matrix 通用 log 判定 + quantile 兜底）
  condition_col_hint = c(
    "title",                      # GSE114007 分组在 title（Normal_Cart_x/OA_x）；source_name 全是 cartilage 无分组信息
    "characteristics_ch1",        # GSE55235 分组在 disease state 列（另有 10 例 RA 自动剔除）
    "characteristics_ch1.2",      # GSE55457 分组在 clinical status 列（另有 13 例 RA 自动剔除）
    "characteristics_ch1.2",      # GSE12021 分组在 disease 列（仅用 GPL96 口径，见 PREFERRED_GPL）
    "characteristics_ch1",        # GSE51588 胫骨平台软骨，供体两部位结构
    "characteristics_ch1",        # GSE169077 池化样本，分组在 disease state 列
    "characteristics_ch1.1",      # GSE168505 分组在 disease state 列（osteoarthritis/non-arthritic）
    "characteristics_ch1.3",      # GSE57218 RAAK 分组在 disease state 列（OA/Preserved/Healthy）
    "characteristics_ch1",        # GSE206848（另有 2 例 RA 须剔除）
    "source_name_ch1",            # GSE48556 GARP: PBMC_OA / PBMC_Control（OA 血液验证）
    "description",                # GSE1428 无 characteristics 列，分组在 description（Young/Older ... years old）
    "characteristics_ch1.1",      # GSE25941 分组在 age: Old/Young（纯衰老代理，勿标肌少症）
    "characteristics_ch1.1",      # GSE111016 分组在 sarcopenia status 列（锚点 20/20）
    "characteristics_ch1.1",      # GSE111010 分组在 sarcopenia status 列（肌少症仅 9 例）
    "characteristics_ch1.1",      # GSE111006 分组在 sarcopenia status 列（肌少症仅 4 例）
    "title",                      # GSE136344 三组在 title（Muscle_young/Old/Old SX），disease state 列会把 young 混入 healthy
    "characteristics_ch1.9",      # GSE117525 基线分组在 ch1.9（Healthy older/Frail/Young baseline）；随访行该列为 NA → other 自动剔除
    "characteristics_ch1.3",      # GSE144304 分组在 ch1.3（group: Fit/Young/Frail；Fit=健康老年对照）
    NA, NA, NA),
  stringsAsFactors = FALSE
)

# 多平台芯片数据集的人工选定口径（load_dataset_eset 按此挑选 series matrix）
PREFERRED_GPL <- c(GSE12021 = "GPL96")   # GPL96 口径 9 正常/10 OA；GPL97 弃用

# 设计警示（写进样本表 design_flag 列；模块 B 读取后走对应设计公式）
DESIGN_FLAGS <- list(
  GSE114007 = "known_batch_platform: VST 后 ComBat(platform=GPL11154/GPL18573)",
  GSE169077 = "POOLED_SAMPLES: 合并样本无生物学重复，仅方向参考，禁止进 RRA",
  GSE48556  = "family_relatedness: GARP 同胞对家系相关（129 女/10 男 OA 血液队列）；GEO 元数据无 family_id，模块 E 若无法从原文补充表获取则退回非配对并在论文 Limitations 声明",
  GSE57218  = "paired_cartilage: 受累/保留配对（block=patient）；与 7 例健康比较为非配对",
  GSE51588  = "donor_correlation: 20 OA 供体×2 部位 + 5 对照供体×2 部位；模块 B 宜 duplicateCorrelation(donor)，需 pheno 补 donor_id 列",
  GSE55235  = "has_RA: 仅取正常+OA（另有 10 例 RA 自动标 other 剔除）",
  GSE55457  = "has_RA: 仅取正常+OA（另有 13 例 RA 自动标 other 剔除）",
  GSE206848 = "has_RA: 仅取正常+OA（另有 2 例 RA 自动标 other 剔除）",
  GSE12021  = "donor_overlap: 与 GSE55457 供体重叠，仅训练侧 RRA；口径 GPL96（9 正常/10 OA）",
  GSE168505 = "small_n: n=7（3 正常/4 OA），仅验证侧",
  GSE1428   = "aging_proxy: 全男性 young vs old（衰老+轻度肌少症代理），case=aging/ctrl=young",
  GSE25941  = "aging_proxy: 纯衰老代理表型（勿标肌少症），case=aging/ctrl=young",
  GSE111010 = "few_cases: 肌少症仅 9 例；与 GSE111016 同系列注意批次",
  GSE111006 = "few_cases_very: 肌少症仅 4 例，统计效力弱",
  GSE136344 = "three_group_all_male: 年轻/健康老年/老年+代谢综合征；主对比=老年+MetS vs 健康老年",
  GSE117525 = "has_followup: 总 259 含 6 个月运动随访样本，主分析仅取基线——pData 复核时剔除随访",
  GSE144304 = "has_young: 另含 26 例年轻人；主对比=(pre-)frail vs 健康老年",
  GSE111016 = "anchor: 确诊肌少症 20/20，肌肉侧表型锚点"
)

# ----------------------------- 工具函数 -------------------------------------
log_msg <- function(...) {
  txt <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = ""))
  cat(txt, "\n")
  cat(txt, "\n", file = CONFIG$LOG_FILE, append = TRUE)
}

# 找到数据集目录中已下载的表达对象（芯片数据集：单平台单 eSet）
load_dataset_eset <- function(gse) {
  d <- file.path(CONFIG$DATA_DIR, gse)
  rds <- file.path(d, paste0(gse, "_eset.rds"))
  if (file.exists(rds)) return(readRDS(rds))
  sm <- list.files(d, pattern = "_series_matrix\\.txt(\\.gz)?$", full.names = TRUE)
  if (length(sm) > 1) {
    pref <- if (gse %in% names(PREFERRED_GPL)) PREFERRED_GPL[[gse]] else NA_character_
    hit  <- if (!is.na(pref)) sm[grepl(pref, basename(sm), fixed = TRUE)] else character(0)
    if (length(hit) == 1) {
      log_msg(gse, ": 多平台数据集（", length(sm), " 个 series matrix），按 PREFERRED_GPL=",
              pref, " 选用 ", basename(hit))
      return(getGEO(filename = hit, getGPL = TRUE))
    }
    stop(gse, ": 存在 ", length(sm), " 个 series matrix（多平台）且 PREFERRED_GPL 未指定/未命中——",
         "请在本脚本 PREFERRED_GPL 中登记选定平台后重跑")
  }
  if (length(sm) > 0) return(getGEO(filename = sm[1], getGPL = TRUE))
  stop(gse, ": 未找到 eSet（请先运行 01_GEO_download_check.R）")
}

# 多平台数据集（GSE114007）的表型合并：读取全部 series matrix 的 pData 并 rbind
# （seq 分支只需要表型；各平台样本互斥，行合并安全）
load_dataset_pheno_all <- function(gse, hint) {
  d <- file.path(CONFIG$DATA_DIR, gse)
  rds <- file.path(d, paste0(gse, "_eset.rds"))
  sm <- list.files(d, pattern = "_series_matrix\\.txt(\\.gz)?$", full.names = TRUE)
  if (file.exists(rds) && length(sm) == 0) {
    return(extract_pheno(readRDS(rds), gse, hint))
  }
  if (length(sm) == 0) return(NULL)
  pd_list <- lapply(sm, function(f) pData(getGEO(filename = f, getGPL = FALSE)))
  pd_all <- do.call(rbind, pd_list)   # 列不一致时 rbind 报错 → 转人工，属显眼错误
  log_msg(gse, ": 合并 ", length(sm), " 个平台的 pData（共 ", nrow(pd_all), " 行）")
  extract_pheno_from_pd(pd_all, gse, hint)
}

# 探针 → 基因 symbol（GPL 注释；多探针取行均值最大）
map_probes_to_symbol <- function(expr, feature_tbl, gse) {
  # Brainarray ENTREZG 重注释平台（如 GPL20880 / GSE117525）：行名 = "<EntrezID>_at"，
  # GPL 表只有 ID/Description/SPOT_ID、无 symbol 列 —— 改走 org.Hs.eg.db 映射
  if (mean(grepl("^[0-9]+_at$", rownames(expr))) > 0.5) {
    log_msg(gse, ": 检测到 ENTREZG 口径（<Entrez>_at），改走 org.Hs.eg.db ENTREZID→SYMBOL")
    ez <- sub("_at$", "", rownames(expr))
    mapped <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(ez),
                                    keytype = "ENTREZID", columns = "SYMBOL")
    mapped <- mapped[!duplicated(mapped$ENTREZID), ]
    sym <- mapped$SYMBOL[match(ez, mapped$ENTREZID)]
    keep <- !is.na(sym) & sym != ""
    expr2 <- expr[keep, , drop = FALSE]; sym <- sym[keep]
    ord <- order(rowMeans(expr2, na.rm = TRUE), decreasing = TRUE)
    expr2 <- expr2[ord, , drop = FALSE]; sym <- sym[ord]
    keep_row <- !duplicated(sym)
    collapsed <- expr2[keep_row, , drop = FALSE]
    rownames(collapsed) <- sym[keep_row]
    log_msg(gse, ": 探针 ", nrow(expr), " → 基因 ", nrow(collapsed), "（ENTREZG 口径）")
    return(collapsed)
  }
  # feature_tbl: GPL 注释 data.frame；需在其中定位 probe id 列与 symbol 列
  probe_col <- colnames(feature_tbl)[1]  # GPL 表首列恒为 ID
  sym_col <- grep("gene.?symbol", colnames(feature_tbl), ignore.case = TRUE, value = TRUE)
  if (length(sym_col) == 0) sym_col <- grep("symbol|gene_assignment", colnames(feature_tbl),
                                            ignore.case = TRUE, value = TRUE)
  if (length(sym_col) == 0) stop(gse, ": GPL 表中找不到 gene symbol 列，请人工检查")
  sym_col <- sym_col[1]
  mp <- data.frame(probe = as.character(feature_tbl[[probe_col]]),
                   symbol = as.character(feature_tbl[[sym_col]]), stringsAsFactors = FALSE)
  # 清洗：取 "///" 分隔的第一个 symbol；去掉空/NA/--- 与 Affy 控制探针
  mp$symbol <- sub("\\s*///.*$", "", mp$symbol)
  mp <- mp[!is.na(mp$symbol) & mp$symbol != "" & mp$symbol != "---", ]
  mp <- mp[!grepl("^AFFX", mp$probe), ]
  common <- intersect(rownames(expr), mp$probe)
  if (length(common) < 100) stop(gse, ": 探针与 GPL 注释匹配数 < 100，映射异常，请人工检查")
  expr2 <- expr[common, , drop = FALSE]
  mp <- mp[match(common, mp$probe), ]
  # 按行均值降序后取每个 symbol 的第一行 = 均值最大探针
  # （不用 avereps——它对同 ID 行取平均，与本设计"最大均值探针"声明不符）
  ord <- order(rowMeans(expr2, na.rm = TRUE), decreasing = TRUE)
  expr2 <- expr2[ord, , drop = FALSE]; mp <- mp[ord, ]
  keep_row <- !duplicated(mp$symbol)
  collapsed <- expr2[keep_row, , drop = FALSE]
  rownames(collapsed) <- mp$symbol[keep_row]
  log_msg(gse, ": 探针 ", nrow(expr), " → 基因 ", nrow(collapsed),
          "（注释列 '", sym_col, "'）")
  collapsed
}

# Ensembl → symbol（RNA-seq 计数；自动剥版本号）
ensembl_to_symbol <- function(ids, gse) {
  clean <- sub("\\..*$", "", ids)
  is_ens <- grepl("^ENSG", clean)
  out <- data.frame(id = ids, symbol = NA_character_, stringsAsFactors = FALSE)
  if (any(is_ens)) {
    mapped <- AnnotationDbi::select(org.Hs.eg.db, keys = clean[is_ens],
                                    keytype = "ENSEMBL", columns = "SYMBOL")
    mapped <- mapped[!duplicated(mapped$ENSEMBL), ]
    out$symbol[is_ens] <- mapped$SYMBOL[match(clean[is_ens], mapped$ENSEMBL)]
  }
  out$symbol[!is_ens] <- ids[!is_ens]  # 已是 symbol 的原样保留
  out
}

# ----------------------------- 平台标准化 -----------------------------------
# Affymetrix：直接对 eSet 表达矩阵做 RMA 等价处理需 CEL；GEO series matrix
# 通常已 RMA/MAS5。判定标准：若 max < 100 视为已 log，跳过；否则 log2(x+1)。
normalize_affy_matrix <- function(expr, gse) {
  mx <- max(expr, na.rm = TRUE)
  if (mx < 100) {
    log_msg(gse, ": 矩阵已在 log 尺度（max=", round(mx, 1), "），不再变换")
    return(expr)
  }
  log_msg(gse, ": 线性尺度（max=", round(mx, 1), "），执行 log2(x+1) + quantile 标准化")
  expr <- log2(expr + 1)
  normalizeBetweenArrays(expr, method = "quantile")
}

# Illumina 表达芯片：limma::neqc 需要原始强度+对照探针；series matrix 一般
# 已 quantile。做 log 判定 + quantile 兜底。
normalize_ilmn_matrix <- function(expr, gse) {
  mx <- max(expr, na.rm = TRUE)
  if (mx < 100) {
    log_msg(gse, ": 矩阵已在 log 尺度，执行 quantile 兜底标准化")
    return(normalizeBetweenArrays(expr, method = "quantile"))
  }
  log_msg(gse, ": 线性尺度，log2(x+1) + quantile")
  normalizeBetweenArrays(log2(expr + 1), method = "quantile")
}

# RNA-seq：读取原始计数 → DESeq2 VST（返回 log 尺度矩阵）
normalize_seq_counts <- function(gse) {
  d <- file.path(CONFIG$DATA_DIR, gse)
  f <- list.files(d, pattern = paste0(CONFIG$SEQ_COUNT_FILENAME, ".*\\.(txt|csv|tsv)(\\.gz)?$"),
                  full.names = TRUE, ignore.case = TRUE)
  if (length(f) == 0) {
    stop(gse, ": 未找到原始计数文件。请从 GEO supplementary 下载计数矩阵，命名为 ",
         gse, "_raw_counts.txt 放入 ", d, "/（", gse, " 的 series matrix 不含可用计数）")
  }
  d_raw <- fread(f[1], data.table = FALSE)
  counts <- as.matrix(d_raw[, -1])
  rownames(counts) <- make.unique(as.character(d_raw[[1]]))  # as.matrix 不认 rownames 参数
  storage.mode(counts) <- "integer"
  log_msg(gse, ": 读入计数 ", nrow(counts), " 基因 × ", ncol(counts), " 样本")
  dds <- DESeqDataSetFromMatrix(counts, colData = data.frame(row.names = colnames(counts)),
                                design = ~ 1)
  dds <- dds[rowSums(counts(dds)) >= 10, ]   # 低表达粗过滤（正式过滤在模块 B）
  vsd <- vst(dds, blind = TRUE)
  assay(vsd)
}

# GSE114007 特例：VST 后 ComBat（已知批次 = 测序平台 GPL11154 vs GPL18573）
combat_gse114007 <- function(vst_mat, pheno) {
  batch <- ifelse(grepl("GPL11154", pheno$platform_id), "GPL11154", "GPL18573")
  if (length(unique(batch)) < 2) {
    log_msg("GSE114007: ⚠ WARNING 仅检测到单一平台，跳过 ComBat——与计划书'已知批次必须校正'",
            "矛盾，请人工核对平台标注（pheno$platform_id）后再继续")
    return(vst_mat)
  }
  log_msg("GSE114007: ComBat 校正平台批次（", paste(table(batch), collapse = " / "), "）")
  ComBat(vst_mat, batch = batch, par.prior = TRUE)
}

# ----------------------------- QC 绘图 --------------------------------------
qc_pdf <- function(expr_pre, expr_post, pheno, gse) {
  if (!CONFIG$MAKE_QC_PDF) return(invisible(NULL))
  out <- file.path(CONFIG$OUT_DIR, gse, paste0(gse, "_QC.pdf"))
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  cond <- if (!is.null(pheno$condition)) as.factor(pheno$condition) else factor("NA")
  pal <- c("#7A9B76", "#C17C5B", "#8B7D6B", "#4A5D43", "#B09A7A", "#6B7F8C")
  pdf(out, width = 9, height = 6)
  # 1) 标准化前后箱线图
  par(mfrow = c(1, 2), mar = c(8, 4, 3, 1))
  boxplot(expr_pre, las = 2, outline = FALSE, main = paste0(gse, " 标准化前"),
          col = pal[as.integer(cond) %% length(pal) + 1], cex.axis = 0.6)
  boxplot(expr_post, las = 2, outline = FALSE, main = paste0(gse, " 标准化后"),
          col = pal[as.integer(cond) %% length(pal) + 1], cex.axis = 0.6)
  # 2) 密度曲线
  par(mfrow = c(1, 1))
  d_max <- max(vapply(seq_len(min(ncol(expr_post), 60)),
                      function(j) max(density(na.omit(expr_post[, j]))$y), numeric(1)))
  plot(density(na.omit(expr_post[, 1])), main = paste0(gse, " 表达密度"), lwd = 1,
       xlab = "expression", col = pal[1], ylim = c(0, d_max * 1.05))
  for (j in 2:min(ncol(expr_post), 60)) {
    lines(density(na.omit(expr_post[, j])), col = pal[as.integer(cond[j]) %% length(pal) + 1])
  }
  legend("topright", legend = levels(cond), col = pal[seq_along(levels(cond))], lty = 1, cex = 0.8)
  # 3) PCA（按 condition 着色；剔除零方差行防止 scale 除零）
  v <- apply(expr_post, 1, var, na.rm = TRUE)
  expr_pca <- expr_post[v > 0, , drop = FALSE]
  pc <- prcomp(t(expr_pca), scale. = TRUE)
  var_exp <- round(100 * summary(pc)$importance[2, 1:2], 1)
  plot(pc$x[, 1:2], col = pal[as.integer(cond) %% length(pal) + 1], pch = 19,
       xlab = paste0("PC1 (", var_exp[1], "%)"), ylab = paste0("PC2 (", var_exp[2], "%)"),
       main = paste0(gse, " PCA"))
  text(pc$x[, 1:2], labels = colnames(expr_post), pos = 3, cex = 0.5)
  legend("topright", legend = levels(cond), col = pal[seq_along(levels(cond))], pch = 19, cex = 0.8)
  dev.off()
  log_msg(gse, ": QC 图 → ", out)
}

# ----------------------------- 表型整理 -------------------------------------
# 从 pData 提取 condition 列：优先使用 hint 列，随后正则兜底；
# 输出统一 pheno 表：sample / condition_raw / condition / platform_id / design_flag
extract_pheno <- function(eset, gse, hint) extract_pheno_from_pd(pData(eset), gse, hint)

extract_pheno_from_pd <- function(pd, gse, hint) {
  cand <- c(hint, grep("characteristics|source_name|title", colnames(pd), value = TRUE))
  cand <- cand[!is.na(cand) & cand %in% colnames(pd)]
  hit <- NULL
  for (cc in cand) {
    vals <- tolower(as.character(pd[[cc]]))
    # 注：OA 词界用 (^|[^a-z])oa($|[^a-z]) 而非 \\boa\\b——"pbmc_oa" 中下划线是词字符，\\b 不触发
    if (any(grepl("osteoarthriti|(^|[^a-z])oa($|[^a-z])|sarcopen|control|healthy|normal|preserved|affected|young|\\bold\\b|older|elderly|frail|metabolic", vals))) {
      hit <- cc; break
    }
  }
  if (is.null(hit)) {
    log_msg(gse, ": ⚠ 未能自动定位表型列，condition 全部置 NA，请人工填写 pheno 表")
    cond_raw <- rep(NA_character_, nrow(pd))
  } else {
    cond_raw <- as.character(pd[[hit]])
    log_msg(gse, ": 表型取自列 '", hit, "'")
  }
  cond <- tolower(cond_raw)
  cond[is.na(cond)] <- "other"   # 如 GSE117525 随访行在 baseline 分组列为 NA → 一律 other（自动剔除待复核）
  condition <- standardize_condition(cond)
  platform_id <- if ("platform_id" %in% colnames(pd)) as.character(pd$platform_id) else NA_character_
  out_df <- data.frame(
    sample = rownames(pd),
    condition_raw = cond_raw,
    condition = condition,
    platform_id = platform_id,
    design_flag = ifelse(gse %in% names(DESIGN_FLAGS), DESIGN_FLAGS[[gse]], ""),
    stringsAsFactors = FALSE
  )
  # GSE51588 供体编号：title 形如 "Normal-LT-5"/"OA-MT-12"，尾号即供体 ID（2026-08-22 实测定）
  if (gse == "GSE51588" && "title" %in% colnames(pd)) {
    did <- sub("^.*-([0-9]+)$", "\\1", as.character(pd$title))
    if (all(grepl("^[0-9]+$", did))) {
      out_df$donor_id <- did
      log_msg(gse, ": 已从 title 提取 donor_id（", length(unique(did)), " 供体 × 2 部位）")
    } else log_msg(gse, ": ⚠ donor_id 提取异常，请人工核对 title 列")
  }
  out_df
}

# 条件标签统一化（保守匹配；不能确定的一律 other，留人工核对）
standardize_condition <- function(cond_lower) {
  out <- rep("other", length(cond_lower))
  # 先判对照（"non-oa" 含子串 "oa"、"non-frail" 含 "frail"，必须先于疾病规则判定）
  # "non-arthritic"（GSE168505/GSE169454 的对照写法）也要判对照
  out[grepl("non-?oa|non-?arthritic|non-?frail|robust|healthy|normal|control", cond_lower) &
        !grepl("osteoarthriti|sarcopen", cond_lower)] <- "control"
  out[grepl("osteoarthriti", cond_lower) |
        (grepl("(^|[^a-z])oa($|[^a-z])", cond_lower) & !grepl("non-?oa|control|healthy", cond_lower))] <- "OA"
  out[grepl("preserved", cond_lower)] <- "OA_preserved"
  out[grepl("affected", cond_lower) & grepl("cartilage", cond_lower)] <- "OA_affected"
  out[grepl("sarcopen", cond_lower)] <- "sarcopenia"
  # 否定式后置覆盖："sarcopenia status: no"（GSE111006/111010/111016 的对照写法）必须判对照
  out[grepl("sarcopen", cond_lower) & grepl(": ?no\\b", cond_lower)] <- "control"
  out[grepl("pre-?sarcopen", cond_lower)] <- "pre_sarcopenia"
  # 虚弱队列（GSE117525/GSE144304，验证侧）：(pre-)frail 合并标 frailty
  out[grepl("pre-?frail|frail", cond_lower) & !grepl("non-?frail", cond_lower)] <- "frailty"
  # 代谢综合征老年（GSE136344：字面量 + 标题缩写 "Muscle_Old SX_..."）
  out[grepl("metabolic syndrome", cond_lower)] <- "metabolic_syndrome"
  out[out == "other" & grepl("(^|[^a-z])sx($|[^a-z])", cond_lower)] <- "metabolic_syndrome"
  # GSE144304 的 "group: Fit" = 健康老年对照
  out[out == "other" & grepl("(^|[^a-z])fit($|[^a-z])", cond_lower)] <- "control"
  # 衰老代理队列（GSE1428/GSE25941/GSE136344）：仅在未归类时标 young / aging
  # 注意①必须先判 young——"18-25 years old" 的 "old" 会误中 aging 规则；
  # 注意②词界用 [^a-z] 而非 \\b——下划线是词字符，"muscle_old_biological" 中 \\b 不触发
  out[out == "other" & grepl("(^|[^a-z])young($|[^a-z])", cond_lower)] <- "young"
  out[out == "other" & grepl("(^|[^a-z])old($|[^a-z])|(^|[^a-z])older($|[^a-z])|elderly|aged", cond_lower)] <- "aging"
  # 随访/干预后样本一律置 other，强制人工复核（GSE117525 含 6 个月运动随访）
  out[grepl("follow-?up|post-?training|after training", cond_lower)] <- "other"
  out
}

# ----------------------------- 单数据集主流程 -------------------------------
process_one <- function(gse, row) {
  log_msg("==== ", gse, " （", row$platform_type, "）====")
  outdir <- file.path(CONFIG$OUT_DIR, gse)
  done_marker <- file.path(outdir, ".moduleA_done")
  if (file.exists(done_marker)) { log_msg(gse, ": 已处理，跳过（删除 ", done_marker, " 可重跑）"); return(invisible(TRUE)) }

  if (row$platform_type == "seq" && row$tissue_group == "singlecell") {
    log_msg(gse, ": 单细胞数据集，模块 A 不处理（移交模块 G）；仅登记")
    return(invisible(TRUE))
  }

  # ---- 1. 读取表达矩阵（预处理前）----
  if (row$platform_type == "seq") {
    expr_pre <- normalize_seq_counts(gse)          # 已是 VST（log 尺度）
    pheno <- load_dataset_pheno_all(gse, row$condition_col_hint)
    if (is.null(pheno)) {
      log_msg(gse, ": ⚠ 无 series matrix 可提供表型，condition 全部置 other，需人工填写")
      pheno <- data.frame(sample = colnames(expr_pre), condition_raw = NA,
                          condition = "other", platform_id = NA, design_flag = "",
                          stringsAsFactors = FALSE)
    }
    # 计数矩阵列名 ↔ pheno$sample 必须显式对齐（列名多为样本标题而非 GSM 时
    # 需要人工在 pheno CSV 中补 sample 映射；此处硬失败，杜绝静默错位）
    pheno <- pheno[match(colnames(expr_pre), pheno$sample), ]
    if (any(is.na(pheno$sample))) {
      miss <- colnames(expr_pre)[is.na(pheno$sample)]
      stop(gse, ": ", length(miss), " 个计数矩阵列名无法匹配 pheno 样本（前 5 个：",
           paste(head(miss, 5), collapse = ", "),
           "）。请人工核对计数文件列名与 GSM/样本标题映射后重跑")
    }
    expr_post <- expr_pre
    if (gse == "GSE114007") expr_post <- combat_gse114007(expr_pre, pheno)
    # Ensembl → symbol（同 symbol 多 Ensembl 时取均值最大者，与芯片口径一致）
    mp <- ensembl_to_symbol(rownames(expr_post), gse)
    keep <- !is.na(mp$symbol) & mp$symbol != ""
    expr_post <- expr_post[keep, , drop = FALSE]
    sym <- mp$symbol[keep]
    ord <- order(rowMeans(expr_post), decreasing = TRUE)
    expr_post <- expr_post[ord, , drop = FALSE]; sym <- sym[ord]
    first <- !duplicated(sym)
    expr_post <- expr_post[first, , drop = FALSE]
    rownames(expr_post) <- sym[first]
  } else {
    eset <- load_dataset_eset(gse)
    expr_raw <- exprs(eset)
    pheno <- extract_pheno(eset, gse, row$condition_col_hint)
    expr_pre <- if (row$platform_type == "affy") normalize_affy_matrix(expr_raw, gse)
                else normalize_ilmn_matrix(expr_raw, gse)
    gpl <- annotation(eset)
    feature_tbl <- if (gpl != "" ) {
      tbl <- tryCatch(Table(getGEO(gpl)), error = function(e) NULL)
      if (is.null(tbl)) fData(eset) else tbl
    } else fData(eset)
    expr_post <- map_probes_to_symbol(expr_pre, feature_tbl, gse)
    expr_pre <- expr_pre[rownames(expr_pre) %in% rownames(expr_raw), ]
  }

  # ---- 2. QC ----
  qc_pdf(if (row$platform_type == "seq") expr_pre else expr_pre[1:min(20000, nrow(expr_pre)), ],
         expr_post[1:min(20000, nrow(expr_post)), ], pheno, gse)

  # ---- 3. 导出 ----
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(expr_post, file.path(outdir, paste0(gse, "_expr_gene.rds")))
  write.csv(pheno, file.path(outdir, paste0(gse, "_pheno.csv")), row.names = FALSE)
  writeLines(c(sprintf("processed_at: %s", Sys.time()),
               sprintf("genes: %d", nrow(expr_post)),
               sprintf("samples: %d", ncol(expr_post)),
               sprintf("platform_type: %s", row$platform_type),
               sprintf("design_flag: %s", pheno$design_flag[1])),
             done_marker)
  log_msg(gse, ": ✅ 完成（", nrow(expr_post), " 基因 × ", ncol(expr_post), " 样本）")
  invisible(TRUE)
}

# ----------------------------- 主循环 ---------------------------------------
## 冒烟测试开关：在 source() 本脚本【之前】先运行
##   SMOKE_GSES <- c("GSE55235", "GSE48556")   # 只跑这两个
## 即可进入冒烟模式；不设（或设为 NULL）则全量 21 个数据集
if (!exists("SMOKE_GSES")) SMOKE_GSES <- NULL
DATASETS_RUN <- if (is.null(SMOKE_GSES)) DATASETS else {
  miss <- setdiff(SMOKE_GSES, DATASETS$gse)
  if (length(miss) > 0) stop("SMOKE_GSES 含未知数据集：", paste(miss, collapse = ", "))
  log_msg("⚠ 冒烟测试模式：仅处理 ", paste(SMOKE_GSES, collapse = ", "))
  DATASETS[DATASETS$gse %in% SMOKE_GSES, ]
}

dir.create(CONFIG$OUT_DIR, showWarnings = FALSE)
log_msg("模块 A 启动；待处理数据集：", nrow(DATASETS_RUN))

results <- data.frame(gse = DATASETS_RUN$gse, status = "pending", stringsAsFactors = FALSE)
for (i in seq_len(nrow(DATASETS_RUN))) {
  gse <- DATASETS_RUN$gse[i]
  st <- tryCatch({ process_one(gse, DATASETS_RUN[i, ]); "OK" },
                 error = function(e) { log_msg(gse, ": ❌ ", conditionMessage(e)); paste0("ERROR: ", conditionMessage(e)) })
  results$status[results$gse == gse] <- st
}

write.csv(results, file.path(CONFIG$OUT_DIR, sprintf("moduleA_status_%s.csv", Sys.Date())), row.names = FALSE)
log_msg("模块 A 全部结束。状态汇总：")
print(table(substr(results$status, 1, 5)))
cat("\n下一步：检查 ", CONFIG$OUT_DIR, "/moduleA_status_", Sys.Date(),
    ".csv 与日志 ", CONFIG$LOG_FILE, "；\n",
    "对 condition=other 或表型列定位失败的样本，人工核对 <GSE>_pheno.csv 后再进入模块 B。\n", sep = "")
