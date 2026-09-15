# ============================================================================
# 07_moduleF_immune_infiltration_fix3_20260827.R
# 模块 F｜免疫浸润分析（xCell + ssGSEA 28 细胞 + CIBERSORT 22 细胞三方法共识，
# 队列内统计 + 跨队列方向汇总）
# 版本：v1.1（fix3，2026-08-27）
# fix3：CIBERSORT 由"授权占位"升级为实跑——Stanford 官网申请长期 pending，
#   改用完全合规的双通道方案（不下载/不使用 Stanford 分发的 CIBERSORT.R 与
#   官方 LM22.txt，无任何许可风险）：
#   ① 算法：CRAN 包 dicepro（v1.0.2，2026-06 维护）内置 CIBERSORT ν-SVR 完整
#      开源实现——e1071 nu-regression、三档 nu（0.25/0.50/0.75）按 RMSE 择优、
#      QN 开关、anti-log（max<50 → 2^Y）、置换 P 值，与 Newman et al. 2015
#      Nat Methods 算法描述逐项对应（沙箱逐行源码核验存档）；
#   ② 签名矩阵：LM22 即 Newman 2015 原文 Supplementary Table 1（PMC4739640，
#      NIH 官方托管的已发表补充材料）——浏览器下载后放项目根（txt 或 xlsx
#      均可，本脚本自动识别+结构校验 547×22）；此渠道为学界千篇论文标准
#      实践，不触 Stanford 许可墙。
#   方法学建议写法：Deconvolution was performed with the CIBERSORT nu-SVR
#   algorithm (Newman et al., 2015) as implemented in the CRAN package
#   dicepro v1.0.2, using the LM22 signature matrix from the original
#   publication's Supplementary Table 1 (perm=100).
#   QN 口径与官方一致：芯片 QN=TRUE；RNA-seq QN=FALSE（按 cf$rnaseq）。
#   依赖：install.packages("dicepro")；BiocManager::install("preprocessCore")
# ----------------------------------------------------------------------------
# 口径与设计：
#   1) 四层队列：
#      - OA_training（OA vs control）：GSE114007/GSE55235/GSE55457/GSE12021/
#        GSE51588；GSE51588 为双部位 donor 结构，组间比较优先
#        limma::duplicateCorrelation(donor_id)，失败降级 donor 去重 + 普通 limma；
#        相关分析一律用 donor 去重样本，防伪重复。
#      - Muscle_sarcopenia（sarcopenia vs control）：GSE111006/GSE111010/
#        GSE111016；GSE111006（case=4）与 GSE111010（case=9）自动 low_power 标记。
#      - Blood_exploratory（OA vs control）：GSE48556，EXPLORATORY 口径，
#        不做跨队列硬合并。
#      - Aging_proxy（aging vs young）：GSE1428/GSE25941，仅作 aging proxy，
#        不得写作 sarcopenia 证据。
#   2) 方法：xCell::xCellAnalysis(as.matrix(expr))（tryCatch 降级）；
#      ssGSEA 28 免疫细胞（Charoentong et al. 2017 28 immune-cell metagene
#      集合的开放复现表，782 行 × 28 细胞类型，列 Metagene/Cell type/Immunity，
#      默认读 moduleF_resources/ssGSEA_28immune_Charoentong2017.csv）；
#      GSVA 新旧 API 双分支（>= 1.50 用 ssgseaParam() + gsva()，否则旧
#      gsva(..., method="ssgsea")；新 API 失败自动 fallback 旧 API）。
#   3) 统计：组间比较只在队列内做（limma moderated t，design=~condition），
#      BH 按 队列×方法 内部校正；跨队列仅汇总方向（禁止 pooled p）；
#      相关分析为队列内 Spearman（10 标志物 × 各细胞类型），BH 按 队列×方法。
#   4) 不做跨队列原始表达/得分合并统计；aging proxy 不写作 sarcopenia。
# ----------------------------------------------------------------------------
# 输入：
#   processed/<GSE>/<GSE>_expr_gene.rds 与 <GSE>_pheno.csv（02_moduleA 产物）
#   results/moduleC/moduleC_gene_tracks.csv（04_moduleC 产物）
#   results/moduleD_v2/moduleD_v2_biomarkers.rds（$biomarkers，05b 产物）
#   moduleF_resources/ssGSEA_28immune_Charoentong2017.csv（见上来源）
# CIBERSORT（fix3 实跑）：LM22.txt 或 LM22.xlsx 放项目根即可（获取指引见头部）；
#   缺文件则跳过 CIBERSORT 并写 moduleF_cibersort_status.csv（xCell+ssGSEA 不受影响）
# 运行：source("07_moduleF_immune_infiltration.R")；输出 results/moduleF/
# ============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(GSVA)
  library(ggplot2)
  library(pheatmap)
  library(data.table)
})
## xCell 按需 requireNamespace（缺失时给出安装提示而非静默失败）
set.seed(20260823)

# ----------------------------- CONFIG --------------------------------------
CONFIG <- list(
  IN_DIR  = "processed",
  OUT_DIR = "results/moduleF",
  FIG_DIR = "results/moduleF/figs",
  LOG_FILE = file.path("results/moduleF", sprintf("moduleF_log_%s.txt", Sys.Date())),
  TRACK_FILE = "results/moduleC/moduleC_gene_tracks.csv",
  V2_RDS     = "results/moduleD_v2/moduleD_v2_biomarkers.rds",
  SSGSEA_CSV = "moduleF_resources/ssGSEA_28immune_Charoentong2017.csv",
  CIBERSORT_LM22_TXT  = "LM22.txt",    # 优先：tab 分隔（首列基因名+22 细胞列）
  CIBERSORT_LM22_XLSX = "LM22.xlsx",   # 备选：PMC4739640 Supplementary Table 1 原样放入，自动转换
  CIBERSORT_PERM = 100,                # 置换次数（Newman 2015 默认；0=跳过 P 值提速）
  CIBERSORT_P_CUTOFF = 0.05,           # 样本级全局 P 值质控标记阈值
  LOW_POWER_MIN_CASE = 10,   # n_case < 10 → low_power（GSE111006=4, GSE111010=9 命中）
  COR_MIN_N = 8,             # Spearman 最少样本数
  TOP_N_BOXPLOT = 6,         # 每层 top 差异细胞 boxplot 个数
  XCELL_PARALLEL_SZ = 1,     # xCell 底层 GSVA 线程数；Windows/SOCK 下 1 最稳
  ## 队列配置（tier / case-control 候选标签；标签实测匹配，防口径漂移）
  ## rnaseq 口径与 02_moduleA_preprocess.R platform_type 对齐：seq=TRUE；affy/ilmn/agilent=FALSE
  COHORTS = list(
    list(gse = "GSE114007", tier = "OA_training", rnaseq = TRUE,
         case_cand = c("OA"), ctrl_cand = c("control", "healthy", "normal")),
    list(gse = "GSE55235",  tier = "OA_training", rnaseq = FALSE,
         case_cand = c("OA"), ctrl_cand = c("control", "healthy", "normal")),
    list(gse = "GSE55457",  tier = "OA_training", rnaseq = FALSE,
         case_cand = c("OA"), ctrl_cand = c("control", "healthy", "normal")),
    list(gse = "GSE12021",  tier = "OA_training", rnaseq = FALSE,
         case_cand = c("OA"), ctrl_cand = c("control", "healthy", "normal")),
    list(gse = "GSE51588",  tier = "OA_training", rnaseq = FALSE,   # 双部位 donor 结构
         case_cand = c("OA"), ctrl_cand = c("control", "healthy", "normal"),
         expect_donor = TRUE),   # pheno 必须带 donor_id，缺失即 stop（口径破坏，禁止静默退化）
    list(gse = "GSE111006", tier = "Muscle_sarcopenia", rnaseq = TRUE,
         case_cand = c("sarcopenia"), ctrl_cand = c("control", "healthy", "non_sarcopenic")),
    list(gse = "GSE111010", tier = "Muscle_sarcopenia", rnaseq = TRUE,
         case_cand = c("sarcopenia"), ctrl_cand = c("control", "healthy", "non_sarcopenic")),
    list(gse = "GSE111016", tier = "Muscle_sarcopenia", rnaseq = TRUE,
         case_cand = c("sarcopenia"), ctrl_cand = c("control", "healthy", "non_sarcopenic")),
    list(gse = "GSE48556",  tier = "Blood_exploratory", rnaseq = FALSE,   # EXPLORATORY：不跨队列硬合并
         case_cand = c("OA"), ctrl_cand = c("control", "healthy")),
    list(gse = "GSE1428",   tier = "Aging_proxy", rnaseq = FALSE,   # aging proxy，非 sarcopenia 证据
         case_cand = c("aging", "aged", "old"), ctrl_cand = c("young")),
    list(gse = "GSE25941",  tier = "Aging_proxy", rnaseq = FALSE,
         case_cand = c("aging", "aged", "old"), ctrl_cand = c("young"))
  ),
  TIER_NOTES = c(
    OA_training        = "OA 训练层",
    Muscle_sarcopenia  = "肌肉 sarcopenia 层",
    Blood_exploratory  = "EXPLORATORY：血液探索层，不做跨队列硬合并",
    Aging_proxy        = "aging proxy（衰老代理，不得写作 sarcopenia 证据）"
  )
)
dir.create(CONFIG$OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CONFIG$FIG_DIR, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  txt <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = ""))
  cat(txt, "\n"); cat(txt, "\n", file = CONFIG$LOG_FILE, append = TRUE)
}
log_msg("==== 模块 F 免疫浸润分析启动（", Sys.Date(), "）====")

# ----------------------------- 硬校验：全部输入存在性 -----------------------
input_files <- c(CONFIG$TRACK_FILE, CONFIG$V2_RDS, CONFIG$SSGSEA_CSV,
                 unlist(lapply(CONFIG$COHORTS, function(cf)
                   file.path(CONFIG$IN_DIR, cf$gse,
                             paste0(cf$gse, c("_expr_gene.rds", "_pheno.csv"))))))
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0) {
  if (CONFIG$SSGSEA_CSV %in% missing_files)
    log_msg("⚠ 缺少 ", CONFIG$SSGSEA_CSV,
            "——请下载 Charoentong et al. 2017 28 免疫细胞 metagene 开放复现表",
            "（782 行 × 28 细胞类型，列 Metagene/Cell type/Immunity）放入该路径后重跑")
  stop("缺少输入文件：\n  ", paste(missing_files, collapse = "\n  "),
       "\n——请先运行上游模块（02_moduleA / 04_moduleC / 05b_moduleD_v2）并备好 ssGSEA CSV")
}

# ----------------------------- 上游输入读取与口径校验 ------------------------
track_tab <- read.csv(CONFIG$TRACK_FILE, stringsAsFactors = FALSE)
if (!all(c("gene", "track") %in% names(track_tab)))
  stop("moduleC_gene_tracks.csv 列名口径漂移，请核对模块 C 输出")
v2 <- readRDS(CONFIG$V2_RDS)
biomarkers <- v2$biomarkers
if (is.null(biomarkers) || length(biomarkers) == 0)
  stop("moduleD_v2_biomarkers.rds 缺少 $biomarkers——请回查 05b 输出")
log_msg("模块 D v2 标志物（", length(biomarkers), " 个）：",
        paste(biomarkers, collapse = ", "))

## ssGSEA 28 免疫细胞基因集（Charoentong 2017 开放复现表）
## 兼容 Windows/Excel 可能加入的 UTF-8 BOM 或列名首尾空格；失配时打印实际列名便于排查
ssgsea_tab <- tryCatch(
  data.table::fread(CONFIG$SSGSEA_CSV, data.table = FALSE,
                    encoding = "UTF-8", check.names = FALSE),
  error = function(e)
    read.csv(CONFIG$SSGSEA_CSV, stringsAsFactors = FALSE,
             fileEncoding = "UTF-8-BOM", check.names = FALSE)
)
## 用字节级 UTF-8 BOM 常量清洗列名，避免不同 R 正则方言对 \u 转义解释不一致
utf8_bom <- rawToChar(as.raw(c(0xEF, 0xBB, 0xBF)))
names(ssgsea_tab) <- trimws(sub(paste0("^", utf8_bom), "", enc2utf8(names(ssgsea_tab))))
if (!all(c("Metagene", "Cell type", "Immunity") %in% names(ssgsea_tab)))
  stop("ssGSEA CSV 列名口径漂移（预期 Metagene/Cell type/Immunity）：", CONFIG$SSGSEA_CSV,
       "\n实际读入列名为：", paste(names(ssgsea_tab), collapse = " | "),
       "\n请确认下载的是随 07 脚本一起交付的 ssGSEA_28immune_Charoentong2017.csv，且未被 Excel 改写或重命名覆盖")
if (anyNA(ssgsea_tab$Metagene) || anyNA(ssgsea_tab[["Cell type"]]))
  stop("ssGSEA CSV 含缺失 Metagene/Cell type——请回查资源文件")
gene_sets <- split(ssgsea_tab$Metagene, ssgsea_tab[["Cell type"]])
immunity_map <- ssgsea_tab$Immunity[match(names(gene_sets), ssgsea_tab[["Cell type"]])]
names(immunity_map) <- names(gene_sets)   # 未命名时 immunity_map[rownames(ss)] 按位置错配
if (length(gene_sets) != 28)
  log_msg("⚠ ssGSEA 细胞类型数 = ", length(gene_sets), "（预期 28，请核对资源版本）")
log_msg("ssGSEA 基因集：", length(gene_sets), " 细胞类型，",
        length(unique(ssgsea_tab$Metagene)), " 唯一 metagene（Charoentong 2017）")

# ----------------------------- 通用函数 -------------------------------------
## case/control 标签实测匹配（候选制，复用 06_moduleE resolve_case 思路）
resolve_label <- function(pheno, candidates, gse, what = "case") {
  avail <- unique(pheno$condition)
  hit <- intersect(candidates, avail)
  if (length(hit) == 0)
    stop(gse, ": ", what, " 候选（", paste(candidates, collapse = "/"),
         "）均不在 condition 实测取值（", paste(avail, collapse = ","), "）中——请人工核对")
  if (hit[1] != candidates[1] || length(hit) > 1)
    log_msg("  ", gse, ": ", what, " 标签实测匹配为 '", hit[1], "'（候选优先级自动选择）")
  hit[1]
}

## 单队列装配：pheno/expr 样本硬一致 + 标签实测匹配 + 仅保留 case/control
assemble_cohort <- function(cf) {
  gse <- cf$gse
  expr  <- readRDS(file.path(CONFIG$IN_DIR, gse, paste0(gse, "_expr_gene.rds")))
  pheno <- read.csv(file.path(CONFIG$IN_DIR, gse, paste0(gse, "_pheno.csv")),
                    stringsAsFactors = FALSE)
  if (!all(c("sample", "condition") %in% names(pheno)))
    stop(gse, ": pheno 缺少 sample/condition 列——请回查模块 A 输出")
  if (!setequal(pheno$sample, colnames(expr)))
    stop(gse, ": pheno 样本与 expr 列不完全一致（pheno ",
         length(setdiff(pheno$sample, colnames(expr))), " 个多出 / expr ",
         length(setdiff(colnames(expr), pheno$sample)), " 个多出）——请回查模块 A 输出")
  pheno <- pheno[match(colnames(expr), pheno$sample), ]
  case_lab <- resolve_label(pheno, cf$case_cand, gse, "case")
  ctrl_lab <- resolve_label(pheno, cf$ctrl_cand, gse, "control")
  keep <- pheno$condition %in% c(ctrl_lab, case_lab)
  expr <- expr[, keep, drop = FALSE]; pheno <- pheno[keep, ]
  if (sum(pheno$condition == case_lab) == 0 || sum(pheno$condition == ctrl_lab) == 0)
    stop(gse, ": 病例/对照有一侧为 0")
  list(expr = as.matrix(expr), pheno = pheno, gse = gse, tier = cf$tier,
       case_lab = case_lab, ctrl_lab = ctrl_lab)
}

## xCell（缺失包时给安装提示；显式加载 xCell.data 并按队列平台选择 RNA-seq/芯片 spillover）
run_xcell <- function(expr, rnaseq = TRUE) {
  if (!requireNamespace("xCell", quietly = TRUE))
    stop("缺少 xCell 包——本脚本不自动安装，请运行\n",
         "  remotes::install_github(\"dviraran/xCell\")   # 或 devtools::install_github(\"dviraran/xCell\")\n",
         "后重跑")
  e <- new.env()
  data("xCell.data", package = "xCell", envir = e)
  xcd <- e$xCell.data
  if (is.null(xcd) || is.null(xcd$signatures) || is.null(xcd$genes))
    stop("xCell 包内 xCell.data 加载失败或结构不完整——请重装 xCell 后重跑")
  spill <- if (isTRUE(rnaseq)) xcd$spill else xcd$spill.array
  if (is.null(spill))
    stop("xCell.data 缺少 ", ifelse(isTRUE(rnaseq), "spill", "spill.array"),
         "——请检查 xCell 安装完整性")
  scores <- xCell::xCellAnalysis(
    expr, signatures = xcd$signatures, genes = xcd$genes, spill = spill,
    rnaseq = isTRUE(rnaseq), parallel.sz = CONFIG$XCELL_PARALLEL_SZ,
    parallel.type = "SOCK"
  )
  if (!is.matrix(scores) && !is.data.frame(scores))
    stop("xCellAnalysis 未返回得分矩阵（可能提示 not enough genes）——请检查输入基因名/覆盖度")
  scores <- as.matrix(scores)
  sig_cover <- tryCatch({
    sig_lists <- lapply(xcd$signatures, function(s)
      if (is.character(s)) s else s@genes)   # xCell signatures 实际为命名字符向量列表
    sig_genes <- unique(unlist(sig_lists))
    sum(sig_genes %in% rownames(expr)) / length(sig_genes)
  }, error = function(err) NA_real_)
  list(scores = scores, sig_cover = sig_cover, n_input_genes = nrow(expr))
}

## ssGSEA：GSVA 新旧 API 双分支（>=1.50 新 API，失败 fallback 旧 API）
run_ssgsea <- function(expr, gene_sets) {
  m <- as.matrix(expr)
  if (packageVersion("GSVA") >= "1.50") {
    tryCatch({
      par <- GSVA::ssgseaParam(m, gene_sets)
      GSVA::gsva(par, verbose = FALSE)
    }, error = function(e) {
      log_msg("    ⚠ GSVA 新 API（ssgseaParam）失败：", conditionMessage(e),
              "——fallback 旧 API")
      GSVA::gsva(m, gene_sets, method = "ssgsea", verbose = FALSE)
    })
  } else {
    GSVA::gsva(m, gene_sets, method = "ssgsea", verbose = FALSE)
  }
}

## ---- CIBERSORT：LM22 加载（fix3 官方补充材料渠道）----------------------------
## 优先 LM22.txt（tab 分隔，首列基因名，22 列细胞类型）；备选 LM22.xlsx
## （Newman 2015 PMC4739640 Supplementary Table 1 原样放项目根，readxl 自动
## 转换）。结构校验：400–700 基因 × 22 细胞类型 + 典型细胞名命中 ≥10/14。
load_lm22 <- function() {
  f_txt  <- CONFIG$CIBERSORT_LM22_TXT
  f_xlsx <- CONFIG$CIBERSORT_LM22_XLSX
  sig <- NULL; src <- NA_character_
  if (file.exists(f_txt)) {
    sig <- utils::read.delim(f_txt, row.names = 1, check.names = FALSE,
                             stringsAsFactors = FALSE, fileEncoding = "UTF-8")
    src <- f_txt
  } else if (file.exists(f_xlsx)) {
    if (!requireNamespace("readxl", quietly = TRUE))
      stop("检测到 LM22.xlsx 但缺 readxl 包——install.packages(\"readxl\") 后重跑")
    x <- as.data.frame(readxl::read_excel(f_xlsx))
    rownames(x) <- make.unique(trimws(as.character(x[[1]])))
    sig <- x[, -1, drop = FALSE]
    src <- paste0(f_xlsx, "（readxl 自动转换）")
  } else return(NULL)
  colnames(sig) <- trimws(gsub('"', '', colnames(sig)))
  sig <- as.data.frame(lapply(sig, function(v) as.numeric(as.character(v))),
                       row.names = rownames(sig), check.names = FALSE)
  ## 剔除说明行/空行（数值转换后整行 NA 或基因名为空）
  keep <- rowSums(is.na(sig)) < ncol(sig) * 0.5 &
          !is.na(rownames(sig)) & rownames(sig) != "" & rownames(sig) != "NA"
  sig <- sig[keep, , drop = FALSE]
  ok_dim <- nrow(sig) >= 400 && nrow(sig) <= 700 && ncol(sig) == 22
  typical <- c("B cells naive", "B cells memory", "Plasma cells", "T cells CD8",
               "T cells CD4 naive", "NK cells resting", "Monocytes",
               "Macrophages M0", "Macrophages M1", "Macrophages M2",
               "Dendritic cells resting", "Mast cells resting",
               "Eosinophils", "Neutrophils")
  n_hit <- sum(typical %in% colnames(sig))
  if (!ok_dim || n_hit < 10)
    stop("LM22 结构校验失败：", nrow(sig), " 基因 × ", ncol(sig),
         " 列，典型细胞名命中 ", n_hit, "/14（期望约 547×22）——请核对文件是否为 ",
         "Newman 2015 Supplementary Table 1（PMC4739640 补充材料）")
  log_msg("LM22 加载成功（", src, "）：", nrow(sig), " 基因 × ", ncol(sig),
          " 细胞类型；典型细胞名命中 ", n_hit, "/14")
  as.matrix(sig)
}

## ---- CIBERSORT ν-SVR 去卷积（dicepro 开源实现；Newman et al. 2015 算法）----
## QN 口径与官方惯例一致：芯片 QN=TRUE；RNA-seq QN=FALSE。
## perm>0 返回样本级全局 P 值（置换经验零分布）。返回 22 细胞 × 样本得分矩阵
## + 样本级 P/Correlation/RMSE + LM22 基因覆盖率。
run_cibersort_dicepro <- function(expr, sig, rnaseq, perm = 100) {
  for (pkg in c("dicepro", "e1071", "preprocessCore"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("缺包 ", pkg, "——安装：install.packages(\"dicepro\")；",
           "BiocManager::install(\"preprocessCore\")")
  m <- as.matrix(expr)
  common <- intersect(rownames(m), rownames(sig))
  if (length(common) < 100)
    stop("LM22 与表达矩阵共同基因不足（", length(common),
         " < 100）——检查表达矩阵基因名口径（需 gene symbol）")
  core <- utils::getFromNamespace(".run_cibersort_core", "dicepro")
  raw <- core(sig_matrix = sig[common, , drop = FALSE],
              mixture   = m[common, , drop = FALSE],
              perm = perm, QN = !rnaseq,
              absolute = FALSE, abs_method = "sig.score")
  raw <- as.data.frame(raw)
  ct <- colnames(sig)
  list(scores = t(as.matrix(raw[, ct, drop = FALSE])),
       pval = unname(raw$P.value), cor = unname(raw$Correlation),
       rmse = unname(raw$RMSE), coverage = length(common) / nrow(sig))
}

## donor 去重（固定规则）：每 donor 优先保留对照样本，无对照的 donor 才保留病例；
## 返回保留索引与实际比较用 case/ctrl 数，并日志留痕
dedup_donor <- function(pheno, ctrl_lab, gse, method) {
  ord <- order(ifelse(pheno$condition == ctrl_lab, 0L, 1L))
  keep_idx <- sort(ord[!duplicated(pheno$donor_id[ord])])
  n_ctrl_used <- sum(pheno$condition[keep_idx] == ctrl_lab)
  n_case_used <- length(keep_idx) - n_ctrl_used
  log_msg("  ", gse, "/", method, ": donor 去重（优先保留对照，无对照的 donor 保留病例）后 ",
          length(keep_idx), " 样本（case=", n_case_used, " / ctrl=", n_ctrl_used, "）")
  list(keep_idx = keep_idx, n_case_used = n_case_used, n_ctrl_used = n_ctrl_used)
}

## 队列内组间比较：limma moderated t（design=~condition）；
## donor 结构（GSE51588）优先 duplicateCorrelation，失败降级 donor 去重
limma_compare <- function(score_mat, pheno, case_lab, ctrl_lab, gse, tier, method,
                          expect_donor = FALSE) {
  n_case <- sum(pheno$condition == case_lab)
  n_ctrl <- sum(pheno$condition == ctrl_lab)
  n_case_used <- n_case; n_ctrl_used <- n_ctrl   # donor 去重后更新为实际比较用样本数
  low_power <- n_case < CONFIG$LOW_POWER_MIN_CASE
  cond <- factor(pheno$condition, levels = c(ctrl_lab, case_lab))
  design <- model.matrix(~ cond)
  paired_method <- "limma_independent"
  fit <- NULL
  if (expect_donor && !"donor_id" %in% names(pheno))
    stop(gse, "/", method, ": 队列配置 expect_donor=TRUE 但 pheno 无 donor_id 列",
         "——donor 结构口径破坏，请回查模块 A 输出，禁止静默退化为独立样本分析")
  if ("donor_id" %in% names(pheno) && anyDuplicated(pheno$donor_id)) {
    if (anyNA(pheno$donor_id)) stop(gse, ": donor_id 含 NA，请人工核对")
    fit <- tryCatch({
      dup <- limma::duplicateCorrelation(score_mat, design, block = pheno$donor_id)
      limma::lmFit(score_mat, design, block = pheno$donor_id,
                   correlation = dup$consensus.correlation)
    }, error = function(e) {
      log_msg("  ⚠ ", gse, "/", method, ": duplicateCorrelation 失败（",
              conditionMessage(e), "）——降级 donor 去重 + 普通 limma")
      NULL
    })
    if (!is.null(fit)) {
      paired_method <- "duplicateCorrelation"
    } else {
      dd <- dedup_donor(pheno, ctrl_lab, gse, method)
      score_mat <- score_mat[, dd$keep_idx, drop = FALSE]
      pheno <- pheno[dd$keep_idx, ]
      n_case_used <- dd$n_case_used; n_ctrl_used <- dd$n_ctrl_used
      cond <- factor(pheno$condition, levels = c(ctrl_lab, case_lab))
      design <- model.matrix(~ cond)
      fit <- limma::lmFit(score_mat, design)
      paired_method <- "donor_dedup_fallback"
    }
  } else {
    fit <- limma::lmFit(score_mat, design)
  }
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
  cts <- rownames(tt)
  is_case <- pheno$condition == case_lab
  mean_case <- rowMeans(score_mat[cts, is_case, drop = FALSE])
  mean_ctrl <- rowMeans(score_mat[cts, !is_case, drop = FALSE])
  delta <- mean_case - mean_ctrl
  data.frame(
    gse = gse, tier = tier, method = method, cell_type = cts,
    n_case = n_case, n_ctrl = n_ctrl,
    n_case_used = n_case_used, n_ctrl_used = n_ctrl_used,   # 实际进入比较的样本数
    mean_case = round(mean_case, 6), mean_ctrl = round(mean_ctrl, 6),
    delta = round(delta, 6), t = tt$t, P = tt$P.Value,
    adjP = p.adjust(tt$P.Value, method = "BH"),   # BH 按队列×方法内部校正
    direction = ifelse(delta > 0, "up_in_case", "down_in_case"),
    low_power = low_power, paired_method = paired_method,
    stringsAsFactors = FALSE
  )
}

## 队列内 Spearman 相关（10 标志物 × 各细胞类型）；donor 结构先去重防伪重复
cor_biomarker_celltype <- function(score_mat, pheno, expr, biomarkers,
                                   gse, tier, method, case_lab, ctrl_lab,
                                   expect_donor = FALSE) {
  if (expect_donor && !"donor_id" %in% names(pheno))
    stop(gse, "/", method, ": 队列配置 expect_donor=TRUE 但 pheno 无 donor_id 列",
         "——donor 结构口径破坏，请回查模块 A 输出，禁止静默退化为独立样本分析")
  if ("donor_id" %in% names(pheno) && anyDuplicated(pheno$donor_id)) {
    if (anyNA(pheno$donor_id)) stop(gse, ": donor_id 含 NA，请人工核对")
    dd <- dedup_donor(pheno, ctrl_lab, gse, method)   # 相关分析用 donor 去重样本防伪重复
    score_mat <- score_mat[, dd$keep_idx, drop = FALSE]
    expr <- expr[, dd$keep_idx, drop = FALSE]
  }
  genes <- intersect(biomarkers, rownames(expr))
  if (length(genes) < length(biomarkers))
    log_msg("  ⚠ ", gse, ": ", length(biomarkers) - length(genes),
            " 个标志物不在表达矩阵（", paste(setdiff(biomarkers, genes), collapse = ","), "）")
  if (length(genes) == 0) {
    log_msg("  ⚠ ", gse, "/", method, ": 无任何标志物落在表达矩阵——相关分析跳过（return NULL）")
    return(NULL)
  }
  rows <- do.call(rbind, lapply(genes, function(g) {
    x <- as.numeric(expr[g, ])
    do.call(rbind, lapply(rownames(score_mat), function(ct) {
      y <- as.numeric(score_mat[ct, ])
      n <- sum(complete.cases(x, y))
      rho <- NA_real_; p <- NA_real_
      if (n >= CONFIG$COR_MIN_N && sd(x, na.rm = TRUE) > 0 && sd(y, na.rm = TRUE) > 0) {
        ct_res <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
        rho <- unname(ct_res$estimate); p <- ct_res$p.value
      }   # n<8 或零方差 → 记 NA（不做无信息推断）
      data.frame(gse = gse, tier = tier, method = method, gene = g, cell_type = ct,
                 n = n, rho = round(rho, 4), P = p, stringsAsFactors = FALSE)
    }))
  }))
  rows$adjP <- p.adjust(rows$P, method = "BH")   # BH 按队列×方法（NA 保留 NA）
  rows
}

# ----------------------------- CIBERSORT：LM22 实载（fix3） -------------------
## 渠道：Newman 2015 PMC4739640 Supplementary Table 1（官方补充材料），
## 算法为 dicepro 开源 ν-SVR——不使用 Stanford 分发的 CIBERSORT.R/官方 LM22.txt
lm22_sig <- tryCatch(load_lm22(), error = function(e) {
  log_msg("⚠ LM22 加载失败：", conditionMessage(e),
          "——CIBERSORT 跳过（xCell + ssGSEA 不受影响）")
  NULL
})
cib_status <- data.frame(
  check_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  lm22_txt = CONFIG$CIBERSORT_LM22_TXT, lm22_xlsx = CONFIG$CIBERSORT_LM22_XLSX,
  lm22_loaded = !is.null(lm22_sig),
  lm22_dim = if (!is.null(lm22_sig))
    paste0(nrow(lm22_sig), "x", ncol(lm22_sig)) else NA_character_,
  algorithm = "dicepro v1.0.2 内置 CIBERSORT ν-SVR 开源实现（Newman 2015）",
  perm = CONFIG$CIBERSORT_PERM,
  status = if (!is.null(lm22_sig))
    "ready（LM22 已载，逐队列实跑 CIBERSORT）"
  else
    "skipped_no_lm22（未检测到 LM22.txt/LM22.xlsx；获取指引见脚本头部 fix3 块）",
  note = "许可合规：LM22 取自原文 Supplementary Table 1（PMC4739640），算法为 CRAN 开源实现；未使用 Stanford 分发文件",
  stringsAsFactors = FALSE
)
write.csv(cib_status, file.path(CONFIG$OUT_DIR, "moduleF_cibersort_status.csv"),
          row.names = FALSE)
log_msg("CIBERSORT 状态：", cib_status$status)

# ----------------------------- 逐队列主循环（单队列降级不拖垮全局）-----------
scores_long_all <- list()   # 得分长表（xCell/ssGSEA/CIBERSORT 分开存）
cmp_all         <- list()   # 组间比较
cor_all         <- list()   # 相关分析
run_rows        <- list()   # 运行摘要
cibersort_qc_all <- list()  # CIBERSORT 样本级 QC（P/Cor/RMSE）

for (cf in CONFIG$COHORTS) {
  gse <- cf$gse; tier <- cf$tier
  log_msg("—— ", gse, "（", tier, "：", CONFIG$TIER_NOTES[[tier]], "）——")
  run_row <- data.frame(
    gse = gse, tier = tier, tier_note = CONFIG$TIER_NOTES[[tier]],
    n_case = NA_integer_, n_ctrl = NA_integer_, low_power = NA,
    xcell_status = "not_run", xcell_rnaseq = isTRUE(cf$rnaseq),
    xcell_sig_coverage = NA_real_,
    xcell_n_input_genes = NA_integer_, xcell_n_celltypes = NA_integer_,
    xcell_paired_method = NA_character_,
    ssgsea_status = "not_run", ssgsea_mean_gene_coverage = NA_real_,
    ssgsea_paired_method = NA_character_,
    cibersort_status = "not_run", cibersort_gene_coverage = NA_real_,
    cibersort_n_p_pass = NA_integer_, cibersort_paired_method = NA_character_,
    error = NA_character_, stringsAsFactors = FALSE
  )
  tryCatch({
    co <- assemble_cohort(cf)
    run_row$n_case <- sum(co$pheno$condition == co$case_lab)
    run_row$n_ctrl <- sum(co$pheno$condition == co$ctrl_lab)
    run_row$low_power <- run_row$n_case < CONFIG$LOW_POWER_MIN_CASE
    if (run_row$low_power)
      log_msg("  ⚠ ", gse, ": n_case=", run_row$n_case, " < ",
              CONFIG$LOW_POWER_MIN_CASE, "，标记 low_power（结果仅作参考）")
    log_msg("  样本：case=", run_row$n_case, " / control=", run_row$n_ctrl,
            "；基因数=", nrow(co$expr))

    ## ---- xCell ----
    xc <- tryCatch(run_xcell(co$expr, rnaseq = isTRUE(cf$rnaseq)), error = function(e) {
      log_msg("  ⚠ ", gse, " xCell 失败：", conditionMessage(e)); NULL
    })
    if (!is.null(xc)) {
      run_row$xcell_status <- "success"
      run_row$xcell_sig_coverage <- round(xc$sig_cover, 4)
      run_row$xcell_n_input_genes <- xc$n_input_genes
      run_row$xcell_n_celltypes <- nrow(xc$scores)
      log_msg("  xCell 完成：", nrow(xc$scores), " 细胞类型；签名覆盖率=",
              ifelse(is.na(xc$sig_cover), "NA（包内对象不可取）",
                     sprintf("%.1f%%", 100 * xc$sig_cover)),
              "；输入基因=", xc$n_input_genes)
      scores_long_all[[paste(gse, "xCell")]] <- data.frame(
        gse = gse, tier = tier, method = "xCell",
        sample = rep(colnames(xc$scores), each = nrow(xc$scores)),
        condition = rep(co$pheno$condition[match(colnames(xc$scores), co$pheno$sample)],
                        each = nrow(xc$scores)),
        cell_type = rep(rownames(xc$scores), times = ncol(xc$scores)),
        immunity = NA_character_,   # 与 ssGSEA 长表列对齐（xCell 无 immunity 口径），防 safe_rbind 列数不一致
        score = as.numeric(xc$scores), stringsAsFactors = FALSE)
      cmp <- tryCatch(
        limma_compare(xc$scores, co$pheno, co$case_lab, co$ctrl_lab, gse, tier, "xCell",
                      expect_donor = isTRUE(cf$expect_donor)),
        error = function(e) { log_msg("  ⚠ ", gse, " xCell 组间比较失败：",
                                      conditionMessage(e))
          if (grepl("expect_donor", conditionMessage(e)))
            run_row$error <<- conditionMessage(e)   # donor 口径破坏在 run_summary 留痕
          NULL })
      if (!is.null(cmp)) {
        cmp_all[[paste(gse, "xCell")]] <- cmp
        run_row$xcell_paired_method <- unique(cmp$paired_method)
      }
      ## 相关分析与组间比较解耦：cmp 失败不阻止 cor 执行
      cor_all[[paste(gse, "xCell")]] <- tryCatch(
        cor_biomarker_celltype(xc$scores, co$pheno, co$expr, biomarkers, gse, tier, "xCell",
                               co$case_lab, co$ctrl_lab,
                               expect_donor = isTRUE(cf$expect_donor)),
        error = function(e) { log_msg("  ⚠ ", gse, " xCell 相关分析失败：",
                                      conditionMessage(e))
          if (grepl("expect_donor", conditionMessage(e)))
            run_row$error <<- conditionMessage(e)
          NULL })
    } else run_row$xcell_status <- "failed（见日志）"

    ## ---- ssGSEA 28 细胞 ----
    ss <- tryCatch(run_ssgsea(co$expr, gene_sets), error = function(e) {
      log_msg("  ⚠ ", gse, " ssGSEA 失败：", conditionMessage(e)); NULL
    })
    if (!is.null(ss)) {
      ss <- as.matrix(ss)
      run_row$ssgsea_status <- "success"
      run_row$ssgsea_mean_gene_coverage <- round(mean(sapply(
        gene_sets[rownames(ss)], function(gs) mean(gs %in% rownames(co$expr)))), 4)
      log_msg("  ssGSEA 完成：", nrow(ss), " 细胞类型；平均 metagene 覆盖率=",
              sprintf("%.1f%%", 100 * run_row$ssgsea_mean_gene_coverage))
      scores_long_all[[paste(gse, "ssGSEA")]] <- data.frame(
        gse = gse, tier = tier, method = "ssGSEA",
        sample = rep(colnames(ss), each = nrow(ss)),
        condition = rep(co$pheno$condition[match(colnames(ss), co$pheno$sample)],
                        each = nrow(ss)),
        cell_type = rep(rownames(ss), times = ncol(ss)),
        immunity = rep(unname(immunity_map[rownames(ss)]), times = ncol(ss)),
        score = as.numeric(ss), stringsAsFactors = FALSE)
      cmp <- tryCatch(
        limma_compare(ss, co$pheno, co$case_lab, co$ctrl_lab, gse, tier, "ssGSEA",
                      expect_donor = isTRUE(cf$expect_donor)),
        error = function(e) { log_msg("  ⚠ ", gse, " ssGSEA 组间比较失败：",
                                      conditionMessage(e))
          if (grepl("expect_donor", conditionMessage(e)))
            run_row$error <<- conditionMessage(e)   # donor 口径破坏在 run_summary 留痕
          NULL })
      if (!is.null(cmp)) {
        cmp_all[[paste(gse, "ssGSEA")]] <- cmp
        run_row$ssgsea_paired_method <- unique(cmp$paired_method)
      }
      ## 相关分析与组间比较解耦：cmp 失败不阻止 cor 执行
      cor_all[[paste(gse, "ssGSEA")]] <- tryCatch(
        cor_biomarker_celltype(ss, co$pheno, co$expr, biomarkers, gse, tier, "ssGSEA",
                               co$case_lab, co$ctrl_lab,
                               expect_donor = isTRUE(cf$expect_donor)),
        error = function(e) { log_msg("  ⚠ ", gse, " ssGSEA 相关分析失败：",
                                      conditionMessage(e))
          if (grepl("expect_donor", conditionMessage(e)))
            run_row$error <<- conditionMessage(e)
          NULL })
    } else run_row$ssgsea_status <- "failed（见日志）"

    ## ---- CIBERSORT 22 细胞（dicepro ν-SVR + LM22；fix3 实跑）----
    if (!is.null(lm22_sig)) {
      cb <- tryCatch(
        run_cibersort_dicepro(co$expr, lm22_sig,
                              rnaseq = isTRUE(cf$rnaseq),
                              perm = CONFIG$CIBERSORT_PERM),
        error = function(e) {
          log_msg("  ⚠ ", gse, " CIBERSORT 失败：", conditionMessage(e)); NULL
        })
      if (!is.null(cb)) {
        run_row$cibersort_status <- "success"
        run_row$cibersort_gene_coverage <- round(cb$coverage, 4)
        run_row$cibersort_n_p_pass <- sum(cb$pval < CONFIG$CIBERSORT_P_CUTOFF,
                                          na.rm = TRUE)
        log_msg("  CIBERSORT 完成：", nrow(cb$scores), " 细胞类型；LM22 覆盖=",
                sprintf("%.1f%%", 100 * cb$coverage), "；P<",
                CONFIG$CIBERSORT_P_CUTOFF, " 样本 ",
                run_row$cibersort_n_p_pass, "/", ncol(cb$scores))
        scores_long_all[[paste(gse, "CIBERSORT")]] <- data.frame(
          gse = gse, tier = tier, method = "CIBERSORT",
          sample = rep(colnames(cb$scores), each = nrow(cb$scores)),
          condition = rep(co$pheno$condition[match(colnames(cb$scores), co$pheno$sample)],
                          each = nrow(cb$scores)),
          cell_type = rep(rownames(cb$scores), times = ncol(cb$scores)),
          immunity = NA_character_,   # 与 xCell 同口径（LM22 无 immunity 分组列）
          score = as.numeric(cb$scores), stringsAsFactors = FALSE)
        cibersort_qc_all[[gse]] <- data.frame(
          gse = gse, tier = tier,
          sample = colnames(cb$scores),
          condition = co$pheno$condition[match(colnames(cb$scores), co$pheno$sample)],
          cibersort_p = cb$pval, cibersort_cor = cb$cor,
          cibersort_rmse = cb$rmse,
          p_pass = cb$pval < CONFIG$CIBERSORT_P_CUTOFF,
          stringsAsFactors = FALSE)
        cmp <- tryCatch(
          limma_compare(cb$scores, co$pheno, co$case_lab, co$ctrl_lab, gse, tier,
                        "CIBERSORT", expect_donor = isTRUE(cf$expect_donor)),
          error = function(e) { log_msg("  ⚠ ", gse, " CIBERSORT 组间比较失败：",
                                        conditionMessage(e))
            if (grepl("expect_donor", conditionMessage(e)))
              run_row$error <<- conditionMessage(e)
            NULL })
        if (!is.null(cmp)) {
          cmp_all[[paste(gse, "CIBERSORT")]] <- cmp
          run_row$cibersort_paired_method <- unique(cmp$paired_method)
        }
        cor_all[[paste(gse, "CIBERSORT")]] <- tryCatch(
          cor_biomarker_celltype(cb$scores, co$pheno, co$expr, biomarkers, gse,
                                 tier, "CIBERSORT", co$case_lab, co$ctrl_lab,
                                 expect_donor = isTRUE(cf$expect_donor)),
          error = function(e) { log_msg("  ⚠ ", gse, " CIBERSORT 相关分析失败：",
                                        conditionMessage(e))
            if (grepl("expect_donor", conditionMessage(e)))
              run_row$error <<- conditionMessage(e)
            NULL })
      } else run_row$cibersort_status <- "failed（见日志）"
    } else run_row$cibersort_status <- "skipped_no_lm22"
  }, error = function(e) {
    log_msg("  ⚠⚠ ", gse, " 整队列降级跳过：", conditionMessage(e))
    run_row$error <<- conditionMessage(e)
  })
  run_rows[[gse]] <- run_row
}

# ----------------------------- 空表守卫 + 汇总写出 ---------------------------
safe_rbind <- function(lst) {
  lst <- lst[!sapply(lst, is.null)]
  if (length(lst) == 0) return(NULL)
  do.call(rbind, lst)
}
scores_long <- safe_rbind(scores_long_all)
cmp_tab     <- safe_rbind(cmp_all)
cor_tab     <- safe_rbind(cor_all)
run_summary <- safe_rbind(run_rows)
if (is.null(run_summary)) stop("全部队列失败——请回查日志与输入")
write.csv(run_summary, file.path(CONFIG$OUT_DIR, "moduleF_run_summary.csv"),
          row.names = FALSE)
log_msg("运行摘要已写出：", nrow(run_summary), " 队列；成功 xCell=",
        sum(run_summary$xcell_status == "success"), "，成功 ssGSEA=",
        sum(run_summary$ssgsea_status == "success"), "，成功 CIBERSORT=",
        sum(run_summary$cibersort_status == "success"))

## CIBERSORT 样本级 QC（P 值/相关/RMSE；fix3）
cib_qc <- safe_rbind(cibersort_qc_all)
if (!is.null(cib_qc)) {
  write.csv(cib_qc, file.path(CONFIG$OUT_DIR, "moduleF_cibersort_sample_qc.csv"),
            row.names = FALSE)
  log_msg("CIBERSORT 样本 QC：", nrow(cib_qc), " 样本；P<",
          CONFIG$CIBERSORT_P_CUTOFF, " 占比 ",
          sprintf("%.1f%%", 100 * mean(cib_qc$p_pass, na.rm = TRUE)))
}

## 得分长表（幂等覆盖）
if (!is.null(scores_long)) {
  write.csv(scores_long[scores_long$method == "xCell",
            c("gse", "tier", "sample", "condition", "cell_type", "score")],
            file.path(CONFIG$OUT_DIR, "moduleF_xcell_scores_long.csv"), row.names = FALSE)
  write.csv(scores_long[scores_long$method == "ssGSEA",
            c("gse", "tier", "sample", "condition", "cell_type", "immunity", "score")],
            file.path(CONFIG$OUT_DIR, "moduleF_ssgsea_scores_long.csv"), row.names = FALSE)
  write.csv(scores_long[scores_long$method == "CIBERSORT",
            c("gse", "tier", "sample", "condition", "cell_type", "score")],
            file.path(CONFIG$OUT_DIR, "moduleF_cibersort_scores_long.csv"), row.names = FALSE)
} else log_msg("⚠ 无任何队列产出得分——得分长表跳过（空表守卫）")

## 组间比较（按方法拆分写出）
if (!is.null(cmp_tab)) {
  write.csv(cmp_tab[cmp_tab$method == "xCell", ],
            file.path(CONFIG$OUT_DIR, "moduleF_xcell_group_compare.csv"), row.names = FALSE)
  write.csv(cmp_tab[cmp_tab$method == "ssGSEA", ],
            file.path(CONFIG$OUT_DIR, "moduleF_ssgsea_group_compare.csv"), row.names = FALSE)
  write.csv(cmp_tab[cmp_tab$method == "CIBERSORT", ],
            file.path(CONFIG$OUT_DIR, "moduleF_cibersort_group_compare.csv"), row.names = FALSE)
  log_msg("组间比较：", nrow(cmp_tab), " 行（显著 adjP<0.05：",
          sum(cmp_tab$adjP < 0.05, na.rm = TRUE), " 行）")
} else log_msg("⚠ 无任何组间比较结果——相关输出跳过（空表守卫）")

# ----------------------------- 跨队列方向汇总（禁止 pooled p）----------------
## 只在 tier × method × cell_type 层面汇总：n_tested / n_sig / n_same_dir /
## verdict（同方向显著队列数 >= 2 且反方向显著数 = 0 → replicated；
## 双向显著 → conflicting；否则 exploratory / none）
if (!is.null(cmp_tab)) {
  dir_sum <- do.call(rbind, lapply(
    split(cmp_tab, list(cmp_tab$tier, cmp_tab$method, cmp_tab$cell_type), drop = TRUE),
    function(df) {
      sig <- df[!is.na(df$adjP) & df$adjP < 0.05, ]
      n_up <- sum(sig$direction == "up_in_case")
      n_dn <- sum(sig$direction == "down_in_case")
      n_same <- max(n_up, n_dn)
      verdict <- if (n_same >= 2 && min(n_up, n_dn) == 0)
        paste0("replicated_", ifelse(n_up >= n_dn, "up_in_case", "down_in_case"))
      else if (n_up >= 1 && n_dn >= 1) "conflicting"
      else if (nrow(sig) >= 1) "exploratory" else "none"
      data.frame(
        tier = df$tier[1], method = df$method[1], cell_type = df$cell_type[1],
        n_tested = nrow(df), n_sig = nrow(sig),
        n_sig_up = n_up, n_sig_down = n_dn, n_same_dir = n_same,
        verdict = verdict,
        note = "仅方向汇总，禁止跨队列 pooled p；BH 已在各队列内部完成",
        stringsAsFactors = FALSE)
    }))
  write.csv(dir_sum, file.path(CONFIG$OUT_DIR,
            "moduleF_crosscohort_direction_summary.csv"), row.names = FALSE)
  log_msg("跨队列方向汇总：", nrow(dir_sum), " 行；replicated=",
          sum(grepl("^replicated", dir_sum$verdict)), "，conflicting=",
          sum(dir_sum$verdict == "conflicting"), " 个 tier×method×cell_type")
}

## 相关分析长表
if (!is.null(cor_tab)) {
  write.csv(cor_tab, file.path(CONFIG$OUT_DIR,
            "moduleF_biomarker_celltype_spearman.csv"), row.names = FALSE)
  log_msg("相关分析：", nrow(cor_tab), " 行（可算 ", sum(!is.na(cor_tab$rho)),
          " 行，显著 adjP<0.05：", sum(cor_tab$adjP < 0.05, na.rm = TRUE), " 行）")
} else log_msg("⚠ 无相关分析结果——输出跳过（空表守卫）")

# ----------------------------- 图形（tryCatch 不阻断 CSV）-------------------
## 1) 每方法：队列 × 细胞类型 signed heatmap（sign(delta) * -log10(adjP)，NA 容忍）
plot_signed_heat <- function(cmp_df, method) {
  df <- cmp_df[cmp_df$method == method, ]
  if (nrow(df) == 0) return(invisible(NULL))
  df$signed <- sign(df$delta) * -log10(pmax(df$adjP, 1e-300))
  cells <- sort(unique(df$cell_type)); gses <- unique(df$gse)
  mat <- matrix(NA_real_, nrow = length(cells), ncol = length(gses),
                dimnames = list(cells, gses))
  for (i in seq_len(nrow(df)))
    mat[df$cell_type[i], df$gse[i]] <- df$signed[i]
  main_ttl <- paste0("Module F ", method, "：signed -log10(adjP)（正=case 上调）")
  cl_rows <- !anyNA(mat) && nrow(mat) > 1   # signed 含 NA（或单行）时关闭行聚类，防 dist/hclust 崩
  pheatmap(mat, cluster_rows = cl_rows, cluster_cols = FALSE, na_col = "grey90",
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           main = main_ttl, fontsize = 7,
           filename = file.path(CONFIG$FIG_DIR,
                                paste0("moduleF_", tolower(method), "_signed_heatmap.pdf")),
           width = 9, height = max(5, 0.22 * nrow(mat) + 2))
  pheatmap(mat, cluster_rows = cl_rows, cluster_cols = FALSE, na_col = "grey90",
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           main = main_ttl, fontsize = 7,
           filename = file.path(CONFIG$FIG_DIR,
                                paste0("moduleF_", tolower(method), "_signed_heatmap.png")),
           width = 9, height = max(5, 0.22 * nrow(mat) + 2))
}

## 2) 每层 top 差异细胞 boxplot（按该层最小 adjP 取 top N）
plot_tier_boxplots <- function(cmp_df, scores_long, method, tier) {
  df <- cmp_df[cmp_df$method == method & cmp_df$tier == tier, ]
  if (nrow(df) == 0) return(invisible(NULL))
  ct_minp <- tapply(df$adjP, df$cell_type,
                    function(p) if (all(is.na(p))) Inf else min(p, na.rm = TRUE))
  top_ct <- names(sort(ct_minp))[seq_len(min(CONFIG$TOP_N_BOXPLOT, length(ct_minp)))]
  sd_df <- scores_long[scores_long$method == method & scores_long$tier == tier &
                         scores_long$cell_type %in% top_ct, ]
  if (nrow(sd_df) == 0) return(invisible(NULL))
  sd_df$cell_type <- factor(sd_df$cell_type, levels = top_ct)
  ttl <- paste0("Module F ", method, " top 差异细胞｜", tier,
                "（", CONFIG$TIER_NOTES[[tier]], "）")
  p <- ggplot(sd_df, aes(x = condition, y = score, fill = condition)) +
    geom_boxplot(outlier.size = 0.6, alpha = 0.8) +
    facet_grid(cell_type ~ gse, scales = "free_y") +
    labs(title = ttl, x = NULL, y = "浸润得分") +
    theme_bw(base_size = 9) +
    theme(legend.position = "none", strip.text.y = element_text(size = 7))
  ggsave(file.path(CONFIG$FIG_DIR, paste0("moduleF_", tolower(method), "_boxplot_",
                                          gsub("[^A-Za-z0-9]", "_", tier), ".pdf")),
         p, width = 2.2 * length(unique(sd_df$gse)) + 2,
         height = 1.8 * length(top_ct) + 1.5, limitsize = FALSE)
  ggsave(file.path(CONFIG$FIG_DIR, paste0("moduleF_", tolower(method), "_boxplot_",
                                          gsub("[^A-Za-z0-9]", "_", tier), ".png")),
         p, width = 2.2 * length(unique(sd_df$gse)) + 2,
         height = 1.8 * length(top_ct) + 1.5, dpi = 150, limitsize = FALSE)
}

## 3) 标志物 × 细胞类型相关热图（tier 内跨队列 median rho）
plot_corr_heat <- function(cor_df, method, tier) {
  df <- cor_df[cor_df$method == method & cor_df$tier == tier, ]
  if (nrow(df) == 0) return(invisible(NULL))
  med <- tapply(df$rho, list(df$gene, df$cell_type),
                function(r) if (all(is.na(r))) NA_real_ else median(r, na.rm = TRUE))
  if (is.null(dim(med))) return(invisible(NULL))
  cl_rows <- !anyNA(med) && nrow(med) > 1   # 含 NA 或单行/单列时关闭聚类，防 dist/hclust 报错
  cl_cols <- !anyNA(med) && ncol(med) > 1
  ttl <- paste0("Module F ", method, "｜", tier,
                "：标志物 × 细胞类型 Spearman median rho（", CONFIG$TIER_NOTES[[tier]], "）")
  fn_base <- file.path(CONFIG$FIG_DIR, paste0("moduleF_", tolower(method),
                       "_corr_heatmap_", gsub("[^A-Za-z0-9]", "_", tier)))
  pheatmap(med, cluster_rows = cl_rows, cluster_cols = cl_cols, na_col = "grey90",
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           breaks = seq(-1, 1, length.out = 101), main = ttl, fontsize = 7,
           filename = paste0(fn_base, ".pdf"),
           width = max(8, 0.28 * ncol(med) + 3), height = 2 + 0.35 * nrow(med))
  pheatmap(med, cluster_rows = cl_rows, cluster_cols = cl_cols, na_col = "grey90",
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           breaks = seq(-1, 1, length.out = 101), main = ttl, fontsize = 7,
           filename = paste0(fn_base, ".png"),
           width = max(8, 0.28 * ncol(med) + 3), height = 2 + 0.35 * nrow(med))
}

## 三类图分别 tryCatch：一张图失败不连累其他图
tryCatch({
  if (!is.null(cmp_tab)) {
    for (mth in unique(cmp_tab$method))
      tryCatch(plot_signed_heat(cmp_tab, mth), error = function(e)
        log_msg("⚠ ", mth, " signed heatmap 失败：", conditionMessage(e)))
    if (!is.null(scores_long))
      for (mth in unique(cmp_tab$method))
        for (tr in unique(cmp_tab$tier[cmp_tab$method == mth]))
          tryCatch(plot_tier_boxplots(cmp_tab, scores_long, mth, tr), error = function(e)
            log_msg("⚠ ", mth, "/", tr, " 分层 boxplot 失败：", conditionMessage(e)))
  }
  if (!is.null(cor_tab))
    for (mth in unique(cor_tab$method))
      for (tr in unique(cor_tab$tier[cor_tab$method == mth]))
        tryCatch(plot_corr_heat(cor_tab, mth, tr), error = function(e)
          log_msg("⚠ ", mth, "/", tr, " 相关热图失败：", conditionMessage(e)))
  log_msg("图形输出完成（", length(list.files(CONFIG$FIG_DIR)), " 个文件）")
}, error = function(e) log_msg("⚠ 图形绘制失败（CSV 不受影响）：", conditionMessage(e)))

# ----------------------------- 收尾 -----------------------------------------
cat("\n模块 F 完成。产物：\n",
    "  ", CONFIG$OUT_DIR, "/moduleF_xcell_scores_long.csv / moduleF_ssgsea_scores_long.csv / moduleF_cibersort_scores_long.csv —— 得分长表\n",
    "  ", CONFIG$OUT_DIR, "/moduleF_xcell_group_compare.csv / moduleF_ssgsea_group_compare.csv / moduleF_cibersort_group_compare.csv —— 队列内组间比较（BH 按队列×方法）\n",
    "  ", CONFIG$OUT_DIR, "/moduleF_cibersort_sample_qc.csv —— CIBERSORT 样本级 P/Cor/RMSE 质控\n",
    "  ", CONFIG$OUT_DIR, "/moduleF_crosscohort_direction_summary.csv —— tier×method×cell_type 方向汇总（无 pooled p）\n",
    "  ", CONFIG$OUT_DIR, "/moduleF_biomarker_celltype_spearman.csv —— 10 标志物 × 细胞类型 Spearman 长表\n",
    "  ", CONFIG$OUT_DIR, "/moduleF_cibersort_status.csv —— CIBERSORT（dicepro+LM22）运行状态\n",
    "  ", CONFIG$OUT_DIR, "/moduleF_run_summary.csv —— 逐队列运行/覆盖率/low_power 摘要\n",
    "  ", CONFIG$FIG_DIR, "/ —— signed heatmap / 分层 boxplot / 相关热图（PDF+PNG）\n",
    "注意：统计只在队列内进行；跨队列仅方向汇总；Aging_proxy 层为衰老代理，",
    "不得写作 sarcopenia 证据；Blood_exploratory 层为 EXPLORATORY。\n", sep = "")
log_msg("==== 模块 F 结束 ====")
