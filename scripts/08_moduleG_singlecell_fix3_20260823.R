#!/usr/bin/env Rscript
# ============================================================================
# 08_moduleG_singlecell_fix3_20260823.R
# ----------------------------------------------------------------------------
# 模块 G / G1 脚本：三个 GEO 单细胞/单核数据集的 Seurat 对象构建与 QC。
# 编码: UTF-8（中文 Windows 用户请确认 RStudio 以 UTF-8 打开）
#
# 版本日期: 2026-08-23 fix3 — 实机三跑修复 1 个问题（GSE152805/GSE169454
#   已全绿，仅 GSE167186 报错）：
#   问题(GSE167186): fix2 诊断显示 H5 barcode 非 NULL、无重复，却仍报
#     "No cell overlap between new meta data and Seurat object"。沙箱解剖
#     GSM5098737 H5 确认为标准 10x v2、barcode 干净；查 SeuratObject 5.x
#     源码确认该校验的真实触发条件是「【命名向量】长度!=细胞数且名字与
#     细胞名零重叠」。根因：process_gse167186 中
#     donors <- vapply(gsms, ...) —— vapply 对字符型 X 的结果自动带名
#     （名=GSM 号），donors[i] 为命名长度1向量；create_one_seurat 中
#     seu$donor <- ifelse(is.na(donor), "unknown", donor) 的 ifelse 保留
#     test 的名字属性 -> 命名标量触发零重叠校验。scRNA 两数据集的
#     parse_meta_fn 返回未命名标量故无恙。
#     修复：create_one_seurat 内 5 个标量赋值统一 unname()（chokepoint
#     防御，覆盖 condition/tissue/donor 一切上游命名路径）；对已全绿的
#     两数据集幂等（值不变）。
# fix2 — 实机二跑修复 2 个问题（GSE169454 已全绿）：
#   问题1(GSE152805): Seurat 5.5 的 Read10X 只找 features.tsv.gz，不认 v2 的
#     genes.tsv.gz -> prepare_read10x_dir 中 v2 基因文件目标名也改为
#     features.tsv.gz（2 列格式兼容；旧目录残留 genes.tsv.gz 无害）。
#   问题2(GSE167186): unname 修复后仍报 "No cell overlap between new meta
#     data and Seurat object"，肇事点缩小到 create_one_seurat（H5 barcode
#     列为空/重复）-> create_one_seurat 开头加 colnames NULL/重复守卫
#     （三数据集统一生效）；GSE167186 第二遍循环加逐步诊断日志 + tryCatch
#     （报错消息带 [GSM] step= 前缀精确定位）。
# fix1 — 实机首跑修复 2 个 bug：
#   Bug A(致命): run_dataset 命名实参 gse= 精确匹配其形参导致位置实参顺延、
#     fn 接到字符串 -> "could not find function \"fn\""；改为 run_dataset
#     透传 gse（fn(inv, gse = gse, ...)），调用点删除 gse= 实参，
#     process_gse167186 签名加 gse 形参。
#   Bug B(严重): GSE167186 报 "No cell overlap between new meta data and
#     Seurat object"；add_qc_metrics 改 unname() 按序赋值绕开 Seurat v5
#     [[<- 名称重叠校验；两个处理器在 create_one_seurat 后加 0 细胞守卫。
# v1 — 依据 G0 闭环产物（moduleG_file_inventory.csv,
#   66 行）作为权威文件索引；supp_extracted 为平铺目录，按 GSM 重组为
#   Read10X 标准结构后建对象。
#
# 【红线（必须遵守，下游引用时同样适用）】
#   1. GSE152805 全部为 OA 样本（SY/oLT/MT x 供体 113/116/118），无健康对照，
#      只能做供体内组织间定位，严禁做病例-对照比较。
#   2. GSE167186 分组为 Old vs Young（aging proxy，衰老代理），
#      严禁表述为 sarcopenia（肌少症）证据。
#   3. G1 不做任何组间统计推断（无差异表达、无富集、无比例检验、无整合）；
#      仅建对象、QC、LogNormalize、FindVariableFeatures。
#      ScaleData / PCA / 聚类 / 整合 / 标志物定位均属于 G2。
#
# 【输入（G0 已闭环，用户实机确认）】
#   - results/moduleG/moduleG_file_inventory.csv  （权威文件索引，66 行）
#   - data/GEO/<GSE>/supp_extracted/               （平铺目录）
#       GSE152805: 27 文件 = 9 GSM x (barcodes/genes/matrix)，10x v2；
#                  实际文件名形如 GSM4626763_SY_113.barcodes.tsv.gz（点为分隔，
#                  SY 样本无 OA_ 前缀），正则 (SY|oLT|MT)_(113|116|118) 覆盖全部。
#       GSE169454: 21 文件 = 7 GSM x filtered 三件套（features/barcodes/matrix），
#                  10x v3；condition(normal/oa) 从文件名解析（normal1-3, oa1-4）。
#       GSE167186: 17 个 GSMxxxxxxx_HM<n>_filtered_feature_bc_matrix.h5；
#                  HM 编号为 HM1-13,15,21-23（非连续，必须正则解析不得假设 1-17）。
#   - data/GEO/GSE167186/GSE167186_SimplifiedMetadataSheet.xlsx（Old/Young 分组，
#     列名未知 -> 先打印列名再模糊匹配；失败回退本地 pData，再失败 stop 留痕）
#
# 【输出】
#   results/moduleG/<GSE>_seurat_G1.rds             （过滤后 + LogNormalize + VST2000）
#   results/moduleG/moduleG_G1_qc_summary.csv       （数据集 x 样本：过滤前后指标与阈值）
#   results/moduleG/moduleG_G1_status.csv           （每数据集 OK/FAIL + 备注）
#   results/moduleG/QC/*.png|pdf                    （过滤前/后 violin + scatter）
#   results/moduleG/moduleG_G1_log_<Sys.Date()>.txt
#
# 【环境】Windows + R 4.6.0，工作目录含中文（C:\Users\凉月\Desktop\ENG论文）。
#   本脚本不调用任何 system2/shell/untar，全部为纯 R 函数，规避中文路径编码问题。
# ============================================================================

# ---- 路径设置（假定从项目根目录运行，与 08a 完全一致） -------------------------
proj_root <- getwd()
geo_dir   <- file.path(proj_root, "data", "GEO")
res_dir   <- file.path(proj_root, "results", "moduleG")
qc_dir    <- file.path(res_dir, "QC")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(qc_dir,  showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(res_dir, paste0("moduleG_G1_log_", Sys.Date(), ".txt"))

log_msg <- function(..., .sep = "") {
  txt <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ",
                paste(..., collapse = .sep))
  cat(txt, "\n")
  cat(txt, "\n", file = log_file, append = TRUE)
  invisible(txt)
}

log_msg("=== 模块 G G1 启动：Seurat 对象构建与 QC ===")
log_msg("proj_root = ", proj_root)

# ---- 依赖检查（缺包即 stop 并给出安装命令） ------------------------------------
check_deps <- function() {
  cran_pkgs <- c("Seurat", "Matrix", "readxl", "data.table", "ggplot2")
  missing <- character(0)
  for (p in cran_pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) missing <- c(missing, p)
  }
  if (length(missing) > 0) {
    log_msg("ERROR: 缺少 R 包: ", paste(missing, collapse = ", "))
    stop(paste0(
      "请先安装缺失包再运行 G1：\n  install.packages(c(",
      paste(sprintf('"%s"', missing), collapse = ", "),
      "))\n（Seurat 依赖较多，建议在稳定的 CRAN 镜像下安装。）"),
      call. = FALSE)
  }
  sv <- as.character(utils::packageVersion("Seurat"))
  if (utils::compareVersion(sv, "5.1.0") < 0) {
    # 5.0.x 早期版本多层（multi-layer）对象有已知 bug，要求 >= 5.1.0
    log_msg("ERROR: Seurat 版本 ", sv, " < 5.1.0（本脚本按 Seurat v5 API 编写，",
            "5.0.x 多层对象有已知 bug）。")
    stop("请升级 Seurat：install.packages(\"Seurat\") 后重跑。", call. = FALSE)
  }
  gv <- as.character(utils::packageVersion("ggplot2"))
  if (utils::compareVersion(gv, "3.4.0") < 0) {
    # QC 图使用 linewidth 美学（ggplot2 >= 3.4.0 引入）
    log_msg("ERROR: ggplot2 版本 ", gv, " < 3.4.0（linewidth 美学需要 >= 3.4.0）。")
    stop("请升级 ggplot2：install.packages(\"ggplot2\") 后重跑。", call. = FALSE)
  }
  log_msg("依赖检查通过: Seurat ", sv, " | ggplot2 ", gv,
          " | Matrix ", as.character(utils::packageVersion("Matrix")),
          " | readxl ", as.character(utils::packageVersion("readxl")),
          " | data.table ", as.character(utils::packageVersion("data.table")))
  invisible(TRUE)
}
check_deps()

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

# ---- 常量 ----------------------------------------------------------------------
SNRNA_CELL_THRESHOLD <- 60000   # GSE167186 合并前细胞数上限，超过走逐样本过滤路径
SCRNA_MT_MAX         <- 15      # scRNA percent.mt 上限 (%)
SNRNA_MT_MAX         <- 5       # snRNA percent.mt 上限 (%)（细胞核线粒体预期低）
SCRNA_NFEATURE_HARD  <- 6000    # scRNA nFeature 硬上限兜底
NFEATURE_LOWER       <- 200     # nFeature 下限（scRNA/snRNA 相同）

# ---- 通用工具函数 ---------------------------------------------------------------

# 读取 G0 权威文件清单；缺失/为空即 stop（G1 不以 list.files 重新猜测文件）
load_inventory <- function() {
  f <- file.path(res_dir, "moduleG_file_inventory.csv")
  if (!file.exists(f)) {
    log_msg("ERROR: 未找到 G0 文件清单: ", f)
    stop("缺少 results/moduleG/moduleG_file_inventory.csv，请先成功运行 08a (G0)。",
         call. = FALSE)
  }
  inv <- tryCatch(
    as.data.frame(data.table::fread(f, stringsAsFactors = FALSE,
                                    na.strings = c("", "NA"))),
    error = function(e) {
      log_msg("ERROR: 读取文件清单失败: ", conditionMessage(e))
      stop(paste0("moduleG_file_inventory.csv 读取失败: ", conditionMessage(e)),
           call. = FALSE)
    })
  need_cols <- c("gse", "gsm", "file_name", "condition_hint", "tissue_site",
                 "donor_hint", "extracted_path")
  miss_cols <- setdiff(need_cols, names(inv))
  if (length(miss_cols) > 0) {
    log_msg("ERROR: 文件清单缺列: ", paste(miss_cols, collapse = ", "))
    stop("moduleG_file_inventory.csv 列结构与预期不符（缺 ",
         paste(miss_cols, collapse = ", "), "）。", call. = FALSE)
  }
  log_msg("文件清单载入: ", nrow(inv), " 行 (", f, ")")
  inv
}

# 三件套/H5 角色判定（与 08a 口径一致：先判 .h5，再判关键词）
role_of_file <- function(fname) {
  if (grepl("\\.h5$", fname, ignore.case = TRUE)) return("h5")
  if (grepl("barcodes", fname, ignore.case = TRUE)) return("barcodes")
  if (grepl("genes", fname, ignore.case = TRUE)) return("genes")
  if (grepl("features", fname, ignore.case = TRUE)) return("features")
  if (grepl("matrix", fname, ignore.case = TRUE)) return("matrix")
  "other"
}

# 解析清单行的实际文件路径：优先 extracted_path，失效时回退到
# data/GEO/<GSE>/supp_extracted/<file_name>（相对路径/机器迁移时的保险）
resolve_file_path <- function(gse, file_name, extracted_path) {
  if (!is.na(extracted_path) && nzchar(extracted_path) &&
      file.exists(extracted_path)) {
    return(extracted_path)
  }
  alt <- file.path(geo_dir, gse, "supp_extracted", file_name)
  if (file.exists(alt)) return(alt)
  NA_character_
}

# 按 GSM 重组 Read10X 标准目录（平铺 -> data/GEO/<GSE>/seurat_input/<GSM>/）。
# 复制并改名为 Read10X 标准名：barcodes.tsv.gz / (genes|features).tsv.gz /
# matrix.mtx.gz。幂等：目标已存在且字节数与源一致则跳过复制；
# 复制后做字节数校验，失败即报错（不静默）。
prepare_read10x_dir <- function(gse, gsm, files_df, format = c("v2", "v3")) {
  format <- match.arg(format)
  gene_role <- if (format == "v2") "genes" else "features"   # 源文件角色判定用
  # 键名随 format 联动（修复：v3 时键名必须为 "features"，否则 std_names["features"]
  # 返回 NA -> file.copy(src, NA) 必崩）
  std_names <- c(barcodes = "barcodes.tsv.gz", matrix = "matrix.mtx.gz")
  # fix2：Seurat 5.5 的 Read10X 只找 features.tsv.gz，不认 v2 的 genes.tsv.gz
  # （报错 "Gene name or features file missing. Expecting features.tsv.gz"）。
  # v2 的 2 列 genes 文件复制后统一改名为 features.tsv.gz（2 列格式 Read10X 兼容）。
  # 幂等不受影响：旧 seurat_input 目录残留的 genes.tsv.gz 无害，新目标名
  # features.tsv.gz 不存在会触发重新复制。
  std_names[gene_role] <- "features.tsv.gz"
  need_roles <- c("barcodes", gene_role, "matrix")

  roles <- vapply(files_df$file_name, role_of_file, character(1))
  missing_roles <- setdiff(need_roles, roles)
  if (length(missing_roles) > 0) {
    stop(paste0(gsm, " 三件套角色缺失: ", paste(missing_roles, collapse = ", "),
                "（现有文件: ", paste(files_df$file_name, collapse = ", "), "）"),
         call. = FALSE)
  }

  out_dir <- file.path(geo_dir, gse, "seurat_input", gsm)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  for (rl in need_roles) {
    src_row <- files_df[roles == rl, , drop = FALSE][1, ]
    src <- resolve_file_path(gse, src_row$file_name, src_row$extracted_path)
    if (is.na(src)) {
      stop(paste0(gsm, " 源文件不存在（清单路径与 supp_extracted 回退均失败）: ",
                  src_row$file_name), call. = FALSE)
    }
    dst <- file.path(out_dir, std_names[rl])
    src_sz <- file.size(src)
    if (file.exists(dst) && !is.na(file.size(dst)) &&
        file.size(dst) == src_sz) {
      next  # 幂等跳过：已存在且字节数匹配
    }
    ok <- file.copy(src, dst, overwrite = TRUE)
    if (!ok || is.na(file.size(dst)) || file.size(dst) != src_sz) {
      stop(paste0(gsm, " 复制/校验失败: ", basename(src), " -> ",
                  basename(dst)), call. = FALSE)
    }
  }
  out_dir
}

# 读取单个 10x H5 为稀疏矩阵（多 assay 时取 Gene Expression）
read_h5_counts <- function(h5_path) {
  m <- Read10X_h5(h5_path)
  if (is.list(m)) {
    if ("Gene Expression" %in% names(m)) {
      m <- m[["Gene Expression"]]
    } else {
      log_msg("WARN: ", basename(h5_path),
              " H5 含多个 assay 但无 'Gene Expression'，取第一个: ",
              names(m)[1])
      m <- m[[1]]
    }
  }
  m
}

# 计算 QC 指标：nFeature/nCount 由 CreateSeuratObject 自带；此处补 percent.mt/ribo。
# 修复 Bug B：PercentageFeatureSet 返回以细胞名命名的向量，Seurat v5 的
# [[<- 会对命名向量做名称重叠校验，多层/merge 场景可能误报
# "No cell overlap between new meta data and Seurat object"。
# 其返回值顺序 = Cells(object) 顺序，unname() 后按序赋值安全，绕开该校验。
add_qc_metrics <- function(seu) {
  seu$percent.mt   <- unname(PercentageFeatureSet(seu, pattern = "^MT-"))
  seu$percent.ribo <- unname(PercentageFeatureSet(seu, pattern = "^RP[SL]"))
  seu
}

# QC 阈值计算（在过滤前对象的 nFeature_RNA 上计算；全程写日志）
# mode: "scRNA" -> 硬上限 6000 兜底 + mt<=15; "snRNA" -> mt<=5
compute_qc_thresholds <- function(seu, gse, mode = c("scRNA", "snRNA"),
                                  context = "数据集合并对象") {
  mode <- match.arg(mode)
  nf  <- seu$nFeature_RNA
  med <- stats::median(nf, na.rm = TRUE)
  md  <- stats::mad(nf, na.rm = TRUE)   # 默认常数 1.4826 缩放后的 MAD
  upper <- med + 3 * md
  hard_cap <- NA_real_
  if (mode == "scRNA") {
    hard_cap <- SCRNA_NFEATURE_HARD
    upper <- min(upper, hard_cap)
  }
  mt_max <- if (mode == "scRNA") SCRNA_MT_MAX else SNRNA_MT_MAX
  # 防御：MAD=0 或异常导致 upper 不合法时退回 99% 分位数，并显式留痕
  if (!is.finite(upper) || upper <= NFEATURE_LOWER) {
    upper <- as.numeric(stats::quantile(nf, 0.99, na.rm = TRUE))
    log_msg("WARN: ", gse, " (", context, ") median+3*MAD 异常，",
            "nFeature 上限退回 99% 分位数 = ", round(upper, 1))
  }
  log_msg(gse, " (", context, ") QC 阈值: nFeature >= ", NFEATURE_LOWER,
          " 且 <= median(", round(med, 1), ") + 3*MAD(", round(md, 1),
          ") = ", round(med + 3 * md, 1),
          if (mode == "scRNA") paste0("，硬上限 ", hard_cap, " 兜底后取 ",
                                      round(upper, 1)) else "",
          "；percent.mt <= ", mt_max, "%")
  list(lower = NFEATURE_LOWER, upper = upper, median = med, mad = md,
       hard_cap = hard_cap, mt_max = mt_max, mode = mode)
}

# 应用 QC 过滤
apply_qc <- function(seu, thr) {
  subset(seu,
         subset = nFeature_RNA >= thr$lower &
                  nFeature_RNA <= thr$upper &
                  percent.mt   <= thr$mt_max)
}

# ---- QC 绘图（逐图 tryCatch，PNG + PDF 双格式） --------------------------------

# 保存单图（png+pdf），内部不再 tryCatch——由调用方逐图保护
save_plot_both <- function(p, base_name, width = 12, height = 7) {
  ggplot2::ggsave(paste0(base_name, ".png"), p, width = width, height = height,
                  dpi = 300, bg = "white", limitsize = FALSE)
  ggplot2::ggsave(paste0(base_name, ".pdf"), p, width = width, height = height,
                  limitsize = FALSE)
  invisible(TRUE)
}

# 带逐图保护的包装：任一图失败只记日志，不中断流程
safe_plot <- function(tag, plot_expr, base_name, width = 12, height = 7) {
  tryCatch({
    p <- plot_expr
    save_plot_both(p, base_name, width, height)
    log_msg("    图已写出: ", basename(base_name), ".png/.pdf")
  }, error = function(e) {
    log_msg("WARN: 绘图失败 [", tag, "]（不中断）: ", conditionMessage(e))
  })
}

# violin: 单指标按样本分面
qc_violin <- function(md, metric, title) {
  ggplot(md, aes(x = gsm, y = .data[[metric]], fill = gsm)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
    geom_boxplot(width = 0.12, outlier.size = 0.2, alpha = 0.6, linewidth = 0.2) +
    facet_wrap(~ gsm, scales = "free_x", nrow = 1) +
    labs(title = title, x = NULL, y = metric) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          legend.position = "none",
          plot.title = element_text(size = 11))
}

# scatter: nCount vs nFeature 按 percent.mt 着色，按样本分面
qc_scatter <- function(md, title) {
  ggplot(md, aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
    geom_point(size = 0.25, alpha = 0.5) +
    scale_color_gradientn(colors = c("#2166AC", "#FEE08B", "#B2182B"),
                          name = "percent.mt") +
    facet_wrap(~ gsm, scales = "free") +
    labs(title = title, x = "nCount_RNA", y = "nFeature_RNA") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(size = 11),
          axis.text = element_text(size = 7))
}

# 输出一个数据集过滤前/后的全套 QC 图（3 violin + 1 scatter，各两套）
make_qc_plots <- function(seu_before, seu_after, gse, thr) {
  md_pre  <- seu_before@meta.data
  md_post <- seu_after@meta.data
  n_samples <- length(unique(md_pre$gsm))
  w <- max(10, min(26, 1.6 * n_samples))   # 随样本数加宽，17 样本约 26in
  metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
  for (stage in c("prefilter", "postfilter")) {
    md <- if (stage == "prefilter") md_pre else md_post
    for (m in metrics) {
      safe_plot(
        paste(gse, stage, m, "violin"),
        qc_violin(md, m, paste0(gse, " ", stage, " QC violin: ", m,
                                " (thr: ", NFEATURE_LOWER, " <= nFeature <= ",
                                round(thr$upper, 1), ", percent.mt <= ",
                                thr$mt_max, "%)")),
        file.path(qc_dir, paste0(gse, "_QC_violin_", m, "_", stage)),
        width = w, height = 5.5)
    }
    safe_plot(
      paste(gse, stage, "scatter"),
      qc_scatter(md, paste0(gse, " ", stage,
                            " nCount vs nFeature (colored by percent.mt)")),
      file.path(qc_dir, paste0(gse, "_QC_scatter_nCount_nFeature_", stage)),
      width = w, height = max(5, ceiling(n_samples / 4) * 3.2))
  }
}

# ---- 样本级元数据解析 -----------------------------------------------------------

# GSE152805：从文件名解析 (SY|oLT|MT)_(113|116|118)；清单列兜底
parse_meta_gse152805 <- function(file_name, inv_row) {
  tissue <- inv_row$tissue_site
  donor  <- inv_row$donor_hint
  m <- regmatches(file_name, regexpr("(SY|oLT|MT)_(113|116|118)", file_name,
                                     perl = TRUE))
  if (length(m) == 1L) {
    parts <- strsplit(m, "_", fixed = TRUE)[[1]]
    if (is.na(tissue)) tissue <- if (parts[1] == "SY") "synovium" else parts[1]
    if (is.na(donor))  donor  <- parts[2]
  }
  list(condition = "OA",                 # 全部为 OA，无健康对照（红线 1）
       tissue = ifelse(is.na(tissue), NA_character_, tissue),
       donor  = ifelse(is.na(donor),  NA_character_, donor))
}

# GSE169454：condition 与供体编号从文件名 (normal|oa)(1-4) 解析
parse_meta_gse169454 <- function(file_name, inv_row) {
  cond <- inv_row$condition_hint
  m <- regmatches(file_name, regexpr("(normal|oa)([0-9]+)", file_name,
                                     ignore.case = TRUE, perl = TRUE))
  donor <- NA_character_
  if (length(m) == 1L) {
    donor <- m
    if (is.na(cond)) {
      cond <- if (grepl("^normal", m, ignore.case = TRUE)) "normal" else "OA"
    }
  }
  if (!is.na(cond) && tolower(cond) == "normal") cond <- "normal"
  if (!is.na(cond) && toupper(cond) == "OA")     cond <- "OA"
  list(condition = ifelse(is.na(cond), NA_character_, cond),
       tissue = "cartilage",
       donor = donor)
}

# GSE167186：donor = HM 编号（HM1-13,15,21-23，非连续，正则解析）
parse_meta_gse167186 <- function(file_name) {
  m <- regmatches(file_name, regexpr("HM[0-9]+", file_name, perl = TRUE))
  list(tissue = "skeletal_muscle",
       donor = if (length(m) == 1L) m else NA_character_)
}

# ---- pData 兜底工具（沿用 08a 实现思路，读取本地缓存的 *_pData.csv） ------------

read_pdata_all <- function(gse) {
  gse_dir <- file.path(geo_dir, gse)
  pfiles <- list.files(gse_dir, pattern = "_pData\\.csv$", full.names = TRUE)
  if (length(pfiles) == 0) {
    log_msg("WARN: ", gse, " 未找到 *_pData.csv（GEOquery pData 本地缓存缺失）。")
    return(NULL)
  }
  dfs <- lapply(pfiles, function(f) {
    tryCatch(utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE),
             error = function(e) {
               log_msg("WARN: 读取 pData 失败 ", basename(f), " : ",
                       conditionMessage(e))
               NULL
             })
  })
  dfs <- Filter(Negate(is.null), dfs)
  if (length(dfs) == 0) return(NULL)
  all_cols <- unique(unlist(lapply(dfs, names)))
  dfs <- lapply(dfs, function(d) {
    miss <- setdiff(all_cols, names(d))
    for (mc in miss) d[[mc]] <- NA
    d[all_cols]
  })
  do.call(rbind, dfs)
}

build_gsm_map <- function(pd) {
  if (is.null(pd)) return(NULL)
  gsm_col <- names(pd)[grepl("geo_accession|gsm", names(pd), ignore.case = TRUE)]
  if (length(gsm_col) == 0) {
    cand <- vapply(pd, function(x) any(grepl("^GSM[0-9]+", as.character(x))),
                   logical(1))
    gsm_col <- names(pd)[cand][1]
  }
  if (length(gsm_col) == 0 || is.na(gsm_col[1])) return(NULL)
  gsm <- as.character(pd[[gsm_col[1]]])
  title_col <- names(pd)[grepl("^title$|source_name|characteristics", names(pd),
                               ignore.case = TRUE)]
  hint <- if (length(title_col) > 0) {
    apply(pd[, title_col, drop = FALSE], 1,
          function(r) paste(na.omit(as.character(r)), collapse = " | "))
  } else rep(NA_character_, length(gsm))
  data.frame(gsm = gsm, pdata_hint = hint, stringsAsFactors = FALSE)
}

# ---- GSE167186 Old/Young 分组映射 -------------------------------------------------
# 优先级：1) SimplifiedMetadataSheet.xlsx（打印列名 + 模糊匹配 HM/group/age 列）
#        2) 本地 pData 提示文本中 old/young 关键词
#        3) 全部失败 -> stop（由调用方 tryCatch 捕获，状态记 FAIL 留痕）
# 返回: data.frame(gsm, group) 或 NULL

# 把任意文本规范化为 Old/Young/NA（aging proxy，严禁写成 sarcopenia）。
# 用词边界 \b 防 "threshold"/"smolder" 等子串误命中。
norm_old_young <- function(x) {
  x <- as.character(x)
  out <- rep(NA_character_, length(x))
  out[grepl("\\bold\\b",   x, ignore.case = TRUE, perl = TRUE)] <- "Old"
  out[grepl("\\byoung\\b", x, ignore.case = TRUE, perl = TRUE)] <- "Young"
  out
}

map_gse167186_from_xlsx <- function(gsms, donors) {
  xlsx <- file.path(geo_dir, "GSE167186", "GSE167186_SimplifiedMetadataSheet.xlsx")
  if (!file.exists(xlsx)) {
    log_msg("WARN: 未找到 ", xlsx, "，跳过 xlsx 分组通道。")
    return(NULL)
  }
  md <- tryCatch(as.data.frame(readxl::read_excel(xlsx)),
                 error = function(e) {
                   log_msg("WARN: readxl 读取失败: ", conditionMessage(e))
                   NULL
                 })
  if (is.null(md) || nrow(md) == 0) return(NULL)
  log_msg("GSE167186 xlsx 列名: ", paste(names(md), collapse = " | "))

  # 1) 找 HM/样本列：列值与 donor (HM 编号) 重合度最高者
  donors_up <- toupper(trimws(donors))
  best_col <- NULL; best_score <- 0L
  for (nm in names(md)) {
    vals <- toupper(trimws(as.character(md[[nm]])))
    score <- sum(donors_up %in% vals)
    if (score > best_score) { best_score <- score; best_col <- nm }
  }
  if (is.null(best_col) || best_score < ceiling(length(donors) / 2)) {
    log_msg("WARN: xlsx 中未找到与 HM 编号匹配的样本列（最佳匹配 ",
            best_score, "/", length(donors), "）。")
    return(NULL)
  }
  log_msg("GSE167186 xlsx 样本列判定: '", best_col, "' (匹配 ", best_score,
          "/", length(donors), " 个 HM 编号)")
  hm_vals <- toupper(trimws(as.character(md[[best_col]])))

  # 2) 找分组列：任一列 old/young 覆盖率 >= 50%
  group_col <- NULL
  for (nm in setdiff(names(md), best_col)) {
    g <- norm_old_young(md[[nm]])
    if (mean(!is.na(g)) >= 0.5) { group_col <- nm; break }
  }
  if (!is.null(group_col)) {
    grp <- norm_old_young(md[[group_col]])
    log_msg("GSE167186 xlsx 分组列判定: '", group_col, "' (Old/Young 文本)")
  } else {
    # 3) 退化通道：age 数值列中位数二分（显式 WARN 留痕）
    age_col <- names(md)[grepl("age", names(md), ignore.case = TRUE)]
    age_col <- setdiff(age_col, c(best_col, group_col))[1]
    grp <- NULL
    if (!is.na(age_col)) {
      av <- suppressWarnings(as.numeric(md[[age_col]]))
      if (sum(!is.na(av)) >= length(av) / 2) {
        cut <- stats::median(av, na.rm = TRUE)
        grp <- ifelse(av >= cut, "Old", "Young")
        log_msg("WARN: xlsx 无 Old/Young 文本列，改用年龄列 '", age_col,
                "' 中位数 (", cut, ") 二分推导 Old/Young（aging proxy 推导值，",
                "请在论文方法学部分如实说明）。")
      }
    }
    if (is.null(grp)) {
      log_msg("WARN: xlsx 中既无 Old/Young 文本列也无可用 age 数值列。")
      return(NULL)
    }
  }

  lut <- data.frame(donor = hm_vals, group = grp, stringsAsFactors = FALSE)
  lut <- lut[!is.na(lut$donor) & !is.na(lut$group), , drop = FALSE]
  hit <- match(donors_up, lut$donor)
  if (any(is.na(hit))) {
    log_msg("WARN: xlsx 分组映射不完整，未匹配 HM: ",
            paste(donors[is.na(hit)], collapse = ", "))
    return(NULL)
  }
  data.frame(gsm = gsms, group = lut$group[hit], stringsAsFactors = FALSE)
}

map_gse167186_from_pdata <- function(gsms) {
  pd <- read_pdata_all("GSE167186")
  gm <- build_gsm_map(pd)
  if (is.null(gm)) {
    log_msg("WARN: pData 兜底通道不可用（无 pData 或无法定位 GSM 列）。")
    return(NULL)
  }
  idx <- match(gsms, gm$gsm)
  if (any(is.na(idx))) {
    log_msg("WARN: pData 中缺少 GSM: ",
            paste(gsms[is.na(idx)], collapse = ", "))
    return(NULL)
  }
  grp <- norm_old_young(gm$pdata_hint[idx])
  if (any(is.na(grp))) {
    log_msg("WARN: pData 提示文本未能解析出全部 Old/Young 分组（缺失 ",
            sum(is.na(grp)), " 个 GSM）。")
    return(NULL)
  }
  log_msg("GSE167186 分组来自 pData 兜底通道（GEOquery 本地缓存）。")
  data.frame(gsm = gsms, group = grp, stringsAsFactors = FALSE)
}

map_gse167186_groups <- function(gsms, donors) {
  mp <- map_gse167186_from_xlsx(gsms, donors)
  if (is.null(mp)) mp <- map_gse167186_from_pdata(gsms)
  if (is.null(mp)) {
    stop(paste0(
      "GSE167186 Old/Young 分组映射失败（xlsx 与 pData 双通道均不可用/不完整）。",
      "请人工检查 data/GEO/GSE167186/GSE167186_SimplifiedMetadataSheet.xlsx ",
      "的列名与本日志中打印的列名，修正后再跑；严禁臆造分组。"), call. = FALSE)
  }
  log_msg("GSE167186 分组映射完成: Old=", sum(mp$group == "Old"),
          " / Young=", sum(mp$group == "Young"),
          "（aging proxy，不得表述为 sarcopenia）")
  mp
}

# ---- qc_summary 行构建 -----------------------------------------------------------
# 对每个 GSM 汇总过滤前后细胞数与中位指标；阈值来自 thr（可能为数据集级或样本级）
qc_summary_rows <- function(seu_before, seu_after, gse, thr, filter_note) {
  md_pre  <- seu_before@meta.data
  md_post <- seu_after@meta.data
  gsms <- sort(unique(md_pre$gsm))
  med_pre <- function(v, g) stats::median(v[md_pre$gsm == g], na.rm = TRUE)
  med_post <- function(v, g) {
    x <- v[md_post$gsm == g]
    if (length(x) == 0) return(NA_real_)
    stats::median(x, na.rm = TRUE)
  }
  rows <- lapply(gsms, function(g) {
    meta <- md_pre[md_pre$gsm == g, , drop = FALSE][1, ]
    data.frame(
      gse           = gse,
      gsm           = g,
      condition     = as.character(meta$condition),
      tissue        = as.character(meta$tissue),
      donor         = as.character(meta$donor),
      cells_before  = sum(md_pre$gsm == g),
      cells_after   = sum(md_post$gsm == g),
      median_nFeature_before = med_pre(md_pre$nFeature_RNA, g),
      median_nCount_before   = med_pre(md_pre$nCount_RNA, g),
      median_pctMT_before    = med_pre(md_pre$percent.mt, g),
      median_pctRibo_before  = med_pre(md_pre$percent.ribo, g),
      median_nFeature_after  = med_post(md_post$nFeature_RNA, g),
      median_nCount_after    = med_post(md_post$nCount_RNA, g),
      median_pctMT_after     = med_post(md_post$percent.mt, g),
      median_pctRibo_after   = med_post(md_post$percent.ribo, g),
      thr_nFeature_lower = thr$lower,
      thr_nFeature_upper = round(thr$upper, 1),
      thr_nFeature_median = round(thr$median, 1),
      thr_nFeature_mad    = round(thr$mad, 1),
      thr_nFeature_hard_cap = ifelse(is.na(thr$hard_cap), NA, thr$hard_cap),
      thr_pct_mt_max = thr$mt_max,
      thr_mode       = thr$mode,
      filter_note    = filter_note,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# ---- 建对象与合并 -----------------------------------------------------------------

# 每个 GSM 一个 CreateSeuratObject（min.cells=3, min.features=200），
# 附 orig.ident/gsm/condition/tissue/donor 元数据
create_one_seurat <- function(counts, gse, gsm, condition, tissue, donor) {
  # fix2 防御性守卫（对三个数据集统一生效）：H5/三件套读出 counts 的 barcode
  # 列名为 NULL 或重复时，CreateSeuratObject 内部 object[[]] <- meta 的行名
  # 重叠校验会报 "No cell overlap between new meta data and Seurat object"。
  # dgCMatrix 直接改 colnames 合法（Matrix 包支持）。
  if (is.null(colnames(counts))) {
    log_msg("WARN: [", gsm, "] counts 无 barcode 列名，按 GSM 合成: ",
            gsm, "_cell_<n>")
    colnames(counts) <- paste0(gsm, "_cell_", seq_len(ncol(counts)))
  }
  if (anyDuplicated(colnames(counts))) {
    log_msg("WARN: [", gsm, "] barcode 重复 ",
            sum(duplicated(colnames(counts))), " 个，make.unique 去重")
    colnames(counts) <- make.unique(colnames(counts))
  }
  seu <- CreateSeuratObject(counts = counts, project = gse,
                            min.cells = 3, min.features = 200)
  # fix3：5 个标量赋值统一 unname()。上游可能传入命名长度1向量（实证：
  # process_gse167186 的 donors <- vapply(gsms, ...) 结果自动带名），
  # ifelse 会保留 test 的名字属性；Seurat v5 的 [[<- 对【命名向量】做
  # 名称重叠校验，名字（GSM 号）与细胞 barcode 零重叠即报
  # "No cell overlap between new meta data and Seurat object"。
  seu$orig.ident <- unname(gsm)
  seu$gsm        <- unname(gsm)
  seu$condition  <- unname(ifelse(is.na(condition), "unknown", condition))
  seu$tissue     <- unname(ifelse(is.na(tissue),    "unknown", tissue))
  seu$donor      <- unname(ifelse(is.na(donor),     "unknown", donor))
  seu
}

merge_seurat_list <- function(objs, gsms) {
  if (length(objs) == 1L) return(objs[[1]])
  merge(x = objs[[1]], y = objs[-1], add.cell.ids = gsms)
}

# ---- scRNA 数据集通用处理（GSE152805 v2 / GSE169454 v3） --------------------------
process_scrna_dataset <- function(inv, gse, format = c("v2", "v3"),
                                  parse_meta_fn) {
  format <- match.arg(format)
  log_msg("---- 数据集: ", gse, " (10x ", format, ", scRNA) ----")
  sub <- inv[inv$gse == gse & !is.na(inv$gsm), , drop = FALSE]
  if (nrow(sub) == 0) stop(gse, " 在文件清单中无任何记录。", call. = FALSE)
  gsms <- sort(unique(sub$gsm))
  log_msg(gse, ": 清单内 ", length(gsms), " 个 GSM，", nrow(sub), " 个文件。")

  objs <- vector("list", length(gsms))
  names(objs) <- gsms
  for (gsm in gsms) {
    gdf <- sub[sub$gsm == gsm, , drop = FALSE]
    meta <- parse_meta_fn(gdf$file_name[1], gdf[1, ])
    log_msg("  [", gsm, "] tissue=", meta$tissue, " donor=", meta$donor,
            " condition=", meta$condition)
    mtx_dir <- prepare_read10x_dir(gse, gsm, gdf, format = format)
    counts <- Read10X(data.dir = mtx_dir)
    if (is.list(counts)) {
      # 多 assay list：优先 Gene Expression，缺失时取第一个并 WARN（与
      # read_h5_counts 的回退逻辑一致，避免 [[NULL]] 静默丢数据）
      if ("Gene Expression" %in% names(counts)) {
        counts <- counts[["Gene Expression"]]
      } else {
        log_msg("WARN: [", gsm, "] Read10X 返回多 assay 但无 'Gene Expression'，",
                "取第一个: ", names(counts)[1])
        counts <- counts[[1]]
      }
    }
    objs[[gsm]] <- create_one_seurat(counts, gse, gsm, meta$condition,
                                     meta$tissue, meta$donor)
    # 0 细胞守卫（Bug B 双保险之二）：min.features=200 可能把退化样本过滤空，
    # 空对象进 merge 会触发难以定位的元数据重叠错误；scRNA 路径直接 stop 留痕
    if (ncol(objs[[gsm]]) == 0) {
      stop(gsm, " 建对象后细胞数为 0（min.cells=3/min.features=200 过滤后为空），",
           "请检查该样本三件套是否完整/对应。", call. = FALSE)
    }
    log_msg("  [", gsm, "] 建对象: ", ncol(objs[[gsm]]), " 细胞 x ",
            nrow(objs[[gsm]]), " 基因")
    rm(counts)
  }

  log_msg(gse, ": 合并 ", length(objs), " 个样本 (add.cell.ids=GSM) ...")
  seu <- merge_seurat_list(objs, gsms)
  rm(objs); gc(verbose = FALSE)
  log_msg(gse, ": 合并后 ", ncol(seu), " 细胞 x ", nrow(seu), " 基因。")

  seu <- add_qc_metrics(seu)
  thr <- compute_qc_thresholds(seu, gse, mode = "scRNA",
                               context = "数据集合并对象(过滤前)")
  seu_f <- apply_qc(seu, thr)
  log_msg(gse, ": QC 过滤 ", ncol(seu), " -> ", ncol(seu_f), " 细胞 (保留 ",
          round(100 * ncol(seu_f) / max(1, ncol(seu)), 1), "%)")
  if (ncol(seu_f) == 0) stop(gse, " QC 后细胞数为 0，请检查阈值与数据。",
                             call. = FALSE)

  make_qc_plots(seu, seu_f, gse, thr)

  log_msg(gse, ": LogNormalize + FindVariableFeatures(vst, 2000) ...")
  seu_f <- NormalizeData(seu_f, normalization.method = "LogNormalize",
                         scale.factor = 10000, verbose = FALSE)
  seu_f <- FindVariableFeatures(seu_f, selection.method = "vst",
                                nfeatures = 2000, verbose = FALSE)
  log_msg(gse, ": 标准化完成，VST 高变基因 ",
          length(VariableFeatures(seu_f)), " 个。（G1 不做 ScaleData/PCA/聚类/整合）")

  rds_path <- file.path(res_dir, paste0(gse, "_seurat_G1.rds"))
  saveRDS(seu_f, rds_path)
  log_msg(gse, ": RDS 已保存: ", rds_path, " (",
          round(file.size(rds_path) / 1e6, 1), " MB)")

  qs <- qc_summary_rows(seu, seu_f, gse, thr,
                        filter_note = "dataset-level thresholds on merged object")
  cells_before <- ncol(seu)
  rm(seu, seu_f); gc(verbose = FALSE)
  list(qc_rows = qs, n_samples = length(gsms), cells_before = cells_before,
       cells_after = sum(qs$cells_after), rds_path = rds_path,
       note = "OK: merged -> dataset-level QC -> LogNormalize+VST2000")
}

# ---- GSE167186 snRNA 处理（H5 + 内存预案 + Old/Young 分组） -----------------------
# 修复 Bug A：签名加 gse 形参以兼容 run_dataset 的 fn(inv, gse = gse, ...) 透传
process_gse167186 <- function(inv, gse = "GSE167186") {
  log_msg("---- 数据集: ", gse, " (10x H5, snRNA; Old vs Young aging proxy) ----")
  sub <- inv[inv$gse == gse & !is.na(inv$gsm), , drop = FALSE]
  sub <- sub[vapply(sub$file_name, role_of_file, character(1)) == "h5", ,
             drop = FALSE]
  if (nrow(sub) == 0) stop(gse, " 文件清单中无任何 H5 记录。", call. = FALSE)
  gsms <- sort(unique(sub$gsm))
  donors <- vapply(gsms, function(g) {
    parse_meta_gse167186(sub$file_name[sub$gsm == g][1])$donor
  }, character(1))
  if (any(is.na(donors))) {
    # 修复：donors 与 gsms 对齐（每 GSM 一个），索引必须用 gsms 侧
    stop(gse, " 无法从文件名解析 HM 编号: ",
         paste(gsms[is.na(donors)], collapse = ", "), call. = FALSE)
  }
  log_msg(gse, ": 清单内 ", length(gsms), " 个 GSM (",
          paste(donors, collapse = ", "), ")。注意 HM 编号非连续，均按正则解析。")

  # Old/Young 分组映射（xlsx -> pData -> stop 留痕）
  grp_map <- map_gse167186_groups(gsms, donors)

  # ---- 第一遍：只读维度估算总细胞数 -------------------------------------------------
  # 修复（内存预案必须保护峰值）：路径决策不得发生在 17 个 H5 全部驻留内存之后。
  # 每个 H5 读入后只取维度并立即 rm + gc，此阶段峰值任何时候只有 1 个样本矩阵。
  total_cells <- 0
  for (i in seq_along(gsms)) {
    gsm <- gsms[i]
    row <- sub[sub$gsm == gsm, , drop = FALSE][1, ]
    h5 <- resolve_file_path(gse, row$file_name, row$extracted_path)
    if (is.na(h5)) stop(gsm, " H5 文件不存在: ", row$file_name, call. = FALSE)
    mat_i <- read_h5_counts(h5)
    nc <- ncol(mat_i)
    total_cells <- total_cells + nc
    log_msg("  [", gsm, " / ", donors[i], "] 维度探测: ", nc, " 核 x ",
            nrow(mat_i), " 基因（矩阵已立即释放）")
    rm(mat_i)
    gc(verbose = FALSE)
  }
  log_msg(gse, ": 估算合并前总细胞数 = ", total_cells,
          "（内存预案阈值 ", SNRNA_CELL_THRESHOLD, "）")

  # 内存预案：BPCells 可用性仅记录留痕；为避免磁盘矩阵 API 不兼容风险，
  # 统一走内存路径——超阈值时逐样本「读入->建对象->过滤->释放」后再 merge。
  bpcells_ok <- requireNamespace("BPCells", quietly = TRUE)
  log_msg(gse, ": BPCells 包", if (bpcells_ok) "可用" else "不可用",
          "；按既定策略选择内存路径（规避 on-disk 矩阵 API 兼容风险）。")

  pre_filter <- total_cells > SNRNA_CELL_THRESHOLD
  if (pre_filter) {
    log_msg(gse, ": 总细胞数 ", total_cells, " > ", SNRNA_CELL_THRESHOLD,
            "，走【逐样本 QC 过滤后再合并】路径（阈值按样本独立计算，留痕于 qc_summary）。")
  } else {
    log_msg(gse, ": 总细胞数 ", total_cells, " <= ", SNRNA_CELL_THRESHOLD,
            "，走【先合并后数据集级 QC】路径。")
  }

  # ---- 第二遍：逐样本读入 -> 建对象 ->（超阈值时）样本级 QC -> 释放 ------------------
  objs <- vector("list", length(gsms))
  names(objs) <- gsms
  sample_thr <- vector("list", length(gsms))
  names(sample_thr) <- gsms
  dropped <- character(0)   # QC 后 0 细胞被丢弃的样本（日志 + status 留痕）
  for (i in seq_along(gsms)) {
    gsm <- gsms[i]
    row <- sub[sub$gsm == gsm, , drop = FALSE][1, ]
    h5 <- resolve_file_path(gse, row$file_name, row$extracted_path)
    if (is.na(h5)) stop(gsm, " H5 文件不存在: ", row$file_name, call. = FALSE)
    mat_i <- read_h5_counts(h5)
    # fix2 诊断日志：逐步留痕，定位 "No cell overlap" 根因
    log_msg("  [", gsm, "] 诊断: H5 dim=", nrow(mat_i), " x ", ncol(mat_i),
            " | colnames 为 NULL: ", is.null(colnames(mat_i)),
            " | 重复 barcode 数: ",
            if (is.null(colnames(mat_i))) NA else sum(duplicated(colnames(mat_i))))
    grp <- grp_map$group[grp_map$gsm == gsm][1]
    # fix2：create_one_seurat + QC 四步包一层 tryCatch，报错消息带
    # [GSM] step= 前缀——即使防御性守卫没盖住房因，下次报错也能精确定位
    # 到样本与步骤。tryCatch 的表达式块不另建作用域，step/thr 赋值外层可见。
    step <- "create_one_seurat"
    thr <- NULL
    seu <- tryCatch({
      s <- create_one_seurat(mat_i, gse, gsm, grp, "skeletal_muscle", donors[i])
      log_msg("  [", gsm, "] 诊断: 建对象完成 ncol=", ncol(s))
      if (pre_filter) {
        step <- "add_qc_metrics"
        s <- add_qc_metrics(s)
        step <- "compute_qc_thresholds"
        thr <- compute_qc_thresholds(s, gse, mode = "snRNA",
                                     context = paste0("样本 ", gsm))
        step <- "apply_qc"
        n_before <- ncol(s)
        s <- apply_qc(s, thr)
        log_msg("  [", gsm, "] 样本级 QC: ", n_before, " -> ", ncol(s), " 核")
      }
      s
    }, error = function(e) {
      stop(paste0("[", gsm, "] step=", step, ": ", conditionMessage(e)),
           call. = FALSE)
    })
    rm(mat_i); gc(verbose = FALSE)   # 峰值 = 1 个样本矩阵 + 已建对象
    if (ncol(seu) == 0) {
      # 0 细胞守卫（合并 fix1 的两处守卫：建对象后为空，或样本级 QC 后为空，
      # 均不允许进 merge）：同通道丢样留痕——本样本不进入 objs/sample_thr，
      # kept_gsms/sample_thr 对齐逻辑天然兼容
      log_msg("ERROR: [", gsm, "] 建对象/QC 后细胞数为 0，丢弃该样本",
              "（不中断其余样本；status 与日志留痕）。")
      dropped <- c(dropped, gsm)
      objs[[gsm]] <- NULL
      rm(seu); gc(verbose = FALSE)
      next
    }
    if (pre_filter) sample_thr[[gsm]] <- thr
    objs[[gsm]] <- seu
    rm(seu); gc(verbose = FALSE)
  }
  kept_gsms <- names(objs)   # 注意：objs[[gsm]] <- NULL 会删除该元素
  if (length(kept_gsms) == 0) {
    stop(gse, " 全部样本 QC 后均为 0 细胞，无法合并。", call. = FALSE)
  }
  if (length(dropped) > 0) {
    log_msg("WARN: ", gse, " 共丢弃 ", length(dropped), " 个 0 细胞样本: ",
            paste(dropped, collapse = ", "), "；实际合并 ", length(kept_gsms),
            " 个样本: ", paste(kept_gsms, collapse = ", "))
    sample_thr <- sample_thr[kept_gsms]   # 与实际进入 merge 的 GSM 对齐
  }

  log_msg(gse, ": 合并 ", length(objs), " 个样本 (add.cell.ids=GSM) ...")
  seu <- merge_seurat_list(objs, kept_gsms)
  rm(objs); gc(verbose = FALSE)
  log_msg(gse, ": 合并后 ", ncol(seu), " 核 x ", nrow(seu), " 基因。")

  if (pre_filter) {
    # 逐样本已过滤：seu 即过滤后对象。为出过滤前/后对比图与 qc_summary，
    # 需要一个"过滤前"参照——此处无法重读全部原始数据（内存预案的初衷），
    # 故 QC 图仅画过滤后一套；qc_summary 的 before 列由样本级日志值回填 NA。
    log_msg(gse, ": 逐样本过滤路径下不重建过滤前大对象；QC 图仅输出 postfilter，",
            "qc_summary 过滤前统计记为 NA（样本级 before/after 细胞数见上方日志）。")
    seu_f <- seu
    make_qc_plots_postfilter_only(seu_f, gse, sample_thr)
    qs <- qc_summary_rows_postfilter_only(seu_f, gse, sample_thr,
      filter_note = paste0("per-sample QC before merge (total estimated cells=",
                           total_cells, ">", SNRNA_CELL_THRESHOLD, ")"))
  } else {
    seu <- add_qc_metrics(seu)
    thr <- compute_qc_thresholds(seu, gse, mode = "snRNA",
                                 context = "数据集合并对象(过滤前)")
    seu_f <- apply_qc(seu, thr)
    log_msg(gse, ": QC 过滤 ", ncol(seu), " -> ", ncol(seu_f), " 核 (保留 ",
            round(100 * ncol(seu_f) / max(1, ncol(seu)), 1), "%)")
    if (ncol(seu_f) == 0) stop(gse, " QC 后细胞数为 0，请检查阈值与数据。",
                               call. = FALSE)
    make_qc_plots(seu, seu_f, gse, thr)
    qs <- qc_summary_rows(seu, seu_f, gse, thr,
                          filter_note = "dataset-level thresholds on merged object")
  }

  log_msg(gse, ": LogNormalize + FindVariableFeatures(vst, 2000) ...")
  seu_f <- NormalizeData(seu_f, normalization.method = "LogNormalize",
                         scale.factor = 10000, verbose = FALSE)
  seu_f <- FindVariableFeatures(seu_f, selection.method = "vst",
                                nfeatures = 2000, verbose = FALSE)
  log_msg(gse, ": 标准化完成，VST 高变基因 ",
          length(VariableFeatures(seu_f)), " 个。（G1 不做 ScaleData/PCA/聚类/整合）")

  rds_path <- file.path(res_dir, paste0(gse, "_seurat_G1.rds"))
  saveRDS(seu_f, rds_path)
  log_msg(gse, ": RDS 已保存: ", rds_path, " (",
          round(file.size(rds_path) / 1e6, 1), " MB)")

  cells_after <- ncol(seu_f)
  rm(seu_f); if (exists("seu")) rm(seu); gc(verbose = FALSE)
  note <- if (pre_filter)
    "OK: per-sample QC -> merge -> LogNormalize+VST2000 (memory-safe path)"
  else
    "OK: merged -> dataset-level QC -> LogNormalize+VST2000"
  if (length(dropped) > 0) {
    note <- paste0(note, "; dropped ", length(dropped), " zero-cell sample(s): ",
                   paste(dropped, collapse = ", "))
  }
  list(qc_rows = qs, n_samples = length(kept_gsms),
       cells_before = if (pre_filter) NA_integer_ else sum(qs$cells_before),
       cells_after = cells_after, rds_path = rds_path,
       note = note)
}

# 逐样本过滤路径专用：仅 postfilter 图（对象已过滤，无过滤前参照）。
# 图注写入样本级阈值的范围（各样本阈值独立计算，见 qc_summary 与日志）。
make_qc_plots_postfilter_only <- function(seu_after, gse, sample_thr) {
  md_post <- seu_after@meta.data
  n_samples <- length(unique(md_post$gsm))
  w <- max(10, min(26, 1.6 * n_samples))
  thr_ok <- Filter(Negate(is.null), sample_thr)
  thr_note <- "per-sample QC before merge"
  if (length(thr_ok) > 0) {
    uppers <- vapply(thr_ok, function(x) x$upper, numeric(1))
    mt_maxs <- vapply(thr_ok, function(x) x$mt_max, numeric(1))
    thr_note <- paste0("per-sample thr: nFeature ", NFEATURE_LOWER, "-[",
                       round(min(uppers), 1), ", ", round(max(uppers), 1),
                       "], percent.mt <= ", max(mt_maxs), "%")
  }
  for (m in c("nFeature_RNA", "nCount_RNA", "percent.mt")) {
    safe_plot(
      paste(gse, "postfilter", m, "violin"),
      qc_violin(md_post, m, paste0(gse, " postfilter QC violin: ", m,
                                   " (", thr_note, ")")),
      file.path(qc_dir, paste0(gse, "_QC_violin_", m, "_postfilter")),
      width = w, height = 5.5)
  }
  safe_plot(
    paste(gse, "postfilter scatter"),
    qc_scatter(md_post, paste0(gse, " postfilter nCount vs nFeature",
                               " (colored by percent.mt; ", thr_note, ")")),
    file.path(qc_dir, paste0(gse, "_QC_scatter_nCount_nFeature_postfilter")),
    width = w, height = max(5, ceiling(n_samples / 4) * 3.2))
}

# 逐样本过滤路径专用 qc_summary：before 列 NA，阈值列为各样本自身阈值
qc_summary_rows_postfilter_only <- function(seu_after, gse, sample_thr,
                                            filter_note) {
  md_post <- seu_after@meta.data
  gsms <- sort(unique(md_post$gsm))
  rows <- lapply(gsms, function(g) {
    thr <- sample_thr[[g]]
    meta <- md_post[md_post$gsm == g, , drop = FALSE]
    sub_md <- meta
    data.frame(
      gse           = gse,
      gsm           = g,
      condition     = as.character(meta$condition[1]),
      tissue        = as.character(meta$tissue[1]),
      donor         = as.character(meta$donor[1]),
      cells_before  = NA_integer_,   # 未保留过滤前对象；样本级 before 见日志
      cells_after   = nrow(sub_md),
      median_nFeature_before = NA_real_,
      median_nCount_before   = NA_real_,
      median_pctMT_before    = NA_real_,
      median_pctRibo_before  = NA_real_,
      median_nFeature_after  = stats::median(sub_md$nFeature_RNA, na.rm = TRUE),
      median_nCount_after    = stats::median(sub_md$nCount_RNA, na.rm = TRUE),
      median_pctMT_after     = stats::median(sub_md$percent.mt, na.rm = TRUE),
      median_pctRibo_after   = stats::median(sub_md$percent.ribo, na.rm = TRUE),
      thr_nFeature_lower = thr$lower,
      thr_nFeature_upper = round(thr$upper, 1),
      thr_nFeature_median = round(thr$median, 1),
      thr_nFeature_mad    = round(thr$mad, 1),
      thr_nFeature_hard_cap = ifelse(is.na(thr$hard_cap), NA, thr$hard_cap),
      thr_pct_mt_max = thr$mt_max,
      thr_mode       = thr$mode,
      filter_note    = filter_note,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# ---- 主流程：逐数据集独立处理，tryCatch 隔离 --------------------------------------
inv <- load_inventory()

status_rows <- list()
qc_rows_all <- list()

# 修复 Bug A：gse 由 run_dataset 透传（fn(inv, gse = gse, ...)），调用点不得
# 再传 gse= 实参——否则命名实参精确匹配 run_dataset 形参 gse，位置实参顺延，
# fn 会接到字符串而报 could not find function "fn"（实机已验证）。
run_dataset <- function(gse, fn, ...) {
  tryCatch(
    fn(inv, gse = gse, ...),
    error = function(e) {
      log_msg("ERROR: ", gse, " 处理失败（不影响其他数据集）: ",
              conditionMessage(e))
      list(qc_rows = NULL, n_samples = NA_integer_, cells_before = NA_integer_,
           cells_after = NA_integer_, rds_path = NA_character_,
           note = paste0("FAIL: ", conditionMessage(e)))
    })
}

# 1) GSE152805：OA 膝 scRNA（10x v2），仅供体内组织间定位，无健康对照（红线 1）
res_152805 <- run_dataset("GSE152805", process_scrna_dataset,
                          format = "v2",
                          parse_meta_fn = parse_meta_gse152805)

# 2) GSE169454：软骨细胞 scRNA（10x v3 filtered 三件套）
res_169454 <- run_dataset("GSE169454", process_scrna_dataset,
                          format = "v3",
                          parse_meta_fn = parse_meta_gse169454)

# 3) GSE167186：肌肉 snRNA（H5），Old vs Young aging proxy（红线 2）
res_167186 <- run_dataset("GSE167186", process_gse167186)

results_by_gse <- list(GSE152805 = res_152805, GSE169454 = res_169454,
                       GSE167186 = res_167186)

# ---- 汇总输出 ---------------------------------------------------------------------
for (gse in names(results_by_gse)) {
  r <- results_by_gse[[gse]]
  status_rows[[length(status_rows) + 1L]] <- data.frame(
    gse          = gse,
    status       = if (grepl("^FAIL", r$note)) "FAIL" else "OK",
    n_samples    = r$n_samples,
    cells_before = r$cells_before,
    cells_after  = r$cells_after,
    rds_path     = r$rds_path,
    note         = r$note,
    stringsAsFactors = FALSE
  )
  if (!is.null(r$qc_rows)) qc_rows_all[[length(qc_rows_all) + 1L]] <- r$qc_rows
}

status_df <- do.call(rbind, status_rows)
status_csv <- file.path(res_dir, "moduleG_G1_status.csv")
tryCatch({
  utils::write.csv(status_df, status_csv, row.names = FALSE)
  log_msg("写出: moduleG_G1_status.csv (", nrow(status_df), " 行)")
}, error = function(e) log_msg("ERROR: 写 status CSV 失败: ", conditionMessage(e)))

qc_csv <- file.path(res_dir, "moduleG_G1_qc_summary.csv")
tryCatch({
  if (length(qc_rows_all) > 0) {
    qc_df <- do.call(rbind, qc_rows_all)
    utils::write.csv(qc_df, qc_csv, row.names = FALSE)
    log_msg("写出: moduleG_G1_qc_summary.csv (", nrow(qc_df), " 行)")
  } else {
    log_msg("WARN: 无任何 QC 汇总行可写出（所有数据集均失败）。")
  }
}, error = function(e) log_msg("ERROR: 写 qc_summary CSV 失败: ", conditionMessage(e)))

# ---- 完成提示 ---------------------------------------------------------------------
log_msg("==== G1 汇总 ====")
for (i in seq_len(nrow(status_df))) {
  log_msg(status_df$gse[i], ": ", status_df$status[i],
          " | samples=", status_df$n_samples[i],
          " | cells_before=", status_df$cells_before[i],
          " | cells_after=", status_df$cells_after[i],
          " | ", status_df$note[i])
}

if (all(status_df$status == "OK")) {
  msg <- paste0("模块 G G1 完成。红线提醒：GSE152805 无健康对照仅供体内组织间定位；",
                "GSE167186 为 Old vs Young aging proxy 不得写成 sarcopenia；",
                "组间统计推断与整合均在 G2 及以后。")
  cat(msg, "\n"); log_msg(msg)
} else {
  bad <- status_df[status_df$status != "OK", ]
  msg <- paste0("模块 G G1 存在失败数据集：\n",
                paste0("  - ", bad$gse, " | ", bad$note, collapse = "\n"),
                "\n修复后直接重跑本脚本即可（复制/读入均幂等）。")
  cat(msg, "\n"); log_msg(msg)
  # 不 stop：已成功数据集的 RDS/CSV 均已落盘，status CSV 已留痕
}
