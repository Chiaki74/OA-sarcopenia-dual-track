#!/usr/bin/env Rscript
# ============================================================================
# 09_moduleG_downstream_20260823.R
# ----------------------------------------------------------------------------
# 模块 G / G2 脚本：标志物定位、三轨 UCell 评分与 donor 级 pseudobulk 统计。
# 编码: UTF-8（中文 Windows 用户请确认 RStudio 以 UTF-8 打开）
#
# 版本日期: 2026-08-24 fix3 — 实机三跑审查注释质量，修归一化 argmax 两漏洞：
#   ① 归一化得分并列 1.0 时按原始得分决胜（fix2 按列序会错标：
#   GSE167186 cl8 内皮 0.65 被 Immune 0.21 抢占）；② 原始得分下限 0.05，
#   剔除近零表达被归一化放大的噪声标签（cl11 Adipocyte 0.0057），全部低于
#   下限标 unknown；改标 cluster 写日志留痕。
# 2026-08-24 fix2 — 实机二跑三数据集全绿，但 GSE167186 自动
#   注释塌缩（12 cluster -> 1 类）：原始平均表达打分被高丰度肌核标志物
#   （TTN/NEB/MYH）淹没。修复：cluster×celltype 得分矩阵按列 max 归一化
#   到 [0,1] 后再 argmax（"相对自身最高"胜出）；CSV 同时留原始+归一化
#   双口径得分矩阵供人工复核。对 GSE152805/GSE169454 可能微调边界
#   cluster 归属，属预期内的更合理结果；经典 DotPlot（cluster 级）不受影响。
# 2026-08-24 fix1 — 实机首跑三数据集全 FAIL 于 UCell 步骤：
#   "Unrecognized input format." 系 UCell ≥2.10 的 ScoreSignatures_UCell
#   移除 Seurat 对象分支所致（只接受 SCE/matrix/dgCMatrix/data.frame），
#   改用官方 Seurat 接口 AddModuleScore_UCell(assay="RNA", slot="data",
#   name="")。该错误发生在自动注释之后，此前产物（UMAP/DotPlot 等）完好。
# 2026-08-23 v1.1 — 独立审查 10 条建议级问题全折叠修订：
#   1) 三轨表剔除 gene/track NA 行 + 轨归属改 match 取行；
#   2) aggregate_counts 聚合单位 anyNA 即 stop（报数据集与 NA 数）；
#   3) log2CPM 显式 as.matrix 后计算，不依赖 Matrix 隐式回收语义；
#   4) 依赖版本下限：Seurat>=5.3.0 / harmony>=1.1 / UCell>=2.4 / patchwork 存在性；
#   5) 新增 <GSE>_per_unit_values.csv（10 标志物 log2CPM + 三轨 per-unit 中位数）；
#   6) 删除 exclude_groups 死配置字段；7) cluster 按数值排序；
#   8) 全部 write.csv 加 fileEncoding="UTF-8"（GBK Windows 中文 note 安全）；
#   9) wilcox p 非有限（NaN/Inf）归一为 NA 并注明 all-zero differences；
#   10) JoinLayers 后 VariableFeatures<1000 时 WARN（VST2000 存活检查）。
# v1 — 从 G1（08_moduleG_singlecell_fix3，已全绿闭环：
#   GSE152805=35503 / GSE169454=65042 / GSE167186=52269 细胞入库）的 RDS 继续。
#
# 【红线（与 G1 一致，运行日志再打印一次）】
#   1. GSE152805 无健康对照：仅供体内 oLT vs MT 配对对比（n=3 对，仅描述性）
#      与组织定位，严禁病例-对照表述；SY 样本不参与统计比较。
#   2. GSE167186 分组为 Old vs Young（aging proxy），严禁表述为 sarcopenia
#      （肌少症）证据——日志与 status note 双留痕。
#   3. 一切差异统计以 donor/样本级 pseudobulk 为底线，严禁细胞级 p 值
#      （细胞是伪重复）；Harmony 仅用于聚类/UMAP，严禁作为差异或评分输入。
#
# 【输入】
#   results/moduleG/<GSE>_seurat_G1.rds  (Seurat v5；counts/data 按 GSM 分层；
#     meta: orig.ident/gsm/condition/tissue/donor/nFeature_RNA/nCount_RNA/
#     percent.mt/percent.ribo；已 LogNormalize+VST2000；未 ScaleData/PCA/聚类)
#   results/moduleC/moduleC_gene_tracks.csv
#     (列: gene, track, direction_OA, direction_muscle, in_WGCNA_key；
#      track 仅三种: concordant(27) / mirror_OAup_muscleDown(37) /
#      mirror_OAdown_muscleUp(35))
#
# 【流程（逐数据集独立处理，tryCatch 隔离，每数据集结束 gc()）】
#   读 RDS -> JoinLayers -> ScaleData(VST2000) -> PCA(50) -> Harmony
#   (group.by: 152805/167186=donor, 169454=gsm) -> Neighbors(harmony,1:30)
#   -> Clusters(res=0.6) -> UMAP(harmony,1:30)；set.seed(20260823)。
#   -> 经典标志物打分法自动注释（不跑 FindAllMarkers；得分矩阵写 CSV +
#      经典标志物 x cluster DotPlot 供人工复核）
#   -> 10 个模块 D v2 标志物 FeaturePlot/DotPlot/VlnPlot
#   -> 三轨 UCell 评分 + FeaturePlot/VlnPlot + per-donor 中位数
#   -> donor/样本级 pseudobulk（counts 层列分组求和 -> log2CPM）统计
#
# 【细胞类型标志物来源（学界公认图谱，自动注释仅供参考，以人工复核为准）】
#   软骨: Ji et al. 2019 Ann Rheum Dis (ARD)；
#   滑膜: Zhang et al. 2023 Nat Immunol (AMP RA/OA 图谱)；
#   骨骼肌: Lai et al. 2024 Nature (人类骨骼肌衰老图谱)。
#
# 【10 个模块 D v2 锁定标志物（来源：模块 D v2 多队列机器学习筛选，见
#   05b_moduleD_v2_multicohort.R 与备忘录）】
#   NELL1, STEAP1, GADD45A, FBLN5, HEMK1, GDE1, RNF14, BNIP3, EGFR, CBR3
#
# 【输出（全部落 results/moduleG/G2/，逐数据集子目录；幂等覆盖写）】
#   G2/<GSE>/  UMAP x3 / 标志物 FeaturePlot/DotPlot/VlnPlot / 三轨 UCell 图 /
#              经典标志物 x cluster DotPlot（一律 PDF+PNG）/
#              <GSE>_cluster_celltype_scores.csv / <GSE>_seurat_G2.rds /
#              <GSE>_pseudobulk_stats.csv
#   G2/marker_localization_summary.csv        （三数据集合并，顶层）
#   G2/moduleG_G2_pseudobulk_stats.csv        （三数据集合并，顶层）
#   G2/moduleG_G2_status.csv                  （每数据集一行）
#   results/moduleG/moduleG_G2_log_<date>.txt （沿用 G1 log_msg 风格）
#
# 【工程约束】无 system2/无交互；向 Seurat 对象赋元数据一律 unname()
#   （G1 fix3 教训：vapply/ifelse 命名污染触发 v5 零重叠校验）；
#   内存估算：GSE167186 约 5.2 万核 JoinLayers 后单对象 <4GB。
# ============================================================================

# ---- 路径设置（假定从项目根目录运行，与 G1 一致） --------------------------------
proj_root <- getwd()
res_dir   <- file.path(proj_root, "results", "moduleG")
g2_dir    <- file.path(res_dir, "G2")
dir.create(g2_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(res_dir, paste0("moduleG_G2_log_", Sys.Date(), ".txt"))

log_msg <- function(..., .sep = "") {
  txt <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ",
                paste(..., collapse = .sep))
  cat(txt, "\n")
  cat(txt, "\n", file = log_file, append = TRUE)
  invisible(txt)
}

log_msg("=== 模块 G G2 启动：标志物定位 / 三轨 UCell / pseudobulk 统计 ===")
log_msg("红线 1: GSE152805 无健康对照，仅供体内 oLT vs MT 配对对比（描述性）。")
log_msg("红线 2: GSE167186 为 Old vs Young aging proxy，严禁 sarcopenia 表述。")
log_msg("红线 3: 差异统计以 donor/样本级 pseudobulk 为底线；Harmony 仅用于聚类/UMAP。")
log_msg("proj_root = ", proj_root)

# ---- 依赖检查（缺包即 stop 并打印确切安装命令） -----------------------------------
check_deps <- function() {
  # (包名, 安装命令, 最低版本或 NA)
  reqs <- list(
    # v1.1 建议4：Seurat 下限抬到 5.3.0（配 ggplot2 4.x 兼容链）；
    # harmony >= 1.1；UCell >= 2.4；patchwork 存在性检查
    list(pkg = "Seurat",     min = "5.3.0",
         cmd = 'install.packages("Seurat")'),
    list(pkg = "Matrix",     min = NA,
         cmd = 'install.packages("Matrix")'),
    list(pkg = "data.table", min = NA,
         cmd = 'install.packages("data.table")'),
    list(pkg = "ggplot2",    min = "3.4.0",
         cmd = 'install.packages("ggplot2")'),
    list(pkg = "patchwork",  min = NA,
         cmd = 'install.packages("patchwork")'),
    list(pkg = "harmony",    min = "1.1",
         cmd = 'install.packages("harmony")'),
    list(pkg = "UCell",      min = "2.4",
         cmd = 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager"); BiocManager::install("UCell")')
  )
  problems <- character(0)
  for (r in reqs) {
    if (!requireNamespace(r$pkg, quietly = TRUE)) {
      problems <- c(problems, paste0("  缺包 ", r$pkg, " -> 安装: ", r$cmd))
      next
    }
    if (!is.na(r$min)) {
      v <- as.character(utils::packageVersion(r$pkg))
      if (utils::compareVersion(v, r$min) < 0) {
        problems <- c(problems, paste0("  ", r$pkg, " 版本 ", v, " < ", r$min,
                                       " -> 升级: ", r$cmd))
      }
    }
  }
  if (length(problems) > 0) {
    log_msg("ERROR: 依赖不满足:\n", paste(problems, collapse = "\n"))
    stop(paste0("G2 依赖不满足，请先执行：\n",
                paste(problems, collapse = "\n")), call. = FALSE)
  }
  log_msg("依赖检查通过: ",
          paste(vapply(reqs, function(r) {
            paste0(r$pkg, " ", as.character(utils::packageVersion(r$pkg)))
          }, character(1)), collapse = " | "))
  invisible(TRUE)
}
check_deps()

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

set.seed(20260823)

# ---- 10 个模块 D v2 锁定标志物（硬编码；来源见头部注释） ---------------------------
MARKERS10 <- c("NELL1", "STEAP1", "GADD45A", "FBLN5", "HEMK1",
               "GDE1", "RNF14", "BNIP3", "EGFR", "CBR3")

# ---- 细胞类型标志物配置块（来源见头部注释：Ji 2019 ARD / Zhang 2023 Nat Immunol /
#      Lai 2024 Nature；自动注释仅供初筛，以人工复核 DotPlot 为准） ------------------
MARKERS_UNIVERSAL <- list(
  Immune     = c("PTPRC"),
  Endothelial = c("PECAM1", "VWF")
)
MARKERS_GSE152805 <- c(MARKERS_UNIVERSAL, list(
  # 滑膜成纤维细胞亚群（lining/sublining 通用）
  SynFibro       = c("PRG4", "THY1", "CD34", "DPP4", "COL1A1"),
  Macrophage     = c("CD68", "CD14", "MRC1"),
  Tcell          = c("CD3D", "CD3E"),
  Bcell          = c("CD79A", "MS4A1"),
  Chondrocyte    = c("COL2A1", "ACAN", "SOX9"),
  HyperChondro   = c("COL10A1", "IHH", "MMP13"),
  FibroChondro   = c("COL1A1", "COL3A1")
))
MARKERS_GSE169454 <- c(MARKERS_UNIVERSAL, list(
  Chondrocyte    = c("COL2A1", "ACAN", "SOX9"),
  HyperChondro   = c("COL10A1", "IHH", "MMP13"),
  FibroChondro   = c("COL1A1", "COL3A1"),
  HomeoChondro   = c("PRG4", "FRZB")
))
MARKERS_GSE167186 <- c(MARKERS_UNIVERSAL, list(
  Myonuclei      = c("TTN", "NEB"),
  Myonuc_IIx     = c("MYH1"),
  Myonuc_IIa     = c("MYH2"),
  Myonuc_I       = c("MYH7"),
  Satellite      = c("PAX7"),
  FAPs           = c("PDGFRA", "LUM"),
  Pericyte_SMC   = c("RGS5", "ACTG2", "MYH11"),
  Adipocyte      = c("ADIPOQ")
))

# ---- 数据集配置 -------------------------------------------------------------------
# harmony_var: Harmony 分组变量（仅用于聚类/UMAP）；grp_var: 比较/展示分组；
# unit_var: pseudobulk 聚合单位（meta 列名，可多个）；comp: 统计比较定义
DS_CFG <- list(
  GSE152805 = list(
    rds = file.path(res_dir, "GSE152805_seurat_G1.rds"),
    harmony_var = "donor",
    grp_var     = "tissue",
    unit_vars   = c("donor", "tissue"),   # donor x tissue（=9 个 GSM）
    markers     = MARKERS_GSE152805,
    comp = list(groupA = "oLT", groupB = "MT", paired = TRUE, pair_var = "donor",
                bh = FALSE,
                # v1.1 建议6：删除 exclude_groups 死字段——SY 的排除已由
                # groupA/groupB 精确匹配隐式实现（stat_one_feature 只取两组）
                name = "oLT_vs_MT_paired",
                note = "paired n=3, descriptive only; SY excluded; no healthy control in GSE152805")
  ),
  GSE169454 = list(
    rds = file.path(res_dir, "GSE169454_seurat_G1.rds"),
    harmony_var = "gsm",
    grp_var     = "condition",
    unit_vars   = c("gsm"),               # 每 GSM 一个聚合单位（7 个）
    markers     = MARKERS_GSE169454,
    comp = list(groupA = "OA", groupB = "normal", paired = FALSE, pair_var = NULL,
                bh = FALSE,
                name = "OA_vs_normal",
                note = "unpaired n=7 (OA4/normal3), small-sample exploratory")
  ),
  GSE167186 = list(
    rds = file.path(res_dir, "GSE167186_seurat_G1.rds"),
    harmony_var = "donor",
    grp_var     = "condition",
    unit_vars   = c("donor"),             # donor=HM 编号（17 个）
    markers     = MARKERS_GSE167186,
    comp = list(groupA = "Old", groupB = "Young", paired = FALSE, pair_var = NULL,
                bh = TRUE,
                name = "Old_vs_Young",
                note = "unpaired n=17 (Old11/Young6); aging proxy - NOT sarcopenia evidence; BH within feature_type")
  )
)

# ---- 载入模块 C 三轨表 + 10 标志物轨归属（缺基因即 stop） ---------------------------
load_gene_tracks <- function() {
  f <- file.path(proj_root, "results", "moduleC", "moduleC_gene_tracks.csv")
  if (!file.exists(f)) {
    log_msg("ERROR: 未找到三轨表: ", f)
    stop("缺少 results/moduleC/moduleC_gene_tracks.csv，请先确认模块 C 产物。",
         call. = FALSE)
  }
  tr <- tryCatch(
    as.data.frame(data.table::fread(f, stringsAsFactors = FALSE)),
    error = function(e) {
      stop(paste0("moduleC_gene_tracks.csv 读取失败: ", conditionMessage(e)),
           call. = FALSE)
    })
  need <- c("gene", "track", "direction_OA", "direction_muscle")
  miss <- setdiff(need, names(tr))
  if (length(miss) > 0) {
    stop("moduleC_gene_tracks.csv 缺列: ", paste(miss, collapse = ", "),
         call. = FALSE)
  }
  bad_track <- setdiff(unique(tr$track),
                       c("concordant", "mirror_OAup_muscleDown",
                         "mirror_OAdown_muscleUp"))
  if (length(bad_track) > 0) {
    log_msg("WARN: 三轨表出现预期外 track 取值: ",
            paste(bad_track, collapse = ", "), "（这些行不参与轨评分）")
  }
  # v1.1 建议1：剔除 gene/track 为 NA 的行，杜绝 NA 行静默污染下游查表
  n_before_na <- nrow(tr)
  tr <- tr[!is.na(tr$track) & !is.na(tr$gene), ]
  if (nrow(tr) < n_before_na) {
    log_msg("WARN: 三轨表剔除 gene/track 含 NA 的行 ", n_before_na - nrow(tr),
            " 行（剩余 ", nrow(tr), " 行）。")
  }
  # 10 标志物轨归属必须全部查到，否则 stop 并报缺哪些
  miss10 <- setdiff(MARKERS10, tr$gene)
  if (length(miss10) > 0) {
    log_msg("ERROR: 以下模块 D v2 标志物在三轨表中查无轨归属: ",
            paste(miss10, collapse = ", "))
    stop(paste0("三轨表缺少 ", length(miss10), " 个锁定标志物: ",
                paste(miss10, collapse = ", "),
                "；请核对 moduleC_gene_tracks.csv 后重跑。"), call. = FALSE)
  }
  log_msg("三轨表载入: ", nrow(tr), " 行；10 个锁定标志物轨归属全部查到。")
  tr
}
TRACKS_DF <- load_gene_tracks()
TRACKS <- c("concordant", "mirror_OAup_muscleDown", "mirror_OAdown_muscleUp")
TRACK_GENES <- lapply(TRACKS, function(tk) TRACKS_DF$gene[TRACKS_DF$track == tk])
names(TRACK_GENES) <- TRACKS
log_msg("三轨基因数: ",
        paste(sprintf("%s=%d", TRACKS, lengths(TRACK_GENES)), collapse = " | "))

# ---- 绘图工具（逐图 tryCatch；PDF+PNG 双格式） -------------------------------------
save_plot_both <- function(p, base_name, width = 10, height = 7) {
  ggplot2::ggsave(paste0(base_name, ".pdf"), p, width = width, height = height,
                  limitsize = FALSE)
  ggplot2::ggsave(paste0(base_name, ".png"), p, width = width, height = height,
                  dpi = 300, bg = "white", limitsize = FALSE)
  invisible(TRUE)
}

safe_plot <- function(tag, plot_expr, base_name, width = 10, height = 7) {
  tryCatch({
    p <- plot_expr
    save_plot_both(p, base_name, width, height)
    log_msg("    图已写出: ", basename(base_name), ".pdf/.png")
    TRUE
  }, error = function(e) {
    log_msg("WARN: 绘图失败 [", tag, "]（不中断）: ", conditionMessage(e))
    FALSE
  })
}

# ---- pseudobulk 工具 ----------------------------------------------------------------
# counts 层（JoinLayers 后）按列分组求和：genes x units 稀疏矩阵
aggregate_counts <- function(seu, unit_vec, gse = "?") {
  counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
  unit <- as.character(unit_vec)
  # v1.1 建议2：聚合单位含 NA 的细胞会被 sparse.model.matrix 静默丢弃，
  # 必须显式 stop 并报数据集与 NA 数
  if (anyNA(unit)) {
    stop(gse, ": pseudobulk 聚合单位向量含 ", sum(is.na(unit)),
         " 个 NA（对应细胞会被静默丢弃）；请检查 meta 列完整性。", call. = FALSE)
  }
  units <- sort(unique(unit))
  mm <- Matrix::sparse.model.matrix(~ 0 + factor(unit, levels = units))
  colnames(mm) <- units
  pb <- counts %*% mm          # dgCMatrix: genes x units
  colnames(pb) <- units
  pb
}

# log2CPM（library size 归一 + log2(x+1)）
# v1.1 建议3：先转稠密矩阵再显式 t(t(.)/lib)，不依赖 Matrix 稀疏对象的
# 隐式向量回收语义（pseudobulk 规模 = 基因数 x 单位数，稠密化代价可忽略）
to_log2cpm <- function(pb) {
  lib <- colSums(pb)
  if (any(lib <= 0)) {
    log_msg("WARN: pseudobulk 存在 library size=0 的聚合单位: ",
            paste(colnames(pb)[lib <= 0], collapse = ", "),
            "（其 log2CPM 将全 0，该单位比较时慎用）")
  }
  pb2 <- as.matrix(pb)
  log2(t(t(pb2) / lib) * 1e6 + 1)
}

# 单特征两组 Wilcoxon（paired 时按 pair_var 对齐）；全部异常退回 p=NA 并注明
stat_one_feature <- function(values, groups, donors, comp) {
  idxA <- which(groups == comp$groupA)
  idxB <- which(groups == comp$groupB)
  vA <- values[idxA]; vB <- values[idxB]
  nA <- length(vA);    nB <- length(vB)
  medA <- if (nA > 0) stats::median(vA, na.rm = TRUE) else NA_real_
  medB <- if (nB > 0) stats::median(vB, na.rm = TRUE) else NA_real_
  eff  <- medA - medB           # marker: log2CPM 差（≈log2FC）；score: 中位数差
  p <- NA_real_; wnote <- ""
  ok <- nA >= 2 && nB >= 2 &&
        length(unique(c(vA, vB))) > 1 &&
        !any(is.na(c(vA, vB)))
  if (ok && comp$paired) {
    ok <- nA == nB              # 配对要求两组单位数相等
    if (!ok) wnote <- "paired groups size mismatch; p=NA"
  }
  if (ok) {
    p <- tryCatch({
      if (comp$paired) {
        oA <- order(donors[idxA]); oB <- order(donors[idxB])
        if (!all(donors[idxA][oA] == donors[idxB][oB])) {
          wnote <- "pair donors mismatch; p=NA"; NA_real_
        } else {
          stats::wilcox.test(vA[oA], vB[oB], paired = TRUE, exact = FALSE)$p.value
        }
      } else {
        stats::wilcox.test(vA, vB, paired = FALSE, exact = FALSE)$p.value
      }
    }, error = function(e) {
      wnote <<- paste0("wilcox error: ", conditionMessage(e), "; p=NA")
      NA_real_
    })
    # v1.1 建议9：wilcox 可能返回 NaN（如配对差值全 0），统一归一为 NA 并注明
    if (!is.finite(p)) {
      p <- NA_real_
      wnote <- trimws(paste(wnote, "; all-zero differences, p undefined"))
    }
  } else if (!nzchar(wnote)) {
    wnote <- "insufficient/degenerate values; p=NA"
  }
  list(nA = nA, nB = nB, medA = medA, medB = medB, effect = eff,
       p = p, note = wnote)
}

# ---- 细胞类型自动注释（经典标志物打分法；不跑 FindAllMarkers） ----------------------
# 对每套候选细胞类型取其标志物（∩ 实际检出基因）在 data 层的细胞级平均表达，
# 再按 cluster 取均值 -> cluster x celltype 得分矩阵（必须写 CSV 留痕）；
# cluster 标记得分最高者，写回 meta 列 cell_type（赋 meta 一律 unname()）。
annotate_clusters <- function(seu, markers_list, gse, out_dir) {
  data_mat <- GetAssayData(seu, assay = "RNA", layer = "data")
  clusters <- as.character(seu$seurat_clusters)
  # v1.1 建议7：cluster 按数值排序（0,1,2,...,10），不用字典序（0,1,10,2...）
  cl_levels <- as.character(sort(as.numeric(unique(clusters))))
  detected <- rownames(data_mat)
  score_mat <- matrix(NA_real_, nrow = length(cl_levels),
                      ncol = length(markers_list),
                      dimnames = list(cl_levels, names(markers_list)))
  for (ct in names(markers_list)) {
    g <- intersect(markers_list[[ct]], detected)
    if (length(g) == 0) {
      log_msg("WARN: ", gse, " 细胞类型 ", ct, " 的标志物全部未检出，",
              "得分列记 NA（cluster 归属时该类型不可胜出）。")
      next
    }
    if (length(g) < length(markers_list[[ct]])) {
      log_msg("WARN: ", gse, " ", ct, " 标志物部分检出 ", length(g), "/",
              length(markers_list[[ct]]), "（未检出: ",
              paste(setdiff(markers_list[[ct]], g), collapse = ","), "）")
    }
    cell_score <- colMeans(data_mat[g, , drop = FALSE])
    score_mat[, ct] <- tapply(cell_score, factor(clusters, levels = cl_levels),
                              mean)
  }
  # fix2：列归一化后再 argmax。原始得分为 cluster 内标志物平均表达，
  # 高丰度标志物（如骨骼肌 TTN/NEB/MYH）会淹没少数群体细胞类型
  # （实机：GSE167186 12 个 cluster 被判成 1 类）。按列（细胞类型）除以
  # 该类在所有 cluster 的最大值压到 [0,1]，使"相对自身最高表达的
  # cluster"胜出，而非绝对表达量最高的类型通吃。
  # 边界：列全 NA 或 max<=0 的列保持原值，永不胜出（归一化列最大值恒为 1）。
  score_norm <- score_mat
  for (j in seq_len(ncol(score_norm))) {
    mx <- max(score_norm[, j], na.rm = TRUE)
    if (is.finite(mx) && mx > 0) score_norm[, j] <- score_norm[, j] / mx
  }
  # fix3：归一化 argmax 的两个补丁（实机 fix2 暴露）——
  # ① 并列决胜按原始得分：max 归一化下多列常并列 1.0，列序决胜会错标
  #    （GSE167186 cl8 内皮 raw 0.65 被 Immune raw 0.21 抢占）；
  # ② 原始得分下限 0.05：近零表达类型会被归一化放大成噪声标签
  #    （GSE167186 cl11 Adipocyte raw 0.0057 自立门户），先剔除再决胜；
  #    全部低于下限时标 unknown（诚实留空，不硬贴标签）。
  RAW_FLOOR <- 0.05
  assigned <- vapply(seq_len(nrow(score_norm)), function(i) {
    r  <- score_norm[i, ]
    rr <- score_mat[i, ]
    ok <- which(!is.na(r) & !is.na(rr) & rr >= RAW_FLOOR)
    if (length(ok) == 0) return("unknown")
    mx <- max(r[ok])
    tied <- ok[r[ok] >= mx - 1e-9]
    if (length(tied) == 1L) return(colnames(score_norm)[tied])
    colnames(score_norm)[tied[which.max(rr[tied])]]
  }, character(1))
  # fix3 留痕：与朴素归一化 argmax（fix2 口径）不同的 cluster 写日志
  naive <- apply(score_norm, 1, function(r) {
    if (all(is.na(r))) return("unknown")
    colnames(score_norm)[which.max(r)]
  })
  diff_idx <- which(assigned != naive)
  if (length(diff_idx) > 0) {
    log_msg(gse, ": fix3 并列决胜/下限改标 ", length(diff_idx), " 个 cluster: ",
            paste(sprintf("%s %s->%s", cl_levels[diff_idx], naive[diff_idx],
                          assigned[diff_idx]), collapse = "; "))
  }
  ct_map <- stats::setNames(assigned, cl_levels)
  seu$cell_type <- unname(ct_map[clusters])   # unname() 防 v5 零重叠校验
  log_msg(gse, ": 自动注释完成，", length(cl_levels), " 个 cluster -> ",
          length(unique(assigned)), " 种细胞类型 (",
          paste(sprintf("%s=%d", names(table(assigned)),
                        as.integer(table(assigned))), collapse = ", "), ")")

  # 得分矩阵写 CSV 留痕（原始得分 + fix2 归一化得分双口径，供人工复核）
  norm_mat <- score_norm
  colnames(norm_mat) <- paste0(colnames(score_norm), "_norm")
  score_df <- data.frame(cluster = cl_levels, score_mat, norm_mat,
                         assigned_cell_type = unname(assigned),
                         check.names = FALSE, stringsAsFactors = FALSE)
  utils::write.csv(score_df,
                   file.path(out_dir, paste0(gse, "_cluster_celltype_scores.csv")),
                   row.names = FALSE, fileEncoding = "UTF-8")
  log_msg("    写出: ", gse, "_cluster_celltype_scores.csv (",
          nrow(score_df), " 行)")

  # 经典标志物 x cluster DotPlot（供人工复核）
  classic <- unique(unlist(markers_list, use.names = FALSE))
  classic <- intersect(classic, detected)
  if (length(classic) > 0) {
    safe_plot(paste(gse, "classic DotPlot"),
              DotPlot(seu, features = classic, group.by = "seurat_clusters") +
                theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                                 size = 7)) +
                labs(title = paste0(gse, " classic markers x cluster",
                                    " (annotation reference)")),
              file.path(out_dir, paste0(gse, "_classic_markers_DotPlot_cluster")),
              width = max(8, 0.35 * length(classic) + 3),
              height = max(5, 0.35 * length(cl_levels) + 2.5))
  }
  seu
}

# ---- 三轨 UCell 评分 -----------------------------------------------------------------
# 输入 RNA assay data 层（JoinLayers 后）；逐轨记录检出覆盖度；检出 <3 的轨跳过
# 并 WARN。返回 list(seu, track_cols(成功评分的轨->meta列名), coverage(字符向量))
run_ucell <- function(seu, gse, out_dir, grp_var) {
  detected <- rownames(seu)
  sigs <- list()
  coverage <- character(0)
  for (tk in TRACKS) {
    det <- intersect(TRACK_GENES[[tk]], detected)
    coverage[tk] <- paste0(length(det), "/", length(TRACK_GENES[[tk]]))
    log_msg(gse, ": 轨 ", tk, " 检出覆盖度 ", coverage[tk])
    if (length(det) < 3) {
      log_msg("WARN: ", gse, " 轨 ", tk, " 检出 ", length(det),
              " < 3 个基因，跳过该轨评分。")
      next
    }
    sigs[[tk]] <- det
  }
  track_cols <- character(0)
  if (length(sigs) > 0) {
    # fix1：UCell 新版（≥2.10，实机 2.16）的 ScoreSignatures_UCell 已移除
    # Seurat 对象分支（只接受 SCE/matrix/dgCMatrix/data.frame，否则报
    # "Unrecognized input format."）；Seurat 接口为 AddModuleScore_UCell，
    # 其内部按 v5 layer 处理（Layers(search=slot)）并 AddMetaData 回对象。
    # UCell 基于细胞内秩次，counts 与 data 层秩次一致（细胞内经单调变换），
    # 用 data 层（G2 规格）。
    # name="" 使 meta 列名 = 轨名；赋 meta 由 UCell 内部完成
    seu <- UCell::AddModuleScore_UCell(seu, features = sigs, assay = "RNA",
                                       slot = "data", name = "")
    for (tk in names(sigs)) {
      col <- if (tk %in% colnames(seu@meta.data)) tk else paste0(tk, "_UCell")
      if (!col %in% colnames(seu@meta.data)) {
        log_msg("WARN: ", gse, " 轨 ", tk, " 评分列未找到（试过 ", tk, " 与 ",
                tk, "_UCell），该轨视为失败。")
        next
      }
      track_cols[tk] <- col
      # FeaturePlot + 按分组 VlnPlot
      safe_plot(paste(gse, tk, "UCell FeaturePlot"),
                FeaturePlot(seu, features = col) +
                  labs(title = paste0(gse, " UCell: ", tk,
                                      " (coverage ", coverage[tk], ")")),
                file.path(out_dir, paste0(gse, "_UCell_", tk, "_FeaturePlot")),
                width = 7.5, height = 6.5)
      safe_plot(paste(gse, tk, "UCell VlnPlot"),
                VlnPlot(seu, features = col, group.by = grp_var, pt.size = 0) +
                  labs(title = paste0(gse, " UCell: ", tk, " by ", grp_var)),
                file.path(out_dir,
                          paste0(gse, "_UCell_", tk, "_VlnPlot_", grp_var)),
                width = 7, height = 5.5)
    }
  }
  log_msg(gse, ": UCell 完成，成功评分轨 ",
          paste(sprintf("%s(%s)", names(track_cols), coverage[names(track_cols)]),
                collapse = ", "))
  list(seu = seu, track_cols = track_cols, coverage = coverage)
}

# ---- 10 标志物定位汇总行（marker_localization_summary） ------------------------------
# data 层逐 cell_type：n_cells_in_type / pct_expr(%>0) / avg_expr(均值)；
# 未检出基因写一行 pct_expr=NA + note="not detected in dataset"
localization_rows <- function(seu, gse) {
  data_mat <- GetAssayData(seu, assay = "RNA", layer = "data")
  ct <- as.character(seu$cell_type)
  ct_levels <- sort(unique(ct))
  detected <- rownames(data_mat)
  rows <- lapply(MARKERS10, function(gene) {
    # v1.1 建议1：match 取行（全局预检已保证 10 基因全部命中，match 不会 NA；
    # 仍显式防御，杜绝逻辑子集遇到 NA 时的静默污染）
    hit <- match(gene, TRACKS_DF$gene)
    if (is.na(hit)) {
      stop(gse, " 标志物 ", gene, " 在三轨表查无行（预检后不应发生）。",
           call. = FALSE)
    }
    tr <- TRACKS_DF[hit, , drop = FALSE]
    base <- data.frame(gse = gse, gene = gene,
                       track = as.character(tr$track),
                       direction_OA = as.character(tr$direction_OA),
                       direction_muscle = as.character(tr$direction_muscle),
                       stringsAsFactors = FALSE)
    if (!gene %in% detected) {
      return(cbind(base, data.frame(
        cell_type = NA_character_, n_cells_in_type = NA_integer_,
        pct_expr = NA_real_, avg_expr = NA_real_,
        note = "not detected in dataset", stringsAsFactors = FALSE)))
    }
    expr <- as.numeric(data_mat[gene, ])
    per_ct <- lapply(ct_levels, function(c) {
      idx <- ct == c
      cbind(base, data.frame(
        cell_type = c,
        n_cells_in_type = sum(idx),
        pct_expr = round(100 * mean(expr[idx] > 0), 2),
        avg_expr = round(mean(expr[idx]), 4),
        note = "", stringsAsFactors = FALSE))
    })
    do.call(rbind, per_ct)
  })
  do.call(rbind, rows)
}

# ---- donor/样本级 pseudobulk 统计（底线：严禁细胞级 p 值） ----------------------------
# counts 层（JoinLayers 后）按聚合单位列求和 -> log2CPM；
# 10 标志物逐基因 + 三轨 UCell per-unit 中位数做两组 Wilcoxon；
# BH 仅对 cfg$comp$bh=TRUE 的数据集在 feature_type 内分别校正。
run_pseudobulk_stats <- function(seu, gse, cfg, track_cols, out_dir) {
  md  <- seu@meta.data
  comp <- cfg$comp
  unit <- apply(md[, cfg$unit_vars, drop = FALSE], 1, paste, collapse = "__")
  log_msg(gse, ": pseudobulk 聚合单位 = ", paste(cfg$unit_vars, collapse = "x"),
          "，共 ", length(unique(unit)), " 个单位。")
  pb   <- aggregate_counts(seu, unit, gse)
  lcpm <- to_log2cpm(pb)
  units <- colnames(pb)

  # 每个聚合单位的分组与 donor（取该单位首行 meta）
  umeta <- do.call(rbind, lapply(units, function(u) {
    r <- md[unit == u, , drop = FALSE][1, ]
    data.frame(unit  = u,
               group = as.character(r[[cfg$grp_var]]),
               donor = as.character(r$donor),
               stringsAsFactors = FALSE)
  }))
  groups <- umeta$group
  donors <- umeta$donor
  log_msg(gse, ": 单位分组分布: ",
          paste(sprintf("%s=%d", names(table(groups)), as.integer(table(groups))),
                collapse = ", "))

  mk_row <- function(feature, ftype, values) {
    st <- stat_one_feature(values, groups, donors, comp)
    data.frame(gse = gse, comparison = comp$name,
               feature = feature, feature_type = ftype,
               n_groupA = st$nA, n_groupB = st$nB,
               median_A = round(st$medA, 4), median_B = round(st$medB, 4),
               log2FC_or_delta = round(st$effect, 4),
               p_value = st$p, p_adj_BH = NA_real_,
               note = trimws(paste(comp$note, st$note)),
               stringsAsFactors = FALSE)
  }

  rows <- list()
  # 10 标志物（log2CPM 值）
  for (gene in MARKERS10) {
    if (!gene %in% rownames(lcpm)) {
      r <- mk_row(gene, "marker", rep(NA_real_, length(units)))
      r$note <- trimws(paste(r$note, "; not detected in dataset"))
      rows[[length(rows) + 1L]] <- r
      next
    }
    rows[[length(rows) + 1L]] <-
      mk_row(gene, "marker", as.numeric(lcpm[gene, units]))
  }
  # 三轨 UCell per-unit 中位数（同时缓存供 per-unit 值表落盘）
  track_medians <- list()
  for (tk in names(track_cols)) {
    sc <- seu@meta.data[[track_cols[tk]]]
    med_per_unit <- as.numeric(tapply(sc, unit, stats::median, na.rm = TRUE)[units])
    track_medians[[tk]] <- med_per_unit
    rows[[length(rows) + 1L]] <- mk_row(tk, "track_score", med_per_unit)
  }
  df <- do.call(rbind, rows)

  # v1.1 建议5：per-unit 值表落盘（10 标志物 log2CPM + 三轨 UCell per-unit
  # 中位数），供复核复算
  val_df <- data.frame(unit = units, group = groups, donor = donors,
                       stringsAsFactors = FALSE)
  for (gene in MARKERS10) {
    val_df[[paste0("log2CPM_", gene)]] <-
      if (gene %in% rownames(lcpm)) {
        round(as.numeric(lcpm[gene, units]), 4)
      } else NA_real_
  }
  for (tk in names(track_medians)) {
    val_df[[paste0("UCellMedian_", tk)]] <- round(track_medians[[tk]], 4)
  }
  utils::write.csv(val_df, file.path(out_dir, paste0(gse, "_per_unit_values.csv")),
                   row.names = FALSE, fileEncoding = "UTF-8")
  log_msg("    写出: ", gse, "_per_unit_values.csv (", nrow(val_df), " 行；",
          "供复核复算)")

  # BH 校正（仅 cfg$comp$bh=TRUE；标志物与轨评分分别进行）
  if (isTRUE(comp$bh)) {
    for (ft in unique(df$feature_type)) {
      idx <- which(df$feature_type == ft & !is.na(df$p_value))
      if (length(idx) > 0) {
        df$p_adj_BH[idx] <- stats::p.adjust(df$p_value[idx], method = "BH")
      }
    }
    log_msg(gse, ": BH 校正完成（feature_type 内分别进行）。")
  }
  utils::write.csv(df, file.path(out_dir, paste0(gse, "_pseudobulk_stats.csv")),
                   row.names = FALSE, fileEncoding = "UTF-8")
  log_msg("    写出: ", gse, "_pseudobulk_stats.csv (", nrow(df), " 行；",
          "单位级 Wilcoxon，无细胞级 p 值)")
  df
}

# ---- 单数据集 G2 主管线 ---------------------------------------------------------------
process_gse_g2 <- function(gse, cfg) {
  log_msg("---- 数据集: ", gse, " ----")
  out_dir <- file.path(g2_dir, gse)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  if (!file.exists(cfg$rds)) {
    stop(gse, " G1 RDS 不存在: ", cfg$rds, "；请先成功运行 G1。", call. = FALSE)
  }
  seu <- readRDS(cfg$rds)
  log_msg(gse, ": RDS 载入 ", ncol(seu), " 细胞 x ", nrow(seu), " 基因。")

  # 1) JoinLayers（merge 后 counts/data 按 GSM 分层，合并为单一层）
  seu <- JoinLayers(seu)
  log_msg(gse, ": JoinLayers 完成，layers=",
          paste(SeuratObject::Layers(seu), collapse = ","))

  # 2) 归一化 -> PCA -> Harmony -> 聚类 -> UMAP
  vf <- VariableFeatures(seu)
  if (length(vf) == 0) stop(gse, " G1 对象无 VariableFeatures，无法 ScaleData/PCA。",
                            call. = FALSE)
  # v1.1 建议10：VST2000 存活检查——JoinLayers 后高变基因异常少则预警
  # （可能是 G1 RDS 损坏或 assay 结构异常），<1000 时 WARN（不 stop）
  if (length(vf) < 1000) {
    log_msg("WARN: ", gse, " JoinLayers 后 VariableFeatures 仅 ", length(vf),
            " 个（预期 ~2000），VST2000 存活异常，请检查 G1 RDS 完整性。")
  }
  seu <- ScaleData(seu, features = vf, verbose = FALSE)
  seu <- RunPCA(seu, npcs = 50, seed.use = 20260823, verbose = FALSE)
  if (!cfg$harmony_var %in% colnames(seu@meta.data)) {
    stop(gse, " meta 缺 Harmony 分组列: ", cfg$harmony_var, call. = FALSE)
  }
  # Harmony 仅用于聚类/UMAP，严禁作为差异或评分输入（红线 3）
  seu <- harmony::RunHarmony(seu, group.by.vars = cfg$harmony_var,
                             reduction = "pca", dims.use = 1:30,
                             reduction.save = "harmony", verbose = FALSE)
  seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:30, verbose = FALSE)
  seu <- FindClusters(seu, resolution = 0.6, verbose = FALSE)
  Idents(seu) <- "seurat_clusters"
  seu <- RunUMAP(seu, reduction = "harmony", dims = 1:30,
                 seed.use = 20260823, verbose = FALSE)
  n_clusters <- length(unique(seu$seurat_clusters))
  log_msg(gse, ": 聚类完成，", n_clusters, " 个 cluster (Harmony by ",
          cfg$harmony_var, ")。")

  # UMAP 三张：cluster / 比较分组 / donor
  for (gv in c("seurat_clusters", cfg$grp_var, "donor")) {
    safe_plot(paste(gse, "UMAP", gv),
              DimPlot(seu, group.by = gv, label = (gv == "seurat_clusters")) +
                labs(title = paste0(gse, " UMAP by ", gv,
                                    " (harmony, dims 1:30)")),
              file.path(out_dir, paste0(gse, "_UMAP_", gv)),
              width = 7.5, height = 6.5)
  }

  # 3) 细胞类型自动注释（得分矩阵 CSV + 经典 DotPlot 在函数内输出）
  seu <- annotate_clusters(seu, cfg$markers, gse, out_dir)

  # 4) 10 个模块 D v2 标志物
  det10  <- intersect(MARKERS10, rownames(seu))
  miss10 <- setdiff(MARKERS10, rownames(seu))
  if (length(miss10) > 0) {
    log_msg("WARN: ", gse, " 未检出标志物 ", length(miss10), " 个: ",
            paste(miss10, collapse = ", "),
            "（summary 中写 not detected 行）")
  }
  log_msg(gse, ": 标志物检出 ", length(det10), "/10。")
  if (length(det10) > 0) {
    safe_plot(paste(gse, "markers FeaturePlot"),
              # 多子图 patchwork：labs() 只会作用于最后一个子图，
              # 总标题必须用 patchwork::plot_annotation（patchwork 随 Seurat 安装）
              FeaturePlot(seu, features = det10, ncol = 3) +
                patchwork::plot_annotation(
                  title = paste0(gse, " moduleD-v2 markers (detected ",
                                 length(det10), "/10)")),
              file.path(out_dir, paste0(gse, "_markers_FeaturePlot")),
              width = 13, height = 3.6 * ceiling(length(det10) / 3) + 1)
    safe_plot(paste(gse, "markers DotPlot"),
              DotPlot(seu, features = det10, group.by = "cell_type") +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                labs(title = paste0(gse, " moduleD-v2 markers x cell_type")),
              file.path(out_dir, paste0(gse, "_markers_DotPlot_celltype")),
              width = max(7, 0.5 * length(det10) + 3),
              height = max(4.5, 0.4 * length(unique(seu$cell_type)) + 2))
    # VlnPlot 按比较分组（152805: tissue 含 SY 仅展示，统计不含 SY）
    safe_plot(paste(gse, "markers VlnPlot"),
              VlnPlot(seu, features = det10, group.by = cfg$grp_var,
                      pt.size = 0, ncol = 2) +
                patchwork::plot_annotation(
                  title = paste0(gse, " moduleD-v2 markers by ", cfg$grp_var)),
              file.path(out_dir,
                        paste0(gse, "_markers_VlnPlot_", cfg$grp_var)),
              width = 11, height = 3.2 * ceiling(length(det10) / 2) + 1)
  }

  # 5) 三轨 UCell
  uc  <- run_ucell(seu, gse, out_dir, cfg$grp_var)
  seu <- uc$seu

  # 6) donor/样本级 pseudobulk 统计
  stat_df <- run_pseudobulk_stats(seu, gse, cfg, uc$track_cols, out_dir)

  # 7) 定位汇总行
  loc_df <- localization_rows(seu, gse)

  # 8) 存 G2 RDS
  rds_g2 <- file.path(out_dir, paste0(gse, "_seurat_G2.rds"))
  saveRDS(seu, rds_g2)
  log_msg(gse, ": G2 RDS 已保存: ", rds_g2, " (",
          round(file.size(rds_g2) / 1e6, 1), " MB)")

  list(
    stat = stat_df,
    loc  = loc_df,
    status = data.frame(
      gse = gse, status = "OK",
      n_cells = ncol(seu),
      n_clusters = n_clusters,
      n_celltypes = length(unique(seu$cell_type)),
      markers_detected = paste0(length(det10), "/10"),
      ucell_coverage = paste(sprintf("%s=%s", names(uc$coverage), uc$coverage),
                             collapse = "; "),
      rds_g2 = rds_g2, out_dir = out_dir,
      note = paste0("OK; harmony_by=", cfg$harmony_var,
                    "; comparison=", cfg$comp$name,
                    if (gse == "GSE167186")
                      "; aging proxy - NOT sarcopenia" else ""),
      stringsAsFactors = FALSE)
  )
}

# ---- 主流程：逐数据集独立处理，tryCatch 隔离，每数据集结束 gc() -----------------------
status_rows <- list()
stat_all    <- list()
loc_all     <- list()

run_gse <- function(gse) {
  cfg <- DS_CFG[[gse]]
  tryCatch(
    process_gse_g2(gse, cfg),
    error = function(e) {
      log_msg("ERROR: ", gse, " G2 处理失败（不影响其他数据集）: ",
              conditionMessage(e))
      list(stat = NULL, loc = NULL,
           status = data.frame(
             gse = gse, status = "FAIL",
             n_cells = NA_integer_, n_clusters = NA_integer_,
             n_celltypes = NA_integer_, markers_detected = NA_character_,
             ucell_coverage = NA_character_, rds_g2 = NA_character_,
             out_dir = file.path(g2_dir, gse),
             note = paste0("FAIL: ", conditionMessage(e)),
             stringsAsFactors = FALSE))
    })
}

for (gse in names(DS_CFG)) {
  res <- run_gse(gse)
  status_rows[[length(status_rows) + 1L]] <- res$status
  if (!is.null(res$stat)) stat_all[[length(stat_all) + 1L]] <- res$stat
  if (!is.null(res$loc))  loc_all[[length(loc_all) + 1L]]  <- res$loc
  rm(res)
  gc(verbose = FALSE)
  log_msg("---- ", gse, " 处理完毕，已 gc() ----")
}

# ---- 汇总输出 -------------------------------------------------------------------------
status_df <- do.call(rbind, status_rows)
tryCatch({
  utils::write.csv(status_df, file.path(g2_dir, "moduleG_G2_status.csv"),
                   row.names = FALSE, fileEncoding = "UTF-8")
  log_msg("写出: G2/moduleG_G2_status.csv (", nrow(status_df), " 行)")
}, error = function(e) log_msg("ERROR: 写 status CSV 失败: ", conditionMessage(e)))

tryCatch({
  if (length(stat_all) > 0) {
    stat_df <- do.call(rbind, stat_all)
    utils::write.csv(stat_df,
                     file.path(g2_dir, "moduleG_G2_pseudobulk_stats.csv"),
                     row.names = FALSE, fileEncoding = "UTF-8")
    log_msg("写出: G2/moduleG_G2_pseudobulk_stats.csv (", nrow(stat_df), " 行)")
  } else {
    log_msg("WARN: 无任何 pseudobulk 统计行可写出（所有数据集均失败）。")
  }
}, error = function(e) log_msg("ERROR: 写 pseudobulk 汇总失败: ", conditionMessage(e)))

tryCatch({
  if (length(loc_all) > 0) {
    loc_df <- do.call(rbind, loc_all)
    utils::write.csv(loc_df,
                     file.path(g2_dir, "marker_localization_summary.csv"),
                     row.names = FALSE, fileEncoding = "UTF-8")
    log_msg("写出: G2/marker_localization_summary.csv (", nrow(loc_df), " 行；",
            "三数据集合并)")
  } else {
    log_msg("WARN: 无任何定位汇总行可写出（所有数据集均失败）。")
  }
}, error = function(e) log_msg("ERROR: 写定位汇总失败: ", conditionMessage(e)))

# ---- 完成提示 + sessionInfo 留痕 --------------------------------------------------------
log_msg("==== G2 汇总 ====")
for (i in seq_len(nrow(status_df))) {
  log_msg(status_df$gse[i], ": ", status_df$status[i],
          " | cells=", status_df$n_cells[i],
          " | clusters=", status_df$n_clusters[i],
          " | celltypes=", status_df$n_celltypes[i],
          " | markers=", status_df$markers_detected[i],
          " | ucell=", status_df$ucell_coverage[i],
          " | ", status_df$note[i])
}

log_msg("==== sessionInfo ====")
for (ln in capture.output(utils::sessionInfo())) log_msg(ln)

if (all(status_df$status == "OK")) {
  msg <- paste0("模块 G G2 完成。红线提醒：GSE152805 仅供体内 oLT vs MT 描述性配对对比；",
                "GSE167186 为 Old vs Young aging proxy 不得写成 sarcopenia；",
                "一切统计为 donor/样本级 pseudobulk，自动注释需人工复核 DotPlot。")
  cat(msg, "\n"); log_msg(msg)
} else {
  bad <- status_df[status_df$status != "OK", ]
  msg <- paste0("模块 G G2 存在失败数据集：\n",
                paste0("  - ", bad$gse, " | ", bad$note, collapse = "\n"),
                "\n修复后重跑本脚本即可（全部产物幂等覆盖写）。")
  cat(msg, "\n"); log_msg(msg)
  # 不 stop：已成功数据集的产物均已落盘，status CSV 已留痕
}
