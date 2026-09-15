# ============================================================================
# 14c_moduleK_CellChat_分组标签与presto修正_20260914.R（完整替代 14/14b，直接重跑即可）
# 修正（相对 14b）：①GSE169454 对照标签实为 normal（非 control）；②identifyOverExpressedGenes 加 do.fast=FALSE 退回标准 wilcoxon（presto 非必需）。// 沿用 14b：G1/G2 RDS 自动鉴别——按文件名优先级（G2/annot/celltype）排序后逐一加载
#   验证元数据（须同时存在细胞类型列与含配置分组值的分组列），不合格自动试下一个；
#   ②分组列按「实际含有 CONFIG 配置分组值」选取，不再死板取第一候选。
# 其余口径与 14 完全一致（三数据集分层 / CellChat 2.x 标准管线 / 种子显式 / 红线复述）。
# ============================================================================

suppressPackageStartupMessages({ library(data.table) })

# ----------------------------- CONFIG --------------------------------------
CONFIG <- list(
  G2_DIR   = "results/moduleG",
  OUT_DIR  = "results/moduleK",
  LOG_FILE = file.path("results/moduleK", sprintf("moduleK_log_%s.txt", Sys.Date())),
  BIOMARKERS = c("NELL1","STEAP1","GADD45A","FBLN5","HEMK1",
                 "GDE1","RNF14","BNIP3","EGFR","CBR3"),
  DATASETS = list(
    GSE169454 = list(groups = c("normal", "OA"),   layer = "主分析（病例-对照；对照标签实为 normal）"),
    GSE152805 = list(groups = c("MT", "oLT"),     layer = "描述性（供体内病变对比）"),
    GSE167186 = list(groups = c("Young", "Old"),  layer = "探索（aging proxy）")
  ),
  CAND_CELLTYPE = c("celltype", "cell_type", "CellType", "celltype_fix3", "annotation",
                    "celltype_final", "cell_type_fix3"),
  CAND_GROUP    = c("condition", "group", "tissue", "group2", "disease", "age",
                    "old_young", "cohort"),
  MIN_CELLS  = 10,
  SEED_MAIN  = 20260914
)
dir.create(CONFIG$OUT_DIR, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  txt <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = ""))
  cat(txt, "\n"); cat(txt, "\n", file = CONFIG$LOG_FILE, append = TRUE)
}

log_msg("—— 模块 K（14b 修正版）：CellChat 细胞通讯，三数据集分层口径 ——")

for (pkg in c("Seurat", "CellChat", "future"))
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("缺少 ", pkg, " 包——请回查模块 G 环境")
suppressPackageStartupMessages({ library(Seurat); library(CellChat) })
future::plan(future::sequential)
options(future.globals.maxSize = 8 * 1024^3)
log_msg("CellChat ", as.character(packageVersion("CellChat")),
        " | Seurat ", as.character(packageVersion("Seurat")))

## 对象发现：文件名优先级 + 加载验证（细胞类型列 ∈ 候选；分组列含 ≥1 配置分组值）
discover_obj <- function(gse, cfg) {
  pool <- list.files(CONFIG$G2_DIR, pattern = paste0(gse, ".*\\.rds$"),
                     full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  if (length(pool) == 0)
    stop(gse, "：", CONFIG$G2_DIR, " 下无任何 ", gse, " 的 RDS——请确认模块 G 已跑",
         "或把实际路径补进 CONFIG$G2_DIR")
  pri <- grepl("G2|annot|celltype", basename(pool), ignore.case = TRUE)
  pool <- c(pool[pri], pool[!pri])
  log_msg("  候选 RDS ", length(pool), " 个：", paste(basename(pool), collapse = ", "))
  for (f in pool) {
    log_msg("  试载 ", basename(f), " …")
    obj <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(obj)) { log_msg("    读取失败，试下一个"); next }
    meta <- obj@meta.data
    ct <- intersect(CONFIG$CAND_CELLTYPE, colnames(meta))
    gr <- CONFIG$CAND_GROUP[CONFIG$CAND_GROUP %in% colnames(meta)]
    gr_hit <- gr[vapply(gr, function(cl) any(meta[[cl]] %in% cfg$groups), logical(1))]
    if (length(ct) > 0 && length(gr_hit) > 0) {
      ## 分组列择优：含配置分组值最多者
      n_hit <- vapply(gr_hit, function(cl) sum(cfg$groups %in% unique(meta[[cl]])), integer(1))
      gr_col <- gr_hit[which.max(n_hit)]
      log_msg("    ✓ 采用（细胞类型列=", ct[1], "；分组列=", gr_col, "；",
              ncol(obj), " 细胞）")
      return(list(obj = obj, ct_col = ct[1], gr_col = gr_col, file = f))
    }
    log_msg("    元数据不合格（细胞类型候选命中 ", length(ct), "；分组候选命中 ",
            length(gr_hit), "；实际列：", paste(colnames(meta), collapse = ","), "）——试下一个")
    rm(obj); gc(verbose = FALSE)
  }
  stop(gse, "：全部候选 RDS 均无合格元数据——请把实际细胞类型列/分组列名补进 CONFIG$CAND")
}

run_cc <- function(data_norm, meta, celltype_col, gse, grp) {
  cc <- createCellChat(object = data_norm, meta = meta, group.by = celltype_col)
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc, do.fast = FALSE)   # 无 presto 退回标准 wilcoxon
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean", population.size = TRUE,
                          seed.use = CONFIG$SEED_MAIN)
  cc <- filterCommunication(cc, min.cells = CONFIG$MIN_CELLS)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  log_msg("    ", gse, " [", grp, "]：通路数 ", length(cc@netP$pathways))
  cc
}

extract_bmk <- function(cc, gse, grp) {
  df <- subsetCommunication(cc)
  if (is.null(df) || nrow(df) == 0) {
    log_msg("    ", gse, " [", grp, "]：无显著 L-R 互作——诚实阴性记录"); return(NULL) }
  hit <- mapply(function(lig, rec) {
    recs <- strsplit(rec, "_", fixed = TRUE)[[1]]
    lig %in% CONFIG$BIOMARKERS || any(recs %in% CONFIG$BIOMARKERS)
  }, df$ligand, df$receptor)
  sub <- df[hit, , drop = FALSE]
  if (nrow(sub) == 0) {
    log_msg("    ", gse, " [", grp, "]：显著互作 ", nrow(df),
            " 条，无一涉及 10 标志物——诚实阴性记录"); return(NULL) }
  sub$dataset <- gse; sub$group <- grp
  log_msg("    ", gse, " [", grp, "]：显著互作 ", nrow(df), " 条，标志物相关 ",
          nrow(sub), " 条")
  sub
}

# ============================ 主循环 ============================
all_bmk <- list()
for (gse in names(CONFIG$DATASETS)) {
  cfg <- CONFIG$DATASETS[[gse]]
  log_msg("▶ ", gse, "（", cfg$layer, "）")
  dis <- discover_obj(gse, cfg)
  obj <- dis$obj; meta <- obj@meta.data
  ct_col <- dis$ct_col; gr_col <- dis$gr_col

  ct_tab <- table(meta[[ct_col]])
  log_msg("  细胞类型 ", length(ct_tab), " 类：",
          paste(sprintf("%s=%d", names(ct_tab), ct_tab), collapse = ", "))
  if (length(ct_tab) < 3) stop(gse, "：细胞类型 <3 类，请回查 G2 注释")
  if (any(ct_tab < CONFIG$MIN_CELLS))
    log_msg("  ⚠ 有细胞类型 <", CONFIG$MIN_CELLS, " 细胞：",
            paste(names(ct_tab)[ct_tab < CONFIG$MIN_CELLS], collapse = ","), "——留痕")
  log_msg("  分组分布：", paste(sprintf("%s=%d", names(table(meta[[gr_col]])),
                                   table(meta[[gr_col]])), collapse = ", "))

  avail_g <- intersect(CONFIG$BIOMARKERS, rownames(obj))
  log_msg("  标志物可检：", length(avail_g), "/10",
          ifelse(length(avail_g) < 10,
                 paste0("（缺：", paste(setdiff(CONFIG$BIOMARKERS, rownames(obj)), collapse = ","), "）"), ""))
  if (length(avail_g) < 5) stop(gse, "：标志物可检 <5/10——对象或基因名异常，请回查")

  obj <- tryCatch(JoinLayers(obj), error = function(e) obj)
  data_norm <- GetAssayData(obj, assay = "RNA", layer = "data")

  cc_list <- list()
  for (grp in cfg$groups) {
    cells <- rownames(meta)[meta[[gr_col]] == grp]
    if (length(cells) < 100) { log_msg("  ⚠ ", grp, " 组仅 ", length(cells),
                                       " 细胞（<100）——跳过并留痕"); next }
    cc_list[[grp]] <- run_cc(data_norm[, cells], meta[cells, , drop = FALSE],
                             ct_col, gse, grp)
    sub <- extract_bmk(cc_list[[grp]], gse, grp)
    if (!is.null(sub)) all_bmk[[paste(gse, grp)]] <- sub
  }
  saveRDS(cc_list, file.path(CONFIG$OUT_DIR, paste0("moduleK_", gse, "_cellchat.rds")))
  rm(obj, data_norm); gc(verbose = FALSE)

  if (length(cc_list) == 2) {
    merged <- mergeCellChat(cc_list, add.names = names(cc_list))
    cmp <- data.frame(
      dataset = gse,
      group1 = names(cc_list)[1], group2 = names(cc_list)[2],
      n_interaction_g1 = sum(merged@net[[1]]$count), n_interaction_g2 = sum(merged@net[[2]]$count),
      weight_g1 = round(sum(merged@net[[1]]$weight), 3),
      weight_g2 = round(sum(merged@net[[2]]$weight), 3))
    fwrite(cmp, file.path(CONFIG$OUT_DIR, paste0("moduleK_", gse, "_compare_summary.csv")))
    pdf(file.path(CONFIG$OUT_DIR, paste0("moduleK_", gse, "_diffInteraction.pdf")),
        width = 12, height = 6)
    par(mfrow = c(1, 2))
    netVisual_diffInteraction(merged, weight.scale = TRUE)
    netVisual_diffInteraction(merged, weight.scale = TRUE, measure = "weight")
    dev.off()
    log_msg("  组间比较：", names(cc_list)[1], " vs ", names(cc_list)[2],
            "（互作数 ", cmp$n_interaction_g1, " vs ", cmp$n_interaction_g2, "）")
  }
}

bmk_tab <- if (length(all_bmk) > 0) rbindlist(all_bmk, fill = TRUE) else data.frame()
if (nrow(bmk_tab) > 0) {
  fwrite(bmk_tab, file.path(CONFIG$OUT_DIR, "moduleK_biomarker_LR_pairs.csv"))
  log_msg("标志物相关 L-R 总表：", nrow(bmk_tab), " 行 → moduleK_biomarker_LR_pairs.csv")
  if (!any(bmk_tab$dataset == "GSE169454"))
    log_msg("⚠ 主分析层 GSE169454 无标志物相关 L-R——E5 整体定位下调为纯描述性，写稿收窄")
} else {
  log_msg("⚠ 三数据集均无标志物相关显著 L-R——诚实阴性：E5 仅保留通讯格局描述层")
}

log_msg("模块 K 完成。红线复述：GSE152805 为供体内病变对比（无健康对照）；",
        "GSE167186 为 aging proxy，不得写成肌少症证据；全部结果入补充材料。")
