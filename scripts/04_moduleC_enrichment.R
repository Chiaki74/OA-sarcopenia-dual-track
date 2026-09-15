# ============================================================================
# 04_moduleC_enrichment.R
# 模块 C｜双轨核心基因的分轨功能富集分析（计划书 v1.4a §3.3.3；备忘录 §八-24/25）
# ----------------------------------------------------------------------------
# 设计依据（2026-08-22 双轨框架）：
#   核心基因集不再单一，按方向一致性拆为三条分析轨：
#     轨 1  concordant           —— 同向轨 27 个（OA 与肌肉同向失调）
#     轨 2  mirror_OAup_muscleDown —— 镜像轨·OA↑/肌↓（线粒体 OXPHOS 簇所在轨）
#     轨 3  mirror_OAdown_muscleUp —— 镜像轨·OA↓/肌↑（CBLB 萎缩轴所在轨）
#   方法：ORA（clusterProfiler: GO-BP / KEGG；ReactomePA: Reactome），
#         背景宇宙 = 双侧 RRA 共同检测基因（intersect(layer2_RRA_OA, layer2_RRA_muscle)），
#         BH 校正，p.adjust < 0.05 显著；compareCluster 做三轨横向对比。
#   生物学检查点（软核对，仅记录日志、不构成硬断言）：
#     ① 镜像 OA↑/肌↓ 轨应出现线粒体/OXPHOS 相关条目（呼应 PMID 31862890、
#        PMID 38504132 Nature Metab 2024 的肌肉线粒体下调 + OA 代偿性上调）；
#     ② 同向轨应出现补体/凝血与 ECM 组织相关条目（C1QA/C1QB/TIMP1 驱动）；
#     ③ 镜像 OA↓/肌↑ 轨关注蛋白稳态/泛素-蛋白酶体（CBLB E3 连接酶轴）。
# ----------------------------------------------------------------------------
# 前置：已运行 03_moduleB_shared_genes.R（results/moduleB/ 下有
#       CORE_concordant_genes.csv / CORE_mirror_genes.csv /
#       layer2_RRA_OA.csv / layer2_RRA_muscle.csv）
# 运行：source("04_moduleC_enrichment.R")；结果输出 results/moduleC/
# ============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(enrichplot)
  library(ggplot2)
  library(data.table)
})

# ----------------------------- CONFIG --------------------------------------
CONFIG <- list(
  IN_DIR   = "results/moduleB",
  OUT_DIR  = "results/moduleC",
  LOG_FILE = sprintf("moduleC_log_%s.txt", Sys.Date()),
  PADJ_CUT   = 0.05,    # BH 校正显著性阈值
  MIN_GS     = 3,       # 最小基因集大小（小核心集适配）
  MAX_GS     = 500,
  TOP_N_SHOW = 15,      # 每轨每库展示条目数
  FIG_W = 9, FIG_H = 7  # 图尺寸（英寸）
)
dir.create(CONFIG$OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------- 工具函数 -------------------------------------
log_msg <- function(...) {
  txt <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = ""))
  cat(txt, "\n"); cat(txt, "\n", file = CONFIG$LOG_FILE, append = TRUE)
}

# ORA 三库执行器：空结果不崩溃，返回 NULL 并记日志
run_ora <- function(entrez_vec, universe_entrez, track_name) {
  res <- list()
  res$GO_BP <- tryCatch(
    enrichGO(gene = entrez_vec, universe = universe_entrez,
             OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
             pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
             minGSSize = CONFIG$MIN_GS, maxGSSize = CONFIG$MAX_GS,
             readable = TRUE),
    error = function(e) { log_msg("  ⚠ ", track_name, " GO-BP 失败：", conditionMessage(e)); NULL })
  res$KEGG <- tryCatch(
    enrichKEGG(gene = entrez_vec, universe = universe_entrez, organism = "hsa",
               pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
               minGSSize = CONFIG$MIN_GS, maxGSSize = CONFIG$MAX_GS),
    error = function(e) { log_msg("  ⚠ ", track_name, " KEGG 失败：", conditionMessage(e)); NULL })
  res$Reactome <- tryCatch(
    enrichPathway(gene = entrez_vec, universe = universe_entrez,
                  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                  minGSSize = CONFIG$MIN_GS, maxGSSize = CONFIG$MAX_GS,
                  readable = TRUE),
    error = function(e) { log_msg("  ⚠ ", track_name, " Reactome 失败：", conditionMessage(e)); NULL })
  res
}

# 提取显著条目表（统一列口径；无显著则返回 NULL）
sig_table <- function(er, db_name, track_name) {
  if (is.null(er)) return(NULL)
  df <- as.data.frame(er)
  if (nrow(df) == 0) return(NULL)
  df <- df[df$p.adjust < CONFIG$PADJ_CUT, ]
  if (nrow(df) == 0) return(NULL)
  df$database <- db_name; df$track <- track_name
  df[order(df$p.adjust), ]
}

# ============================================================================
# 1. 读取模块 B 核心集 + 构建三轨
# ============================================================================
f_conc  <- file.path(CONFIG$IN_DIR, "CORE_concordant_genes.csv")
f_mirr  <- file.path(CONFIG$IN_DIR, "CORE_mirror_genes.csv")
f_rra_o <- file.path(CONFIG$IN_DIR, "layer2_RRA_OA.csv")
f_rra_m <- file.path(CONFIG$IN_DIR, "layer2_RRA_muscle.csv")
for (f in c(f_conc, f_mirr, f_rra_o, f_rra_m))
  if (!file.exists(f)) stop("缺少模块 B 输出：", f, "——请先运行 03_moduleB_shared_genes.R")

core_conc <- fread(f_conc)
core_mirr <- fread(f_mirr)
rra_oa <- fread(f_rra_o)
rra_mu <- fread(f_rra_m)

# 方向标签硬校验（防 03 口径漂移造成静默错分；独立审查 S2）
valid_dir <- c("up", "down")
if (!all(core_conc$direction_OA %in% valid_dir) || !all(core_conc$direction_muscle %in% valid_dir) ||
    !all(core_mirr$direction_OA %in% valid_dir) || !all(core_mirr$direction_muscle %in% valid_dir))
  stop("direction 列存在非 up/down 取值，请核对 03 脚本输出口径")
if (any(core_mirr$direction_OA == core_mirr$direction_muscle))
  stop("CORE_mirror 中含同向行——镜像轨口径漂移，请核对 03 脚本")
if (any(core_conc$direction_OA != core_conc$direction_muscle))
  stop("CORE_concordant 中含反向行——同向轨口径漂移，请核对 03 脚本")

tracks <- list(
  concordant = core_conc$gene,
  mirror_OAup_muscleDown = core_mirr$gene[core_mirr$direction_OA == "up"],
  mirror_OAdown_muscleUp = core_mirr$gene[core_mirr$direction_OA == "down"]
)
log_msg("三轨基因数：concordant ", length(tracks$concordant),
        " / mirror_OAup_muscleDown ", length(tracks$mirror_OAup_muscleDown),
        " / mirror_OAdown_muscleUp ", length(tracks$mirror_OAdown_muscleUp))
if (any(sapply(tracks, length) < 3)) log_msg("⚠ 存在基因数 <3 的轨，其富集结果可能为空")

# 带轨标签的核心基因总表（下游模块 D/E 输入）
# 逐轨构建 + rbindlist：0 长轨不产生"长度1字面量回收"崩溃（独立审查 F2）
build_track_tab <- function(genes, track_name, src_df, mask) {
  if (length(genes) == 0) return(NULL)
  data.frame(gene = genes, track = track_name,
             direction_OA = src_df$direction_OA[mask],
             direction_muscle = src_df$direction_muscle[mask],
             in_WGCNA_key = src_df$in_WGCNA_key[mask])
}
mask_mup   <- core_mirr$direction_OA == "up"
mask_mdown <- core_mirr$direction_OA == "down"
track_tab <- rbindlist(list(
  build_track_tab(tracks$concordant, "concordant", core_conc, rep(TRUE, nrow(core_conc))),
  build_track_tab(tracks$mirror_OAup_muscleDown, "mirror_OAup_muscleDown", core_mirr, mask_mup),
  build_track_tab(tracks$mirror_OAdown_muscleUp, "mirror_OAdown_muscleUp", core_mirr, mask_mdown)
), fill = TRUE)
write.csv(track_tab, file.path(CONFIG$OUT_DIR, "moduleC_gene_tracks.csv"), row.names = FALSE)

# ============================================================================
# 2. 背景宇宙：双侧 RRA 共同检测基因 → Entrez
# ============================================================================
universe_symbol <- intersect(rra_oa$gene, rra_mu$gene)
log_msg("背景宇宙（双侧共同检测基因）：", length(universe_symbol), " 个")

sym2entrez <- function(symbols, tag) {
  # bitr 对空输入/全落空键在部分版本会报错，兜底返回 character(0)（独立审查 A3）
  mp <- tryCatch(bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
                 error = function(e) { log_msg("  ⚠ ", tag, " bitr 失败：", conditionMessage(e))
                                       data.frame(SYMBOL = character(0), ENTREZID = character(0)) })
  lost <- setdiff(symbols, mp$SYMBOL)
  if (length(lost) > 0)
    log_msg("  ℹ ", tag, "：", length(lost), " 个 symbol 未能映射 Entrez（",
            paste(head(lost, 5), collapse = ","), if (length(lost) > 5) " ..." else "", "）")
  unique(mp$ENTREZID)
}
universe_entrez <- sym2entrez(universe_symbol, "universe")

# ============================================================================
# 3. 逐轨 ORA 三库
# ============================================================================
all_sig <- list()   # 汇总显著条目
ora_results <- list()
for (tr in names(tracks)) {
  log_msg("—— 轨 ", tr, "（", length(tracks[[tr]]), " 基因）——")
  ez <- sym2entrez(tracks[[tr]], tr)
  if (length(ez) < 3) { log_msg("  ⚠ Entrez 基因数 <3，跳过"); next }
  ora <- run_ora(ez, universe_entrez, tr)
  ora_results[[tr]] <- ora
  for (db in names(ora)) {
    st <- sig_table(ora[[db]], db, tr)
    n_sig <- if (is.null(st)) 0 else nrow(st)
    log_msg("  ", db, " 显著条目：", n_sig)
    out_csv <- file.path(CONFIG$OUT_DIR, sprintf("enrich_%s_%s.csv", db, tr))
    if (!is.null(st)) {
      write.csv(st, out_csv, row.names = FALSE)
      all_sig[[paste(tr, db, sep = "|")]] <- st
    } else {
      # 空结果也落盘占位，防"文件缺失"与"无显著"混淆
      write.csv(data.frame(track = tr, database = db, note = "no significant term at p.adjust<0.05"),
                out_csv, row.names = FALSE)
    }
  }
}
if (length(all_sig) > 0) {
  # KEGG 结果多 category/subcategory 两列，base::rbind 必崩——rbindlist(fill=TRUE)（独立审查 F1）
  sig_all <- rbindlist(all_sig, fill = TRUE)
  write.csv(sig_all, file.path(CONFIG$OUT_DIR, "enrich_ALL_tracks_significant.csv"), row.names = FALSE)
  log_msg("全部轨×库显著条目汇总 ", nrow(sig_all), " 行 → enrich_ALL_tracks_significant.csv")
}

# ============================================================================
# 4. 三轨横向对比（compareCluster：GO-BP 与 Reactome）
# ============================================================================
track_entrez <- lapply(names(tracks), function(tr) sym2entrez(tracks[[tr]], paste0("cc-", tr)))
names(track_entrez) <- names(tracks)
track_entrez <- track_entrez[sapply(track_entrez, length) >= 3]

if (length(track_entrez) >= 2) {
  # 背景宇宙必须与逐轨 ORA 一致（双侧 RRA 共同检测基因），否则口径静默漂移（独立审查 S1）
  cc_go <- tryCatch(
    compareCluster(geneClusters = track_entrez, fun = "enrichGO",
                   OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
                   universe = universe_entrez,
                   pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                   minGSSize = CONFIG$MIN_GS, maxGSSize = CONFIG$MAX_GS,
                   readable = TRUE),
    error = function(e) { log_msg("⚠ compareCluster GO 失败：", conditionMessage(e)); NULL })
  cc_re <- tryCatch(
    compareCluster(geneClusters = track_entrez, fun = "enrichPathway",
                   universe = universe_entrez,
                   pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                   minGSSize = CONFIG$MIN_GS, maxGSSize = CONFIG$MAX_GS,
                   readable = TRUE),
    error = function(e) { log_msg("⚠ compareCluster Reactome 失败：", conditionMessage(e)); NULL })
  saveRDS(list(GO_BP = cc_go, Reactome = cc_re),
          file.path(CONFIG$OUT_DIR, "compareCluster_objects.rds"))

  for (cc_name in c("GO_BP", "Reactome")) {
    cc <- if (cc_name == "GO_BP") cc_go else cc_re
    if (is.null(cc)) next
    cc_df <- as.data.frame(cc)
    write.csv(cc_df, file.path(CONFIG$OUT_DIR, sprintf("compareCluster_%s_full.csv", cc_name)),
              row.names = FALSE)
    cc_sig <- cc_df[cc_df$p.adjust < CONFIG$PADJ_CUT, ]
    if (nrow(cc_sig) > 0) {
      # 用副本改写 slot，保留原对象（RDS 已存未过滤版；独立审查 A1）
      cc_plot <- cc
      cc_plot@compareClusterResult <- cc_sig   # 图只画显著条目
      p <- tryCatch(
        dotplot(cc_plot, showCategory = CONFIG$TOP_N_SHOW, by = "geneRatio",
                title = sprintf("三轨功能富集对比（%s, BH<0.05）", cc_name)) +
          theme(axis.text.x = element_text(angle = 30, hjust = 1)),
        error = function(e) { log_msg("⚠ ", cc_name, " dotplot 失败：", conditionMessage(e)); NULL })
      if (!is.null(p)) {
        # Windows 下中文图题需 cairo_pdf，否则 PDF 字形乱码（独立审查 A2）
        ggsave(file.path(CONFIG$OUT_DIR, sprintf("compareCluster_%s_dotplot.pdf", cc_name)),
               p, width = CONFIG$FIG_W + 1, height = CONFIG$FIG_H + 1, device = cairo_pdf)
        ggsave(file.path(CONFIG$OUT_DIR, sprintf("compareCluster_%s_dotplot.png", cc_name)),
               p, width = CONFIG$FIG_W + 1, height = CONFIG$FIG_H + 1, dpi = 300)
        log_msg("  ", cc_name, " 对比 dotplot 已输出（显著条目 ", nrow(cc_sig), "）")
      }
    } else log_msg("  ", cc_name, " 三轨均无显著条目，跳过 dotplot")
  }
} else log_msg("⚠ 有效轨 <2，跳过 compareCluster")

# ============================================================================
# 5. 生物学检查点（软核对，备忘录 §八-24 预期模式）
# ============================================================================
log_msg("—— 生物学检查点（软核对）——")
check_terms <- function(db_csv_pattern, keywords, label) {
  csvs <- list.files(CONFIG$OUT_DIR, pattern = db_csv_pattern, full.names = TRUE)
  hits <- character(0)
  for (f in csvs) {
    d <- tryCatch(fread(f), error = function(e) NULL)
    if (is.null(d) || !"Description" %in% names(d)) next
    m <- Reduce(`|`, lapply(keywords, function(k) grepl(k, d$Description, ignore.case = TRUE)))
    m[is.na(m)] <- FALSE   # Description 含 NA 时防索引带出 NA 行（独立审查 A5）
    hit <- d$Description[m]
    hits <- c(hits, paste0(basename(f), ": ", hit))
  }
  if (length(hits) > 0) log_msg("  ✅ ", label, " 命中：\n     ", paste(hits, collapse = "\n     "))
  else log_msg("  ⚠ ", label, " 未命中（写作时如实报告，不强行解读）")
}
check_terms("^enrich_.*mirror_OAup_muscleDown\\.csv$",
            c("oxidative phosphorylation", "mitochondri", "respiratory chain", "ATP synthesis"),
            "检查点① 镜像OA↑肌↓轨线粒体/OXPHOS")
check_terms("^enrich_.*concordant\\.csv$",
            c("complement", "coagulation", "extracellular matrix", "collagen"),
            "检查点② 同向轨补体/ECM")
check_terms("^enrich_.*mirror_OAdown_muscleUp\\.csv$",
            c("ubiquitin", "proteasome", "proteolysis", "muscle atrophy", "protein catabolic"),
            "检查点③ 镜像OA↓肌↑轨泛素-蛋白酶体/萎缩")

cat("\n模块 C 完成（分轨富集）。核心产物：\n",
    "  ", CONFIG$OUT_DIR, "/moduleC_gene_tracks.csv            —— 三轨基因总表（下游模块 D/E 输入）\n",
    "  ", CONFIG$OUT_DIR, "/enrich_<DB>_<track>.csv            —— 逐轨逐库显著条目\n",
    "  ", CONFIG$OUT_DIR, "/enrich_ALL_tracks_significant.csv  —— 全部显著条目汇总\n",
    "  ", CONFIG$OUT_DIR, "/compareCluster_<DB>_dotplot.pdf/png —— 三轨横向对比图\n",
    "  ", CONFIG$OUT_DIR, "/compareCluster_objects.rds          —— 绘图对象（改图备用）\n",
    "检查日志 ", CONFIG$LOG_FILE, " 中的生物学检查点核对结果后再继续。\n", sep = "")
