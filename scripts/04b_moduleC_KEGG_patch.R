# ============================================================================
# 04b_moduleC_KEGG_patch.R
# 模块 C 补丁｜KEGG 富集专项补跑（2026-08-22）
# ----------------------------------------------------------------------------
# 背景：模块 C 首跑时 rest.kegg.jp 直连失败（SSL 错误 + 超时），三轨 KEGG
#       全部落占位文件；在线重试 5/5 再败。2026-08-23 改版为**零联网方案**：
#       ① 助手侧预下载 KEGG hsa 数据并清洗为两个 TSV（口径与
#          clusterProfiler::download_KEGG 清洗后完全一致，快照日期 2026-08-23），
#          用户只需把 TSV 放入 results/moduleC/；
#       ② 用 enricher() 纯离线富集，口径与 04 主脚本完全一致
#          （背景宇宙 = 双侧 RRA 共同检测基因，BH<0.05，minGSSize=3）；
#       ③ 覆盖写入 enrich_KEGG_<track>.csv，并把显著条目（剔除旧 KEGG 行后）
#          追加进 enrich_ALL_tracks_significant.csv；
#       ④ TSV 缺失时才回退在线下载（带 5 次重试）。
# 前置：已运行 04_moduleC_enrichment.R（results/moduleC/ 已有占位文件）；
#       两个 KEGG TSV 已放入 results/moduleC/
# 运行：source("04b_moduleC_KEGG_patch.R")（全程零联网）
# Methods 报告口径：KEGG 数据快照日期 2026-08-23（rest.kegg.jp）
# ============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(data.table)
})

CONFIG <- list(
  MODB_DIR = "results/moduleB",
  OUT_DIR  = "results/moduleC",
  LOG_FILE = file.path("results/moduleC", sprintf("moduleC_KEGG_patch_log_%s.txt", Sys.Date())),  # 审查 A5：日志入 OUT_DIR
  KEGG_CACHE = "results/moduleC/KEGG_hsa_cache.rds",
  PADJ_CUT = 0.05, MIN_GS = 3, MAX_GS = 500,
  MAX_RETRY = 5
)

log_msg <- function(...) {
  txt <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = ""))
  cat(txt, "\n"); cat(txt, "\n", file = CONFIG$LOG_FILE, append = TRUE)
}

# ---- 0. 读取三轨基因表 + 背景宇宙（与 04 主脚本同口径）----
f_track <- file.path(CONFIG$OUT_DIR, "moduleC_gene_tracks.csv")
f_rra_o <- file.path(CONFIG$MODB_DIR, "layer2_RRA_OA.csv")
f_rra_m <- file.path(CONFIG$MODB_DIR, "layer2_RRA_muscle.csv")
for (f in c(f_track, f_rra_o, f_rra_m))
  if (!file.exists(f)) stop("缺少文件：", f, "——请先按顺序运行 03/04 脚本")

track_tab <- fread(f_track)
rra_oa <- fread(f_rra_o); rra_mu <- fread(f_rra_m)
# 列名硬校验（审查 A3/A4：缺列/空 track 静默丢基因）
stopifnot("gene" %in% names(rra_oa), "gene" %in% names(rra_mu),
          all(c("gene", "track") %in% names(track_tab)),
          !any(is.na(track_tab$track)), !any(is.na(track_tab$gene)))
tracks <- split(track_tab$gene, track_tab$track)
log_msg("读取三轨：", paste(sprintf("%s=%d", names(tracks), sapply(tracks, length)), collapse = " / "))

universe_symbol <- intersect(rra_oa$gene, rra_mu$gene)

sym2entrez_map <- function(symbols) {
  tryCatch(bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
           error = function(e) { log_msg("⚠ bitr 失败：", conditionMessage(e))
                                 data.frame(SYMBOL = character(0), ENTREZID = character(0)) })
}
mp_uni  <- sym2entrez_map(universe_symbol)
universe_entrez <- unique(mp_uni$ENTREZID)
# 审查 S3：universe 映射为空时级联成全假阴性，必须硬中止
if (length(universe_entrez) < 50)
  stop("背景宇宙 Entrez 映射异常（", length(universe_entrez), " 个）——请检查 org.Hs.eg.db 是否完好")
log_msg("背景宇宙：", length(universe_symbol), " symbol → ", length(universe_entrez), " Entrez")

# ---- 1. 加载 KEGG 数据：优先本地 TSV（零联网），缺失才回退在线下载 ----
# 2026-08-23 改版：rest.kegg.jp 在用户网络下持续超时（5/5 失败），改由助手侧
# 预下载并清洗为 TSV（口径与 clusterProfiler::download_KEGG 清洗后完全一致：
# 剥 path:/hsa: 前缀、去物种名后缀）。用户只需把两个 TSV 放入 results/moduleC/。
f_t2g <- file.path(CONFIG$OUT_DIR, "KEGG_hsa_pathway2gene_20260823.tsv")
f_t2n <- file.path(CONFIG$OUT_DIR, "KEGG_hsa_pathway_names_20260823.tsv")

if (file.exists(f_t2g) && file.exists(f_t2n)) {
  t2g <- fread(f_t2g)   # 列：from=通路ID(hsa00010) / to=Entrez
  t2n <- fread(f_t2n)   # 列：from=通路ID / to=通路名称
  log_msg("使用本地 KEGG TSV（快照日期 2026-08-23，零联网）：",
          length(unique(t2g$from)), " 通路 / ", length(unique(t2g$to)), " Entrez 基因")
} else {
  log_msg("本地 TSV 缺失，回退在线下载（带重试）...")
  kegg_data <- NULL
  if (file.exists(CONFIG$KEGG_CACHE)) {
    kegg_data <- readRDS(CONFIG$KEGG_CACHE)
    log_msg("使用本地缓存 ", CONFIG$KEGG_CACHE, "（零联网）")
  } else {
    for (i in seq_len(CONFIG$MAX_RETRY)) {
      log_msg("下载 KEGG hsa 通路数据（第 ", i, "/", CONFIG$MAX_RETRY, " 次）...")
      kegg_data <- tryCatch(download_KEGG("hsa"),
                            error = function(e) { log_msg("  失败：", conditionMessage(e)); NULL })
      if (!is.null(kegg_data)) break
      if (i < CONFIG$MAX_RETRY) { log_msg("  ", 30 * i, " 秒后重试"); Sys.sleep(30 * i) }
    }
    if (is.null(kegg_data))
      stop("本地 TSV 缺失且在线下载连续 ", CONFIG$MAX_RETRY, " 次失败——",
           "请将助手提供的 KEGG_hsa_pathway2gene_20260823.tsv 与 KEGG_hsa_pathway_names_20260823.tsv ",
           "放入 results/moduleC/ 后重跑")
    saveRDS(kegg_data, CONFIG$KEGG_CACHE)
  }
  t2g <- kegg_data$KEGGPATHID2EXTID
  t2n <- kegg_data$KEGGPATHID2NAME
  log_msg("在线数据快照日期：", Sys.Date(), "（Methods 中 KEGG 版本按此日期报告）")
}
# 审查 A2：结构校验
if (is.null(t2g) || is.null(t2n) || ncol(t2g) != 2 || nrow(t2g) == 0)
  stop("KEGG 数据结构异常——请检查 TSV 文件或删除缓存后重跑")
# enricher 要求字符型键（fread 可能把纯数字 to 列读成整数，防御性转换）
t2g <- as.data.frame(t2g); t2n <- as.data.frame(t2n)
t2g$to <- as.character(t2g$to)
log_msg("KEGG 通路数：", length(unique(t2g[, 1])))

# ---- 2. 逐轨离线富集（enricher，口径同 04 主脚本）----
# 审查 A1：symbol 查找表按 Entrez 去重，避免一对多取到别名
mp_dedup <- mp_uni[!duplicated(mp_uni$ENTREZID), ]
sym_lookup <- setNames(mp_dedup$SYMBOL, mp_dedup$ENTREZID)

write_placeholder <- function(out_csv, tr, note_text) {
  write.csv(data.frame(track = tr, database = "KEGG", note = note_text),
            out_csv, row.names = FALSE)
}

sig_new <- list()
kegg_status <- list()   # 三轨状态记录（审查 A6：日志口径一致）
for (tr in names(tracks)) {
  out_csv <- file.path(CONFIG$OUT_DIR, sprintf("enrich_KEGG_%s.csv", tr))
  ez <- unique(sym2entrez_map(tracks[[tr]])$ENTREZID)
  if (length(ez) < 3) {   # 审查 S4：跳过分支也必须落占位，防陈旧文件误导
    write_placeholder(out_csv, tr, "Entrez<3, skipped")
    kegg_status[[tr]] <- "skipped(Entrez<3)"
    log_msg("轨 ", tr, " Entrez <3，已落占位（skipped）"); next
  }
  # 审查 S2：区分 enricher 技术故障与真实阴性
  er <- NULL; er_err <- NULL
  tryCatch(
    er <<- enricher(ez, TERM2GENE = t2g, TERM2NAME = t2n,
                    universe = universe_entrez,
                    pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                    minGSSize = CONFIG$MIN_GS, maxGSSize = CONFIG$MAX_GS),
    error = function(e) er_err <<- conditionMessage(e))
  if (!is.null(er_err)) {
    write_placeholder(out_csv, tr, paste0("enricher error: ", er_err))
    kegg_status[[tr]] <- "ERROR"
    log_msg("⚠ 轨 ", tr, " enricher 报错（已落 ERROR 占位，勿作阴性解读）：", er_err)
    next
  }
  st <- NULL
  if (!is.null(er)) {
    df <- as.data.frame(er)
    df <- df[df$p.adjust < CONFIG$PADJ_CUT, ]
    if (nrow(df) > 0) {
      # geneID 由 Entrez 转 symbol，与 GO/Reactome 输出口径一致
      df$geneID <- vapply(strsplit(df$geneID, "/"),
                          function(x) paste(na.omit(sym_lookup[x]), collapse = "/"),
                          character(1))
      if (any(df$geneID == "")) log_msg("  ⚠ ", tr, " 有条目 geneID 转换后为空，请抽查")
      df$database <- "KEGG"; df$track <- tr
      st <- df[order(df$p.adjust), ]
      sig_new[[tr]] <- st
    }
  }
  if (!is.null(st)) {
    write.csv(st, out_csv, row.names = FALSE)
    kegg_status[[tr]] <- sprintf("significant(%d)", nrow(st))
    log_msg("轨 ", tr, "：KEGG 显著条目 ", nrow(st), " → ", basename(out_csv))
  } else {
    write_placeholder(out_csv, tr, "no significant term at p.adjust<0.05")
    kegg_status[[tr]] <- "true negative"
    log_msg("轨 ", tr, "：KEGG 无显著条目（真实阴性，enricher 正常返回）")
  }
}
log_msg("KEGG 补跑状态总览：", paste(sprintf("%s=%s", names(kegg_status), kegg_status), collapse = " / "))

# ---- 3. 追加进汇总表（审查 S1：先剔除旧 KEGG 行，重跑不产生重复行）----
f_all <- file.path(CONFIG$OUT_DIR, "enrich_ALL_tracks_significant.csv")
if (length(sig_new) > 0) {
  old <- if (file.exists(f_all)) fread(f_all) else NULL
  if (!is.null(old) && nrow(old) > 0 && "database" %in% names(old))
    old <- old[old$database != "KEGG", ]
  sig_all <- rbindlist(c(list(old), sig_new), fill = TRUE)
  write.csv(sig_all, f_all, row.names = FALSE)
  log_msg("汇总表已更新（旧 KEGG 行已剔除再并入）：", nrow(sig_all), " 行 → ", basename(f_all))
}

cat("\nKEGG 补跑完成。请将本日志贴回核对：\n",
    "  ", CONFIG$LOG_FILE, "\n",
    "注意：enricher 离线结果无 category/subcategory 两列（在线 enrichKEGG 才有），\n",
    "属正常差异，不影响显著性结论。\n", sep = "")
