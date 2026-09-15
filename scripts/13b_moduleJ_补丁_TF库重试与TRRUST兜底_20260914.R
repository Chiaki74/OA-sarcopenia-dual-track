# ===== 补丁 13b 20260914：TF 库 HTTP/2 报错修复 + TRRUST 直连兜底 =====
# 现象：Enrichr 的 TRRUST/ChEA_2022 库报 "HTTP/2 stream ... INTERNAL_ERROR"
#   （maayanlab 服务器端已知问题），miRTarBase 正常 → TF 层缺位。
# 修复双通道：A) httr 强制 HTTP/1.1 后逐库重试 ×3；
#   B) 仍失败则下载 TRRUST 官方原始表（grnpedia.org，文献 curated TF→靶基因），
#      以「数据库直接比对」口径建边（adjP=NA，library 标 TRRUST_direct）。
# 前置：与 13 同会话最佳（CONFIG/log_msg 已在）；独立运行亦可（自动重建最小配置）。

suppressPackageStartupMessages({ library(data.table); library(enrichR); library(httr) })

if (!exists("CONFIG"))
  CONFIG <- list(OUT_DIR = "results/moduleJ",
                 LOG_FILE = file.path("results/moduleJ", sprintf("moduleJ_log_%s.txt", Sys.Date())),
                 BIOMARKERS = c("NELL1","STEAP1","GADD45A","FBLN5","HEMK1",
                                "GDE1","RNF14","BNIP3","EGFR","CBR3"),
                 TF_LIBS = c("TRRUST_Transcription_Factors_2019", "ChEA_2022"),
                 MIR_LIBS = c("miRTarBase_2017"),
                 HUB_MIN = 3, PADJ_CUT = 0.05, SEED_MAIN = 20260914)
dir.create(CONFIG$OUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (!exists("log_msg"))
  log_msg <- function(...) {
    txt <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = ""))
    cat(txt, "\n"); cat(txt, "\n", file = CONFIG$LOG_FILE, append = TRUE)
  }

log_msg("—— 补丁 13b：TF 层修复重跑 ——")

# ---- 通道 A：强制 HTTP/1.1 + 重试 ----
httr::set_config(httr::config(http_version = 2L))   # 2 = CURL_HTTP_VERSION_1_1
query_retry <- function(lib, tries = 3) {
  for (i in seq_len(tries)) {
    res <- tryCatch(enrichr(CONFIG$BIOMARKERS, lib)[[lib]],
                    error = function(e) { log_msg("  ", lib, " 第 ", i, " 次失败：",
                                                  conditionMessage(e)); NULL })
    if (!is.null(res) && nrow(res) > 0) {
      log_msg("  ", lib, " 第 ", i, " 次成功：", nrow(res), " 条")
      return(as.data.frame(res)) }
    Sys.sleep(5)
  }
  NULL
}
log_msg("通道 A：HTTP/1.1 重试 Enrichr TF 库")
tf_tabs <- lapply(CONFIG$TF_LIBS, query_retry)
names(tf_tabs) <- CONFIG$TF_LIBS
mir_tab <- query_retry(CONFIG$MIR_LIBS[1], tries = 1)   # miRNA 层顺带重查，保持同次口径

tf_ok <- !vapply(tf_tabs, is.null, logical(1))
log_msg("通道 A 结果：", paste(sprintf("%s=%s", names(tf_ok), ifelse(tf_ok, "OK", "FAIL")),
                               collapse = " "))

# ---- 通道 B：TRRUST 直连兜底（任一 TF 库失败即触发）----
trrust_direct <- data.frame()
if (!all(tf_ok)) {
  log_msg("通道 B：下载 TRRUST 官方原始表（grnpedia.org）")
  raw <- NULL
  for (u in c("https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv",
              "http://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv")) {
    raw <- tryCatch(fread(u, header = FALSE), error = function(e) NULL)
    if (!is.null(raw) && nrow(raw) > 0) break
  }
  if (is.null(raw) || nrow(raw) == 0) {
    log_msg("  ⚠ TRRUST 直连也失败——TF 层本批如实记「数据库暂不可达，未执行」，",
            "不阻塞模块 K；改日重跑本补丁即可补齐")
  } else {
    setnames(raw, c("regulator", "target", "mode", "pmid")[seq_len(min(4, ncol(raw)))])
    hit <- raw[target %in% CONFIG$BIOMARKERS]
    trrust_direct <- unique(data.frame(
      regulator = hit$regulator, target = hit$target, kind = "TF",
      library = "TRRUST_direct", adjP = NA_real_, nominal_only = FALSE,
      stringsAsFactors = FALSE))
    log_msg("  TRRUST 直连：", nrow(raw), " 条全库记录中，靶向 10 标志物的边 ",
            nrow(trrust_direct), " 条（", length(unique(trrust_direct$regulator)),
            " 个 TF；覆盖 ", length(unique(trrust_direct$target)), "/10 标志物）")
  }
}

# ---- 重建边表（Enrichr 显著边 + TRRUST 直连边 + miRNA 边）----
to_edges <- function(tab, kind, lib) {
  if (is.null(tab)) return(data.frame())
  tab$library <- lib
  sig <- tab[tab$Adjusted.P.value < CONFIG$PADJ_CUT, , drop = FALSE]
  if (nrow(sig) == 0) {
    log_msg("  ", lib, "：无 adj.P<", CONFIG$PADJ_CUT, " 条目——取名义显著前 10 并标记")
    sig <- head(tab[order(tab$Adjusted.P.value), ], 10); sig$nominal_only <- TRUE
  } else sig$nominal_only <- FALSE
  out <- do.call(rbind, lapply(seq_len(nrow(sig)), function(i) {
    gs <- intersect(strsplit(sig$Genes[i], ";", fixed = TRUE)[[1]], CONFIG$BIOMARKERS)
    if (length(gs) == 0) return(NULL)
    data.frame(regulator = sig$Term[i], target = gs, kind = kind, library = lib,
               adjP = sig$Adjusted.P.value[i], nominal_only = sig$nominal_only[i],
               stringsAsFactors = FALSE)
  }))
  if (is.null(out)) data.frame() else unique(out)
}
edges <- unique(rbind(
  do.call(rbind, Map(function(tb, lb) to_edges(tb, "TF", lb),
                     tf_tabs[tf_ok], names(tf_tabs)[tf_ok])),
  trrust_direct,
  to_edges(mir_tab, "miRNA", CONFIG$MIR_LIBS[1])))
if (nrow(edges) == 0)
  stop("13b 后仍零边——Enrichr 与 TRRUST 双通道均不可达，请改日重跑；模块 K 不受影响")

tf_src_n <- tapply(edges$library[edges$kind == "TF"], edges$regulator[edges$kind == "TF"],
                   function(x) length(unique(x)))
edges$tf_dual <- edges$kind == "TF" &
  edges$regulator %in% names(tf_src_n[tf_src_n >= 2])
fwrite(edges, file.path(CONFIG$OUT_DIR, "moduleJ_TFmiRNA_edges.csv"))
log_msg("边表重建：", nrow(edges), " 条（TF ", sum(edges$kind == "TF"),
        " / miRNA ", sum(edges$kind == "miRNA"), "；多源 TF 高置信 ",
        length(unique(edges$regulator[edges$tf_dual])), " 个）→ moduleJ_TFmiRNA_edges.csv")

cov <- sapply(CONFIG$BIOMARKERS, function(g) any(edges$target == g & !edges$nominal_only))
log_msg("标志物覆盖：", sum(cov), "/10",
        ifelse(sum(cov) < 10,
               paste0("（未覆盖：", paste(CONFIG$BIOMARKERS[!cov], collapse = ","), "——如实报告）"), ""))

# ---- hub 与网络图重建 ----
deg <- tapply(edges$target[!edges$nominal_only], edges$regulator[!edges$nominal_only],
              function(x) length(unique(x)))
hub <- data.frame(regulator = names(deg), n_targets = as.integer(deg),
                  kind = edges$kind[match(names(deg), edges$regulator)],
                  tf_dual = names(deg) %in% names(tf_src_n[tf_src_n >= 2]))
hub <- hub[order(-hub$n_targets), ]; hub$is_hub <- hub$n_targets >= CONFIG$HUB_MIN
fwrite(hub, file.path(CONFIG$OUT_DIR, "moduleJ_TFmiRNA_hub.csv"))
log_msg("hub（≥", CONFIG$HUB_MIN, "/10 靶）", sum(hub$is_hub), " 个：",
        paste(head(hub$regulator[hub$is_hub], 15), collapse = ", "))

if (requireNamespace("igraph", quietly = TRUE) && nrow(edges[!edges$nominal_only, ]) > 0) {
  suppressPackageStartupMessages(library(igraph))
  e_sig <- edges[!edges$nominal_only, ]
  g <- graph_from_data_frame(e_sig[, c("regulator", "target")], directed = TRUE)
  V(g)$type <- V(g)$name %in% CONFIG$BIOMARKERS
  set.seed(CONFIG$SEED_MAIN)
  pdf(file.path(CONFIG$OUT_DIR, "moduleJ_TFmiRNA_network.pdf"), width = 12, height = 9)
  vcol <- ifelse(V(g)$type, "#C05746",
                 ifelse(V(g)$name %in% hub$regulator[hub$is_hub], "#3A6EA5", "#8C8C8C"))
  plot(g, layout = layout_with_fr(g), vertex.size = ifelse(V(g)$type, 7, 4),
       vertex.color = vcol, vertex.label.cex = ifelse(V(g)$type, 0.8, 0.55),
       edge.arrow.size = 0.25, edge.color = adjustcolor("grey60", 0.5),
       main = "Predicted upstream TF / miRNA regulators of the 10-biomarker panel")
  legend("bottomleft", bty = "n", cex = 0.8, pch = 19,
         col = c("#C05746", "#3A6EA5", "#8C8C8C"),
         legend = c("biomarker (10)", "hub regulator", "other regulator"))
  dev.off()
  log_msg("网络图已重出：moduleJ_TFmiRNA_network.pdf")
}
log_msg("13b 完成。接下来可跑模块 K（14_moduleK_CellChat细胞通讯_20260914.R）。")
