# ============================================================================
# 13_moduleJ_TFmiRNA调控网络_20260914.R
# 模块 J｜v5 修订新增分析（结构对照拍板项 E4，合并 Wang 借鉴项 C10，备忘录 §66/§69）
# ----------------------------------------------------------------------------
# 产出：10 标志物的上游 TF 与 miRNA 调控网络（预测性、数据库驱动，入补充材料）：
#   Part 1：Enrichr 在线库富集——TRRUST_Transcription_Factors_2019 与 ChEA_2022
#     （TF 双库互证）、miRTarBase_2017（实验验证 miRNA-靶基因）；
#   Part 2：边表整理（TF→基因、miRNA→基因），双库 TF 交集为高置信层；
#   Part 3：hub 调控因子（调控 ≥3/10 标志物者）+ 网络图 PDF（igraph，种子固定）。
# 纪律：全部为公共数据库预测，写稿口径为「上游调控假设生成」，不作机制断言；
#   任一库不可用即如实降级并在日志留痕，不静默跳过。
# 前置：仅需联网（Enrichr API）与 10 标志物清单（06 封板口径，硬编码）。
# 运行：source("13_moduleJ_TFmiRNA调控网络_20260914.R")
# 输出：results/moduleJ/
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ----------------------------- CONFIG --------------------------------------
CONFIG <- list(
  OUT_DIR  = "results/moduleJ",
  LOG_FILE = file.path("results/moduleJ", sprintf("moduleJ_log_%s.txt", Sys.Date())),
  BIOMARKERS = c("NELL1","STEAP1","GADD45A","FBLN5","HEMK1",
                 "GDE1","RNF14","BNIP3","EGFR","CBR3"),   # 06 封板，勿改
  TF_LIBS    = c("TRRUST_Transcription_Factors_2019", "ChEA_2022"),
  MIR_LIBS   = c("miRTarBase_2017"),
  HUB_MIN    = 3,          # hub 判定：调控 ≥3/10 标志物
  PADJ_CUT   = 0.05,
  SEED_MAIN  = 20260914
)
dir.create(CONFIG$OUT_DIR, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  txt <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = ""))
  cat(txt, "\n"); cat(txt, "\n", file = CONFIG$LOG_FILE, append = TRUE)
}

log_msg("—— 模块 J：TF-miRNA 上游调控网络（E4，10 标志物，06 封板口径）——")

## 硬自检 0：enrichR 可用性与联网
if (!requireNamespace("enrichR", quietly = TRUE))
  stop("缺少 enrichR 包——install.packages(\"enrichR\") 后重跑")
suppressPackageStartupMessages(library(enrichR))
avail <- tryCatch(listEnrichrDbs()$libraryName,
                  error = function(e) character(0))
if (length(avail) == 0)
  stop("Enrichr API 不可达（无网络或被拦截）——本模块为在线数据库分析，",
       "请检查网络后重跑；若长期不可达，在备忘录标注本模块降级为「未执行」")
miss_libs <- setdiff(c(CONFIG$TF_LIBS, CONFIG$MIR_LIBS), avail)
if (length(miss_libs) > 0)
  stop("Enrichr 库名变更，找不到：", paste(miss_libs, collapse = ", "),
       "——当前可用库已写入日志，请人工挑选等价库后改 CONFIG")
log_msg("Enrichr 在线，目标库全部就位：",
        paste(c(CONFIG$TF_LIBS, CONFIG$MIR_LIBS), collapse = " | "))

# ---- Part 1：逐库富集 ----
query_one <- function(lib) {
  res <- tryCatch(
    enrichr(CONFIG$BIOMARKERS, lib)[[lib]],
    error = function(e) { log_msg("  ⚠ 库 ", lib, " 查询失败：",
                                  conditionMessage(e)); NULL })
  if (is.null(res) || nrow(res) == 0) { log_msg("  ⚠ 库 ", lib, " 无返回"); return(NULL) }
  res <- as.data.frame(res)
  res$library <- lib
  log_msg("  ", lib, "：", nrow(res), " 条（adj.P<", CONFIG$PADJ_CUT, " 者 ",
          sum(res$Adjusted.P.value < CONFIG$PADJ_CUT), " 条）")
  res
}
log_msg("Part 1：TF 与 miRNA 富集查询")
tf_tab  <- do.call(rbind, Filter(Negate(is.null), lapply(CONFIG$TF_LIBS,  query_one)))
mir_tab <- do.call(rbind, Filter(Negate(is.null), lapply(CONFIG$MIR_LIBS, query_one)))
if (is.null(tf_tab) && is.null(mir_tab))
  stop("TF 与 miRNA 库全部无返回——网络或库名问题，请回查，不产出空结果")

# ---- Part 2：边表整理（解析 Genes 列，仅保留 10 标志物内的靶）----
to_edges <- function(tab, kind) {
  if (is.null(tab)) return(data.frame())
  sig <- tab[tab$Adjusted.P.value < CONFIG$PADJ_CUT, , drop = FALSE]
  if (nrow(sig) == 0) { log_msg("  ", kind, "：无 adj.P<", CONFIG$PADJ_CUT,
                                " 条目——如实阴性，降级输出名义显著前 10"); 
    sig <- head(tab[order(tab$Adjusted.P.value), ], 10); sig$nominal_only <- TRUE
  } else sig$nominal_only <- FALSE
  edges <- do.call(rbind, lapply(seq_len(nrow(sig)), function(i) {
    gs <- intersect(strsplit(sig$Genes[i], ";", fixed = TRUE)[[1]], CONFIG$BIOMARKERS)
    if (length(gs) == 0) return(NULL)
    data.frame(regulator = sig$Term[i], target = gs, kind = kind,
               library = sig$library[i], adjP = sig$Adjusted.P.value[i],
               nominal_only = sig$nominal_only[i], stringsAsFactors = FALSE)
  }))
  if (is.null(edges)) data.frame() else unique(edges)
}
edges <- rbind(to_edges(tf_tab, "TF"), to_edges(mir_tab, "miRNA"))
## 硬自检 1：必须存在边
if (nrow(edges) == 0)
  stop("显著调节因子无一靶向 10 标志物——罕见，请回查富集表原始 Genes 列")
## 双库 TF 交集为高置信层
tf_lib_n <- tapply(edges$library[edges$kind == "TF"],
                   edges$regulator[edges$kind == "TF"],
                   function(x) length(unique(x)))
edges$tf_dual <- edges$kind == "TF" & edges$regulator %in% names(tf_lib_n[tf_lib_n >= 2])
fwrite(edges, file.path(CONFIG$OUT_DIR, "moduleJ_TFmiRNA_edges.csv"))
log_msg("Part 2：边表 ", nrow(edges), " 条（TF ", sum(edges$kind == "TF"),
        " / miRNA ", sum(edges$kind == "miRNA"), "；双库 TF 高置信 ",
        length(unique(edges$regulator[edges$tf_dual])), " 个）→ moduleJ_TFmiRNA_edges.csv")

## 覆盖度（诚实报告，不设硬停）
cov <- sapply(CONFIG$BIOMARKERS, function(g) any(edges$target == g & !edges$nominal_only))
log_msg("标志物覆盖：", sum(cov), "/10 有显著上游调节因子",
        ifelse(sum(cov) < 10,
               paste0("（未覆盖：", paste(CONFIG$BIOMARKERS[!cov], collapse = ","), "——如实报告）"),
               ""))

# ---- Part 3：hub 与网络图 ----
deg <- tapply(edges$target[!edges$nominal_only], edges$regulator[!edges$nominal_only],
              function(x) length(unique(x)))
hub <- data.frame(regulator = names(deg), n_targets = as.integer(deg),
                  kind = edges$kind[match(names(deg), edges$regulator)],
                  tf_dual = names(deg) %in% names(tf_lib_n[tf_lib_n >= 2]))
hub <- hub[order(-hub$n_targets), ]
hub$is_hub <- hub$n_targets >= CONFIG$HUB_MIN
fwrite(hub, file.path(CONFIG$OUT_DIR, "moduleJ_TFmiRNA_hub.csv"))
log_msg("Part 3：hub（≥", CONFIG$HUB_MIN, "/10 靶）", sum(hub$is_hub), " 个：",
        paste(head(hub$regulator[hub$is_hub], 15), collapse = ", "),
        " → moduleJ_TFmiRNA_hub.csv")

if (requireNamespace("igraph", quietly = TRUE)) {
  suppressPackageStartupMessages(library(igraph))
  e_sig <- edges[!edges$nominal_only, ]
  g <- graph_from_data_frame(e_sig[, c("regulator", "target")], directed = TRUE)
  V(g)$type <- V(g)$name %in% CONFIG$BIOMARKERS
  set.seed(CONFIG$SEED_MAIN)
  pdf(file.path(CONFIG$OUT_DIR, "moduleJ_TFmiRNA_network.pdf"), width = 12, height = 9)
  lay <- layout_with_fr(g)
  vcol <- ifelse(V(g)$type, "#C05746", ifelse(V(g)$name %in% hub$regulator[hub$kind == "TF" & hub$is_hub],
                 "#3A6EA5", "#8C8C8C"))
  plot(g, layout = lay, vertex.size = ifelse(V(g)$type, 7, 4),
       vertex.color = vcol, vertex.label.cex = ifelse(V(g)$type, 0.8, 0.55),
       vertex.label.color = "black", edge.arrow.size = 0.25,
       edge.color = adjustcolor("grey60", 0.5), main = "Predicted upstream TF / miRNA regulators of the 10-biomarker panel")
  legend("bottomleft", bty = "n", cex = 0.8, pch = 19,
         col = c("#C05746", "#3A6EA5", "#8C8C8C"),
         legend = c("biomarker (10)", "hub TF (dual-library)", "other regulator"))
  dev.off()
  log_msg("网络图已出：moduleJ_TFmiRNA_network.pdf")
} else log_msg("⚠ 无 igraph 包——网络图跳过（install.packages(\"igraph\") 后可单补），边表/Hub 表不受影响")

log_msg("模块 J 完成。产物：moduleJ_TFmiRNA_edges.csv / moduleJ_TFmiRNA_hub.csv",
        " / moduleJ_TFmiRNA_network.pdf（若有 igraph）")
log_msg("写稿口径：数据库预测性上游调控，入补充材料；hub 与双库 TF 交集可在讨论",
        "机制段以「假设生成」一句带过，不作因果/机制断言。")
