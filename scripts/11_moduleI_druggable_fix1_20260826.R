#!/usr/bin/env Rscript
# ============================================================================
# 11_moduleI_druggable_20260826.R
# ----------------------------------------------------------------------------
# 模块 I：可成药基因组比对 + 候选药物富集（DSigDB ORA）。
# 编码: UTF-8（中文 Windows 用户请确认 RStudio 以 UTF-8 打开）
#
# 版本日期: 2026-08-26 fix1 — 门禁审查修复（模块 H 真实输出结构以
#   10_moduleH_MR_coloc_fix3_20260826.R 为唯一事实依据）：
#   B1(blocker): extract_causal_genes 重构为三档降级链——
#     档1 主路径读 moduleH_known_causal_overlap.csv 的 causal_established_any；
#     档2 跨文件组合回退（MR_summary x coloc_summary 按 gene||outcome 关联，
#     口径 = sig_bonferroni + steiger_p<0.05 + PPH4>0.8，NA 不当 TRUE）；
#     档3 保底 markers_only。三档命中数与启用档位均 log。
#   B2(major): 删除 pos() 数值列 >0.5 语义；逻辑列优先，数值列回退语义为
#     p 值类 <0.05 / PPH4 >0.8；列名匹配不依赖列序。
#   m1: norm_gene 后过滤改 x[!is.na(x) & nzchar(x)]（宇宙构建与 GMT 背景）；
#   m2: 模块 H CSV 读加 fileEncoding="UTF-8"；readLines(gmt, encoding="UTF-8")；
#   m3: log 与 ORA 输出显式声明 "FDR 校正族 = 有重叠(k>=1)的条目"；
#   m4: data.table 从必需依赖移除（脚本实际未使用）。
# v1 — 正式版。范围已收敛：不做调控网络/CMap/分子对接/MD。
#
# 【输入】（一律 file.path 相对项目根目录拼接；用户机器 Windows + R 4.6.0，
#   工作目录 C:\Users\凉月\Desktop\ENG论文）
#   1. 10 个锁定标志物（硬编码兜底）：
#      NELL1/STEAP1/GADD45A/FBLN5/HEMK1/GDE1/RNF14/BNIP3/EGFR/CBR3
#   2. 可选模块 H 结果（存在则读，不存在则降级为仅标志物模式并 log）：
#      results/moduleH/moduleH_MR_summary.csv / moduleH_coloc_summary.csv
#      ——列结构未知：先 read.csv 读列名再适配（causal_established 列优先；
#      否则 Bonferroni+Steiger+PPH4 组合列；再不行 WARN 降级）。
#   3. 可成药宇宙（data/druggable/ 下同名文件，列结构已在起草时实查确认）：
#      - druggable_genome_Finan2017.csv（Finan 2017，4,463 唯一基因；
#        基因列 = hgnc_names，另带 druggability_tier）
#      - 可成药基因组_DGIdb类别全量_20260820.csv（DGIdb 全量 11,665 行；
#        基因列 = gene_symbol，类别列 = categories；
#        【关键口径】须先按 categories 含 "DRUGGABLE GENOME" 过滤，
#        过滤后 5,809 基因——实查验证：union(Finan, DGIdb-DG) = 5,928，
#        intersect = 4,344，与既定口径一致）
#   4. DSigDB 基因集：data/druggable/DSigDB_All.gmt
#      ——起草时代下验证：Enrichr 当前库名为 "DSigDB"（"DSigDB_All" 已 404），
#      下载地址 https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=DSigDB
#      实得 2,991,990 字节、4,026 个药物条目、制表符分隔 GMT，已存
#      /mnt/agents/output/moduleI_data/DSigDB_All.gmt 供交付拷贝。
#
# 【分析逻辑】
#   1. 宇宙合并：union(Finan, DGIdb[DRUGGABLE GENOME])；运行时 log 实际
#      并集/交集数，与口径 5,928/4,344 对照，不符则 WARN（不 stop）；
#      逐基因标注 in_universe/in_Finan/in_DGIdb/finan_tier/dgidb_categories。
#   2. 比对：标志物（+H 因果基因，若有）∩ 宇宙；命中表含 evidence_tier
#      （"H因果+标志物"/"仅H因果"/"仅标志物"）；Finan-only 敏感性口径重复判定并 log。
#   3. DSigDB ORA：手写 Fisher 精确检验（单侧 greater），BH 校正；
#      背景 = 可成药宇宙 ∩ DSigDB 覆盖基因（log 背景大小）；
#      FDR<0.05 为主，不足则报告 top20 原始 p 并在表/日志如实标注。
#   4. 可视化：top15 候选药物条形图 PDF（ggplot2；条形=-log10FDR，标命中基因数）。
#
# 【输出】（results/moduleI/；write.csv 一律 fileEncoding="UTF-8"）
#   moduleI_druggable_overlap.csv / moduleI_dsigdb_enrichment.csv /
#   moduleI_summary.csv / moduleI_drug_candidates_barplot.pdf /
#   moduleI_status.csv / moduleI_log_<date>.txt（控制台同步）
#
# 【工程纪律】无 system2 / setwd / 绝对路径；幂等覆盖写；命名向量赋值一律
#   unname()；网络访问仅限 §0 的可选存在性检查（base R url()，不下载、
#   失败仅 WARN）；sessionInfo 落盘。
# ============================================================================

# ---- 路径设置（假定从项目根目录运行） ---------------------------------------------
proj_root <- getwd()
drug_dir  <- file.path(proj_root, "data", "druggable")
res_dir   <- file.path(proj_root, "results", "moduleI")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(res_dir, paste0("moduleI_log_", Sys.Date(), ".txt"))

log_msg <- function(..., .sep = "") {
  txt <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ",
                paste(..., collapse = .sep))
  cat(txt, "\n")
  cat(txt, "\n", file = log_file, append = TRUE)
  invisible(txt)
}

log_msg("=== 模块 I 启动：可成药基因组比对 + DSigDB 候选药物富集 ===")
log_msg("proj_root = ", proj_root)

# ---- 依赖检查（缺包即 stop 并打印确切安装命令） --------------------------------------
check_deps <- function() {
  # fix1-m4：data.table 从必需依赖移除（本脚本实际未使用 fread/fwrite，
  # 全部读写走 utils::read.csv / write.csv），自检清单与实际保持一致
  reqs <- list(
    list(pkg = "ggplot2",    cmd = 'install.packages("ggplot2")')
  )
  problems <- character(0)
  for (r in reqs) {
    if (!requireNamespace(r$pkg, quietly = TRUE)) {
      problems <- c(problems, paste0("  缺包 ", r$pkg, " -> 安装: ", r$cmd))
    }
  }
  if (length(problems) > 0) {
    log_msg("ERROR: 依赖不满足:\n", paste(problems, collapse = "\n"))
    stop(paste0("模块 I 依赖不满足，请先执行：\n",
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
  library(ggplot2)
})

# ---- 常量 ---------------------------------------------------------------------------
MARKERS10 <- c("NELL1", "STEAP1", "GADD45A", "FBLN5", "HEMK1",
               "GDE1", "RNF14", "BNIP3", "EGFR", "CBR3")
EXPECTED_UNION     <- 5928   # 既定口径（起草时实查验证）
EXPECTED_INTERSECT <- 4344
DSIGDB_URL <- paste0("https://maayanlab.cloud/Enrichr/geneSetLibrary?",
                     "mode=text&libraryName=DSigDB")  # 注意: "DSigDB_All" 已 404

status_rows <- list()
add_status <- function(step, status, note = "") {
  status_rows[[length(status_rows) + 1L]] <<- data.frame(
    step = step, status = status, note = note, stringsAsFactors = FALSE)
  invisible(TRUE)
}

# ---- §0 自检：输入文件精确清单，缺即 stop ---------------------------------------------
check_inputs <- function() {
  need_files <- c(
    file.path(drug_dir, "druggable_genome_Finan2017.csv"),
    file.path(drug_dir, "可成药基因组_DGIdb类别全量_20260820.csv"),
    file.path(drug_dir, "DSigDB_All.gmt")
  )
  missing <- need_files[!file.exists(need_files)]
  if (length(missing) > 0) {
    log_msg("ERROR: 模块 I 缺输入文件:\n",
            paste0("  - ", missing, collapse = "\n"))
    stop(paste0(
      "模块 I 缺输入文件，请将以下文件放到 data/druggable/ 后重跑：\n",
      paste0("  - ", basename(missing), collapse = "\n"),
      "\n（DSigDB_All.gmt 可由 ", DSIGDB_URL, " 下载）"), call. = FALSE)
  }
  for (f in need_files) {
    log_msg("§0 输入在位: ", basename(f), " (",
            round(file.size(f) / 1e6, 2), " MB)")
  }
  # 可选网络存在性检查（不下载；失败仅 WARN，不影响主流程）
  net_note <- tryCatch({
    con <- url(DSIGDB_URL, open = "rb")
    on.exit(close(con), add = TRUE)
    invisible(readBin(con, what = "raw", n = 64L))
    "reachable"
  }, error = function(e) paste0("unreachable: ", conditionMessage(e)))
  log_msg("§0 DSigDB 源 URL 存在性检查（可选，不下载）: ", net_note)
  add_status("S0_input_check", "OK",
             paste0("3 files present; DSigDB url ", net_note))
  invisible(TRUE)
}
check_inputs()

# ---- 1. 可成药宇宙合并 ----------------------------------------------------------------
norm_gene <- function(x) toupper(trimws(as.character(x)))

build_universe <- function() {
  # Finan 2017（列结构实查：基因列 hgnc_names，另带 druggability_tier）
  f_finan <- file.path(drug_dir, "druggable_genome_Finan2017.csv")
  finan <- tryCatch(
    # fileEncoding="UTF-8"：description 列含非 ASCII 字符，GBK Windows 上
    # 不显式声明可能报 invalid multibyte string
    utils::read.csv(f_finan, stringsAsFactors = FALSE, check.names = FALSE,
                    fileEncoding = "UTF-8"),
    error = function(e) stop("Finan2017 读取失败: ", conditionMessage(e),
                             call. = FALSE))
  if (!"hgnc_names" %in% names(finan)) {
    stop("Finan2017 缺基因列 hgnc_names（现有列: ",
         paste(names(finan), collapse = ", "), "）", call. = FALSE)
  }
  finan_genes <- unique(norm_gene(finan$hgnc_names))
  # fix1-m1：nzchar(NA) 产生 NA 索引会静默引入 NA 行，必须显式 !is.na
  finan_genes <- finan_genes[!is.na(finan_genes) & nzchar(finan_genes)]
  finan_tier  <- stats::setNames(finan$druggability_tier,
                                 norm_gene(finan$hgnc_names))
  log_msg("Finan2017: ", nrow(finan), " 行 -> 唯一基因 ", length(finan_genes),
          "（口径预期 4,463）")

  # DGIdb 类别全量（列结构实查：gene_symbol / categories；
  # 【关键口径】只取 categories 含 "DRUGGABLE GENOME" 的行，全量 11,665 行 ->
  # 5,809 基因；用全量会把宇宙错扩到 11,784）
  f_dgidb <- file.path(drug_dir, "可成药基因组_DGIdb类别全量_20260820.csv")
  dgidb <- tryCatch(
    utils::read.csv(f_dgidb, stringsAsFactors = FALSE, check.names = FALSE,
                    fileEncoding = "UTF-8"),
    error = function(e) stop("DGIdb 类别全量读取失败: ", conditionMessage(e),
                             call. = FALSE))
  if (!all(c("gene_symbol", "categories") %in% names(dgidb))) {
    stop("DGIdb 类别全量缺列（需 gene_symbol/categories；现有列: ",
         paste(names(dgidb), collapse = ", "), "）", call. = FALSE)
  }
  dg_rows <- grepl("DRUGGABLE GENOME", dgidb$categories, fixed = TRUE)
  dg_dg <- dgidb[dg_rows & !is.na(dg_rows), , drop = FALSE]
  dgidb_genes <- unique(norm_gene(dg_dg$gene_symbol))
  dgidb_genes <- dgidb_genes[!is.na(dgidb_genes) & nzchar(dgidb_genes)]  # fix1-m1
  # 类别合并（同一基因多行 -> 去重拼接）
  dgidb_cat <- vapply(dgidb_genes, function(g) {
    paste(unique(dg_dg$categories[norm_gene(dg_dg$gene_symbol) == g]),
          collapse = "; ")
  }, character(1))
  log_msg("DGIdb: 全量 ", nrow(dgidb), " 行 -> DRUGGABLE GENOME 类别 ",
          length(dgidb_genes), " 基因（口径预期 5,809）")

  universe <- sort(union(finan_genes, dgidb_genes))
  inter <- length(intersect(finan_genes, dgidb_genes))
  log_msg("宇宙合并: union = ", length(universe), " | intersect = ", inter,
          "（口径预期 ", EXPECTED_UNION, "/", EXPECTED_INTERSECT, "）")
  if (length(universe) != EXPECTED_UNION || inter != EXPECTED_INTERSECT) {
    log_msg("WARN: 宇宙规模与既定口径不符（实际 ", length(universe), "/",
            inter, " vs 预期 ", EXPECTED_UNION, "/", EXPECTED_INTERSECT,
            "）——不 stop，请人工核对源文件版本。")
    add_status("universe_build", "WARN",
               paste0("union=", length(universe), " intersect=", inter,
                      " expected ", EXPECTED_UNION, "/", EXPECTED_INTERSECT))
  } else {
    add_status("universe_build", "OK",
               paste0("union=", length(universe), " intersect=", inter))
  }

  uni_df <- data.frame(
    gene = universe,
    in_universe = TRUE,
    in_Finan = universe %in% finan_genes,
    in_DGIdb = universe %in% dgidb_genes,
    finan_tier = unname(finan_tier[universe]),
    dgidb_categories = unname(dgidb_cat[universe]),
    stringsAsFactors = FALSE
  )
  list(uni_df = uni_df, finan_genes = finan_genes, dgidb_genes = dgidb_genes,
       universe = universe)
}
UNI <- build_universe()

# ---- 2. 可选模块 H 因果基因（fix1-B1/B2：三档降级链；列结构以模块 H fix3 定稿为准） ----
# 模块 H 真实输出（results/moduleH/，UTF-8）：
#   moduleH_known_causal_overlap.csv: gene, ..., causal_established_outcomes,
#     causal_established_any(逻辑)   ——档1 主路径
#   moduleH_MR_summary.csv (长表 基因x结局xeQTL源): gene, outcome, method, pval,
#     p_bonferroni, sig_bonferroni(逻辑), steiger_p(数值), steiger_correct_dir(逻辑)
#   moduleH_coloc_summary.csv: gene, outcome, PPH4(数值), PPH3,
#     pass_relaxed_pph4_gt_0.8(逻辑), pass_strict_pph4_gt_0.9(逻辑), status
#   ——档2 跨文件组合：某基因存在任一结局同时满足 sig_bonferroni(或
#     p_bonferroni<0.05) + steiger_p<0.05 + pass_relaxed_pph4_gt_0.8(或 PPH4>0.8)
extract_causal_genes <- function() {
  h_dir <- file.path(proj_root, "results", "moduleH")
  f_overlap <- file.path(h_dir, "moduleH_known_causal_overlap.csv")
  f_mr      <- file.path(h_dir, "moduleH_MR_summary.csv")
  f_coloc   <- file.path(h_dir, "moduleH_coloc_summary.csv")

  # fix1-m2：模块 H CSV 一律 UTF-8 读取
  read_h <- function(f) {
    tryCatch(utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE,
                             fileEncoding = "UTF-8"),
             error = function(e) {
               log_msg("WARN: 读取模块 H 文件失败 ", basename(f), " : ",
                       conditionMessage(e))
               NULL
             })
  }
  # 逻辑列真值判定：NA 一律不当 TRUE；接受 logical / 数值 1 / 常见字符串
  as_true <- function(v) {
    if (is.logical(v)) return(!is.na(v) & v)
    if (is.numeric(v)) return(!is.na(v) & v == 1)
    !is.na(v) & toupper(trimws(as.character(v))) %in%
      c("TRUE", "YES", "1", "PASS")
  }
  tier3_fallback <- function(reason) {
    log_msg("档3(保底): ", reason, "——模块 H 结果不可用，本次为标志物模式。")
    add_status("moduleH_causal", "DEGRADED",
               paste0("tier3 markers_only; ", reason))
    list(genes = character(0), mode = "markers_only")
  }

  # ---- 档 1（主路径）：known_causal_overlap 的 causal_established_any ----
  if (file.exists(f_overlap)) {
    df <- read_h(f_overlap)
    if (!is.null(df) && nrow(df) > 0 &&
        all(c("gene", "causal_established_any") %in% names(df))) {
      hit <- as_true(df$causal_established_any)
      causal <- unique(norm_gene(df$gene[hit]))
      causal <- causal[!is.na(causal) & nzchar(causal)]
      log_msg("档1(主路径): ", basename(f_overlap),
              " causal_established_any=TRUE 基因 ", length(causal), " 个",
              if (length(causal) > 0)
                paste0(": ", paste(causal, collapse = ", "))
              else "（H 无因果命中，属正常结果）")
      if (length(causal) > 0) {
        add_status("moduleH_causal", "OK",
                   paste0("tier1 known_causal_overlap; n=", length(causal)))
        return(list(genes = causal, mode = "markers_plus_H"))
      }
      # 文件存在但全 FALSE：正常结果，如实 log 后进入档 3（不静默跳过）
      return(tier3_fallback("tier1 file present but causal_established_any all FALSE (0 hits)"))
    }
    log_msg("WARN: ", basename(f_overlap),
            " 存在但读取失败或缺 gene/causal_established_any 列（现有列: ",
            if (is.null(df)) "< unreadable >" else paste(names(df), collapse = ", "),
            "）-> 回退档2。")
  } else {
    log_msg("档1 不可用: ", basename(f_overlap), " 不存在 -> 回退档2。")
  }

  # ---- 档 2（跨文件组合回退）：MR_summary x coloc_summary 按 gene||outcome 关联 ----
  # fix1-B2：一律优先逻辑列；逻辑列缺失才回退数值列，语义 p 值类 <0.05 /
  # PPH4 >0.8；列名匹配一律按名称，不依赖文件内列序
  if (file.exists(f_mr) && file.exists(f_coloc)) {
    mr <- read_h(f_mr)
    co <- read_h(f_coloc)
    need_mr <- c("gene", "outcome", "steiger_p")
    need_co <- c("gene", "outcome")
    mr_sig_ok <- any(c("sig_bonferroni", "p_bonferroni", "pval") %in% names(mr))
    co_sig_ok <- any(c("pass_relaxed_pph4_gt_0.8", "PPH4") %in% names(co))
    if (!is.null(mr) && !is.null(co) && nrow(mr) > 0 && nrow(co) > 0 &&
        all(need_mr %in% names(mr)) && all(need_co %in% names(co)) &&
        mr_sig_ok && co_sig_ok) {
      log_msg("档2: MR 列 = ", paste(names(mr), collapse = " | "))
      log_msg("档2: coloc 列 = ", paste(names(co), collapse = " | "))
      # MR 显著：逻辑列 sig_bonferroni 优先，数值回退 p_bonferroni<0.05 -> pval<0.05
      sig_ok <- if ("sig_bonferroni" %in% names(mr)) {
        as_true(mr$sig_bonferroni)
      } else if ("p_bonferroni" %in% names(mr)) {
        !is.na(mr$p_bonferroni) & mr$p_bonferroni < 0.05
      } else {
        !is.na(mr$pval) & mr$pval < 0.05
      }
      st_ok <- !is.na(mr$steiger_p) & mr$steiger_p < 0.05
      # coloc 显著：逻辑列 pass_relaxed_pph4_gt_0.8 优先，数值回退 PPH4>0.8
      co_ok <- if ("pass_relaxed_pph4_gt_0.8" %in% names(co)) {
        as_true(co$pass_relaxed_pph4_gt_0.8)
      } else {
        !is.na(co$PPH4) & co$PPH4 > 0.8
      }
      key_mr <- paste(norm_gene(mr$gene), mr$outcome, sep = "||")
      key_co <- paste(norm_gene(co$gene), co$outcome, sep = "||")
      both <- intersect(key_mr[sig_ok & st_ok], key_co[co_ok])
      causal <- unique(sub("\\|\\|.*$", "", both))
      causal <- causal[!is.na(causal) & nzchar(causal)]
      log_msg("档2(跨文件组合): MR 三条件满足键 ", sum(sig_ok & st_ok),
              " 个；coloc 满足键 ", sum(co_ok), " 个；gene||outcome 交集 ",
              length(both), " 个 -> 因果基因 ", length(causal), " 个",
              if (length(causal) > 0)
                paste0(": ", paste(causal, collapse = ", ")) else "")
      if (length(causal) > 0) {
        add_status("moduleH_causal", "OK",
                   paste0("tier2 MR+coloc combined; n=", length(causal)))
        return(list(genes = causal, mode = "markers_plus_H"))
      }
      return(tier3_fallback("tier2 combined criteria yielded 0 causal genes"))
    }
    return(tier3_fallback("tier2 inputs unreadable or required columns missing"))
  }
  tier3_fallback("tier1/tier2 files absent")
}
H <- extract_causal_genes()
log_msg("模块 H 因果基因提取档位完成，模式 = ", H$mode, "，因果基因 ",
        length(H$genes), " 个。")

# ---- 3. 标志物(+H 因果基因) ∩ 宇宙 ------------------------------------------------------
overlap_analysis <- function() {
  input_genes <- unique(c(MARKERS10, H$genes))
  tier <- ifelse(input_genes %in% MARKERS10 & input_genes %in% H$genes,
                 "H因果+标志物",
                 ifelse(input_genes %in% H$genes, "仅H因果", "仅标志物"))
  idx <- match(input_genes, UNI$uni_df$gene)   # 未命中 -> NA
  ov <- data.frame(
    gene = input_genes,
    evidence_tier = unname(tier),
    in_universe = !is.na(idx),
    in_Finan  = ifelse(is.na(idx), FALSE, UNI$uni_df$in_Finan[idx]),
    in_DGIdb  = ifelse(is.na(idx), FALSE, UNI$uni_df$in_DGIdb[idx]),
    finan_tier = ifelse(is.na(idx), NA, UNI$uni_df$finan_tier[idx]),
    dgidb_categories = ifelse(is.na(idx), NA, UNI$uni_df$dgidb_categories[idx]),
    stringsAsFactors = FALSE
  )
  n_hit <- sum(ov$in_universe)
  n_hit_finan <- sum(ov$in_Finan)
  log_msg("比对: 输入基因 ", nrow(ov), " 个（标志物 ", length(MARKERS10),
          " + H 因果 ", length(H$genes), "，模式 ", H$mode, "）")
  log_msg("命中宇宙 ", n_hit, "/", nrow(ov), "；Finan-only 敏感性口径命中 ",
          n_hit_finan, "/", nrow(ov))
  log_msg("未命中基因: ",
          paste(ov$gene[!ov$in_universe], collapse = ", "))
  utils::write.csv(ov, file.path(res_dir, "moduleI_druggable_overlap.csv"),
                   row.names = FALSE, fileEncoding = "UTF-8")
  log_msg("写出: moduleI_druggable_overlap.csv (", nrow(ov), " 行)")
  add_status("overlap", "OK",
             paste0("hits=", n_hit, "/", nrow(ov),
                    " finan_only=", n_hit_finan, " mode=", H$mode))
  list(overlap = ov, input_genes = input_genes, n_hit = n_hit,
       n_hit_finan = n_hit_finan)
}
OV <- overlap_analysis()

# ---- 4. DSigDB ORA（手写 Fisher 精确检验 + BH） -----------------------------------------
parse_gmt <- function(gmt_path) {
  lines <- readLines(gmt_path, encoding = "UTF-8", warn = FALSE)   # fix1-m2
  sets <- lapply(lines, function(ln) {
    f <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    f <- f[nzchar(f)]
    if (length(f) < 3L) return(NULL)
    g <- unique(norm_gene(f[-c(1, 2)]))
    g <- g[!is.na(g) & nzchar(g)]                                  # fix1-m1
    if (length(g) == 0L) return(NULL)
    list(term = f[1], genes = g)
  })
  sets <- Filter(Negate(is.null), sets)
  log_msg("GMT 解析: ", basename(gmt_path), " -> ", length(sets), " 个条目")
  sets
}

dsigdb_ora <- function() {
  gmt <- parse_gmt(file.path(drug_dir, "DSigDB_All.gmt"))
  gmt_genes_all <- unique(unlist(lapply(gmt, function(s) s$genes),
                                 use.names = FALSE))
  # 背景 = 可成药宇宙 ∩ DSigDB 覆盖基因
  background <- intersect(UNI$universe, gmt_genes_all)
  query <- intersect(OV$input_genes, background)
  dropped <- setdiff(OV$input_genes, background)
  if (length(dropped) > 0) {
    log_msg("WARN: ", length(dropped), " 个输入基因不在 ORA 背景内，已剔除: ",
            paste(dropped, collapse = ", "))
  }
  N <- length(background); n <- length(query)
  log_msg("ORA 背景 = 宇宙 ∩ DSigDB 覆盖 = ", N, " 基因；查询基因 ", n, " 个。")
  if (n < 2) stop("背景内查询基因 < 2，无法做富集。", call. = FALSE)

  res <- lapply(gmt, function(s) {
    term_bg <- intersect(s$genes, background)
    K <- length(term_bg)
    if (K == 0) return(NULL)
    hit <- intersect(term_bg, query)
    k <- length(hit)
    if (k == 0) return(NULL)
    # 2x2: [k, K-k; n-k, N-K-n+k]，单侧 greater
    a <- k; b <- K - k; cc <- n - k; dd <- N - K - n + k
    if (dd < 0) return(NULL)
    p <- tryCatch(
      stats::fisher.test(matrix(c(a, b, cc, dd), nrow = 2),
                         alternative = "greater")$p.value,
      error = function(e) NA_real_)
    data.frame(drug = s$term, genes_hit = paste(hit, collapse = ";"),
               overlap = k, term_genes_in_bg = K, p_value = p,
               stringsAsFactors = FALSE)
  })
  res <- Filter(Negate(is.null), res)
  if (length(res) == 0) stop("DSigDB 无任何条目与查询基因重叠。", call. = FALSE)
  df <- do.call(rbind, res)
  df <- df[!is.na(df$p_value), , drop = FALSE]
  df$p_adjust_BH <- stats::p.adjust(df$p_value, method = "BH")
  df <- df[order(df$p_value), , drop = FALSE]
  n_sig <- sum(df$p_adjust_BH < 0.05)
  # fix1-m3：显式声明 BH 校正族口径（仅与查询有重叠 k>=1 的条目），
  # 写日志与输出表双留痕
  log_msg("DSigDB ORA: 测了 ", nrow(df), " 个有重叠条目；FDR<0.05 显著 ",
          n_sig, " 个。【FDR 校正族 = 有重叠(k>=1)的条目，共 ", nrow(df),
          " 条】")
  df$fdr_family <- paste0("overlap(k>=1) terms only; family size = ",
                          nrow(df))
  if (n_sig == 0) {
    log_msg("WARN: 无 FDR<0.05 条目；如实报告 top20 原始 p（探索性，",
            "不得作为显著性结论）: top5 = ",
            paste(sprintf("%s(p=%.3g)", df$drug[1:min(5, nrow(df))],
                          df$p_value[1:min(5, nrow(df))]), collapse = ", "))
    df$report_basis <- ifelse(seq_len(nrow(df)) <= 20,
                              "top20 raw p (no FDR<0.05 hit)", "")
  } else {
    df$report_basis <- ifelse(df$p_adjust_BH < 0.05, "FDR<0.05", "")
  }
  utils::write.csv(df, file.path(res_dir, "moduleI_dsigdb_enrichment.csv"),
                   row.names = FALSE, fileEncoding = "UTF-8")
  log_msg("写出: moduleI_dsigdb_enrichment.csv (", nrow(df), " 行)")
  add_status("dsigdb_ora", "OK",
             paste0("tested=", nrow(df), " bg=", N, " sig_fdr005=", n_sig))
  list(df = df, n_sig = n_sig, bg = N)
}
ORA <- dsigdb_ora()

# ---- 5. top15 候选药物条形图（PDF） -------------------------------------------------------
make_barplot <- function() {
  df <- ORA$df
  use_fdr <- ORA$n_sig > 0
  ord <- if (use_fdr) order(df$p_adjust_BH) else order(df$p_value)
  top <- df[ord, , drop = FALSE][1:min(15, nrow(df)), , drop = FALSE]
  top$neglog10 <- if (use_fdr) -log10(top$p_adjust_BH) else -log10(top$p_value)
  ylab <- if (use_fdr) "-log10(FDR)" else "-log10(raw p)  [no FDR<0.05 hit]"
  top$drug <- factor(top$drug, levels = rev(top$drug))
  p <- ggplot(top, aes(x = drug, y = neglog10)) +
    geom_col(fill = "#2C7FB8", width = 0.7) +
    geom_text(aes(label = paste0(overlap, " genes")), hjust = -0.1, size = 3) +
    coord_flip() +
    labs(title = paste0("Module I top", nrow(top), " DSigDB drug candidates"),
         x = NULL, y = ylab) +
    theme_bw(base_size = 11) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  out <- file.path(res_dir, "moduleI_drug_candidates_barplot.pdf")
  tryCatch({
    ggplot2::ggsave(out, p, width = 8.5, height = max(4, 0.42 * nrow(top) + 1.5),
                    limitsize = FALSE)
    log_msg("写出: moduleI_drug_candidates_barplot.pdf (top", nrow(top), ")")
    add_status("barplot", "OK", paste0("top", nrow(top), " by ",
                                       if (use_fdr) "FDR" else "raw p"))
  }, error = function(e) {
    log_msg("WARN: 条形图失败（不中断）: ", conditionMessage(e))
    add_status("barplot", "WARN", conditionMessage(e))
  })
  invisible(TRUE)
}
make_barplot()

# ---- 6. 一行式 summary + status + sessionInfo -------------------------------------------
summary_row <- data.frame(
  mode = H$mode,                                   # markers_only / markers_plus_H
  n_markers = length(MARKERS10),
  n_H_causal_genes = length(H$genes),
  n_input_genes = length(OV$input_genes),
  universe_size = length(UNI$universe),
  finan_size = length(UNI$finan_genes),
  dgidb_DG_size = length(UNI$dgidb_genes),
  universe_intersect = length(intersect(UNI$finan_genes, UNI$dgidb_genes)),
  hits_universe = OV$n_hit,
  hits_finan_only = OV$n_hit_finan,
  dsigdb_background = ORA$bg,
  dsigdb_terms_tested = nrow(ORA$df),
  n_sig_fdr005 = ORA$n_sig,
  top_drug = as.character(ORA$df$drug[1]),
  note = ifelse(ORA$n_sig > 0,
                "FDR<0.05 hits available; see moduleI_dsigdb_enrichment.csv",
                "no FDR<0.05 hit; top20 raw p reported (exploratory)"),
  stringsAsFactors = FALSE
)
utils::write.csv(summary_row, file.path(res_dir, "moduleI_summary.csv"),
                 row.names = FALSE, fileEncoding = "UTF-8")
log_msg("写出: moduleI_summary.csv（模式 ", H$mode, "；命中宇宙 ", OV$n_hit,
        "/", length(OV$input_genes), "；显著条目 ", ORA$n_sig, "）")

status_df <- do.call(rbind, status_rows)
utils::write.csv(status_df, file.path(res_dir, "moduleI_status.csv"),
                 row.names = FALSE, fileEncoding = "UTF-8")
log_msg("写出: moduleI_status.csv (", nrow(status_df), " 行)")

log_msg("==== 汇总 ====")
log_msg("模式: ", H$mode, " | 宇宙: ", length(UNI$universe),
        " | 命中: ", OV$n_hit, " | Finan-only 命中: ", OV$n_hit_finan,
        " | ORA 背景: ", ORA$bg, " | FDR<0.05: ", ORA$n_sig)

log_msg("==== sessionInfo ====")
for (ln in capture.output(utils::sessionInfo())) log_msg(ln)

if (all(status_df$status == "OK")) {
  log_msg("模块 I 完成，全部步骤 OK。")
} else {
  log_msg("模块 I 完成（含 WARN/DEGRADED 步骤，详见 moduleI_status.csv）。")
}
