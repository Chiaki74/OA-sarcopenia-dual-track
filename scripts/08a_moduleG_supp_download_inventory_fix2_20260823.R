#!/usr/bin/env Rscript
# ============================================================================
# 08a_moduleG_supp_download_inventory_fix2_20260823.R
# ----------------------------------------------------------------------------
# 模块 G / G0 脚本：下载并选择性解压三个 GEO 单细胞/单核补充包，
# 生成可追溯文件清单（inventory / tar contents / download status / 日志）。
#
# 版本日期: 2026-08-23 fix2 — 修复 tar -tvf 输出在中文 Windows 下含无效 UTF-8 字节导致 trimws 崩溃（解析块移入 tryCatch + iconv 清洗，大小校验可退回 NA）
#
# 【为什么不走 01 脚本】
#   01 脚本的 `.done` 断点逻辑在 DOWNLOAD_SUPP=TRUE 时仍会因已有 .done 标记
#   而跳过补下大文件；且本脚本刻意不使用 GEOquery::getGEOSuppFiles。
#   这里用 base R 直接下载指定的补充文件，互不干扰，用户无需修改 01 脚本。
#
# 【四个下载资源】（总下载约 2.2GB，需稳定网络与足够磁盘空间）
#   1. GSE152805_RAW.tar            329,338,880 bytes
#   2. GSE169454_RAW.tar          1,684,111,360 bytes
#   3. GSE167186_RAW.tar            292,526,080 bytes
#   4. GSE167186_SimplifiedMetadataSheet.xlsx   10,124 bytes
#   注意：刻意不下载 GSE169454 的 ~1.3G
#   GSE169454_OA_gene_expression_matrix.tsv.gz（走 RAW tar 内 filtered 三件套）。
#
# 【选择性解压口径】（解压到 data/GEO/<GSE>/supp_extracted/）
#   - GSE152805: 全部 GSM* 的 .barcodes.tsv.gz/.genes.tsv.gz/.matrix.mtx.gz
#                （预计 27 个 = 9 样本 x 3 件套）
#   - GSE169454: 仅文件名含 filtered_ 的 GSM* features/barcodes/matrix 三件套
#                （预计 21 个 = 7 样本 x 3），绝不提取 raw_ 文件
#   - GSE167186: 仅 *_filtered_feature_bc_matrix.h5（预计 17 个）；
#                csv.gz 成员只登记在 tar contents，不提取
#   - 先 untar(list=TRUE) 列成员，再 untar(files=selected) 选择性解压；
#     防路径穿越（拒绝含 ".." / 绝对路径的成员）；解压后整理为平铺。
#
# 【幂等】
#   tar 大小与 expected 完全一致 → 重跑跳过下载；tar 不完整 → 删除重下。
#   解压跳过需目标 basename 精确集合逐一通过存在/非零/大小校验（大小取自
#   tar listing），任何缺失或不符都重新选择性解压；额外文件不掩盖目标缺失。
#
# 【输出】
#   results/moduleG/moduleG_file_inventory.csv
#   results/moduleG/moduleG_tar_contents.csv
#   results/moduleG/moduleG_download_status.csv
#   results/moduleG/moduleG_G0_log_<Sys.Date()>.txt
# ============================================================================

options(timeout = 7200)  # 大文件下载需长超时

# ---- 路径设置（假定从项目根目录运行） ---------------------------------------
proj_root  <- getwd()
geo_dir    <- file.path(proj_root, "data", "GEO")
res_dir    <- file.path(proj_root, "results", "moduleG")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(res_dir, paste0("moduleG_G0_log_", Sys.Date(), ".txt"))

log_msg <- function(..., .sep = "") {
  txt <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ",
                paste(..., collapse = .sep))
  cat(txt, "\n")
  cat(txt, "\n", file = log_file, append = TRUE)
  invisible(txt)
}

log_msg("=== 模块 G G0 启动：GEO 单细胞补充包下载与清单 ===")
log_msg("提示：总下载量约 2.2GB，请确保网络稳定、磁盘空间充足。")

# ---- 资源定义 ---------------------------------------------------------------
resources <- list(
  list(
    id  = "GSE152805_RAW",
    gse = "GSE152805",
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE152nnn/GSE152805/suppl/GSE152805_RAW.tar",
    expected_bytes = 329338880,
    kind = "tar",
    extract_pred = function(m) {
      grepl("^GSM[0-9]+", m) &
        grepl("\\.(barcodes\\.tsv\\.gz|genes\\.tsv\\.gz|matrix\\.mtx\\.gz)$", m)
    },
    expected_extract_n = 27
  ),
  list(
    id  = "GSE169454_RAW",
    gse = "GSE169454",
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169454/suppl/GSE169454_RAW.tar",
    expected_bytes = 1684111360,
    kind = "tar",
    extract_pred = function(m) {
      grepl("^GSM[0-9]+", m) & grepl("filtered_", m, fixed = TRUE) &
        grepl("(features|barcodes|matrix)\\.(tsv|mtx)\\.gz$", m)
    },
    expected_extract_n = 21
  ),
  list(
    id  = "GSE167186_RAW",
    gse = "GSE167186",
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167186/suppl/GSE167186_RAW.tar",
    expected_bytes = 292526080,
    kind = "tar",
    extract_pred = function(m) {
      # 严格要求 GSM 前缀 + _filtered_feature_bc_matrix.h5 结尾；
      # 选中数 != 17 时由 extract_tar_selective() 的严格数量检查判 FAIL
      grepl("^GSM[0-9]+", m) & grepl("_filtered_feature_bc_matrix\\.h5$", m)
    },
    expected_extract_n = 17
  ),
  list(
    id  = "GSE167186_metadata_xlsx",
    gse = "GSE167186",
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167186/suppl/GSE167186_SimplifiedMetadataSheet.xlsx",
    expected_bytes = 10124,
    kind = "file",
    extract_pred = NULL,
    expected_extract_n = 0
  )
)

# ---- 通用工具函数 -----------------------------------------------------------

# 路径穿越安全检查：先把反斜杠归一化为 "/"，再拒绝 ".." 段、
# 绝对路径、UNC（\\server\share）与盘符路径（C:/...）
is_safe_member <- function(m) {
  m <- gsub("\\", "/", m, fixed = TRUE)
  if (grepl("^//", m)) return(FALSE)               # UNC: \\server\share -> //server/share
  if (grepl("^(/|~)", m)) return(FALSE)            # 绝对路径 / home 路径
  if (grepl("^[A-Za-z]:", m)) return(FALSE)        # Windows 盘符路径
  if (grepl("(^|/)\\.\\.(/|$)", m)) return(FALSE)  # 任何 ".." 路径段
  TRUE
}

# tar contents 统一列 schema（0 行空表）：所有 extract_tar_selective() 返回
# 路径（含失败路径）都必须使用同一套列，避免下游 do.call(rbind, ...) 崩溃
empty_tar_contents <- function() {
  data.frame(
    gse_tar        = character(),
    member         = character(),
    member_base    = character(),
    member_size    = numeric(),
    extracted      = logical(),
    safe_member    = logical(),
    extracted_path = character(),
    stringsAsFactors = FALSE
  )
}

# 精确校验目标文件集合：逐一检查存在、非零；tar listing 大小可得时逐项比对。
# 不依赖目录内文件总数——额外文件不能掩盖目标缺失。
verify_extracted_set <- function(target_bases, target_sizes, out_dir) {
  if (length(target_bases) == 0) return(FALSE)
  tgt <- file.path(out_dir, target_bases)
  fsz <- file.size(tgt)                      # 不存在时返回 NA
  if (any(is.na(fsz)) || any(fsz <= 0)) return(FALSE)
  known <- !is.na(target_sizes)
  if (any(known) && any(fsz[known] != target_sizes[known])) return(FALSE)
  TRUE
}

# 下载到 .part，大小校验通过后才改名；失败删除 .part。
# 多通道回退（应对 Windows 上 NCBI HTTPS 的 "SSL connect error"）：
#   1) R download.file 多 method 依次尝试（Windows 先 wininet 再 libcurl；
#      非 Windows 先 libcurl 再 internal）；
#   2) R 方法全失败时尝试系统 curl（https URL 及 ftp:// 变体，Windows 加
#      --ssl-no-revoke 绕过 CRL 吊销检查问题）；
#   3) 任何方法产生的 .part 必须 file.size == expected_bytes 才 rename；
#      失败/0 字节/不完整一律删除 .part 再试下一方法。
download_one <- function(url, dest, expected_bytes) {
  tmp <- paste0(dest, ".part")
  if (file.exists(dest) && file.size(dest) == expected_bytes) {
    log_msg("已存在且大小匹配，跳过下载: ", basename(dest))
    return(list(action = "skipped", actual_bytes = file.size(dest),
                ok = TRUE, note = "existing file size matched"))
  }
  if (file.exists(dest) && file.size(dest) != expected_bytes) {
    log_msg("WARN: 已存在文件大小不符 (", file.size(dest),
            " != ", expected_bytes, ")，删除重下: ", basename(dest))
    unlink(dest)
  }
  if (file.exists(tmp)) unlink(tmp)
  log_msg("开始下载: ", url)

  is_win <- .Platform$OS.type == "windows"
  errors <- character(0)          # 每个方法的错误记录
  success_method <- NA_character_

  # .part 是否已精确完整（唯一成功判据）
  tmp_complete <- function() {
    file.exists(tmp) && !is.na(file.size(tmp)) &&
      file.size(tmp) == expected_bytes
  }
  # 某方法失败后记录原因并清理 .part
  record_failure <- function(tag, why) {
    log_msg("ERROR [", tag, "]: ", why, " (", url, ")")
    errors <<- c(errors, paste0(tag, ": ", why))
    unlink(tmp)
  }

  # ---- 通道 1: R download.file 多 method 回退 ----
  r_methods <- if (is_win) c("wininet", "libcurl") else c("libcurl", "internal")
  for (m in r_methods) {
    if (file.exists(tmp)) unlink(tmp)
    log_msg("尝试 download.file(method=\"", m, "\"): ", url)
    ok <- tryCatch({
      download.file(url, tmp, mode = "wb", method = m, quiet = FALSE)
      TRUE
    }, error = function(e) {
      record_failure(m, paste0("download.file error: ", conditionMessage(e)))
      FALSE
    })
    if (ok && tmp_complete()) { success_method <- m; break }
    if (ok) {
      sz <- file.size(tmp)
      record_failure(m, paste0("incomplete/empty file (got ",
                               ifelse(is.na(sz), "NA", sz),
                               " bytes, expected ", expected_bytes, ")"))
    }
  }

  # ---- 通道 2: 系统 curl 回退（R 方法全部失败时） ----
  if (is.na(success_method)) {
    curl_bin <- Sys.which("curl")
    if (nzchar(curl_bin)) {
      curl_variants <- list(list(tag = "curl-https", u = url))
      ftp_url <- sub("^https://", "ftp://", url)
      if (ftp_url != url) {
        curl_variants[[length(curl_variants) + 1L]] <-
          list(tag = "curl-ftp", u = ftp_url)
        log_msg("已准备 ftp 变体 URL 作为 curl 备选: ", ftp_url)
      }
      for (v in curl_variants) {
        if (file.exists(tmp)) unlink(tmp)
        log_msg("尝试系统 curl (", v$tag, "): ", v$u)
        args <- c("-L", "--fail", "--retry", "5", "--retry-delay", "5",
                  "--connect-timeout", "60", "-o", tmp, v$u)
        if (is_win) args <- c("--ssl-no-revoke", args)  # 绕过 Windows CRL 吊销检查失败
        rc <- tryCatch(suppressWarnings(system2(curl_bin, args)),
                       error = function(e) {
                         record_failure(v$tag, paste0("system2 error: ",
                                                      conditionMessage(e)))
                         NA_integer_
                       })
        if (!is.na(rc) && rc == 0L && tmp_complete()) {
          success_method <- v$tag
          break
        }
        if (is.na(rc)) next  # 错误已在 record_failure 中记录
        if (rc != 0L) {
          record_failure(v$tag, paste0("curl exit code ", rc))
        } else {
          sz <- file.size(tmp)
          record_failure(v$tag, paste0("incomplete/empty file (got ",
                                       ifelse(is.na(sz), "NA", sz),
                                       " bytes, expected ", expected_bytes, ")"))
        }
      }
    } else {
      log_msg("WARN: Sys.which(\"curl\") 不可用，跳过系统 curl 回退通道。")
      errors <- c(errors, "system curl: not found on PATH")
    }
  }

  # ---- 全部通道失败：删除 .part，给出浏览器手动下载提示 ----
  if (is.na(success_method)) {
    unlink(tmp)
    note <- paste0(
      "download failed via all methods (",
      paste(errors, collapse = " | "),
      "); you may also manually download ", url,
      " via a web browser to ", dest,
      " (must be exactly ", expected_bytes,
      " bytes) and rerun this script - it will detect the file by byte size and skip downloading")
    log_msg("ERROR: ", note)
    return(list(action = "downloaded", actual_bytes = NA_real_, ok = FALSE,
                note = note))
  }

  # ---- 成功：精确字节数校验已通过，rename 为最终文件 ----
  sz <- file.size(tmp)
  renamed <- file.rename(tmp, dest)
  if (!renamed) {
    unlink(tmp)
    return(list(action = "downloaded", actual_bytes = sz, ok = FALSE,
                note = "failed to rename .part to final file; rerun this script"))
  }
  log_msg("下载完成并校验通过: ", basename(dest), " (", sz, " bytes, via ",
          success_method, ")")
  list(action = "downloaded", actual_bytes = sz, ok = TRUE,
       note = paste0("downloaded and size-verified via ", success_method))
}

# 选择性解压：list 成员 -> 谓词筛选 -> 安全检查 -> untar(files=) -> 平铺
extract_tar_selective <- function(tar_path, extract_pred, out_dir, expected_n) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  members <- tryCatch(untar(tar_path, list = TRUE),
                      error = function(e) {
                        log_msg("ERROR: untar(list=TRUE) 失败 ", basename(tar_path),
                                " : ", conditionMessage(e))
                        NULL
                      })
  if (is.null(members)) {
    return(list(extracted_n = 0L, contents = empty_tar_contents(), ok = FALSE))
  }

  basenames <- basename(members)

  # 尝试用系统 tar -tvf 取得成员大小；不可用时退回 NA
  sizes <- rep(NA_real_, length(members))
  tv <- tryCatch(system2("tar", c("-tvf", shQuote(tar_path)), stdout = TRUE,
                         stderr = FALSE),
                 error = function(e) character(0),
                 warning = function(w) character(0))
  if (length(tv) > 0) {
    # fix2：中文 Windows 下 system2 捕获的 tar 输出可能含无效 UTF-8 字节
    # （GBK 环境），原实现中 trimws 的 perl 正则会因校验失败而崩溃并拖垮
    # 整个解压。成员大小仅用于提取后校验，可退回 NA；故先 iconv 清洗、
    # 再把整个解析块包进 tryCatch。
    tryCatch({
      tv <- iconv(tv, from = "", to = "UTF-8", sub = "")
      tv[is.na(tv)] <- ""
      # 典型行: -rw-r--r-- 0/0  12345 2020-01-01 00:00 path/to/file
      parts <- strsplit(trimws(tv), "[[:space:]]+")
      sz <- suppressWarnings(as.numeric(vapply(parts, function(p) p[3], character(1))))
      nm <- vapply(parts, function(p) p[length(p)], character(1))
      idx <- match(members, nm)
      hit <- !is.na(idx) & !is.na(sz[idx])
      sizes[hit] <- sz[idx[hit]]
    }, error = function(e) {
      log_msg("WARN: tar -tvf 输出解析失败，成员大小记为 NA（不影响解压与校验）: ",
              conditionMessage(e))
    })
  }

  safe <- vapply(members, is_safe_member, logical(1))
  if (any(!safe)) {
    log_msg("WARN: ", sum(!safe), " 个 tar 成员未通过路径穿越安全检查，仅登记不解压。")
  }
  selected <- vapply(basenames, extract_pred, logical(1)) & safe
  log_msg(basename(tar_path), "：tar 成员 ", length(members),
          " 个，选中待解压 ", sum(selected), " 个（预期 ", expected_n, "）")

  # 严格数量检查：选中数不等于预期即 FAIL（不静默继续）
  count_mismatch <- expected_n > 0 && sum(selected) != expected_n
  if (count_mismatch) {
    log_msg("ERROR: ", basename(tar_path), " 选中成员数 ", sum(selected),
            " != 预期 ", expected_n, "，本资源判定 FAIL。")
  }

  # 平铺前检查：selected basename 不得重复，否则平铺后互相覆盖，直接 FAIL
  sel_bases <- basenames[selected]
  dup_bases <- unique(sel_bases[duplicated(sel_bases)])
  if (length(dup_bases) > 0) {
    log_msg("ERROR: ", basename(tar_path), " 选中成员 basename 重复，平铺会互相覆盖: ",
            paste(dup_bases, collapse = ", "))
    contents <- data.frame(
      gse_tar     = basename(tar_path),
      member      = members,
      member_base = basenames,
      member_size = sizes,
      extracted   = selected,
      safe_member = safe,
      extracted_path = NA_character_,
      stringsAsFactors = FALSE
    )
    return(list(extracted_n = 0L, contents = contents, ok = FALSE))
  }

  contents <- data.frame(
    gse_tar     = basename(tar_path),
    member      = members,
    member_base = basenames,
    member_size = sizes,
    extracted   = selected,
    safe_member = safe,
    stringsAsFactors = FALSE
  )

  # 幂等：对目标 basename 精确集合逐一校验（存在 + 非零 + 大小匹配 tar
  # listing），全部通过才跳过；任何缺失/零字节/大小不符都触发重新选择性解压。
  # 目录中的额外文件不会掩盖目标缺失（不比较文件总数）。
  target <- contents[contents$extracted, , drop = FALSE]
  skip_extract <- !count_mismatch && expected_n > 0 &&
    nrow(target) == expected_n &&
    verify_extracted_set(target$member_base, target$member_size, out_dir)
  if (skip_extract) {
    log_msg("目标精确集合 ", nrow(target), " 个文件全部通过存在/非零/大小校验，",
            "跳过解压: ", basename(out_dir))
    contents$extracted_path <- ifelse(contents$extracted,
                                      file.path(out_dir, contents$member_base),
                                      NA_character_)
    return(list(extracted_n = nrow(target), contents = contents, ok = TRUE))
  }
  if (expected_n > 0 && length(list.files(out_dir)) > 0 && !skip_extract) {
    log_msg("已提取内容未通过目标精确集合校验（缺失/零字节/大小不符），重新选择性解压: ",
            basename(out_dir))
  }

  if (any(selected)) {
    ok_ext <- tryCatch({
      untar(tar_path, files = members[selected], exdir = out_dir)
      TRUE
    }, error = function(e) {
      log_msg("ERROR: 选择性解压失败 ", basename(tar_path), " : ", conditionMessage(e))
      FALSE
    })
    if (!ok_ext) {
      contents$extracted_path <- NA_character_
      return(list(extracted_n = 0L, contents = contents, ok = FALSE))
    }

    # 平铺：若 untar 还原出子目录，把目标文件移动到 out_dir 顶层；
    # file.rename 失败不静默，直接 FAIL
    moved <- 0L
    rename_failed <- character(0)
    for (m in members[selected]) {
      src <- file.path(out_dir, m)
      dst <- file.path(out_dir, basename(m))
      if (file.exists(src) && normalizePath(src, mustWork = FALSE) !=
          normalizePath(dst, mustWork = FALSE)) {
        dir.create(dirname(dst), showWarnings = FALSE, recursive = TRUE)
        if (file.rename(src, dst)) {
          moved <- moved + 1L
        } else {
          rename_failed <- c(rename_failed, m)
        }
      }
    }
    if (length(rename_failed) > 0) {
      log_msg("ERROR: 平铺 file.rename 失败 ", length(rename_failed), " 个: ",
              paste(rename_failed, collapse = ", "))
      contents$extracted_path <- NA_character_
      return(list(extracted_n = 0L, contents = contents, ok = FALSE))
    }
    # 清理空子目录
    subdirs <- list.dirs(out_dir, recursive = TRUE, full.names = TRUE)
    subdirs <- subdirs[subdirs != out_dir]
    for (d in rev(subdirs)) {
      if (length(list.files(d, all.files = TRUE, no.. = TRUE)) == 0) unlink(d, recursive = TRUE)
    }
    if (moved > 0) log_msg("平铺整理完成：移动 ", moved, " 个文件到 ", out_dir)
  }

  extracted_files <- list.files(out_dir, full.names = FALSE)
  contents$extracted_path <- ifelse(
    contents$extracted & contents$member_base %in% extracted_files,
    file.path(out_dir, contents$member_base), NA_character_)
  # 最终 ok：数量严格匹配预期 + 目标精确集合通过存在/非零/大小校验
  final_ok <- !count_mismatch &&
    (expected_n == 0 ||
       (nrow(target) == expected_n &&
          verify_extracted_set(target$member_base, target$member_size, out_dir)))
  list(extracted_n = length(extracted_files), contents = contents, ok = final_ok)
}

# 读取 pData：合并 data/GEO/<GSE>/ 下全部 *_pData.csv；缺失不中断
read_pdata_all <- function(gse) {
  gse_dir <- file.path(geo_dir, gse)
  pfiles <- list.files(gse_dir, pattern = "_pData\\.csv$", full.names = TRUE)
  if (length(pfiles) == 0) {
    log_msg("WARN: ", gse, " 未找到 *_pData.csv，GSM 映射列将为 NA（不中断）。")
    return(NULL)
  }
  dfs <- lapply(pfiles, function(f) {
    tryCatch(utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE),
             error = function(e) {
               log_msg("WARN: 读取 pData 失败 ", basename(f), " : ", conditionMessage(e))
               NULL
             })
  })
  dfs <- Filter(Negate(is.null), dfs)
  if (length(dfs) == 0) return(NULL)
  # 列名不齐时做 fill 合并
  all_cols <- unique(unlist(lapply(dfs, names)))
  dfs <- lapply(dfs, function(d) {
    miss <- setdiff(all_cols, names(d))
    for (mc in miss) d[[mc]] <- NA
    d[all_cols]
  })
  do.call(rbind, dfs)
}

# 从 pData 中提取 GSM -> 提示信息 的映射
build_gsm_map <- function(pd) {
  if (is.null(pd)) return(NULL)
  gsm_col <- names(pd)[grepl("geo_accession|gsm", names(pd), ignore.case = TRUE)]
  if (length(gsm_col) == 0) {
    cand <- vapply(pd, function(x) any(grepl("^GSM[0-9]+", as.character(x))), logical(1))
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

# 文件名层面的 condition / tissue / donor 推断
infer_hints <- function(gse, fname) {
  condition_hint <- NA_character_
  tissue_site    <- NA_character_
  donor_hint     <- NA_character_
  if (gse == "GSE152805") {
    condition_hint <- "OA"
    # 通用正则：覆盖 SY/oLT/MT × 113/116/118 全部 9 个样本（如
    # GSMxxx_OA_SY_113_...、GSMxxx_OA_oLT_116_...、GSMxxx_OA_MT_118_...）
    m <- regmatches(fname, regexpr("(SY|oLT|MT)_(113|116|118)", fname,
                                   perl = TRUE))
    if (length(m) == 1L) {
      parts <- strsplit(m, "_", fixed = TRUE)[[1]]
      tissue_site <- if (parts[1] == "SY") "synovium" else parts[1]
      donor_hint  <- parts[2]
    }
  } else if (gse == "GSE169454") {
    if (grepl("normal", fname, ignore.case = TRUE)) condition_hint <- "normal"
    else if (grepl("oa", fname, ignore.case = TRUE)) condition_hint <- "OA"
  } else if (gse == "GSE167186") {
    # 文件名只含 HM 编号；condition 靠 pData/xlsx 映射，此处留 NA
    condition_hint <- NA_character_
  }
  list(condition_hint = condition_hint, tissue_site = tissue_site,
       donor_hint = donor_hint)
}

# 三件套角色判定
role_of_file <- function(fname) {
  # 先判 .h5（如 *_filtered_feature_bc_matrix.h5 同时含 "matrix"），再判三件套
  if (grepl("\\.h5$", fname, ignore.case = TRUE)) return("h5")
  if (grepl("barcodes", fname, ignore.case = TRUE)) return("barcodes")
  if (grepl("genes", fname, ignore.case = TRUE)) return("genes")
  if (grepl("features", fname, ignore.case = TRUE)) return("features")
  if (grepl("matrix", fname, ignore.case = TRUE)) return("matrix")
  if (grepl("\\.xlsx$", fname, ignore.case = TRUE)) return("metadata")
  "other"
}

# ---- 主流程 ------------------------------------------------------------------
status_rows   <- list()
tar_contents  <- list()
inventory_rows <- list()

for (res in resources) {
  log_msg("---- 处理资源: ", res$id, " ----")
  gse_dir <- file.path(geo_dir, res$gse)
  dir.create(gse_dir, showWarnings = FALSE, recursive = TRUE)
  dest <- file.path(gse_dir, basename(res$url))

  dl <- tryCatch(
    download_one(res$url, dest, res$expected_bytes),
    error = function(e) {
      log_msg("ERROR: ", res$id, " 下载异常（不影响其他资源）: ", conditionMessage(e))
      list(action = "downloaded", actual_bytes = NA_real_, ok = FALSE,
           note = paste0("unexpected error: ", conditionMessage(e), "; rerun this script"))
    })

  extracted_n <- 0L
  note <- dl$note
  res_ok <- dl$ok

  if (dl$ok && res$kind == "tar") {
    out_dir <- file.path(gse_dir, "supp_extracted")
    ext <- tryCatch(
      extract_tar_selective(dest, res$extract_pred, out_dir, res$expected_extract_n),
      error = function(e) {
        log_msg("ERROR: ", res$id, " 解压异常（不影响其他资源）: ", conditionMessage(e))
        list(extracted_n = 0L, contents = empty_tar_contents(), ok = FALSE)
      })
    extracted_n <- ext$extracted_n
    if (!is.null(ext$contents)) tar_contents[[res$id]] <- ext$contents
    if (!ext$ok) {
      res_ok <- FALSE
      note <- paste0(note, "; extraction incomplete: got ", extracted_n,
                     ", expected ", res$expected_extract_n)
    }

    # 构建提取文件清单行
    if (!is.null(ext$contents) && any(ext$contents$extracted, na.rm = TRUE)) {
      pd  <- read_pdata_all(res$gse)
      gsm_map <- build_gsm_map(pd)
      sel <- ext$contents[!is.na(ext$contents$extracted) & ext$contents$extracted, ]
      for (i in seq_len(nrow(sel))) {
        fn <- sel$member_base[i]
        gsm <- regmatches(fn, regexpr("GSM[0-9]+", fn))
        if (length(gsm) == 0) gsm <- NA_character_
        hints <- infer_hints(res$gse, fn)
        pdata_cond <- NA_character_
        if (!is.null(gsm_map) && !is.na(gsm)) {
          hit <- gsm_map$pdata_hint[gsm_map$gsm == gsm]
          if (length(hit) > 0) pdata_cond <- hit[1]
        }
        cond <- hints$condition_hint
        if (is.na(cond) && !is.na(pdata_cond) && nzchar(pdata_cond)) cond <- pdata_cond
        fp <- ifelse(is.na(sel$extracted_path[i]), NA,
                     file.path(sel$extracted_path[i]))
        inventory_rows[[length(inventory_rows) + 1L]] <- data.frame(
          gse           = res$gse,
          gsm           = gsm,
          file_name     = fn,
          file_type     = sub("^.*\\.", "", fn),
          role_in_matrix = role_of_file(fn),
          condition_hint = cond,
          tissue_site   = hints$tissue_site,
          donor_hint    = hints$donor_hint,
          extracted_path = ifelse(is.na(fp), NA, fp),
          size_bytes    = ifelse(file.exists(fp), file.size(fp), NA),
          md5_optional  = NA_character_,
          stringsAsFactors = FALSE
        )
      }
    }
  } else if (dl$ok && res$kind == "file") {
    # 元数据 xlsx：仅登记清单，不解压
    inventory_rows[[length(inventory_rows) + 1L]] <- data.frame(
      gse = res$gse, gsm = NA_character_, file_name = basename(dest),
      file_type = "xlsx", role_in_matrix = "metadata",
      condition_hint = NA_character_, tissue_site = NA_character_,
      donor_hint = NA_character_, extracted_path = dest,
      size_bytes = file.size(dest), md5_optional = NA_character_,
      stringsAsFactors = FALSE)
    log_msg("元数据文件已登记清单: ", basename(dest))
  }

  status_rows[[length(status_rows) + 1L]] <- data.frame(
    resource    = res$id,
    gse         = res$gse,
    url         = res$url,
    expected_bytes = res$expected_bytes,
    actual_bytes   = ifelse(is.null(dl$actual_bytes), NA, dl$actual_bytes),
    action      = dl$action,   # downloaded / skipped
    extracted_n = extracted_n,
    expected_extract_n = res$expected_extract_n,
    status      = ifelse(res_ok, "OK", "FAIL"),
    note        = note,
    stringsAsFactors = FALSE
  )
}

# ---- 输出 --------------------------------------------------------------------
status_df <- do.call(rbind, status_rows)
utils::write.csv(status_df, file.path(res_dir, "moduleG_download_status.csv"),
                 row.names = FALSE)
log_msg("写出: moduleG_download_status.csv (", nrow(status_df), " 行)")

if (length(tar_contents) > 0) {
  tar_df <- do.call(rbind, tar_contents)
  utils::write.csv(tar_df, file.path(res_dir, "moduleG_tar_contents.csv"),
                   row.names = FALSE)
  log_msg("写出: moduleG_tar_contents.csv (", nrow(tar_df), " 行)")
} else {
  utils::write.csv(
    data.frame(gse_tar = character(), member = character(),
               member_base = character(), member_size = numeric(),
               extracted = logical(), safe_member = logical(),
               extracted_path = character()),
    file.path(res_dir, "moduleG_tar_contents.csv"), row.names = FALSE)
  log_msg("WARN: 无任何 tar contents 可写出（可能全部解压失败）。")
}

if (length(inventory_rows) > 0) {
  inv_df <- do.call(rbind, inventory_rows)
  utils::write.csv(inv_df, file.path(res_dir, "moduleG_file_inventory.csv"),
                   row.names = FALSE)
  log_msg("写出: moduleG_file_inventory.csv (", nrow(inv_df), " 行)")
} else {
  utils::write.csv(
    data.frame(gse = character(), gsm = character(), file_name = character(),
               file_type = character(), role_in_matrix = character(),
               condition_hint = character(), tissue_site = character(),
               donor_hint = character(), extracted_path = character(),
               size_bytes = numeric(), md5_optional = character()),
    file.path(res_dir, "moduleG_file_inventory.csv"), row.names = FALSE)
  log_msg("WARN: 文件清单为空（可能所有资源均失败）。")
}

# ---- 汇总与完成提示 -----------------------------------------------------------
log_msg("==== 汇总 ====")
for (res in resources[ vapply(resources, function(x) x$kind, character(1)) == "tar" ]) {
  r <- status_df[status_df$resource == res$id, ]
  log_msg(res$gse, ": extracted_n = ", r$extracted_n,
          " / expected = ", r$expected_extract_n)
}

if (all(status_df$status == "OK")) {
  msg <- "模块 G G0 完成，下一步可进入 Seurat 建对象脚本 08_moduleG_singlecell.R"
  cat(msg, "\n"); log_msg(msg)
} else {
  bad <- status_df[status_df$status != "OK", ]
  msg <- paste0("模块 G G0 存在失败资源，请检查以下 URL/文件后重跑：\n",
                paste0("  - ", bad$resource, " | ", bad$url, " | ", bad$note,
                       collapse = "\n"))
  cat(msg, "\n"); log_msg(msg)
  stop("G0 未全部通过，详见 results/moduleG/moduleG_download_status.csv 与日志。")
}
