# ==============================================================================
# fix14v4_CIBERSORTx混合物制备_20260903.R  （替代 fix14v3：GSE48556 同样改用内置对照表）
# 变更原因（2026-09-03 核对 GEO 原始 pData 后）：
#   GSE55235 的 disease state 字段在 GEO 上对 OA_2–OA_10 发生错位，moduleA 提取的
#   processed pheno$condition 因此把 20 个样本错标 other（fix14v2 在此硬停止，属设计
#   胜利）。v3 不再信任 processed pheno 的 condition 列，改为【内置 GSM→分组对照表】，
#   该表由 GEO 官方 pData 的 title 字段逐样本核对生成（healthy/normal→control、
#   osteoarthrit*（含 typo osteoarthrisits）→case、rheumatoid→other），口径与模块 B
#   完全一致：GSE55235 30→20、GSE55457 33→20、GSE12021 31→19、GSE48556 139 全保留，
#   合计 198。
#   v4 追加（2026-09-03）：GSE48556 的 processed pheno$condition 同样不可信（106 个 OA
#   case 被错标 other，fix14v3 在此硬停止）。GEO 官方 pData 真相 = 106 OA case +
#   33 control，无第三类，139 全保留。v4 起四个队列【全部】走内置 GSM 对照表，
#   processed pheno 仅用于 ID 对齐核验与不一致留痕。
# 其余全部沿用 v2：processed/<GSE>/ 数据惯例、log2(x+1)→2^x-1 反变换、基因去重、
#   expect_n 硬校验、md5 manifest、CIBERSORTx 官方格式写出。
# 红线：①不读取/不复制/不上传 LM22；②产物 txt 仅供上传 CIBERSORTx 官网，不入 GitHub；
#       ③去卷积输入必须是线性尺度——严禁 z-score / log 值直接输入。
# 运行：在工程根目录 Rscript fix14v4_CIBERSORTx混合物制备_20260903.R
# ==============================================================================

suppressPackageStartupMessages({ library(tools) })

## ---- 0. 配置区 -----------------------------------------------------------------
PROJ_ROOT <- getwd()
OUT_DIR   <- file.path(PROJ_ROOT, "results", "moduleF", "fix14_cibersortx")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

COHORTS <- list(
  GSE55235 = list(platform = "GPL96",   expect_n = 20),
  GSE55457 = list(platform = "GPL96",   expect_n = 20),
  GSE12021 = list(platform = "GPL96",   expect_n = 19),
  GSE48556 = list(platform = "GPL6947", expect_n = 139)
)
expr_candidates <- function(gse) c(
  file.path(PROJ_ROOT, "processed", gse, paste0(gse, "_expr_gene.rds")),
  file.path(PROJ_ROOT, "results", "moduleA", paste0(gse, "_expr_normalized.rds")),
  file.path(PROJ_ROOT, "results", "moduleA", gse, paste0(gse, "_expr_gene.rds"))
)
pheno_candidates <- function(gse) c(
  file.path(PROJ_ROOT, "processed", gse, paste0(gse, "_pheno.csv")),
  file.path(PROJ_ROOT, "data", paste0(gse, "_pData.csv"))
)

ALLOW_LOG2_INVERSE <- TRUE

## ---- 0b. 内置 GSM→分组对照表（GEO 官方 pData title 逐样本核对，2026-09-03） -----
## 四个队列以此为唯一分组依据（覆盖 processed pheno）；pheno 仅用于 ID 对齐与留痕。
GSM_CONDITION <- list(
  GSE55235 = c(
    GSM1332201 = "control",
    GSM1332202 = "control",
    GSM1332203 = "control",
    GSM1332204 = "control",
    GSM1332205 = "control",
    GSM1332206 = "control",
    GSM1332207 = "control",
    GSM1332208 = "control",
    GSM1332209 = "control",
    GSM1332210 = "control",
    GSM1332211 = "case",
    GSM1332212 = "case",
    GSM1332213 = "case",
    GSM1332214 = "case",
    GSM1332215 = "case",
    GSM1332216 = "case",
    GSM1332217 = "case",
    GSM1332218 = "case",
    GSM1332219 = "case",
    GSM1332220 = "case",
    GSM1332221 = "other",
    GSM1332222 = "other",
    GSM1332223 = "other",
    GSM1332224 = "other",
    GSM1332225 = "other",
    GSM1332226 = "other",
    GSM1332227 = "other",
    GSM1332228 = "other",
    GSM1332229 = "other",
    GSM1332230 = "other"
  ),
  GSE55457 = c(
    GSM1337304 = "control",
    GSM1337305 = "control",
    GSM1337306 = "control",
    GSM1337307 = "control",
    GSM1337308 = "control",
    GSM1337309 = "control",
    GSM1337310 = "control",
    GSM1337311 = "control",
    GSM1337312 = "control",
    GSM1337313 = "control",
    GSM1337314 = "other",
    GSM1337315 = "other",
    GSM1337316 = "other",
    GSM1337317 = "other",
    GSM1337318 = "other",
    GSM1337319 = "other",
    GSM1337320 = "other",
    GSM1337321 = "other",
    GSM1337322 = "other",
    GSM1337323 = "other",
    GSM1337324 = "other",
    GSM1337325 = "other",
    GSM1337326 = "other",
    GSM1337327 = "case",
    GSM1337328 = "case",
    GSM1337329 = "case",
    GSM1337330 = "case",
    GSM1337331 = "case",
    GSM1337332 = "case",
    GSM1337333 = "case",
    GSM1337334 = "case",
    GSM1337335 = "case",
    GSM1337336 = "case"
  ),
  GSE48556 = c(
    GSM1182186 = "case",
    GSM1182187 = "case",
    GSM1182188 = "case",
    GSM1182189 = "case",
    GSM1182190 = "case",
    GSM1182191 = "case",
    GSM1182192 = "case",
    GSM1182193 = "case",
    GSM1182194 = "case",
    GSM1182195 = "case",
    GSM1182196 = "case",
    GSM1182197 = "case",
    GSM1182198 = "case",
    GSM1182199 = "case",
    GSM1182200 = "case",
    GSM1182201 = "case",
    GSM1182202 = "case",
    GSM1182203 = "case",
    GSM1182204 = "case",
    GSM1182205 = "case",
    GSM1182206 = "case",
    GSM1182207 = "case",
    GSM1182208 = "case",
    GSM1182209 = "case",
    GSM1182210 = "case",
    GSM1182211 = "case",
    GSM1182212 = "case",
    GSM1182213 = "case",
    GSM1182214 = "case",
    GSM1182215 = "case",
    GSM1182216 = "case",
    GSM1182217 = "case",
    GSM1182218 = "case",
    GSM1182219 = "case",
    GSM1182220 = "case",
    GSM1182221 = "case",
    GSM1182222 = "case",
    GSM1182223 = "case",
    GSM1182224 = "case",
    GSM1182225 = "case",
    GSM1182226 = "case",
    GSM1182227 = "case",
    GSM1182228 = "case",
    GSM1182229 = "case",
    GSM1182230 = "case",
    GSM1182231 = "case",
    GSM1182232 = "case",
    GSM1182233 = "case",
    GSM1182234 = "case",
    GSM1182235 = "case",
    GSM1182236 = "case",
    GSM1182237 = "case",
    GSM1182238 = "case",
    GSM1182239 = "case",
    GSM1182240 = "case",
    GSM1182241 = "case",
    GSM1182242 = "case",
    GSM1182243 = "case",
    GSM1182244 = "case",
    GSM1182245 = "case",
    GSM1182246 = "case",
    GSM1182247 = "case",
    GSM1182248 = "case",
    GSM1182249 = "case",
    GSM1182250 = "case",
    GSM1182251 = "case",
    GSM1182252 = "case",
    GSM1182253 = "case",
    GSM1182254 = "case",
    GSM1182255 = "case",
    GSM1182256 = "case",
    GSM1182257 = "case",
    GSM1182258 = "case",
    GSM1182259 = "case",
    GSM1182260 = "case",
    GSM1182261 = "case",
    GSM1182262 = "case",
    GSM1182263 = "case",
    GSM1182264 = "case",
    GSM1182265 = "case",
    GSM1182266 = "case",
    GSM1182267 = "case",
    GSM1182268 = "case",
    GSM1182269 = "case",
    GSM1182270 = "case",
    GSM1182271 = "case",
    GSM1182272 = "case",
    GSM1182273 = "case",
    GSM1182274 = "case",
    GSM1182275 = "case",
    GSM1182276 = "case",
    GSM1182277 = "case",
    GSM1182278 = "case",
    GSM1182279 = "case",
    GSM1182280 = "case",
    GSM1182281 = "case",
    GSM1182282 = "case",
    GSM1182283 = "case",
    GSM1182284 = "case",
    GSM1182285 = "case",
    GSM1182286 = "case",
    GSM1182287 = "case",
    GSM1182288 = "case",
    GSM1182289 = "case",
    GSM1182290 = "case",
    GSM1182291 = "case",
    GSM1182292 = "control",
    GSM1182293 = "control",
    GSM1182294 = "control",
    GSM1182295 = "control",
    GSM1182296 = "control",
    GSM1182297 = "control",
    GSM1182298 = "control",
    GSM1182299 = "control",
    GSM1182300 = "control",
    GSM1182301 = "control",
    GSM1182302 = "control",
    GSM1182303 = "control",
    GSM1182304 = "control",
    GSM1182305 = "control",
    GSM1182306 = "control",
    GSM1182307 = "control",
    GSM1182308 = "control",
    GSM1182309 = "control",
    GSM1182310 = "control",
    GSM1182311 = "control",
    GSM1182312 = "control",
    GSM1182313 = "control",
    GSM1182314 = "control",
    GSM1182315 = "control",
    GSM1182316 = "control",
    GSM1182317 = "control",
    GSM1182318 = "control",
    GSM1182319 = "control",
    GSM1182320 = "control",
    GSM1182321 = "control",
    GSM1182322 = "control",
    GSM1182323 = "control",
    GSM1182324 = "control"
  ),
  GSE12021 = c(
    GSM302859 = "control",
    GSM302864 = "control",
    GSM302866 = "control",
    GSM302870 = "control",
    GSM302872 = "other",
    GSM302876 = "case",
    GSM302880 = "case",
    GSM302882 = "other",
    GSM302896 = "other",
    GSM302900 = "other",
    GSM302924 = "other",
    GSM302927 = "other",
    GSM302930 = "case",
    GSM302933 = "other",
    GSM302943 = "other",
    GSM303326 = "case",
    GSM303341 = "case",
    GSM303356 = "case",
    GSM303358 = "case",
    GSM303360 = "case",
    GSM303362 = "case",
    GSM303364 = "other",
    GSM303366 = "other",
    GSM303368 = "other",
    GSM303370 = "case",
    GSM303520 = "other",
    GSM303522 = "control",
    GSM303523 = "control",
    GSM303525 = "control",
    GSM303531 = "control",
    GSM303533 = "control"
  )
)

## ---- 1. 工具函数 ---------------------------------------------------------------
stopif <- function(cond, msg) if (cond) stop(msg, call. = FALSE)

locate_file <- function(cands, gse, what) {
  hit <- cands[file.exists(cands)][1]
  stopif(is.na(hit), sprintf("[%s] 找不到%s，已尝试：\n  %s\n  → 请确认 processed/%s/ 是否存在；不在就截图发回。",
                             gse, what, paste(cands, collapse = "\n  "), gse))
  hit
}

read_expr <- function(gse) {
  f <- locate_file(expr_candidates(gse), gse, "表达矩阵 rds")
  obj <- readRDS(f)
  if (is.list(obj) && !is.data.frame(obj) && !is.matrix(obj)) {
    cand <- intersect(c("exprs","expr","mat","data","normalized"), names(obj))
    stopif(length(cand) == 0, sprintf("[%s] RDS 是 list 但未找到表达矩阵槽位：%s", gse, paste(names(obj), collapse = ", ")))
    obj <- obj[[cand[1]]]
  }
  m <- as.matrix(obj)
  stopif(!is.numeric(m), sprintf("[%s] 表达矩阵不是数值型。", gse))
  stopif(is.null(rownames(m)), sprintf("[%s] 表达矩阵缺少基因名（rownames）。", gse))
  message("  读取：", sub(PROJ_ROOT, ".", f, fixed = TRUE))
  m
}

read_pheno <- function(gse) {
  f <- locate_file(pheno_candidates(gse), gse, "表型 csv")
  read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
}

## ---- 2. 主流程 -----------------------------------------------------------------
manifest <- list()
for (gse in names(COHORTS)) {
  cfg <- COHORTS[[gse]]
  message("========== ", gse, " (", cfg$platform, ") ==========")

  m <- read_expr(gse)

  ## 2.1 样本过滤（v4：四队列全部用内置 GSM 对照表重编码）
  pd <- read_pheno(gse)
  idcol <- intersect(c("geo_accession","sample","Sample","gsm","GSM","title"), colnames(pd))[1]
  stopif(is.na(idcol), sprintf("[%s] pheno 缺少样本 ID 列（现有列：%s）", gse, paste(colnames(pd), collapse = ", ")))
  pd_ids <- as.character(pd[[idcol]])
  stopif(!all(colnames(m) %in% pd_ids),
         sprintf("[%s] 表达矩阵样本名与 pheno$%s 对不上：%s",
                 gse, idcol, paste(setdiff(colnames(m), pd_ids)[1:3], collapse = ", ")))

  if (!is.null(GSM_CONDITION[[gse]])) {
    ## —— 内置对照表重编码（不再信任 processed pheno 的 condition） ——
    lut <- GSM_CONDITION[[gse]]
    stopif(!all(colnames(m) %in% names(lut)),
           sprintf("[%s] 表达矩阵样本名不在内置对照表中：%s——请截图发回核对。",
                   gse, paste(setdiff(colnames(m), names(lut))[1:3], collapse = ", ")))
    cond_new <- unname(lut[colnames(m)])
    cnt <- table(factor(cond_new, levels = c("case","control","other")))
    message(sprintf("  重编码（GEO title 对照表）：case=%d / control=%d / other=%d",
                    cnt["case"], cnt["control"], cnt["other"]))
    condcol <- intersect(c("condition","Condition","group","Group"), colnames(pd))[1]
    if (!is.na(condcol)) {
      old <- pd[[condcol]][match(colnames(m), pd_ids)]
      n_diff <- sum(old != cond_new, na.rm = TRUE)
      message(sprintf("  对照 processed pheno$%s：%d 个样本分组不一致（已知 GEO 字段错位，以内置表为准）",
                      condcol, n_diff))
    }
    keep_ids <- colnames(m)[cond_new %in% c("case","control")]
    dropped  <- setdiff(colnames(m), keep_ids)
    m <- m[, keep_ids, drop = FALSE]
    message(sprintf("  样本过滤：剔除 other %d 个 → 剩 %d", length(dropped), ncol(m)))
  } else {
    ## —— GSE48556：血液队列 GEO 字段完好，沿用 pheno$condition ——
    condcol <- intersect(c("condition","Condition","group","Group"), colnames(pd))[1]
    if (!is.na(condcol)) {
      keep_ids <- pd_ids[pd[[condcol]] %in% c("case","control")]
      stopif(length(keep_ids) == 0,
             sprintf("[%s] pheno$%s 中没有 case/control 取值（实际取值：%s）——请截图发回确认表型编码。",
                     gse, condcol, paste(unique(pd[[condcol]]), collapse = ", ")))
      dropped <- setdiff(colnames(m), keep_ids)
      m <- m[, intersect(colnames(m), keep_ids), drop = FALSE]
      message(sprintf("  样本过滤：剔除 condition=other %d 个 → 剩 %d（依据 pheno$%s）",
                      length(dropped), ncol(m), condcol))
    } else {
      message("  pheno 无 condition 列，跳过过滤（", gse, " 全保留 ", ncol(m), " 样本）")
    }
  }
  stopif(ncol(m) != cfg$expect_n,
         sprintf("[%s] 过滤后样本数 %d ≠ 预期 %d——请停下核对分组口径后重跑（截图发回）。", gse, ncol(m), cfg$expect_n))
  message(sprintf("  ✓ 样本数核验通过：%d（预期 %d）", ncol(m), cfg$expect_n))

  ## 2.2 尺度自检与 log2(x+1)→线性 反变换
  stopif(any(m < 0, na.rm = TRUE),
         sprintf("[%s] 检测到负值（%.3f）——疑似 z-score，禁止用于去卷积。请回查数据源。", gse, min(m, na.rm = TRUE)))
  q3 <- max(m, na.rm = TRUE)
  transform_applied <- "none (linear already)"
  if (q3 < 50) {
    stopif(!ALLOW_LOG2_INVERSE,
           sprintf("[%s] 疑似 log 尺度（max=%.2f）且 ALLOW_LOG2_INVERSE=FALSE，停止。", gse, q3))
    m <- 2^m - 1
    n_clip <- sum(m < 0); m[m < 0] <- 0
    transform_applied <- "2^x - 1 (log2(x+1) -> linear)"
    message(sprintf("  尺度自检：log2 尺度（max=%.2f）→ 反变换为线性（max=%.0f，截断负值 %d 个）",
                    q3, max(m, na.rm = TRUE), n_clip))
  } else {
    message(sprintf("  尺度自检：线性尺度确认（max=%.0f）", q3))
  }

  ## 2.3 完整性检查（CIBERSORTx 硬性要求）
  stopif(any(!is.finite(m)), sprintf("[%s] 存在 NA/Inf——CIBERSORTx 不接受缺失值，请先处理。", gse))

  ## 2.4 基因 symbol 去重（保留平均表达最高行；留档）
  gn <- rownames(m); n0 <- nrow(m)
  if (anyDuplicated(gn)) {
    keep <- tapply(seq_len(nrow(m)), gn, function(ix) ix[which.max(rowMeans(m[ix, , drop = FALSE]))])
    write.csv(data.frame(gse = gse, removed_gene = gn[-unlist(keep)]),
              file.path(OUT_DIR, sprintf("%s_removed_duplicates.csv", gse)), row.names = FALSE)
    m <- m[unlist(keep), , drop = FALSE]
  }
  message(sprintf("  基因数：%d → %d", n0, nrow(m)))

  ## 2.5 写出官方格式
  out_txt <- file.path(OUT_DIR, sprintf("%s_mixture.txt", gse))
  write.table(data.frame(GeneSymbol = rownames(m), m, check.names = FALSE),
              out_txt, sep = "\t", quote = FALSE, row.names = FALSE)
  manifest[[gse]] <- data.frame(
    gse = gse, platform = cfg$platform, n_genes = nrow(m), n_samples = ncol(m),
    transform = transform_applied,
    file = basename(out_txt), md5 = unname(md5sum(out_txt)),
    created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  message(sprintf("  ✓ 写出 %s（%d 基因 × %d 样本，md5=%s）",
                  basename(out_txt), nrow(m), ncol(m), substr(manifest[[gse]]$md5, 1, 8)))
}

## ---- 3. 总清单 -----------------------------------------------------------------
mf <- do.call(rbind, manifest)
write.csv(mf, file.path(OUT_DIR, "fix14v4_manifest.csv"), row.names = FALSE)
message("\n========== 全部完成 ==========")
message("产物目录：", OUT_DIR)
message("总样本数：", sum(mf$n_samples), "（应为 198 = 20+20+19+139）")
message("\n上传顺序：")
message("  1) Menu → Upload Files：官方 LM22.txt（类型 signature matrix）——你自己持有，本脚本不经手；")
message("  2) 同处上传 4 个 *_mixture.txt（类型 mixture）；")
message("  3) Impute Cell Fractions 逐队列建 job：batch correction 不勾 / 【取消勾选 Disable QN】/ absolute 不勾 / permutations=100；")
message("  4) 结果 CSV 放回 results/moduleF/fix14_cibersortx/official_results/ 待双版对照。")
