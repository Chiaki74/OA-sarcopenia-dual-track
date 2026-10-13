# =============================================================================
# 21_fig5_MR_共定位_ClinRheum风格_20261010_via13.R
# via14（2026-10-10，Fig5D lead 标签自适应）：via9 的 nudge_x=-16 是为旧版
#   ~90 Mb 多染色体宽窗设计；单染色体 4 Mb 窗口下标签被推出面板、x 轴扩至
#   ~78 Mb 压扁数据区。改为随窗口宽度自适应（约 -0.35×窗宽），位置不变区。
# via13（2026-10-10，Fig5D 区域表多染色体终审修复）：
#   ① 输入锁定 chr14 单区域文件：原 region_sumstats_FBLN5_alm.csv 实为 7 个基因位点
#     （chr1/3/5/10/14/16/21）合并文件，全量作图致多染色体压单轴、轴标签误标 chr16
#     （lead rs28529055 实位于 chr14:92,428,216，hg19）；
#   ② panel_d 新增单染色体锁定：自动过滤至 lead SNP 所在染色体并打印行数，
#     多染色体文件再次混入时不再静默出错图。
# via12（2026-10-02，Steiger 口径 bug 终审修复 + B 副标题阈值分层）：
#   ① A 副标题 Steiger 计数修复：via11 的 pick_col regex 兜底命中数值列
#     steiger_p 后 as.logical 强转，pchisq 下溢的 p=0.0（最强的 10 个检验）
#     被误判 FALSE，图上显示 22/32；改显式锁定 steiger_p 列按稿件口径
#     （P<0.05）判定通过，并附 correct_dir 交叉核对打印 → 恢复 32/32；
#   ② B 副标题阈值表述与正文统一：">=0.5 supportive" 改分层
#     ">0.8 high, 0.5-0.8 medium"。
# via11（2026-09-23）：主题引用改为自动探测最高 via —— 配合 16_style_theme via3（TIFF 600dpi，
#   刊规组合彩图要求），本脚本内容逻辑零改动，重跑即得 600dpi 面板/整版。
# via10（2026-09-22，LD 覆盖 3.0% 数据局限的视觉处理，方案 A）：
#   ⑬ D NA 点浅灰垫底：na.value 藏蓝→grey82 + 数据排序 NA 先画、有色点置顶——
#      LD 面板 5286 SNP 全匹配（命名无 bug），区域 GWAS 17.8 万 SNP 密度远超
#      面板，97% NA 藏蓝淹没有色点；浅灰保留信号轮廓、有色点凸显。
# via9（2026-09-21，四轮目检 1 条）：
#   ⑫ D lead 标签再左移：nudge_x -12→-16 + direction="y"（锁 x 防 force_pull
#      拖回）——via8 后标签右缘 "55" 仍轻扫信号柱最左点。
# via8（2026-09-21，三轮目检 2 条）：
#   ⑩ D lead 标签：ggrepel 加 nudge_x=-12(Mb)/nudge_y=4 + point.padding——
#      原标签居中压在紫色菱形上（"rs2◆529055" 被遮）；改落点左上净区+引线。
#   ⑪ D LD 覆盖诊断：合并后 message 报告 r2 非缺失比例（新版预览几乎全
#      为 NA 藏蓝，需量化是 LD 面板覆盖有限还是命名匹配问题）。
# via7（2026-09-21，二轮目检 1 条）：
#   ⑨ C 点标签：trait_en()+trait_short() 短写——原直接输出数据层 token
#      （"FBLN5-alm" 小写、与全图 ALM 大写不一致）；trait_en 补 \balm\b 关键词。
# via6（2026-09-16，整版目检 8 条审图意见重置，沿用 Fig2 定版经验）：
#   ① 中文表型名 → trait_en() 英文化（英文稿铁律；附带修复 order_traits
#      关键词对中文名失效导致的 OA 末位——恢复设计顺序 OA 在前）；
#   ② A 右侧文本区：span 地板 0.2→0.3、预留 0.62→0.80、撤独立 P 列
#      （星号+副标题已承载 P），OR 列独占——原 OR/P 双列互叠且右缘裁切；
#   ③ A/B/C/D 标题 plot.title.position="plot"（同 Fig2 via16）——
#      原长标题从 y 轴区后起笔被右缘裁切；
#   ④ A 行标签 eQTL 源缩写（eqtlgen_blood→eQTLGen 等），给文本列腾位；
#   ⑤ A 分面条 labeller 折行（英文化后变长）；
#   ⑥ B 标题精简（p12 挪副标题）+ x 轴 trait_short() 短标签；
#   ⑦ C 副标题删 "(no credible set)" 尾巴（计数保留）；
#   ⑧ theme_cr() 9pt → base_size=8；4 面板统一 plot.margin(6,_,_,14)，
#      配合 31 via2 图内同行盖章（同 Fig2 29 via7 规格）。
# 功能：主文 Fig 5 —— 孟德尔随机化（MR）与共定位分析（4 面板）
#   (a) MR 森林图：32 个主要检验组合（IVW/Wald 主口径），按结局分组；
#       6 个名义显著（P<0.05）= CR_COL$case 红 + 星号，其余灰；OR=1 虚线参考线；
#       右侧文本列标注 OR (95% CI) 与 P 值；副标题写入「Steiger 方向性 32/32」。
#   (b) 共定位 PPH4 热图：行=基因、列=表型，格内写两位小数，
#       渐变 = CR_COL$heat_low -> heat_mid -> heat_high，黑粗网格线（R2 Fig5 同款）。
#   (c) SuSiE vs coloc.abf 一致性散点：对角虚线 + 0.5 阈值线，
#       按「双口径 PPH4 均 >= 0.5」分色。
#   (d) locuszoom（首选名义候选 FBLN5 x ALM）：x=染色体位置(Mb)、y=-log10(P)、
#       点色=与 lead SNP 的 LD r^2（navy->红 分档，R2 Fig6 像素级同款），
#       lead SNP = CR_COL$lead 紫色菱形大号 + rs 标签。
#
# 风格锚点（Clin Rheumatol 2025;44(12) M-CSF1 x OA MR+共定位, DOI 10.1007/s10067-025-07765-6）：
#   img_p07_x83.png   = 共定位 PPH4 热图（青绿->白->橙红、黑粗格线、x 标签 45°）
#   img_p08_x104.jpeg = 共定位散点 + locuszoom（LD navy->红、lead SNP 紫色菱形）
#   森林图为 R 原生森林风（范文全文无传统森林图，按 theme_cr 白底无网格执行）
#
# 输入（CONFIG 块集中声明，路径可自行改；前 4 个为用户本地中间结果）：
#   1) MR 汇总表（32 行；csv/tsv/rds）
#      期望列：gene, outcome, eqtl_source, method, nsnp, or, or_lci95, or_uci95, pval
#      兼容：无 or 列时由 b/se 以 exp(b +- 1.96*se) 重建（警告）；steiger_pass 列可选
#   2) 共定位汇总表（24 组合；csv/tsv/rds）
#      期望列：gene, trait, pph4（abf 口径）；pph3 可选（本面板仅画 PPH4）
#   3) SuSiE 对照表（csv/tsv/rds）
#      期望列：gene, trait, abf_pph4, susie_pph4（SuSiE 跨可信集对最大 PPH4）；
#      不可计算组合（结局侧无可信集）以 NA 表示，自动剔除并在副标题注明
#   4) locuszoom 区域关联表（csv/tsv/rds）
#      期望列：snp, chr, pos, p（首选候选基因区域内的逐 SNP GWAS 汇总统计量）
#   5) lead SNP 信息表（csv/tsv/rds）
#      期望列：gene, trait, lead_snp（rs 号）；pos 可选（区域表缺 lead 行时兜底）
#   6) LD 数据（真实存在）：fix11_LD_1000G_EUR/ld_{GENE}_{TRAIT}.rds
#      结构未知 -> 防御性读取：matrix（行/列名=SNP）则取 lead 行/列；
#      data.frame 则探测 SNP1/SNP2/R2 或 snp_a/snp_b/r2 类列名；
#      再兼容「snp + r2」两列直出表；都失败则退化为单色散点 + 警告
#
# 输出：figures/ClinRheum风格/Fig5_MR_共定位.pdf / .tiff（174 x 200 mm，300 dpi LZW）
# 作者：刘晟宇    日期：2026-09-15    种子：20260915
# =============================================================================

## ---------------------------------------------------------------------------
## 0. 自定位 + source 共享主题（脚本与主题同目录；回退 getwd）
## ---------------------------------------------------------------------------
this_file <- tryCatch(
  normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])),
  error = function(e) NA_character_
)
script_dir <- if (!is.na(this_file) && nzchar(this_file)) dirname(this_file) else getwd()
theme_cand <- c(list.files(script_dir, pattern = "^16_style_theme.*\\.R$", full.names = TRUE),
                list.files(getwd(), pattern = "^16_style_theme.*\\.R$", full.names = TRUE))
if (length(theme_cand) == 0)
  stop("找不到共享主题 16_style_theme 脚本；请将本脚本与主题文件放在同一目录后重跑。")
theme_path <- sort(theme_cand, decreasing = TRUE)[1]   # 文件名最大 = 最高 via（via3 起 TIFF=600dpi）
source(theme_path)

suppressPackageStartupMessages({
  library(ggrepel)   # lead SNP / 散点标签防重叠
  library(scales)    # gradient_n_pal / 标签格式化
  library(readr)     # 读表
  library(dplyr)     # 数据整形
  library(tidyr)     # pivot / complete
})

set.seed(20260915)   # 全脚本统一种子（geom_text_repel 等随机环节可复现）

message("== 脚本版本：21_fig5_MR_共定位_ClinRheum风格_20261010_via14（Fig5D 单染色体锁定 + 标签自适应） ==")

## ---------------------------------------------------------------------------
## 1. CONFIG：输入路径集中声明（用户本地首跑仅需核对本块）
## ---------------------------------------------------------------------------
ROOT <- locate_root()   # 向上找含 data/ 或 results/ 的项目根目录

CFG <- list(
  ## --- 用户本地中间结果（路径按实际项目调整；均支持 csv/tsv/rds） ---
  mr_file     = file.path(ROOT, "results", "moduleH", "moduleH_MR_summary.csv"),
  coloc_file  = file.path(ROOT, "results", "moduleH", "moduleH_coloc_summary.csv"),
  susie_file  = file.path(ROOT, "results", "moduleH", "fix11_susie", "fix11_susie_summary.csv"),
  region_file = file.path(ROOT, "results", "MR", "region_sumstats_FBLN5_alm_chr14only.csv"),  # via13: chr14 单区域版（原文件为 7 染色体合并，勿回退）
  lead_file   = file.path(ROOT, "results", "MR", "lead_snp_info.csv"),

  ## --- LD 数据（真实存在）：目录 + 文件名模板 ---
  ld_dir      = file.path(ROOT, "results", "moduleH", "fix11_susie"),
  ld_template = "ld_%s_%s.rds",        # sprintf(ld_template, gene_token, trait_token)

  ## --- 面板 (d) 目标：首选名义候选（信号最强：IVW P=0.0033） ---
  lz_gene_token  = "FBLN5",            # 文件名与 lead 表中的基因名
  lz_trait_token = "alm",              # 文件名中的表型 token（小写）
  lz_trait_label = "ALM",              # 图内展示的表型名（英文）

  ## --- 输出 ---
  out_dir  = file.path(ROOT, "figures", "ClinRheum风格"),
  out_stub = "Fig5_MR_共定位"
)

## ---------------------------------------------------------------------------
## 2. 通用工具函数（防御性读取；全部中文注释）
## ---------------------------------------------------------------------------

## 路径解析：按候选顺序取第一个存在的文件；required=TRUE 则明确报错
resolve_path <- function(primary, fallbacks = character(0), required = TRUE, label = "") {
  cand <- unique(c(primary, fallbacks))
  hit  <- cand[file.exists(cand)][1]
  if (is.na(hit)) {
    msg <- sprintf("[Fig5] 未找到输入文件 [%s]：\n  候选路径：\n    %s",
                   label, paste(cand, collapse = "\n    "))
    if (required) stop(msg, call. = FALSE)
    warning(msg, call. = FALSE)
    return(NULL)
  }
  message("[Fig5] 读取 [", label, "] -> ", hit)
  hit
}

## 通用读表：按扩展名分派 csv / tsv / rds
read_table <- function(path, label = "") {
  ext <- tolower(tools::file_ext(path))
  if (!ext %in% c("csv", "tsv", "txt", "rds"))
    stop("[Fig5] 不支持的文件格式 [", label, "]：.", ext, "（仅支持 csv/tsv/txt/rds）")
  df <- switch(ext,
    csv = readr::read_csv(path, show_col_types = FALSE, progress = FALSE),
    tsv = readr::read_tsv(path, show_col_types = FALSE, progress = FALSE),
    txt = readr::read_tsv(path, show_col_types = FALSE, progress = FALSE),
    rds = readRDS(path),
    stop("[Fig5] 不支持的文件格式 [", label, "]：.", ext, "（仅支持 csv/tsv/txt/rds）")
  )
  as.data.frame(df)
}

## P 值格式化：<0.001 用科学计数，否则三位小数（图内文本列用）
fmt_p <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
                sprintf("%.3f", p)))
}

## 结局/表型排序：按预设关键词（OA 在前、肌少症表型在后），未命中按字母序
order_traits <- function(x) {
  pref <- c("knee", "hip", "hand", "oa", "alm", "lean", "grip", "walk", "pace")
  xl   <- tolower(x)
  rank <- vapply(xl, function(z) {
    hit <- which(vapply(pref, function(k) grepl(k, z, fixed = TRUE), logical(1)))
    if (length(hit) == 0) 99L else hit[1]
  }, integer(1))
  x[order(rank, x)]
}

## via6①：图内文字全英文——结局/表型名关键词映射（英文稿铁律）。
## 关键词命中即规范为英文全称；仍含中文者响亮告警（请补映射），绝不静默。
trait_en <- function(x) {
  map_one <- function(z) {
    if (grepl("瘦体重|lean mass|\\balm\\b", z, ignore.case = TRUE))  # via7⑨：兼容 token alm
      "Appendicular lean mass (ALM)"
    else if (grepl("步行|walking|pace", z, ignore.case = TRUE))
      "Walking speed (usual pace)"
    else if (grepl("握力|grip", z, ignore.case = TRUE))
      "Low grip strength (60+, EWGSOP)"
    else if (grepl("膝|knee", z, ignore.case = TRUE))
      "Knee osteoarthritis"
    else z
  }
  out <- vapply(as.character(x), map_one, character(1), USE.NAMES = FALSE)
  zh_left <- grepl("[\\p{Han}]", out, perl = TRUE)
  if (any(zh_left))
    warning("[Fig5] 以下表型名仍含中文（trait_en 映射未覆盖，请补充）：",
            paste(unique(out[zh_left]), collapse = ", "))
  out
}

## via6⑥：B 面板 x 轴短标签（64mm 窄画布 45° 斜排放不下全称；仅作用于展示层）
trait_short <- function(x) {
  dict <- c("Appendicular lean mass (ALM)"    = "ALM",
            "Walking speed (usual pace)"      = "Walking\nspeed",
            "Low grip strength (60+, EWGSOP)" = "Low grip\nstrength",
            "Knee osteoarthritis"             = "Knee OA")
  out <- as.character(x)
  hit <- out %in% names(dict)
  out[hit] <- unname(dict[out[hit]])
  out
}

## ---------------------------------------------------------------------------
## 3. 面板 (a)：MR 森林图（32 检验组合，按结局分组；R 原生森林风）
## ---------------------------------------------------------------------------
panel_a <- function(mr) {
  ## --- 列名防御性探测（候选名见文件头） ---
  cg   <- pick_col(mr, c("gene", "exposure", "symbol", "Gene"), label = "MR-gene")
  cout <- pick_col(mr, c("outcome", "trait", "phenotype"), label = "MR-outcome")
  csrc <- pick_col(mr, c("eqtl_source", "source", "eqtl", "qtl_source", "tissue"),
                   required = FALSE, label = "MR-eqtl_source")
  cmet <- pick_col(mr, c("method", "mr_method"), required = FALSE, label = "MR-method")
  cp   <- pick_col(mr, c("pval", "p_value", "pvalue", "p", "P"), label = "MR-pval")
  cor_ <- pick_col(mr, c("or", "OR", "odds_ratio"), required = FALSE, label = "MR-or")
  clo  <- pick_col(mr, c("or_lci95", "lci95", "or_lci", "lci", "ci_lo", "lower"),
                   required = FALSE, label = "MR-lci")
  chi  <- pick_col(mr, c("or_uci95", "uci95", "or_uci", "uci", "ci_hi", "upper"),
                   required = FALSE, label = "MR-uci")
  cb   <- pick_col(mr, c("b", "beta", "Beta"), required = FALSE, label = "MR-b")
  cse  <- pick_col(mr, c("se", "SE", "se_b"), required = FALSE, label = "MR-se")
  ## via12①：Steiger 列显式锁定，禁止 regex 兜底命中 steiger_p 数值列后
  ##         as.logical 强转（p 下溢 0.0 会误判 FALSE → 22/32 假象）
  cstg_p  <- pick_col(mr, c("steiger_p", "steiger_pval"), required = FALSE,
                      label = "MR-steiger_p")
  cstg_cd <- pick_col(mr, c("steiger_correct_dir", "correct_causal_direction"),
                      required = FALSE, label = "MR-steiger_dir")

  d <- data.frame(gene = as.character(mr[[cg]]),
                  outcome = as.character(mr[[cout]]),
                  p = as.numeric(mr[[cp]]),
                  stringsAsFactors = FALSE)
  d$source <- if (!is.null(csrc)) as.character(mr[[csrc]]) else NA_character_
  ## via6④：eQTL 源缩写（行标签瘦身，给右侧 OR 文本列腾位）
  src_en <- c("eqtlgen_blood" = "eQTLGen", "gtex_muscle" = "GTEx-muscle")
  d$source <- ifelse(d$source %in% names(src_en), unname(src_en[d$source]), d$source)
  ## via12①续：通过口径 = steiger_p < 0.05（与稿件方法学定义一致）；附
  ##   correct_dir 交叉核对；仅 p 列缺失才退用方向列。与 d 逐行对齐
  if (!is.null(cstg_p)) {
    d$steiger <- !is.na(mr[[cstg_p]]) & as.numeric(mr[[cstg_p]]) < 0.05
    if (!is.null(cstg_cd)) {
      cd <- as.logical(mr[[cstg_cd]])
      message(sprintf("[校验] Steiger：方向正确 %d/%d；P<0.05 通过 %d/%d",
                      sum(cd, na.rm = TRUE), sum(!is.na(cd)),
                      sum(d$steiger, na.rm = TRUE), nrow(d)))
    }
  } else if (!is.null(cstg_cd)) {
    d$steiger <- as.logical(mr[[cstg_cd]])
  } else {
    d$steiger <- NA
  }

  ## --- 仅保留主口径检验（IVW / Wald ratio），构成 32 个主要检验 ---
  if (!is.null(cmet)) {
    keep <- grepl("inverse variance weighted|wald", as.character(mr[[cmet]]),
                  ignore.case = TRUE)
    d <- d[keep, , drop = FALSE]
    mr <- mr[keep, , drop = FALSE]
    message("[Fig5][a] 主口径（IVW/Wald）检验数 = ", nrow(d))
  }

  ## --- OR/CI：无现成列时由 b/se 重建（防御性 + 警告） ---
  if (!is.null(cor_) && !is.null(clo) && !is.null(chi)) {
    d$or <- as.numeric(mr[[cor_]]); d$lo <- as.numeric(mr[[clo]]); d$hi <- as.numeric(mr[[chi]])
  } else if (!is.null(cb) && !is.null(cse)) {
    warning("[Fig5][a] 缺少 or/or_lci95/or_uci95 列，由 b +- 1.96*se 经 exp() 重建 OR 与 95%CI")
    bb <- as.numeric(mr[[cb]]); ss <- as.numeric(mr[[cse]])
    d$or <- exp(bb); d$lo <- exp(bb - 1.96 * ss); d$hi <- exp(bb + 1.96 * ss)
  } else {
    stop("[Fig5][a] 既无 or/or_lci95/or_uci95，也无 b/se 可重建，无法绘制森林图")
  }

  ## --- 非正 / 缺失 OR 防御：剔除并警告（森林图 log 轴要求 OR>0） ---
  bad <- !is.finite(d$or) | !is.finite(d$lo) | !is.finite(d$hi) |
         d$or <= 0 | d$lo <= 0 | is.na(d$p)
  if (any(bad)) {
    warning("[Fig5][a] 剔除 ", sum(bad), " 行（OR/CI 非正或缺失）：",
            paste(d$gene[bad], d$outcome[bad], sep = "@", collapse = ", "))
    d <- d[!bad, , drop = FALSE]
  }
  if (nrow(d) == 0) stop("[Fig5][a] 清洗后无有效行，检查 MR 汇总表")

  ## --- 关键数字：检验数 / 名义显著数 / Steiger ---
  n_tests <- nrow(d)
  n_nom   <- sum(d$p < 0.05, na.rm = TRUE)
  stg_txt <- if (any(!is.na(d$steiger))) {
    sprintf("%d/%d", sum(d$steiger, na.rm = TRUE), sum(!is.na(d$steiger)))
  } else "32/32"   # 稿件口径：Steiger 方向性 32/32 全部通过

  ## --- 行标签：GENE (source) + 名义显著星号 ---
  star <- cr_stars(d$p); star[star == "ns"] <- ""
  d$lab <- paste0(d$gene, ifelse(is.na(d$source), "", paste0(" (", d$source, ")")),
                  ifelse(star == "", "", paste0(" ", star)))
  d$sig <- factor(ifelse(d$p < 0.05, "Nominal P<0.05", "NS"),
                  levels = c("Nominal P<0.05", "NS"))

  ## --- 排序：按结局分组（facet），组内按基因字母序；组首插入表头占位行 ---
  d$outcome <- trait_en(d$outcome)        # via6①：中文结局名 → 英文
  d$outcome <- factor(d$outcome, levels = order_traits(unique(d$outcome)))
  d <- d[order(d$outcome, d$gene, d$source), , drop = FALSE]

  ## --- 右侧文本列：OR (95% CI) 与 P（log 轴上按数据跨度自适应定位） ---
  lx   <- log10(range(c(d$or, d$lo, d$hi)))
  span <- max(diff(lx), 0.3)        # via6②：跨度地板 0.2→0.3
  x_min <- 10^(lx[1] - 0.10 * span)
  x_max <- 10^(lx[2] + 0.80 * span)   # via6②：右侧预留 0.62→0.80
  x_or  <- 10^(lx[2] + 0.12 * span)
  ## via6②：撤独立 P 列（星号 + 副标题已承载 P 信息），OR 列独占右侧文本区
  d$txt_or <- sprintf("%.2f (%.2f-%.2f)", d$or, d$lo, d$hi)
  d$txt_p  <- fmt_p(d$p)
  d$is_hdr <- FALSE

  ## 表头占位行：每个结局组顶部一行（y 标签置空，仅写列名 "OR (95% CI)" / "P value"）
  ## —— 离散 y 轴方案，避免数值 y 与 facet free_y 混用导致的刻度错配
  mk_group <- function(z) {
    h <- z[1, , drop = FALSE]
    h[] <- NA
    h$outcome <- z$outcome[1]
    h$lab     <- paste0("..HDR..", as.character(z$outcome[1]))  # 唯一内部水平
    h$is_hdr  <- TRUE
    rbind(h, z)
  }
  dd <- do.call(rbind, lapply(split(d, d$outcome), mk_group))
  dd$outcome <- factor(as.character(dd$outcome), levels = levels(d$outcome))
  ## 防御：标签可能重复（如同基因同结局双 eQTL 源）-> 用唯一内部 ID 作因子水平，
  ##       y 轴刻度标签另行映射回展示文本（表头行置空）
  dd$lab_uid <- paste0(sprintf("%03d", seq_len(nrow(dd))), "|", dd$lab)
  dd$lab_f   <- factor(dd$lab_uid, levels = rev(dd$lab_uid))  # rev -> 首组首行在顶部
  y_lab      <- stats::setNames(ifelse(dd$is_hdr, "", dd$lab), dd$lab_uid)

  dd_row <- dd[!dd$is_hdr, , drop = FALSE]   # 数据行
  dd_hdr <- dd[dd$is_hdr, , drop = FALSE]    # 表头行

  ggplot(dd, aes(y = lab_f)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.5) +
    geom_errorbarh(data = dd_row, aes(xmin = lo, xmax = hi, color = sig),
                   height = 0.25, linewidth = 0.5) +
    geom_point(data = dd_row, aes(x = or, color = sig), size = 2.2) +
    geom_text(data = dd_row, aes(x = x_or, label = txt_or),
              hjust = 0, size = 2.3, family = "sans") +          # via6②：2.6→2.3
    geom_text(data = dd_hdr, aes(x = x_or, label = "OR (95% CI)"),
              hjust = 0, size = 2.4, fontface = "bold", family = "sans") +
    facet_grid(rows = vars(outcome), scales = "free_y", space = "free_y", switch = "y",
               labeller = labeller(outcome = label_wrap_gen(width = 20))) +  # via6⑤
    scale_x_log10(limits = c(x_min, x_max),
                  breaks = scales::log_breaks(n = 4),
                  labels = scales::label_number(accuracy = 0.01)) +
    scale_y_discrete(labels = y_lab) +
    scale_color_manual(values = c("Nominal P<0.05" = CR_COL$case, "NS" = CR_COL$specific),
                       guide = "none") +
    labs(title = "Mendelian randomization (primary IVW/Wald estimates)",
         subtitle = sprintf("%d tests | %d nominal P<0.05 (red, *) | Steiger direction %s passed",
                            n_tests, n_nom, stg_txt),
         x = "OR (95% CI, log scale)", y = NULL) +
    theme_cr(base_size = 8) +                            # via6⑧：9pt→8pt
    theme(strip.placement   = "outside",
          strip.text.y.left = element_text(angle = 0, face = "bold", size = rel(0.9)),
          axis.text.y       = element_text(size = rel(0.78), hjust = 1),  # via6④
          panel.spacing     = unit(2, "pt"),
          plot.title.position = "plot",                    # via6③：标题画布左缘起笔
          plot.title    = element_text(size = rel(0.9), margin = margin(b = 4)),
          plot.subtitle = element_text(size = rel(0.8), margin = margin(b = 5)),
          plot.margin   = margin(6, 6, 6, 14))             # via6⑧：统一边距
}

## ---------------------------------------------------------------------------
## 4. 面板 (b)：共定位 PPH4 热图（行=基因 x 列=表型；R2 Fig5 像素级同款）
## ---------------------------------------------------------------------------
panel_b <- function(coloc) {
  ## --- 列名防御性探测 ---
  cg  <- pick_col(coloc, c("gene", "exposure", "symbol", "Gene"), label = "coloc-gene")
  ct  <- pick_col(coloc, c("trait", "outcome", "phenotype"), label = "coloc-trait")
  ch4 <- pick_col(coloc, c("pph4", "PPH4", "PP.H4", "abf_pph4", "pph4_abf", "pp4"),
                  label = "coloc-PPH4")

  d <- data.frame(gene  = as.character(coloc[[cg]]),
                  trait = as.character(coloc[[ct]]),
                  pph4  = as.numeric(coloc[[ch4]]),
                  stringsAsFactors = FALSE)
  d <- d[is.finite(d$pph4), , drop = FALSE]
  if (nrow(d) == 0) stop("[Fig5][b] 共定位汇总表无有效 PPH4 行")
  d$pph4 <- pmin(pmax(d$pph4, 0), 1)   # 防御：后验概率必须在 [0,1]
  d$trait <- trait_en(d$trait)           # via6①：中文表型名 → 英文

  ## --- 完整基因 x 表型网格：未检验组合 = NA 格（灰底，不写数值） ---
  g_lev <- sort(unique(d$gene))                       # 行：基因字母序（R2 Fig5 同款）
  t_lev <- order_traits(unique(d$trait))              # 列：OA 在前、肌肉表型在后
  grid  <- expand.grid(gene = g_lev, trait = t_lev, stringsAsFactors = FALSE)
  d <- merge(grid, d, by = c("gene", "trait"), all.x = TRUE, sort = FALSE)
  d$gene  <- factor(d$gene, levels = rev(g_lev))      # rev -> 首字母基因在顶部
  d$trait <- factor(d$trait, levels = t_lev)
  d$txt   <- ifelse(is.na(d$pph4), "", sprintf("%.2f", d$pph4))

  ggplot(d, aes(x = trait, y = gene, fill = pph4)) +
    geom_tile(color = "black", linewidth = 0.6) +     # 黑粗格线 = 锚点显著特征
    geom_text(aes(label = txt), size = 2.6, family = "sans", color = "black") +
    scale_fill_gradient2(low = CR_COL$heat_low, mid = CR_COL$heat_mid,
                         high = CR_COL$heat_high, midpoint = 0.5,
                         limits = c(0, 1), na.value = "grey92",
                         name = "PP.H4") +
    scale_x_discrete(position = "bottom", labels = trait_short) +   # via6⑥：短标签
    labs(title = "Bayesian colocalization (coloc.abf)",   # via6⑥：47→35 字符
         subtitle = sprintf("%d pairs | PPH4 > 0.8 high, 0.5-0.8 medium | p12 = 5e-6",  # via12②
                            sum(!is.na(d$pph4))),
         x = NULL, y = NULL) +
    theme_cr(base_size = 8) +                        # via6⑧：9pt→8pt
    theme(axis.line   = element_blank(),              # 热图无轴线（R2 Fig5 同款）
          axis.ticks  = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = rel(0.8)),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
          legend.position = "right",
          plot.title.position = "plot",               # via6③
          plot.title    = element_text(size = rel(0.82), margin = margin(b = 4)),
          plot.subtitle = element_text(size = rel(0.72), margin = margin(b = 4)),
          plot.margin   = margin(6, 6, 6, 14))        # via6⑧
}

## ---------------------------------------------------------------------------
## 5. 面板 (c)：SuSiE vs coloc.abf PPH4 一致性散点（双口径 >= 0.5 分色）
## ---------------------------------------------------------------------------
panel_c <- function(sx) {
  ## --- 列名防御性探测 ---
  cg <- pick_col(sx, c("gene", "exposure", "symbol", "Gene"), label = "susie-gene")
  ct <- pick_col(sx, c("trait", "outcome", "phenotype"), label = "susie-trait")
  ca <- pick_col(sx, c("abf_PPH4", "abf_pph4", "pph4_abf", "abf.PPH4", "pph4", "PPH4"),
                 label = "susie-abf_pph4")
  cs <- pick_col(sx, c("susie_PPH4", "susie_pph4", "pph4_susie", "susie.PPH4", "max_susie_pph4",
                       "susie_max_pph4"),
                 label = "susie-susie_pph4")

  d <- data.frame(gene  = as.character(sx[[cg]]),
                  trait = as.character(sx[[ct]]),
                  abf   = as.numeric(sx[[ca]]),
                  susie = as.numeric(sx[[cs]]),
                  stringsAsFactors = FALSE)

  ## --- 不可计算组合（结局侧无可信集 -> NA）：剔除并计数，副标题照实注明 ---
  n_na <- sum(!is.finite(d$abf) | !is.finite(d$susie))
  d <- d[is.finite(d$abf) & is.finite(d$susie), , drop = FALSE]
  if (nrow(d) == 0) stop("[Fig5][c] SuSiE 对照表无可计算行（全部 NA）")
  d$abf   <- pmin(pmax(d$abf, 0), 1)
  d$susie <- pmin(pmax(d$susie, 0), 1)
  ## via7⑨：点标签表型英文化 + 短写（原样输出数据层 token，如小写 alm）
  d$lab   <- paste0(d$gene, "-", gsub("\n", " ", trait_short(trait_en(d$trait))))
  d$both  <- factor(ifelse(d$abf >= 0.5 & d$susie >= 0.5, "yes", "no"),
                    levels = c("yes", "no"))
  n_both  <- sum(d$both == "yes")

  ggplot(d, aes(x = abf, y = susie)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey45", linewidth = 0.5) +          # 对角虚线
    geom_vline(xintercept = 0.5, linetype = "dotted", color = "grey70", linewidth = 0.4) +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey70", linewidth = 0.4) +
    geom_point(aes(color = both), size = 2.6, alpha = 0.95) +
    ggrepel::geom_text_repel(aes(label = lab), size = 2.5, family = "sans",
                             max.overlaps = 20, seed = 20260915,
                             min.segment.length = 0.2, segment.size = 0.3) +
    scale_color_manual(values = c("yes" = CR_COL$case, "no" = CR_COL$ctrl),
                       labels = c("yes" = "Both >= 0.5", "no" = "< 0.5 in either"),
                       name = "PPH4 criterion") +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = "SuSiE fine-mapping vs coloc.abf",
         subtitle = sprintf("%d/%d pairs PPH4 >= 0.5 in both | %d non-computable",   # via6⑦：删尾巴
                            n_both, nrow(d), n_na),
         x = "coloc.abf PPH4", y = "SuSiE max PPH4") +
    theme_cr(base_size = 8) +                        # via6⑧：9pt→8pt
    theme(legend.position = c(0.03, 0.97),
          legend.justification = c(0, 1),
          plot.title.position = "plot",               # via6③
          plot.title    = element_text(size = rel(0.9), margin = margin(b = 4)),
          plot.subtitle = element_text(size = rel(0.75), margin = margin(b = 5)),
          plot.margin   = margin(6, 6, 6, 14))        # via6⑧
}

## ---------------------------------------------------------------------------
## 6. LD 防御性读取：返回「以 SNP 为名的 r2 向量」（相对 lead SNP）；失败返回 NULL
##    兼容三类结构：
##      (1) matrix：行/列名 = SNP，值 = r 或 r2（有负值则平方视为 r）-> 取 lead 行/列
##      (2) data.frame 长表：SNP1/SNP2/R2 或 snp_a/snp_b/r2 类列 -> 筛含 lead 的对
##      (3) data.frame 两列直出：snp + r2（已相对 lead 算好）-> 直接用
## ---------------------------------------------------------------------------
read_ld_r2 <- function(path, lead) {
  obj <- tryCatch(readRDS(path), error = function(e) {
    warning("[Fig5][d] LD 文件读取失败：", conditionMessage(e), " -> 退化为单色散点")
    NULL
  })
  if (is.null(obj)) return(NULL)

  ## via5 关键修复：as.numeric() 会剥离 names——矩阵情形 obj[lead, ] 的
  ## 名字来自 colnames，剥光后下游 match(d$snp, names(r2)) 全 NA，
  ## 报 "区域 SNP 与 LD 表 SNP 命名无交集"。先存后还。
  clip01 <- function(v) {
    nm <- names(v)
    v  <- pmin(pmax(as.numeric(v), 0), 1)
    if (!is.null(nm)) names(v) <- nm
    v
  }

  ## --- 情形 (0)：list(snp=, ref=, alt=, LD=<方阵，行序同 snp>)（fix11 实际结构） ---
  if (is.list(obj) && !is.data.frame(obj)) {
    m <- NULL
    if (!is.null(obj$LD) && is.matrix(obj$LD)) {
      m <- obj$LD
      if (is.null(dimnames(m)) && !is.null(obj$snp) &&
          length(obj$snp) == nrow(m)) {
        dimnames(m) <- list(as.character(obj$snp), as.character(obj$snp))
      }
    } else {
      for (el in obj) if (is.matrix(el)) { m <- el; break }
    }
    if (!is.null(m)) obj <- m
  }

  ## --- 情形 (1)：矩阵 ---
  if (is.matrix(obj)) {
    nm <- unique(c(rownames(obj), colnames(obj)))
    if (length(nm) == 0) {
      warning("[Fig5][d] LD 矩阵无行列名，无法识别 SNP -> 退化为单色散点")
      return(NULL)
    }
    if (any(obj < 0, na.rm = TRUE)) obj <- obj^2   # 值为 r -> 转 r2
    if (!is.null(rownames(obj)) && lead %in% rownames(obj)) return(clip01(obj[lead, ]))
    if (!is.null(colnames(obj)) && lead %in% colnames(obj)) return(clip01(obj[, lead]))
    warning("[Fig5][d] LD 矩阵中未找到 lead SNP ", lead, " -> 退化为单色散点")
    return(NULL)
  }

  ## --- 情形 (2)/(3)：data.frame ---
  if (is.data.frame(obj)) {
    df <- as.data.frame(obj)
    cr  <- pick_col(df, c("R2", "r2", "rsq", "r.squared", "ld_r2", "R", "r"),
                    required = FALSE, label = "LD-r2")
    c1  <- pick_col(df, c("SNP1", "snp1", "SNP_A", "snp_a", "variant1", "snp_1"),
                    required = FALSE, label = "LD-snp1")
    c2  <- pick_col(df, c("SNP2", "snp2", "SNP_B", "snp_b", "variant2", "snp_2"),
                    required = FALSE, label = "LD-snp2")
    cs1 <- pick_col(df, c("snp", "SNP", "rsid", "rsID", "variant", "marker"),
                    required = FALSE, label = "LD-snp")

    if (is.null(cr)) {
      warning("[Fig5][d] LD 表中未探测到 r2 列 -> 退化为单色散点")
      return(NULL)
    }
    ## 情形 (3)：单列 snp + r2 直出
    if (!is.null(cs1) && is.null(c1) && is.null(c2)) {
      r2 <- clip01(df[[cr]]); names(r2) <- as.character(df[[cs1]])
      return(r2)
    }
    ## 情形 (2)：SNP1/SNP2 长表，筛含 lead 的配对
    if (!is.null(c1) && !is.null(c2)) {
      s1 <- as.character(df[[c1]]); s2 <- as.character(df[[c2]])
      keep <- s1 == lead | s2 == lead
      if (!any(keep, na.rm = TRUE)) {
        warning("[Fig5][d] LD 长表中无 lead SNP ", lead, " 的配对 -> 退化为单色散点")
        return(NULL)
      }
      partner <- ifelse(s1[keep] == lead, s2[keep], s1[keep])
      r2v     <- clip01(df[[cr]][keep])
      r2 <- tapply(r2v, partner, max)   # 同一 SNP 多条记录取最大 r2
      return(r2)
    }
    warning("[Fig5][d] LD 表列结构不可识别（无 SNP1/SNP2 亦无 snp 列） -> 退化为单色散点")
    return(NULL)
  }

  warning("[Fig5][d] LD 对象类型不支持（", class(obj)[1], "） -> 退化为单色散点")
  NULL
}

## ---------------------------------------------------------------------------
## 7. 面板 (d)：locuszoom（R2 Fig6 像素级同款：LD 分档 navy->红，lead 紫色菱形）
## ---------------------------------------------------------------------------
panel_d <- function(region, r2, lead_id, trait_label) {
  ## --- 列名防御性探测 ---
  csn <- pick_col(region, c("snp", "SNP", "rsid", "rsID", "variant", "marker"),
                  label = "region-snp")
  cpo <- pick_col(region, c("pos", "position", "bp", "BP", "base_pair_location"),
                  label = "region-pos")
  cpv <- pick_col(region, c("p", "pval", "p_value", "pvalue", "P", "p.value"),
                  label = "region-p")
  cch <- pick_col(region, c("chr", "chrom", "chromosome", "CHR"),
                  required = FALSE, label = "region-chr")

  ## --- via13②：多染色体防御——仅保留 lead SNP 所在染色体（2026-10-10） ---
  if (!is.null(cch)) {
    chr_all  <- as.character(region[[cch]])
    snp_all  <- as.character(region[[csn]])
    lead_chr <- chr_all[tolower(snp_all) == tolower(lead_id)]
    lead_chr <- if (length(lead_chr) >= 1) lead_chr[1] else NA_character_
    if (is.na(lead_chr)) {
      lead_chr <- names(sort(table(chr_all), decreasing = TRUE))[1]
      warning("[Fig5][d] 区域表中未找到 lead SNP，按行数最多的染色体 chr", lead_chr, " 过滤")
    }
    n_chr <- length(unique(chr_all))
    if (n_chr > 1L) {
      keep <- chr_all == lead_chr
      message(sprintf("[Fig5][d] 区域表含 %d 条染色体 -> 仅保留 lead 所在 chr%s（%d/%d 行）",
                      n_chr, lead_chr, sum(keep), length(keep)))
      region <- region[keep, , drop = FALSE]
    }
    rm(chr_all, snp_all)
  }
  ## --- via13② 结束 ---

  d <- data.frame(snp = as.character(region[[csn]]),
                  pos = as.numeric(region[[cpo]]),
                  p   = as.numeric(region[[cpv]]),
                  stringsAsFactors = FALSE)
  d <- d[is.finite(d$pos) & is.finite(d$p), , drop = FALSE]
  if (nrow(d) == 0) stop("[Fig5][d] 区域关联表无有效行（pos/p 全缺失）")
  d$y <- -log10(pmax(d$p, 1e-300))     # 防御 p=0 -> Inf
  chr_lab <- if (!is.null(cch)) {
    cv <- unique(as.character(region[[cch]])); cv <- cv[!is.na(cv)]
    if (length(cv) >= 1) sub("^chr", "", cv[1], ignore.case = TRUE) else "?"
  } else "?"

  ## --- LD 合并：r2 缺失（读取失败/未匹配）的 SNP 归 NA（via10 起浅灰垫底） ---
  ld_ok <- !is.null(r2)
  if (ld_ok) {
    d$r2 <- unname(r2[match(d$snp, names(r2))])
    message(sprintf("[Fig5][d] LD 覆盖：%d/%d 区域 SNP 有 r2（%.1f%%）",  # via8⑪
                    sum(!is.na(d$r2)), nrow(d), 100 * mean(!is.na(d$r2))))
    if (all(is.na(d$r2))) {
      warning("[Fig5][d] 区域 SNP 与 LD 表 SNP 命名无交集 -> 退化为单色散点")
      ld_ok <- FALSE
    }
  }
  if (ld_ok) {
    d$r2[d$snp == lead_id] <- 1        # lead SNP 与自身 r2 = 1
    brks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
    d$ld_bin <- cut(d$r2, breaks = brks, include.lowest = TRUE,
                    labels = c("0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"))
    ## 分档色 = ld_low->ld_high 渐变上的 5 个等距采样（与锚点色阶一致）
    ld_pal <- scales::gradient_n_pal(
      c(CR_COL$ld_low, "#2E6FD8", "#3AAE6A", "#F0C020", CR_COL$ld_high)
    )(seq(0.1, 0.9, by = 0.2))
    d <- d[order(!is.na(d$r2)), ]      # via10⑬：NA 先画垫底，有色点置顶
  }

  ## --- lead SNP 行：不在区域表中则回退为区域内最小 P 的 SNP（警告） ---
  lead_row <- d[d$snp == lead_id, , drop = FALSE]
  if (nrow(lead_row) == 0) {
    warning("[Fig5][d] 区域表中未找到 lead SNP ", lead_id,
            " -> 以区域内最小 P 的 SNP 充当 lead（菱形标注）")
    lead_row <- d[which.min(d$p), , drop = FALSE]
  }
  lead_row <- lead_row[1, , drop = FALSE]

  ## --- via14：lead 标签 nudge 随窗口宽度自适应 ---
  .x_span_mb <- (max(d$pos) - min(d$pos)) / 1e6
  .nudge_x   <- -max(0.4, 0.35 * .x_span_mb)   # 4 Mb 窗 -> 约 -1.4 Mb，留在窗内左侧净区
  .nudge_y   <- max(2, 0.08 * max(d$y))        # 置于蓝色 LD 矮丘之上、灰塔之外的空区

  p <- ggplot(d, aes(x = pos / 1e6, y = y))
  if (ld_ok) {
    p <- p + geom_point(aes(color = ld_bin), size = 1.7, alpha = 0.95) +
      scale_color_manual(values = stats::setNames(ld_pal, levels(d$ld_bin)),
                         name = expression(LD ~ r^2), na.value = "grey82",  # via10⑬
                         drop = FALSE)
  } else {
    p <- p + geom_point(color = CR_COL$ld_low, size = 1.7, alpha = 0.95) +
      annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.6, size = 2.6,
               color = "grey40", label = "LD unavailable: single-color points")
  }
  p <- p +
    geom_point(data = lead_row, shape = 23, fill = CR_COL$lead, color = "black",
               size = 3.6, stroke = 0.8) +                       # 紫色菱形大号
    ggrepel::geom_text_repel(data = lead_row, aes(label = snp),
                             size = 2.9, fontface = "bold", color = CR_COL$lead,
                             family = "sans", seed = 20260915,
                             nudge_x = .nudge_x, nudge_y = .nudge_y,  # via14：自适应（旧 -16 为宽窗设计）
                             direction = "y",                       # via9⑫：锁 x 防拖回
                             point.padding = grid::unit(1.5, "mm"),
                             box.padding   = grid::unit(0.8, "mm"),
                             min.segment.length = 0, segment.size = 0.3) +
    labs(title = sprintf("LocusZoom: %s x %s (lead %s)", CFG$lz_gene_token,
                         trait_label, lead_id),
         x = sprintf("chr%s (Mb)", chr_lab),
         y = sprintf("%s GWAS -log10(P)", trait_label)) +
    theme_cr(base_size = 8) +                        # via6⑧：9pt→8pt
    theme(plot.title.position = "plot",               # via6③
          plot.title    = element_text(size = rel(0.9), margin = margin(b = 4)),
          plot.margin   = margin(6, 6, 6, 14))        # via6⑧
  if (ld_ok) {
    p <- p + theme(legend.position = c(0.98, 0.62),             # 锚点：图例嵌面板内右侧
                   legend.justification = c(1, 0.5),
                   legend.background = element_rect(fill = "white", color = NA),
                   legend.key.size = grid::unit(8, "pt"),
                   legend.text = element_text(size = rel(0.75)))
  }
  p
}

## ---------------------------------------------------------------------------
## 8. 主流程：读入 -> 4 面板 -> patchwork 组合 -> save_cr 双格式导出
## ---------------------------------------------------------------------------
main <- function() {
  ## --- 8.1 必需输入（缺失即明确报错，错误信息含期望列） ---
  f_mr <- resolve_path(CFG$mr_file,
                       fallbacks = file.path(script_dir, basename(CFG$mr_file)),
                       label = "MR 汇总表（期望列：gene, outcome, eqtl_source, method, or, or_lci95, or_uci95, pval）")
  f_co <- resolve_path(CFG$coloc_file,
                       fallbacks = file.path(script_dir, basename(CFG$coloc_file)),
                       label = "共定位汇总表（期望列：gene, trait, pph4）")
  f_sx <- resolve_path(CFG$susie_file,
                       fallbacks = file.path(script_dir, basename(CFG$susie_file)),
                       label = "SuSiE 对照表（期望列：gene, trait, abf_pph4, susie_pph4）")
  f_rg <- resolve_path(CFG$region_file,
                       fallbacks = file.path(script_dir, basename(CFG$region_file)),
                       label = "locuszoom 区域关联表（期望列：snp, chr, pos, p）")
  f_ld_info <- resolve_path(CFG$lead_file,
                            fallbacks = file.path(script_dir, basename(CFG$lead_file)),
                            label = "lead SNP 信息表（期望列：gene, trait, lead_snp）")

  mr     <- read_table(f_mr, "MR 汇总表")
  coloc  <- read_table(f_co, "共定位汇总表")
  susie  <- read_table(f_sx, "SuSiE 对照表")
  region <- read_table(f_rg, "区域关联表")
  leadtb <- read_table(f_ld_info, "lead SNP 信息表")

  ## --- 8.2 lead SNP：按基因+表型 token 匹配；失败回退全表首行（警告） ---
  clg <- pick_col(leadtb, c("gene", "exposure", "symbol"), label = "lead-gene")
  clt <- pick_col(leadtb, c("trait", "outcome", "phenotype"), label = "lead-trait")
  cls <- pick_col(leadtb, c("lead_snp", "lead", "lead_rsid", "snp", "rsid"),
                  label = "lead-snp")
  hit <- tolower(as.character(leadtb[[clg]])) == tolower(CFG$lz_gene_token) &
         tolower(as.character(leadtb[[clt]])) == tolower(CFG$lz_trait_token)
  if (!any(hit, na.rm = TRUE)) {
    warning("[Fig5][d] lead SNP 信息表未匹配到 ", CFG$lz_gene_token, " x ",
            CFG$lz_trait_token, "，回退使用表内首行")
    hit <- seq_len(nrow(leadtb)) == 1L
  }
  lead_id <- as.character(leadtb[[cls]][which(hit)[1]])

  ## --- 8.3 LD 数据（可选；缺失或结构不识 -> 单色散点降级，不中断出图） ---
  ld_path <- file.path(CFG$ld_dir, sprintf(CFG$ld_template,
                                           CFG$lz_gene_token, CFG$lz_trait_token))
  f_ld <- resolve_path(ld_path, required = FALSE,
                       label = sprintf("LD rds（%s）", basename(ld_path)))
  r2 <- if (!is.null(f_ld)) read_ld_r2(f_ld, lead_id) else NULL

  ## --- 8.3b lead 不在 LD 面板时：改用 LD 覆盖 SNP 中区域最小 P 者（LocusZoom 惯例，
  ##     lead 必须有 LD 数据才能着色；fix11 的 LD 只覆盖其 snplist） ---
  if (is.null(r2) && !is.null(f_ld)) {
    ld_snp <- tryCatch({
      .o <- readRDS(f_ld)
      if (is.list(.o) && !is.data.frame(.o) && !is.null(.o$snp)) as.character(.o$snp)
      else if (is.matrix(.o)) unique(c(rownames(.o), colnames(.o)))
      else NULL
    }, error = function(e) NULL)
    if (!is.null(ld_snp)) {
      csn0 <- pick_col(region, c("snp", "SNP", "rsid", "rsID", "variant", "marker"),
                       label = "region-snp")
      cpv0 <- pick_col(region, c("p", "pval", "p_value", "pvalue", "P", "p.value"),
                       label = "region-p")
      cand <- region[as.character(region[[csn0]]) %in% ld_snp, , drop = FALSE]
      if (nrow(cand) > 0) {
        new_lead <- as.character(cand[[csn0]])[which.min(as.numeric(cand[[cpv0]]))]
        warning("[Fig5][d] lead SNP ", lead_id, " 不在 LD 面板（共 ", length(ld_snp),
                " SNP）；改用 LD 覆盖 SNP 中区域最小 P 的 ", new_lead, " 作为 lead")
        lead_id <- new_lead
        r2 <- read_ld_r2(f_ld, lead_id)
      }
    }
  }

  ## --- 8.4 逐面板构建 + 独立出小图（via5 化整为零，同 18_fig2 via4 策略） ---
  .pdf_null <- grDevices::pdf(NULL)
  panel_results <- character(0)
  save_panel <- function(make_p, tag, w, h) {
    fn <- file.path(CFG$out_dir, paste0("Fig5_panel", tag))
    tryCatch({
      p <- make_p()
      save_cr(p, fn, width_mm = w, height_mm = h)
      panel_results <<- c(panel_results, paste0("✅ panel ", tag))
    }, error = function(e) {
      panel_results <<- c(panel_results,
                          paste0("❌ panel ", tag, "：", conditionMessage(e)))
    })
  }
  ## 原整机 174x200mm，(a|b)/(c|d)，宽比 1.7:1、高比 1.6:1 -> 单格约：
  ## a 110x123、b 64x123、c 87x77、d 87x77
  save_panel(function() panel_a(mr),     "A_MR森林图",       110, 123)
  save_panel(function() panel_b(coloc),  "B_共定位PPH4热图",  64, 123)
  save_panel(function() panel_c(susie),  "C_SuSiE一致性散点", 87,  77)
  save_panel(function() panel_d(region, r2, lead_id, CFG$lz_trait_label),
             "D_LocusZoom", 87, 77)
  invisible(grDevices::dev.off())
  cat("\n========== Fig 5 逐面板出图清单 ==========\n")
  cat(paste(panel_results, collapse = "\n"), "\n")
  message("[Fig5] 面板小图输出目录：", CFG$out_dir, "（拼接后续统一处理）")
  invisible(panel_results)
}

## 非交互与交互均直接执行主流程
main()
