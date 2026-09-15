###############################################################################
# ======================================================================
# fix11e（2026-08-29）：收敛审计版——在 fix11d 基础上仅做四处改动：
#   ①三处 susie_rss 调用 max_iter=100 → 1000（fix11d 警告 9 提示至少
#     一个拟合 IBSS 未在默认 100 迭代内收敛，可信集/PIP 严格来说只是
#     "提前停止"产物，必须升级为收敛确认）；
#   ②逐拟合记录 converged/niter 并写入汇总表（d1/d2_converged/niter 四列）；
#   ③每组合的 d1/d2 拟合对象存盘 fix11e_fits_<tag>.rds（供日后复核）；
#   ④汇总输出改写 fix11e_susie_summary_convergence.csv——
#     fix11d 的 fix11_susie_summary.csv 原样保留，两版对照。
# 判读纪律：若某拟合在 1000 迭代内收敛且 CS 数/PPH4/PPH3 与 fix11d 定性
#   一致 → 备忘录记"收敛审计通过"；若结果位移（尤其 FBLN5/HEMK1 出现
#   PPH4 抬头）→ 立即重估 §45 口径，两版差异不掩盖、必报告。
# 用法：与 fix11d 完全相同（FIX11_MODE <- "stage2"，source 整本运行），
#   依赖 results/moduleH/fix11_susie/ 下既有 ld_*.rds 与 fix9 缓存，幂等。
# ======================================================================
# 模块 H：10 个锁定标志物的 TwoSampleMR + 贝叶斯共定位正式版流水线
# 项目：骨关节炎 × 肌少症共享分子特征整合生物信息学研究（计划书 v1.3）
# 版本：v2.9（fix9）　日期：2026-08-27　编码：UTF-8
# fix9：修实跑 coloc 全部组合 check_dataset 全灭（"dataset : s must be
#       between 0 and 1"，日志 24 条同一句）——根因非数据、非 coloc 版本：
#       align_coloc_pair() 内 sub_ds() 的 `if (!is.null(ds$s)) out$s <- ds$s`
#       踩 R `$` 部分匹配陷阱——ds 无 "s" 字段但有 "snp" 时，ds$s 静默部分
#       匹配返回 SNP 名字符向量（warnPartialMatchDollar 默认 FALSE，无任何
#       警告），out$s 被污染为字符向量 → check_dataset 判 s 非 numeric。
#       eQTL 侧 ds（type="quant"，本无 s）每个组合必踩，与结局无关；
#       continuous 结局的 d2（无 s）同样被污染；binary 结局 d2 有真 s
#       幸免但同组合 d1 先死。修复：names() 精确判断 + [[ 精确取值；
#       check_dataset 失败时增打 d1/d2 字段指纹（names/s 的 class/len/值），
#       便于远程诊断。排查链存档：coloc 5.2.3（CRAN 当前版）与 main 分支
#       （GitHub，v6.0.1，2025-11-19）check_dataset 源码逐行核对——s 检查
#       逻辑两版一致（"s" %in% nd && (!is.numeric(d$s) || d$s<=0 || d$s>=1)），
#       用户 coloc 无需升降级；OUTCOMES 四结局 s=0.189/0.062/NULL/NULL 复核
#       合法；沙箱 28 组合数据层模拟全过但系"理想构造"、未经 sub_ds 真实
#       产物——教训：模拟必须复现完整构造链。MR 主线不经过 check_dataset，
#       不受影响；mr__ 全部 32 个缓存保留，重跑仅补 coloc 组合。
# fix3：修用户首跑崩溃（R 4.6.0，STEAP1 @ gtex_muscle 段，"Rle of type
#       'list' is not supported"）——①library(dplyr) 移到所有 Bioc 包之后
#       attach（加载顺序即掩码纪律）；②单 SNP 兜底等热点改 base R；③全脚本
#       裸 tidyverse 动词清零（无 dplyr:: 前缀者一律改写）；④standardize_eqtl
#       入口加 list 列守卫
# fix4：修用户实跑 stop（fetch_eqtlgen "MAF 表缺少 SNP/MAF 列"）——
#       ①MAF 表列名适配（新增 adapt_maf_table() 共享函数：有 MAF 列照旧，
#       否则用 AlleleB_all 按等位基因一致性条件推导 eaf）；②等位基因一致性
#       守卫（AssessedAllele∈{AlleleA,AlleleB} 且 OtherAllele 配对，不符剔除
#       并计数）；③z_to_beta_se 入参与 eaf.exposure 改用真实 eaf（语义修正，
#       harmonise 回文推断的核心收益）；④coloc 侧 ds$MAF=pmin(eaf,1-eaf)，
#       与结局侧同为 minor AF 口径；⑤实证数据来源：主文件头 50,000 行与 MAF
#       表 join 核验，AssessedAllele==AlleleB 仅 73.45%（2026-08-26 沙箱实测）
# 说明：基于 moduleH_TwoSampleMR_coloc_pipeline.R v1.0（2026-08-13）骨架全量实现；
#       v2.5 热修复：adapt_maf_table() 内 maf_tbl[, need] 变量子集在
#       data.table 下被当列名查找（"未找到列名 need"）——改为
#       as.data.frame(maf_tbl)[, need]（L465，修复合并前崩溃）
#       v2.6 热修复：断点缓存文件名 key 含 "|"（paste(sep="|")），为 Windows
#       非法文件名字符 → saveRDS gzfile "Invalid argument" 崩溃。新增
#       safe_key() 将 Windows 非法字符 [< > : " / \\ | ? *] 统一替换为 "-"，
#       仅用于 mr__/coloc__ 缓存文件名；状态留痕与日志仍用原 key（含 |）。
#       RDS 按内容读取（list.files 模式匹配 ^mr__/^coloc__ 不受影响），
#       且 Windows 下旧命名从未成功落盘，无陈旧缓存兼容问题
#       v2.7 热修复：eQTLGen 全量 cis 文件解压后 12.93GiB，fread 内存映射
#       在用户机器失败（连续虚拟内存不足）。改为读取沙箱预过滤的 10 锁定
#       基因子集（49,389 行 / 1.4MB gz，GeneSymbol 白名单精确匹配，与原
#       全量文件逐字节同源）。实证确认 NELL1/STEAP1/EGFR 在全量文件中
#       零记录（结构性缺失，2026-08-26 沙箱全量 13.9GB 单遍扫描核实）——
#       三基因的 eQTLGen 分支 MR/coloc 均为结构性阴性，GTEx 分支不受影响
#       v2.8 修复：coloc 结局区域远程 VCF 服务 gwas.mrcieu.ac.uk/files/ 已关停
#       （301→opengwas.io→404；gwasvcf 报 "vcf appears to be a string..." 即
#       URL 死链），用户侧 coloc 全灭。改为本地镜像主路径：沙箱从官方镜像
#       预提取 4 结局 × 7 个 eQTLGen 有记录锁定基因 ±2Mb 区域数据（用户落位
#       data/coloc_regions/，沙箱原件 /mnt/agents/output/moduleH_data/coloc_regions/；
#       NELL1/STEAP1/EGFR 全量零记录，无区域可提取，见 v2.7 结构性阴性）。
#       四文件 provenance（统一 schema：SNP,chr,pos,effect_allele,other_allele,
#       eaf,beta,se,p；tab 分隔 gzip）：
#         region_GCST90007526.tsv.gz  低握力  91,213 行  GWAS Catalog FTP
#                                       GRCh37 调和版，原生 rsID
#         region_GCST90000025.tsv.gz  ALM     177,822 行 GWAS Catalog FTP
#                                       GRCh37，原生 rsID
#         region_GCST007090.tsv.gz    膝OA    104,252 行 GWAS Catalog
#                                       Tachmazidou 原文件；chr:pos 经 eQTLGen
#                                       MAF 表映射 rsID（Allele1=effect，
#                                       Freq1=eaf，Effect=logOR）
#         region_ukbb924_neale.tsv.gz 步行速度 90,117 行 【替代源】Neale lab v3
#                                       phenotype 924 raw both_sexes（IEU VCF
#                                       关停）；beta 相对 alt 等位基因，eaf 按
#                                       alt==minor 条件推导，chr:pos 经 MAF 表
#                                       映射 rsID
#       每基因与 eQTL 侧交集已复核（eQTL 侧 MAF 适配 + 结局侧有效行
#       0<eaf<1 且 rsID 去重后）：3,213–8,760 SNP（低握力 3,213–8,260；
#       ALM 3,751–8,731；步行速度 3,539–8,667；膝OA 3,838–8,760；
#       均远超 coloc_min_snps=50）。等位基因协调保守规则后预计保留
#       3,002–8,457 SNP（沙箱按 align_coloc_pair() 规则复算）。
#       【替代源透明声明】步行速度 coloc 用 Neale v3 替代源（表型相同、处理
#       管线不同、coloc_n=358,974 ≠ ukb-b-4711 的 459,915）；MR 主线仍用
#       OpenGWAS API 的 ukb-b-4711，仅 coloc 用替代源，投稿需透明声明。
#       v2.1 按独立审查结论修复：B1 Steiger 静默失效（补 samplesize.exposure /
#       r.outcome）、B2 coloc 等位基因协调、M1–M4 及 minor 项（详见各修复点注释）；
#       v2.2 按复审结论修复：N1 build_outcome_dataset 携带结局侧等位基因列
#       （消除 B2 死代码）、N2 回文 SNP 单独分支杜绝误翻转（保守优先）、
#       N3 coloc cacheable 标记透传、N4 汇总增加 allele_harmonised 可机读列、
#       N5 陈旧 b37 注释更正
# 依据：计划书 v1.3 §4.8；备忘录 v1.1 §6.2–6.3
# 方法学对标：Yin et al. 2024 JCSM（单SNP工具兜底、Steiger P<0.05、PPH4>0.9 严标准）；
#             Wallace 2021 PLoS Genet（coloc-SuSiE 混合策略）；coloc v6 官方文档；
#             STROBE-MR（多重校正与敏感性分析报告口径）
#
# ----------------------------------------------------------------------------
# 【数据源核查】核查日期：2026-08-25（本模块定稿前已完成两项核查，结论如下）
#
# 核查 1｜结局 GWAS ID 有效性（OpenGWAS，api.opengwas.io / opengwas.io）
#   1) ebi-a-GCST90007526 — "Low hand grip strength (60 years and older)
#      (EWGSOP)"，n=256,523（48,596 病例 / 207,927 对照），EUR，
#      Jones et al. 2021 Nat Commun（PMID 33510174）。仍可查询，保留。
#      证据：platform.opentargets.org/study/GCST90007526；
#            多篇 MR 文献数据源表（如 Nutr Hosp 2025;42(6), PMID 表 I）。
#   2) ebi-a-GCST90000025 — "Appendicular lean mass"，n=450,243，EUR，
#      Pei et al. 2020（PMID 33097823）。仍可查询，保留，为肌少症 MR
#      文献主流 ALM 数据集。注：存在更新/更大 ALM 汇总集（如 GCST007841，
#      n≈458,000，来源待人工确认），如需可在 OUTCOMES 增列 local 分支。
#      证据：ebi.ac.uk/gwas/labs/studies/GCST90000025；Song et al. 2024
#            Nutr Metab Cardiovasc Dis。
#   3) ukb-b-4711 — "Usual walking pace"，n=459,915，EUR（UK Biobank
#      MRC-IEU 流水线）。仍可查询，保留。
#      证据：opengwas.io/datasets/ukb-b-4711；Liang et al. 2025 等。
#   4) ebi-a-GCST007090 — 膝骨关节炎，Tachmazidou et al. 2019 Nat Genet
#      （UK Biobank + arcOGEN，PMID 30664745），24,955 病例 / 378,169 对照
#      （n=403,124），EUR。仍可查询，保留，与计划书 §3.2.4 一致。
#      注：更大样本 OA 集（GO Consortium / Boer et al. 2021 Cell，全关节
#      n=826,690）不在 OpenGWAS，需本地下载后走 outcome_region_mode=
#      "local" 分支（见 OUTCOMES 注释）。
#      证据：ebi.ac.uk/gwas/labs/studies/GCST007090；BMC Med Genomics
#            2023 表 1；Arch Med Sci 2024 表 I。
#
# 核查 2｜eQTLGen 2019 全量 cis-eQTL 汇总统计可用性
#   结论：【提供全量下载】。除 FDR<0.05 显著配对文件（308 MB，
#   2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz）
#   外，官网目录同时提供全量 cis 汇总统计（4 GB，
#   2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz，
#   注意文件名中无 "0.05"，文献中称为 "Full cis-eQTL summary statistics"）。
#   证据：molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/ 目录清单；
#         eqtlgen.org/cis-eqtls.html；Nat Commun 2024（doi:10.1038/
#         s41467-024-49921-7）Data availability 明确使用全量文件做共定位。
#   落实：coloc 的 eQTL 侧一律使用全量文件（fetch_eqtlgen_full_region()）；
#        MR 工具变量筛选仍使用 FDR0.05 显著配对文件（fetch_eqtlgen()）。
#        全量文件缺失时 coloc 整列跳过并在 status/日志写清原因，
#        严禁拿显著配对文件冒充全区域数据做 coloc。
#
# ----------------------------------------------------------------------------
# 【本脚本输入】
#   - results/moduleC/moduleC_gene_tracks.csv（模块 C/D v2 轨归属表；
#     列：gene, track, direction_OA, direction_muscle, in_WGCNA_key；
#     10 个锁定基因缺一即 stop）
#   - eQTLGen FDR0.05 显著 cis-eQTL 文件 + SNP MAF 表（本地，见 EQTL_CONFIG）
#   - eQTLGen 全量 cis-eQTL 文件（本地，4 GB；coloc 必需，缺失则 coloc 跳过）
#   - GTEx v8 骨骼肌 signif_variant_gene_pairs + lookup 表（本地）
#   - OpenGWAS JWT（.Renviron 的 OPENGWAS_JWT）
#   - （可选）结局 GWAS 本地汇总文件（outcome_region_mode="local" 时）
# 【本脚本输出】（均落盘于 results/moduleH/，write.csv 一律 UTF-8）
#   - moduleH_MR_summary.csv        主方法（IVW/Wald）+ Bonferroni/BH 校正汇总
#   - moduleH_MR_all_methods.csv    全部方法长表（IVW/WM/Egger/Wald + 诊断量）
#   - moduleH_coloc_summary.csv     共定位汇总（分支/PPH4/双阈值/可信集对数/status）
#   - moduleH_known_causal_overlap.csv  双侧已知因果基因对照清单命中比对
#   - moduleH_status.csv            逐 key 执行状态留痕
#   - moduleH_MR_forest_primary.pdf 主效应森林图
#   - moduleH_log_<yyyymmdd>.txt    运行日志
#   - mr__<gene>|<src>|<ocid>.rds / coloc__<gene>|<src>|<ocid>.rds 断点续跑缓存
#
# 【判定口径】
#   因果成立 = Bonferroni 校正后显著（检验族 = 基因 × 结局 × eQTL 源）
#             + Steiger P<0.05（暴露→结局方向为真）
#             + 共定位 PPH4>0.8（宽标准；PPH4>0.9 严标准并列报告，
#               对齐 Yin et al. 2024）
# 【红线提醒】aging proxy / 人群混杂 / 方向倒置三条红线不涉及本模块
#   分析路径，但结果解读前请对照备忘录 v1.1 附录 D"结果对照"清单逐条勾核。
#
# 【使用前必读】
# 1. OpenGWAS 自 2024-05-01 起强制 JWT 认证：
#    a) 登录 https://api.opengwas.io/profile/ 生成 token
#    b) 在 .Renviron 写入：OPENGWAS_JWT=<your_token>（usethis::edit_r_environ()）
#    c) 重启 R 会话后运行 §0 验证（token 值严禁打印/落盘）
# 2. 所有中间结果落盘，断点可续跑（RDS 已存在则跳过重算）；随机种子固定
###############################################################################

# =============================================================================
# 0. 环境准备、依赖前置检查与认证验证
# =============================================================================

## 0.1 依赖前置检查（缺则 stop 并打印确切安装命令）
REQUIRED_PKGS <- c("TwoSampleMR", "ieugwasr", "MRPRESSO", "coloc", "susieR",
                   "gwasvcf", "data.table", "dplyr", "ggplot2",
                   "AnnotationDbi", "org.Hs.eg.db")
INSTALL_CMDS <- list(
  TwoSampleMR  = 'remotes::install_github("MRCIEU/TwoSampleMR")',
  ieugwasr     = 'remotes::install_github("MRCIEU/ieugwasr")',
  MRPRESSO     = 'install.packages("MRPRESSO")',
  coloc        = 'remotes::install_github("chr1swallace/coloc")  # >=5.2 可运行，推荐 v6（含 SuSiE 接口）',
  susieR       = 'install.packages("susieR")',
  gwasvcf      = 'remotes::install_github("mrcieu/gwasvcf")  # VariantAnnotation/SummarizedExperiment 随其自动安装',
  data.table   = 'install.packages("data.table")',
  dplyr        = 'install.packages("dplyr")',
  ggplot2      = 'install.packages("ggplot2")',
  AnnotationDbi = 'BiocManager::install("AnnotationDbi")',
  org.Hs.eg.db = 'BiocManager::install("org.Hs.eg.db")'
)
.missing <- REQUIRED_PKGS[!vapply(REQUIRED_PKGS, requireNamespace, logical(1), quietly = TRUE)]
if (length(.missing) > 0) {
  stop("模块 H 缺少依赖包：", paste(.missing, collapse = ", "),
       "\n请依次执行确切安装命令后重跑：\n  ",
       paste(unname(unlist(INSTALL_CMDS[.missing])), collapse = "\n  "))
}
## 【m1 口径说明】版本门槛 >=5.2：coloc.abf 与 runsusie/coloc.susie 接口自
## 5.2 起稳定可用，故 5.2 可运行；推荐 v6（SuSiE 集成与 check_dataset 校验改进）
.coloc_ver <- as.character(utils::packageVersion("coloc"))
if (utils::compareVersion(.coloc_ver, "5.2") < 0) {
  stop("coloc 版本过低（当前 ", .coloc_ver, "），需 >=5.2，推荐 v6：",
       'remotes::install_github("chr1swallace/coloc")')
}
message("依赖检查通过（coloc ", .coloc_ver, "）")

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
  library(MRPRESSO)
  library(coloc)
  library(susieR)
  library(gwasvcf)
  library(data.table)
  library(ggplot2)
  ## 【fix3 根治】加载顺序即掩码纪律：Bioconductor 栈（AnnotationDbi/
  ## org.Hs.eg.db → S4Vectors/IRanges）会先 attach，dplyr 必须在其【之后】
  ## attach——否则 Bioc 泛型抢占 dplyr 的 slice/filter/select 等 S4 分派，
  ## 裸管道走 S4 方法后报 "Rle of type 'list' is not supported"（用户首跑
  ## 崩溃现场，R 4.6.0，STEAP1 @ gtex_muscle 段）
  library(AnnotationDbi)   # Bioc 栈先 attach
  library(org.Hs.eg.db)
  library(dplyr)           # dplyr 最后 attach，确保其方法在搜索路径最前
})

set.seed(20260825)

## 0.2 日志函数（时间戳 + 同时写 results/moduleH/moduleH_log_<date>.txt）
## 说明：LOG_FILE 在 CONFIG 区 dir.create 之后初始化；此前仅打印到控制台
LOG_FILE <- NULL
log_msg <- function(...) {
  txt <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(...))
  cat(txt, "\n")
  if (!is.null(LOG_FILE)) cat(txt, "\n", file = LOG_FILE, append = TRUE)
  invisible(txt)
}

## 0.3 OpenGWAS JWT 验证（只经 .Renviron 的 OPENGWAS_JWT 读取；严禁打印 token 值）
opengwas_jwt <- Sys.getenv("OPENGWAS_JWT")
if (!nzchar(opengwas_jwt)) {
  stop("未检测到 OPENGWAS_JWT。请在 .Renviron 写入 OPENGWAS_JWT=<token> 后重启 R",
       "（usethis::edit_r_environ()；token 在 https://api.opengwas.io/profile/ 生成）。")
}
message("OpenGWAS 用户验证：")
print(tryCatch(ieugwasr::user(), error = function(e) {
  stop("JWT 无效或已过期，请到 https://api.opengwas.io/profile/ 重新生成。错误信息：",
       conditionMessage(e))
}))

## 0.4 OpenGWAS 连通性探针（最小查询；不通则给出网络排查提示后 stop）
probe_opengwas <- function() {
  log_msg("OpenGWAS 连通性探针：最小查询 gwasinfo('ieu-a-2') ...")
  ok <- tryCatch({
    info <- ieugwasr::gwasinfo("ieu-a-2")   # 最小公开数据集元数据查询
    !is.null(info) && nrow(info) > 0
  }, error = function(e) {
    log_msg("探针异常：", conditionMessage(e))
    FALSE
  })
  if (!ok) {
    stop(paste0(
      "OpenGWAS 连通性探针失败。网络排查提示：\n",
      "  1) 确认可访问 https://api.opengwas.io（公司/校园网代理需设置 ",
      "http_proxy/https_proxy 环境变量）；\n",
      "  2) 确认 OPENGWAS_JWT 未过期（profile 页重新生成后改 .Renviron 并重启 R）；\n",
      "  3) 若 429/503 为临时限流，稍后重跑（本脚本已内置指数退避重试）；\n",
      "  4) 仍不通则改用 outcome_region_mode='local' 并预下载结局汇总文件。"))
  }
  log_msg("OpenGWAS 连通性探针通过")
  invisible(TRUE)
}
probe_opengwas()

# =============================================================================
# CONFIG：分析配置区（锁定标志物已由模块 D v2 定稿，勿改）
# =============================================================================

CONFIG <- list(
  outdir = file.path("results", "moduleH"),      # 输出目录
  ## 10 个锁定标志物（模块 D v2 定稿；运行时从 moduleC_gene_tracks.csv 查轨归属与方向）
  genes = c("NELL1", "STEAP1", "GADD45A", "FBLN5", "HEMK1",
            "GDE1", "RNF14", "BNIP3", "EGFR", "CBR3"),
  gene_track_file = file.path("results", "moduleC", "moduleC_gene_tracks.csv"),
  cis_window_kb = 1000,                          # cis 区域：基因上下游 1 Mb（对齐 Yin et al.）
  iv_p_threshold = 5e-8,                         # 工具变量首选 P 阈值
  clump_r2 = 0.001,                              # LD clumping 首选（ieugwasr 默认严格档）
  clump_kb = 10000,
  min_ivs = 3,                                   # 少于此数触发分级放宽
  relax_r2_steps = c(0.01, 0.1),                 # 分级放宽 r² 台阶；最终兜底为单 SNP Wald ratio
  f_threshold = 10,                              # 弱工具变量阈值
  alpha_bonferroni = 0.05,                       # Bonferroni 主校正
  alpha_fdr = 0.05,                              # FDR(BH) 敏感性校正
  coloc_pph4_strict = 0.9,                       # 严标准（对齐 Yin et al.）
  coloc_pph4_relaxed = 0.8,                      # 宽标准（判定口径主阈值）
  coloc_prior_p12 = 5e-6,                        # coloc 先验 p12（官方推荐默认）
  coloc_min_snps = 50,                           # 区域内对齐后少于此数不做 coloc
  outcome_region_mode = "api",                   # 结局区域数据来源："api"（gwasvcf 查 OpenGWAS）
                                                 # 或 "local"（本地汇总文件，见 OUTCOMES$file）
  region_chunk_snp_cap = 10000,                  # 区域查询单次返回 SNP 上限（超限自动二分）
  region_block_bp = 1000000,                     # 区域查询初始分块长度（bp）
  coloc_region_dir = "data/coloc_regions",       # 【v2.8】coloc 结局区域本地镜像目录
                                                 # （远程 VCF 服务已关停后的主路径；见头部 v2.8 块）
  api_max_tries = 3,                             # API 调用指数退避重试上限（429/503/超时均重试）
  api_base_sleep = 2,                            # 退避基数（秒）：2, 4, 8
  ## 双侧已发表已知因果基因对照清单（备忘录 v1.1 §6.1b 附录 D 口径）
  ## 【M4 注】OA 侧清单保留 LGALS3：沿用 v1.0 骨架附录 D 口径（13 基因）；
  ## 备忘录重构版为 12 基因，二者差 1（LGALS3），待用户对照附录 D 原文核定。
  ## LGALS3 不在本模块 10 锁定基因内，对本流水线运行无影响。
  known_causal = list(
    sarcopenia = c("HP", "HLA-DRA", "MAP3K3", "MFGE8", "COL15A1", "AURKA"),  # Yin et al. 2024 JCSM
    oa = c("MAPK3", "CSK", "IFNGR2", "MVD", "FES", "ITGA2", "GZMK",
           "LGALS3", "ITIH1", "DNAJB12", "USP8", "IL12B", "RGMB")
  )
)

## 结局 GWAS（计划书 §3.2.4；有效性核查见头部"数据源核查"块，2026-08-25）
## s = 病例比例（仅 case-control 结局 coloc 需要）；n 供 coloc N 与报告；
## ncase/ncontrol（仅二分类）供 Steiger 检验用 get_r_from_lor() 预计算 r.outcome
## 【B1② 修复】samplesize.outcome 由主循环按 n 补齐（extract_outcome_data
## 返回缺列时兜底），二分类结局 r.outcome 在主循环显式预计算（见 §7 注释）
## 【v2.8】region_file = coloc 区域本地镜像文件名（位于 CONFIG$coloc_region_dir）；
## coloc_n = 替代源样本量（仅替代源结局设置，build_outcome_dataset 优先采用）
OUTCOMES <- list(
  low_grip  = list(id = "ebi-a-GCST90007526", label = "低握力(60+,EWGSOP)",
                   type = "binary",     n = 256523, s = 48596 / 256523,
                   ncase = 48596, ncontrol = 207927,
                   region_file = "region_GCST90007526.tsv.gz"),
  alm       = list(id = "ebi-a-GCST90000025", label = "四肢瘦体重(ALM)",
                   type = "continuous", n = 450243, s = NULL,
                   ncase = NULL, ncontrol = NULL,
                   region_file = "region_GCST90000025.tsv.gz"),
  walk_pace = list(id = "ukb-b-4711",         label = "步行速度(平常步速)",
                   type = "continuous", n = 459915, s = NULL,
                   ncase = NULL, ncontrol = NULL,
                   ## 【v2.8 替代源声明】coloc 区域数据用 Neale lab v3
                   ## phenotype 924 raw both_sexes 替代源（IEU VCF 服务关停）：
                   ## 表型相同、处理管线不同、n=358,974 ≠ ukb-b-4711 的 459,915。
                   ## MR 主线仍用 API 的 ukb-b-4711，仅 coloc 用替代源，
                   ## 投稿需透明声明（头部 v2.8 块有完整 provenance）
                   region_file = "region_ukbb924_neale.tsv.gz",
                   coloc_n = 358974),
  knee_oa   = list(id = "ebi-a-GCST007090",   label = "膝骨关节炎",
                   type = "binary",     n = 403124, s = 24955 / 403124,
                   ncase = 24955, ncontrol = 378169,
                   region_file = "region_GCST007090.tsv.gz")
  ## GO Consortium 全关节 OA（Boer et al. 2021 Cell，n=826,690）不在 OpenGWAS，
  ## 需另行下载汇总数据后用 local 分支读取，示例：
  ## go_oa = list(id = "local:GO_kneeOA", label = "OA(GO Consortium)", type = "binary",
  ##              n = 826690, s = 177517 / 826690, ncase = 177517, ncontrol = 649173,
  ##              file = "data/GWAS/GO_OA_sumstats.tsv")   # 列：chr,pos,SNP,beta,se,maf
)

## 暴露 eQTL 来源（计划书 §3.2.4；下载指引见 eQTL资源获取指引_20260820.md）
## - eQTLGen 血液 cis-eQTL（n=31,684，Võsa et al. 2021 Nat Genet）：MR + 共定位主线
## - GTEx v8 骨骼肌 cis-eQTL（n=706 EUR）：仅 MR（v8 起全配对文件因隐私不公开，
##   无法做区域共定位——共定位主线在 eQTLGen，与计划书 §4.8 一致）
EQTL_SOURCES <- c("eqtlgen_blood", "gtex_muscle")

## eQTL 本地文件配置
EQTL_CONFIG <- list(
  ## --- MR 工具变量筛选用：eQTLGen FDR<0.05 显著 cis-eQTL（约 308 MB）---
  ##   URL: https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/
  ##        2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz
  ##   列：Pvalue SNP SNPChr SNPPos AssessedAllele OtherAllele Zscore Gene
  ##       GeneSymbol GeneChr GenePos NrCohorts NrSamples FDR BonferroniP
  ##   ※ 无 beta/SE 列，需用 Zscore + MAF 换算（见 z_to_beta_se()）
  eqtlgen_file = "data/eQTL/eqtlgen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz",
  ## --- coloc 用：eQTLGen 全量 cis-eQTL 汇总统计（约 4 GB，核查 2 已确认存在）---
  ##   URL: https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/
  ##        2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz
  ##   （文件名无 "0.05"，即 eqtlgen.org 所称 "Full cis-eQTL summary statistics"）
  ##   ※ 含区域内全部受检 SNP（非仅显著配对），共定位必须使用本文件；
  ##     缺失则 coloc 整列跳过并留痕，严禁用 FDR0.05 文件冒充。
  ##   ※ fread 该文件建议 >=16 GB 内存；内存不足请先人工按染色体/基因拆分
  ##     （沙箱约定不使用 system2，拆分步骤在 R 外完成）。
  eqtlgen_full_file = "data/eQTL/eqtlgen/2019-12-11-cis-eQTLsFDR_10锁定基因子集.txt.gz",
                                                 # 【v2.7】沙箱预过滤 10 锁定基因子集（内存治本）；
                                                 # 原全量文件 12.93GiB 解压体积超出常规机器 fread 映射能力
  ## eQTLGen MAF 表（同站下载；全量与显著文件共用）：
  ## 【fix4 实证表头】SNP/hg19_chr/hg19_pos/AlleleA/AlleleB/allA_total/
  ##   allAB_total/allB_total/AlleleB_all——末列为 AlleleB 全队列频率
  ##   （非 minor AF）；列名适配与 eaf 推导见 adapt_maf_table()
  eqtlgen_maf_file = "data/eQTL/eqtlgen/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz",
  eqtlgen_n = 31684,
  ## GTEx v8 骨骼肌显著配对（自 GTEx tar 包中按 *Muscle_Skeletal* 提取）：
  ##   列：gene_id variant_id tss_distance ma_samples ma_count maf
  ##       pval_nominal slope slope_se 等；variant_id = chr_pos_ref_alt_b38
  ##   【m5 更正】GTEx v8 坐标为 GRCh38（非 hg19/b37）——与 eQTLGen/OpenGWAS
  ##   的 hg19 不同 build；经 rsID 对接可跨 build（同一变异），但坐标列（pos）
  ##   严禁与 hg19 侧混用，GTEx 分支不做 coloc 也规避了坐标混用风险
  gtex_muscle_file = "data/eQTL/GTEx_v8/Muscle_Skeletal.v8.signif_variant_gene_pairs.txt.gz",
  ## variant_id ↔ rsID 映射表（GTEx tar 包内提取；不设则 GTEx 分支无法做结局提取）
  gtex_lookup_file = "data/eQTL/GTEx_v8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz",
  gtex_muscle_n = 706
)

dir.create(CONFIG$outdir, showWarnings = FALSE, recursive = TRUE)
LOG_FILE <- file.path(CONFIG$outdir,
                      sprintf("moduleH_log_%s.txt", format(Sys.Date(), "%Y%m%d")))
log_msg("模块 H v2.8(fix8) 启动；输出目录：", normalizePath(CONFIG$outdir, mustWork = FALSE))

## 【M2 修复】输入文件统一自检（一次性打印精确缺失清单）：
## 必需文件（MR 主线：eQTLGen FDR0.05 + MAF 表 + GTEx signif + GTEx lookup）
## 缺失即 stop；eQTLGen 全量 4GB 文件缺失仅警告降级（coloc 列跳过，见
## fetch_eqtlgen_full_region() 与 M3 修复——缺失类跳出不缓存，文件到位后重跑即启用）
.required_files <- c("eQTLGen FDR0.05显著cis文件" = EQTL_CONFIG$eqtlgen_file,
                     "eQTLGen SNP MAF 表"        = EQTL_CONFIG$eqtlgen_maf_file,
                     "GTEx v8 骨骼肌signif配对"  = EQTL_CONFIG$gtex_muscle_file,
                     "GTEx v8 lookup映射表"      = EQTL_CONFIG$gtex_lookup_file)
.missing_req <- names(.required_files)[!file.exists(.required_files)]
if (length(.missing_req) > 0) {
  stop("模块 H 缺少必需输入文件（MR 主线无法运行）：\n  ",
       paste(sprintf("%s -> %s", .missing_req, unname(.required_files[.missing_req])),
             collapse = "\n  "),
       "\n下载指引见 eQTL资源获取指引_20260820.md；GTEx lookup 表缺失虽只影响 ",
       "GTEx 分支，但按必需口径处理以保证双源完整性。")
}
EQTLGEN_FULL_AVAILABLE <- file.exists(EQTL_CONFIG$eqtlgen_full_file)
if (!EQTLGEN_FULL_AVAILABLE) {
  log_msg("警告：eQTLGen 全量 cis 文件缺失（", EQTL_CONFIG$eqtlgen_full_file,
          "）——共定位列整体降级跳过（不缓存，文件到位后重跑即自动启用）；",
          "核查 2 已确认官网提供该文件（4 GB），严禁用 FDR0.05 文件冒充。")
} else {
  log_msg("输入文件自检通过：必需文件 4/4 在位，eQTLGen 全量文件在位（coloc 启用）")
}

## 【v2.8】coloc 结局区域镜像文件自检（缺则 WARN 不 stop——对应结局 coloc
## 回退 API 历史路径；VCF 服务已关停，预期失败，镜像落位后重跑即自动启用）
.region_miss <- character(0)
for (.ocn in names(OUTCOMES)) {
  .rf <- OUTCOMES[[.ocn]]$region_file
  if (!is.null(.rf) && !file.exists(file.path(CONFIG$coloc_region_dir, .rf)))
    .region_miss <- c(.region_miss, paste0(.ocn, " -> ", .rf))
}
if (length(.region_miss) > 0) {
  log_msg("警告：coloc 结局区域镜像缺失 ", length(.region_miss), " 个：\n  ",
          paste(.region_miss, collapse = "\n  "),
          "\n对应结局 coloc 回退 API 历史路径（VCF 服务已关停，预期失败）")
} else {
  log_msg("coloc 结局区域镜像自检通过：4/4 在位（", CONFIG$coloc_region_dir, "）")
}

# ---- 运行时轨归属查表（模块 C/D 输出；缺一即 stop）--------------------------
read_gene_tracks <- function() {
  f <- CONFIG$gene_track_file
  if (!file.exists(f)) {
    stop("缺少模块 C 轨归属表：", f,
         "。请先运行模块 C/D v2 生成该文件（列：gene, track, direction_OA,",
         " direction_muscle, in_WGCNA_key）。")
  }
  tracks <- utils::read.csv(f, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  need_cols <- c("gene", "track", "direction_OA", "direction_muscle", "in_WGCNA_key")
  missing_cols <- setdiff(need_cols, colnames(tracks))
  if (length(missing_cols) > 0) {
    stop("轨归属表缺列：", paste(missing_cols, collapse = ", "), "（文件：", f,
         "；实际列名：", paste(colnames(tracks), collapse = ", "), "）")
  }
  idx <- match(CONFIG$genes, tracks$gene)
  if (any(is.na(idx))) {
    stop("以下锁定基因在轨归属表中缺失，缺一即停：",
         paste(CONFIG$genes[is.na(idx)], collapse = ", "), "（文件：", f, "）")
  }
  tracks <- tracks[idx, , drop = FALSE]
  valid_tracks <- c("concordant", "mirror_OAup_muscleDown", "mirror_OAdown_muscleUp")
  bad <- setdiff(tracks$track, valid_tracks)
  if (length(bad) > 0) {
    stop("轨归属表存在非法 track 取值：", paste(bad, collapse = ", "),
         "（合法值：", paste(valid_tracks, collapse = "/"), "）")
  }
  rownames(tracks) <- NULL
  log_msg("轨归属查表完成：", nrow(tracks), " 个锁定基因全部命中；track 分布：",
          paste(sprintf("%s=%d", names(table(tracks$track)),
                        unname(as.integer(table(tracks$track)))), collapse = ", "))
  tracks
}
GENE_TRACKS <- read_gene_tracks()

## 取某基因的轨归属行（供主循环逐基因调用）
track_of <- function(gene_symbol) {
  r <- GENE_TRACKS[GENE_TRACKS$gene == gene_symbol, , drop = FALSE]
  if (nrow(r) != 1) stop("轨归属表中基因 ", gene_symbol, " 行数异常（", nrow(r), "）")
  r
}

# =============================================================================
# 通用工具：指数退避重试（429/503/超时均重试；失败留痕不中断全局）
# =============================================================================

## expr_fun: 无参函数；desc: 留痕描述；返回 NULL 表示最终失败（调用方自行跳过）
with_retry <- function(expr_fun, desc,
                       max_tries = CONFIG$api_max_tries,
                       base_sleep = CONFIG$api_base_sleep) {
  last_err <- ""
  for (i in seq_len(max_tries)) {
    res <- tryCatch(list(ok = TRUE, value = expr_fun()),
                    error = function(e) list(ok = FALSE, err = conditionMessage(e)))
    if (res$ok) {
      if (i > 1) log_msg(sprintf("  [重试成功] %s（第 %d 次尝试）", desc, i))
      return(res$value)
    }
    last_err <- res$err
    tag <- if (grepl("429", last_err)) "限流(429)"
           else if (grepl("503", last_err)) "服务不可用(503)"
           else if (grepl("timeout|timed out", last_err, ignore.case = TRUE)) "超时"
           else "错误"
    log_msg(sprintf("  [%s，重试 %d/%d] %s：%s", tag, i, max_tries, desc, last_err))
    if (i < max_tries) Sys.sleep(base_sleep * 2^(i - 1))   # 2, 4, 8 秒指数退避
  }
  log_msg(sprintf("  [放弃] %s 连续 %d 次失败（最后错误：%s），留痕后继续后续分析",
                  desc, max_tries, last_err))
  NULL
}

# =============================================================================
# 1. 工具变量筛选（分级策略：5e-8/r2=0.001 → 放宽 r2 → 单 SNP Wald ratio）
# =============================================================================

## ---- eQTL 读取辅助函数 ------------------------------------------------------

## Zscore → beta/SE 换算（eQTLGen 文件无 beta/SE 列；SMR/eQTLGen 官方 cookbook
## 推荐近似：se = 1/sqrt(2·p·(1-p)·(n+z²))，beta = z·se，p 用 MAF 近似，n = 样本量）
z_to_beta_se <- function(z, maf, n) {
  se <- 1 / sqrt(2 * maf * (1 - maf) * (n + z^2))
  list(beta = z * se, se = se)
}

## 【fix4】MAF 表列名适配 + 效应等位基因频率（eaf）条件推导
## （fetch_eqtlgen 与 fetch_eqtlgen_full_region 共用本函数，保证两路径口径一致）
## 实证（2026-08-26 沙箱核验，用户实跑触发）：MAF 文件真实表头为
##   SNP / hg19_chr / hg19_pos / AlleleA / AlleleB / allA_total / allAB_total /
##   allB_total / AlleleB_all
## 末列 AlleleB_all = AlleleB 的全队列频率（实测值域 0.15–0.87；>0.5 占 26.6%，
## 不是 minor AF）；主文件头 50,000 行 join 核验 AssessedAllele==AlleleB 仅
## 73.45%——即 26.55% 的 SNP 效应等位基因是 AlleleA。因此必须按等位基因
## 一致性条件推导 eaf，严禁把 AlleleB_all 直接当 MAF 或效应频率使用。
## na_policy: "floor" = eaf 缺失以 0.5 兜底（MR 路径沿用原精神）；
##            "drop"  = eaf 缺失剔除并计数（coloc 路径，不允许 0.5 兜底）
## 【v2.6】Windows 合法文件名：非法字符 [< > : " / \ | ? *] → "-"
safe_key <- function(k) gsub('[<>:"/\\\\|?*]', "-", k)

adapt_maf_table <- function(maf_tbl, dat, gene_symbol, context,
                            na_policy = c("floor", "drop")) {
  na_policy <- match.arg(na_policy)
  cols <- colnames(maf_tbl)
  ## 【m3 修复】merge 前断言 SNP 键唯一，防止重复键静默放大行数
  stopifnot(!anyDuplicated(maf_tbl$SNP))
  if ("MAF" %in% cols) {
    ## 旧表头：MAF 列即次等位基因频率，eaf 沿用近似口径（eaf≈MAF）
    dat <- merge(dat, maf_tbl[, c("SNP", "MAF")], by = "SNP", all.x = TRUE)
    dat$eaf <- dat$MAF
    log_msg(sprintf("  [%s] MAF 表为旧表头（MAF 列），eaf 沿用近似口径（%s）",
                    context, gene_symbol))
  } else if ("AlleleB_all" %in% cols) {
    need <- c("SNP", "AlleleA", "AlleleB", "AlleleB_all")
    miss <- setdiff(need, cols)
    if (length(miss) > 0)
      stop("MAF 表含 AlleleB_all 但缺 eaf 推导所需列：", paste(miss, collapse = ", "),
           "（实际列名：", paste(cols, collapse = ", "), "）——请按实际文件头调整 ",
           context)
    ## merge 时把 AlleleA/AlleleB 一并带上，供等位基因一致性推导
    dat <- merge(dat, as.data.frame(maf_tbl)[, need], by = "SNP", all.x = TRUE)
    ## 等位基因一致性守卫（统一 toupper 后比较）：
    ##   AssessedAllele 应 ∈ {AlleleA, AlleleB}；OtherAllele 应为集合中另一个
    aa <- toupper(dat$AssessedAllele); oa <- toupper(dat$OtherAllele)
    aA <- toupper(dat$AlleleA);       aB <- toupper(dat$AlleleB)
    ok_effect <- !is.na(aa) & (aa == aB | aa == aA)
    ok_other  <- !is.na(oa) & ((oa == aA & aa == aB) | (oa == aB & aa == aA))
    n_bad <- sum(!ok_effect | !ok_other, na.rm = TRUE)
    if (n_bad > 0)
      log_msg(sprintf("  [%s] 等位基因一致性守卫：%s 剔除 %d/%d 个 SNP",
                      context, gene_symbol, n_bad, nrow(dat)),
              "（AssessedAllele 不在 {AlleleA,AlleleB} 或 OtherAllele 不配对的位点）")
    ## 条件推导效应等位基因频率：AssessedAllele==AlleleB 用 AlleleB_all；
    ## ==AlleleA 用 1-AlleleB_all；不一致者为 NA（随后按 na_policy 处理）
    dat$eaf <- ifelse(aa == aB, dat$AlleleB_all,
               ifelse(aa == aA, 1 - dat$AlleleB_all, NA_real_))
    dat$eaf[!ok_effect | !ok_other] <- NA_real_
    log_msg(sprintf("  [%s] eaf 已按等位基因一致性条件推导（%s）：AssessedAllele==AlleleB 占 %.2f%%",
                    context, gene_symbol,
                    100 * sum(aa == aB, na.rm = TRUE) / max(1, sum(!is.na(aa)))))
  } else {
    stop("MAF 表既无 MAF 列也无 AlleleB_all 列（实际列名：",
         paste(cols, collapse = ", "), "）——请按实际文件头调整 ", context,
         "（fix4 适配口径见 adapt_maf_table() 注释）")
  }
  ## eaf 缺失处理
  n_na <- sum(is.na(dat$eaf))
  if (n_na > 0) {
    if (na_policy == "drop") {
      n_before <- nrow(dat)
      dat <- dat[!is.na(dat$eaf), ]
      log_msg(sprintf("  [%s] %s 剔除 eaf 缺失/不一致 SNP %d/%d 个（coloc 路径不允许兜底）",
                      context, gene_symbol, n_before - nrow(dat), n_before))
    } else {
      log_msg(sprintf("  警告：[%s] %s 部分 SNP 缺 eaf（%d/%d），以 0.5 兜底换算，解读时注意",
                      context, gene_symbol, n_na, nrow(dat)))
      dat$eaf[is.na(dat$eaf)] <- 0.5
    }
  }
  dat
}

## 基因坐标查询（org.Hs.eg.db 离线注释；按染色体分组取主染色体，避免跨行混合）
gene_coords <- function(gene_symbol) {
  loc <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_symbol,
                               keytype = "SYMBOL",
                               columns = c("CHR", "CHRLOC", "CHRLOCEND"))
  loc <- loc[!is.na(loc$CHR) & loc$CHR %in% c(1:22, "X"), ]
  if (nrow(loc) == 0) stop("基因坐标查询失败（org.Hs.eg.db 无记录）：", gene_symbol)
  main_chr <- names(which.max(table(loc$CHR)))
  loc <- loc[loc$CHR == main_chr, ]
  list(chr = as.character(unname(main_chr)),
       start = min(abs(loc$CHRLOC), na.rm = TRUE),
       end   = max(abs(loc$CHRLOCEND), na.rm = TRUE))
}

## 大文件缓存：eQTLGen 主文件 ~308MB（全量 ~4GB），按（基因×源）循环时避免重复 fread
.eqtl_cache <- new.env()
cached_fread <- function(path, ...) {
  args <- list(...)
  key <- paste0(path, "|", paste(capture.output(str(args)), collapse = "|"))
  if (!exists(key, envir = .eqtl_cache))
    assign(key, fread(path, ...), envir = .eqtl_cache)
  get(key, envir = .eqtl_cache)
}

## 统一输出口径（TwoSampleMR 兼容——harmonise_data 要求 *.exposure 后缀的
## 等位基因列与 exposure/id.exposure 列，否则 harmonise 必报错）：
##   SNP, chr, pos, beta.exposure, se.exposure, eaf.exposure, pval.exposure,
##   effect_allele.exposure, other_allele.exposure, exposure, id.exposure, gene
## 【B1① 修复】新增 samplesize.exposure（按源取：eQTLGen=31684 / GTEx=706）——
## directionality_test 需 samplesize.exposure + pval.exposure 计算 r.exposure，
## 缺列会 return(NULL) 且不抛错（Steiger 静默失效），故必须显式携带
standardize_eqtl <- function(df, source_tag, samplesize_n) {
  ## 【fix3 诊断加固】list 类型列守卫：Bioc 栈掩码或上游读取异常可能产生
  ## list 列（正是 S4 分派崩溃 "Rle of type 'list'" 的静默前兆）——告警列名；
  ## 单元素 list 列尝试 unlist 转原子向量，长度>1 的列表元素保留原样并告警
  bad <- names(df)[vapply(df, is.list, logical(1))]
  if (length(bad) > 0) {
    log_msg("  ★警告：standardize_eqtl(", source_tag, ") 入参含 list 列：",
            paste(bad, collapse = ", "))
    for (cl in bad) {
      lens <- lengths(df[[cl]])
      if (all(lens == 1)) {
        df[[cl]] <- unlist(df[[cl]], use.names = FALSE)   # 【fix4 minor】去名字向量
        log_msg("    列 ", cl, " 为单元素 list，已 unlist 转原子向量")
      } else {
        log_msg("    列 ", cl, " 含长度>1 的列表元素（", sum(lens > 1),
                " 个），保留原样——下游若报 S4 分派错误请优先排查此列")
      }
    }
  }
  df$exposure <- paste0(df$gene, " cis-eQTL (", source_tag, ")")
  df$id.exposure <- paste(source_tag, df$gene, sep = ":")
  df$samplesize.exposure <- samplesize_n
  df[, c("SNP", "chr", "pos", "beta.exposure", "se.exposure", "eaf.exposure",
         "pval.exposure", "effect_allele.exposure", "other_allele.exposure",
         "exposure", "id.exposure", "gene", "samplesize.exposure")]
}

## eQTLGen 读取（MR IV 筛选用，FDR0.05 显著配对文件）：GeneSymbol 过滤 + Zscore 换算
fetch_eqtlgen <- function(gene_symbol) {
  if (!file.exists(EQTL_CONFIG$eqtlgen_file))
    stop("缺少 eQTLGen FDR0.05 文件：", EQTL_CONFIG$eqtlgen_file,
         "（下载指引见 eQTL资源获取指引_20260820.md）")
  dat <- cached_fread(EQTL_CONFIG$eqtlgen_file)
  dat <- dat[dat$GeneSymbol == gene_symbol, ]
  if (nrow(dat) == 0) return(NULL)
  maf_tbl <- cached_fread(EQTL_CONFIG$eqtlgen_maf_file)
  ## 【fix4】MAF 表列名适配 + eaf 条件推导（共用 adapt_maf_table()；
  ## MR 路径 na_policy="floor"，eaf 缺失沿用 0.5 兜底精神）
  dat <- adapt_maf_table(maf_tbl, dat, gene_symbol,
                         context = "fetch_eqtlgen", na_policy = "floor")
  ## 【fix4 要点 3】z_to_beta_se 的 maf 入参改用 eaf——公式对 maf↔1-maf 对称、
  ## 数值不变，但语义必须正确（该参数本义即效应等位基因频率）；
  ## eaf.exposure 用真实 eaf 而非 MAF 近似——harmonise(action=2) 的回文推断
  ## 依赖它，这是本修复的核心收益；稳健性仍由 action=3（剔回文）敏感性验证。
  conv <- z_to_beta_se(dat$Zscore, dat$eaf, EQTL_CONFIG$eqtlgen_n)
  standardize_eqtl(data.frame(
    SNP = dat$SNP, chr = dat$SNPChr, pos = dat$SNPPos,
    beta.exposure = conv$beta, se.exposure = conv$se, eaf.exposure = dat$eaf,
    pval.exposure = dat$Pvalue,
    effect_allele.exposure = dat$AssessedAllele,
    other_allele.exposure = dat$OtherAllele,
    gene = gene_symbol, stringsAsFactors = FALSE),
    "eQTLGen-blood", samplesize_n = EQTL_CONFIG$eqtlgen_n)
}

## GTEx v8 骨骼肌读取：gene_id（Ensembl 带版本号）过滤 + slope/slope_se 直接取用
fetch_gtex <- function(gene_symbol) {
  if (!file.exists(EQTL_CONFIG$gtex_muscle_file))
    stop("缺少 GTEx 文件：", EQTL_CONFIG$gtex_muscle_file,
         "（下载指引见 eQTL资源获取指引_20260820.md）")
  ens <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_symbol,
                               keytype = "SYMBOL", columns = "ENSEMBL")$ENSEMBL
  if (length(ens) == 0 || all(is.na(ens))) stop("Ensembl ID 查询失败：", gene_symbol)
  dat <- cached_fread(EQTL_CONFIG$gtex_muscle_file)
  dat <- dat[sub("\\..*$", "", dat$gene_id) %in% ens, ]
  if (nrow(dat) == 0) return(NULL)
  vid <- do.call(rbind, strsplit(dat$variant_id, "_"))  # chr_pos_ref_alt_b38（GTEx v8 = GRCh38，见 EQTL_CONFIG m5 注）
  ## rsID 映射（关键）：ld_clump 与 extract_outcome_data 都按 rsID 查询——
  ## 没有 rsID，GTEx 分支的 MR 整体不可用（不止 clump 受限）。
  snp_key <- dat$variant_id
  if (!is.null(EQTL_CONFIG$gtex_lookup_file) && file.exists(EQTL_CONFIG$gtex_lookup_file)) {
    lk <- cached_fread(EQTL_CONFIG$gtex_lookup_file)
    hit <- lk$rs_id[match(dat$variant_id, lk$variant_id)]
    mapped <- !is.na(hit)
    snp_key[mapped] <- hit[mapped]
    log_msg(sprintf("  GTEx rsID 映射：%d/%d 个位点成功", sum(mapped), nrow(dat)))
  } else {
    log_msg("  警告：未找到 GTEx lookup 表（EQTL_CONFIG$gtex_lookup_file）——",
            "GTEx 分支无法做结局提取与 LD clumping，MR 主线将仅剩 eQTLGen")
  }
  standardize_eqtl(data.frame(
    SNP = snp_key,
    chr = sub("^chr", "", vid[, 1]), pos = as.integer(vid[, 2]),
    beta.exposure = dat$slope, se.exposure = dat$slope_se,
    ## 【m2 注】eaf.exposure 以 GTEx maf（次等位基因频率）近似效应(alt)等位基因
    ## 频率——与 eQTLGen 侧同口径的近似；harmonise(action=2) 回文推断依赖 eaf，
    ## 稳健性由 action=3（剔回文）敏感性验证（两源一致处理）
    eaf.exposure = dat$maf,
    pval.exposure = dat$pval_nominal,
    effect_allele.exposure = vid[, 4],       # GTEx 效应相对 alt 等位基因
    other_allele.exposure = vid[, 3],
    gene = gene_symbol, stringsAsFactors = FALSE),
    "GTEx-muscle", samplesize_n = EQTL_CONFIG$gtex_muscle_n)
}

## 提取某基因 cis 区域显著 eQTL 作为工具变量
## 说明：两源均为本地文件读取；eQTLGen 按 GeneSymbol 过滤（FDR<0.05 显著配对），
##       GTEx 按 Ensembl gene_id 过滤（signif_variant_gene_pairs 同样已限显著）。
fetch_cis_eqtl <- function(gene_symbol, source) {
  log_msg(sprintf("  提取 %s @ %s 的 cis-eQTL ...", gene_symbol, source))
  out <- switch(source,
    eqtlgen_blood = fetch_eqtlgen(gene_symbol),
    gtex_muscle   = fetch_gtex(gene_symbol),
    stop("未知 eQTL 来源：", source, "（合法值：", paste(EQTL_SOURCES, collapse = "/"), "）"))
  if (is.null(out) || nrow(out) == 0) {
    log_msg(sprintf("  %s 在 %s 中无显著 cis-eQTL 记录", gene_symbol, source))
    return(NULL)
  }
  ## cis 窗口过滤（CONFIG$cis_window_kb 生效点；两源文件虽已限 cis 显著配对，
  ## 仍按基因坐标 ± 窗口显式裁剪，保证口径可审计）
  gc <- gene_coords(gene_symbol)
  w <- CONFIG$cis_window_kb * 1000
  out <- out[out$chr == gc$chr & out$pos >= gc$start - w & out$pos <= gc$end + w, ]
  if (nrow(out) == 0) {
    log_msg(sprintf("  %s 在 ±%d kb cis 窗口内无记录", gene_symbol, CONFIG$cis_window_kb))
    return(NULL)
  }
  out
}

## 分级 clumping：返回该基因最终采用的工具变量集与实际使用参数（留痕）
## ld_clump 走 with_retry（429/503/超时重试 ≤3 次，失败留痕不中断）
## 【m6 修复】通过 .iv_diag 环境区分 NULL 返回的原因：
##   "真无 IV"（无显著 SNP / clump 成功但档内不足且兜底也不满足）
##   vs "clump 服务不可用"（所有已尝试档位 ld_clump API 均失败）——
##   后者不缓存结论，服务恢复后重跑即可，主循环据此写 status 留痕
.iv_diag <- new.env()
.iv_diag$status <- ""
select_instruments <- function(eqtl_df, gene_symbol, source) {
  if (is.null(eqtl_df) || nrow(eqtl_df) == 0) {
    .iv_diag$status <- "真无IV：该源无显著 cis-eQTL 记录"
    return(NULL)
  }
  any_attempt <- FALSE      # 是否存在 P<阈值的 SNP 可送 clump
  api_fail_all <- TRUE      # 已尝试档位是否全部 API 失败
  for (r2 in c(CONFIG$clump_r2, CONFIG$relax_r2_steps)) {
    ## 【fix3】base R 子集（原 dplyr::filter 裸管道，Bioc 掩码风险点清零）
    sig <- eqtl_df[!is.na(eqtl_df$pval.exposure) &
                   eqtl_df$pval.exposure < CONFIG$iv_p_threshold, , drop = FALSE]
    if (nrow(sig) == 0) next
    any_attempt <- TRUE
    clumped <- with_retry(
      function() ieugwasr::ld_clump(
        dplyr::tibble(rsid = sig$SNP, pval = sig$pval.exposure),
        clump_r2 = r2, clump_kb = CONFIG$clump_kb
      ),
      desc = sprintf("ld_clump %s @ %s (r2=%.3f)", gene_symbol, source, r2))
    if (is.null(clumped)) next            # 该档 API 失败，试下一档
    api_fail_all <- FALSE
    if (nrow(clumped) >= CONFIG$min_ivs) {
      iv <- sig[sig$SNP %in% clumped$rsid, , drop = FALSE]   # 【fix3】base R（原裸 filter 管道）
      attr(iv, "params") <- sprintf("P<%.0e, r2=%.3f, kb=%d", CONFIG$iv_p_threshold, r2, CONFIG$clump_kb)
      attr(iv, "tier")   <- ifelse(r2 == CONFIG$clump_r2, "首选参数", "放宽参数")
      .iv_diag$status <- "ok"
      return(iv)
    }
  }
  ## 兜底：单 SNP（cis 区域最显著位点），Wald ratio —— Yin et al. 2024 主要即此策略
  ## 【fix3】改 base R 实现（order 取最小 P 行），彻底摆脱 dplyr/S4 分派风险——
  ## 用户首跑崩溃即发生在此路径（Bioc 泛型抢占 slice 后 list 列触发
  ## "Rle of type 'list' is not supported"）
  top <- eqtl_df[!is.na(eqtl_df$pval.exposure) &
                 eqtl_df$pval.exposure < CONFIG$iv_p_threshold, , drop = FALSE]
  if (nrow(top) > 0) top <- top[order(top$pval.exposure)[1], , drop = FALSE]
  if (nrow(top) > 0) {
    attr(top, "params") <- "单SNP(top cis eQTL), Wald ratio"
    attr(top, "tier")   <- "单SNP兜底"
    log_msg(sprintf("  %s @ %s 工具变量不足 %d 个，启用单 SNP 兜底（%s）%s",
                    gene_symbol, source, CONFIG$min_ivs, top$SNP[1],
                    ifelse(any_attempt && api_fail_all,
                           "【注意：clump 服务不可用，兜底由 API 失败触发，服务恢复后建议重跑】",
                           "")))
    .iv_diag$status <- if (any_attempt && api_fail_all)
      "clump服务不可用：单SNP兜底为降级结果" else "ok"
    return(top)
  }
  .iv_diag$status <- if (any_attempt && api_fail_all) {
    paste0("clump 服务不可用（ld_clump 三档 API 重试均失败）且无可用兜底 SNP",
           "——非'真无 IV'结论，服务恢复后重跑")
  } else {
    sprintf("真无IV：无 P<%.0e 的 cis-eQTL", CONFIG$iv_p_threshold)
  }
  log_msg(sprintf("  %s @ %s %s", gene_symbol, source, .iv_diag$status))
  NULL
}

## 弱工具变量检验：F = (beta/se)^2，逐 SNP 计算并报告最小值/均值
check_f_stat <- function(iv_df) {
  f <- (iv_df$beta.exposure / iv_df$se.exposure)^2
  list(min = min(f), mean = mean(f), pass = all(f > CONFIG$f_threshold))
}

# =============================================================================
# 2. 结局数据提取与 harmonise（action=2 主分析；action=3 剔回文敏感性）
# =============================================================================

## 结局提取：指数退避重试 ≤3 次 + 失败留痕不中断全局
get_outcome <- function(iv_df, outcome_id) {
  with_retry(
    function() TwoSampleMR::extract_outcome_data(snps = iv_df$SNP, outcomes = outcome_id),
    desc = sprintf("结局提取 %s", outcome_id))
}

## harmonise：action=2（回文尽量按 eaf 推断，主分析）；
##             action=3（剔除所有回文 SNP，敏感性——针对 eaf 为 MAF 近似
##             可能致回文位点方向误判的稳健性验证）
harmonise_pair <- function(iv_df, out_df, action = 2) {
  dat <- tryCatch(TwoSampleMR::harmonise_data(iv_df, out_df, action = action),
                  error = function(e) {
                    log_msg("  harmonise(action=", action, ") 失败：", conditionMessage(e))
                    NULL
                  })
  if (is.null(dat) || !"mr_keep" %in% colnames(dat)) return(NULL)
  dat <- dat[dat$mr_keep, ]   # 剔除 harmonise 判定为模糊/不可用的 SNP
  if (nrow(dat) == 0) return(NULL)
  dat
}

# =============================================================================
# 3. MR 主分析 + 敏感性分析三件套（Steiger / MR-PRESSO / 异质性+多效性）
# =============================================================================

run_mr_suite <- function(dat, gene_symbol, outcome_label) {
  nsnp <- length(unique(dat$SNP))
  res <- list(gene = gene_symbol, outcome = outcome_label, nsnp = nsnp)

  ## 3.1 主分析：>=3 SNP 用 IVW 全家桶；2 SNP 仅 IVW；1 SNP 仅 Wald ratio
  methods <- if (nsnp >= 3) c("mr_ivw", "mr_egger_regression", "mr_weighted_median")
             else if (nsnp == 2) "mr_ivw"
             else "mr_wald_ratio"
  res$mr <- tryCatch(TwoSampleMR::mr(dat, method_list = methods),
                     error = function(e) {
                       log_msg("  mr() 失败（", gene_symbol, " × ", outcome_label,
                               "）：", conditionMessage(e))
                       NULL
                     })

  ## 3.2 异质性（Cochran's Q）与水平多效性（MR-Egger 截距）
  res$heterogeneity <- tryCatch(TwoSampleMR::mr_heterogeneity(dat), error = function(e) NULL)
  res$pleiotropy    <- tryCatch(TwoSampleMR::mr_pleiotropy_test(dat), error = function(e) NULL)

  ## 3.3 Steiger 方向性检验（Steiger P<0.05 判定暴露→结局方向成立，Yin et al. 标准）
  ## 【B1③ 修复】directionality_test 在 samplesize.exposure / r.outcome 等列缺失时
  ## return(NULL) 且不抛错——显式检测并告警，杜绝 steiger_p 静默全 NA
  res$steiger <- tryCatch({
    st <- TwoSampleMR::directionality_test(dat)
    if (is.null(st)) {
      log_msg("  ★警告：directionality_test 返回 NULL（", gene_symbol, " × ",
              outcome_label, "）——请检查 samplesize.exposure / samplesize.outcome /",
              " r.outcome 是否随 harmonise 保留；该 key 的 Steiger P 记 NA 并留痕")
    }
    st
  }, error = function(e) {
    log_msg("  ★警告：directionality_test 抛错（", gene_symbol, " × ", outcome_label,
            "）：", conditionMessage(e))
    NULL
  })

  ## 3.4 MR-PRESSO 全局检验与离群值剔除（需 >=4 SNP，否则不可做，留 NA 并注明）
  res$mrpresso <- tryCatch({
    if (nsnp >= 4) {
      MRPRESSO::mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                          SdOutcome = "se.outcome", SdExposure = "se.exposure",
                          OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                          data = as.data.frame(dat), NbDistribution = 1000, SignifThreshold = 0.05)
    } else { log_msg("  SNP<4，MR-PRESSO 不适用（留 NA）"); NA }
  }, error = function(e) { log_msg("  MR-PRESSO 失败：", conditionMessage(e)); NULL })

  ## 3.5 留一法（敏感性）
  res$loo <- tryCatch(TwoSampleMR::mr_leaveoneout(dat), error = function(e) NULL)
  res
}

# =============================================================================
# 4. 主效应提取 + 多重检验校正（检验族 = 基因 × 结局 × eQTL 来源）
# =============================================================================

## 从 run_mr_suite 的 list 提取 IVW/WM/Egger/Wald 主效应行 + nsnp + F 统计 +
## Steiger P + MR-PRESSO 全局 P + 异质性 Q P + Egger 截距 P 到统一 data.frame
## （每方法一行；is_primary 标记该 key 的主方法：nsnp>=2 取 IVW，否则 Wald ratio）
extract_main_effects <- function(res) {
  if (is.null(res) || is.null(res$mr) || nrow(res$mr) == 0) return(NULL)
  mr <- res$mr
  primary_method <- if (res$nsnp >= 2) "Inverse variance weighted" else "Wald ratio"

  ## 异质性 Q P（IVW 与 Egger 分列）
  het_ivw <- if (!is.null(res$heterogeneity))
    res$heterogeneity[res$heterogeneity$method == "Inverse variance weighted", ] else NULL
  het_eg <- if (!is.null(res$heterogeneity))
    res$heterogeneity[res$heterogeneity$method == "MR Egger", ] else NULL
  q_p_ivw <- if (!is.null(het_ivw) && nrow(het_ivw) > 0) unname(het_ivw$Q_pval[1]) else NA_real_
  q_p_eg  <- if (!is.null(het_eg) && nrow(het_eg) > 0) unname(het_eg$Q_pval[1]) else NA_real_

  ## Steiger 方向性
  steiger_p <- NA_real_; steiger_dir <- NA
  if (!is.null(res$steiger) && nrow(res$steiger) > 0 &&
      "steiger_pval" %in% colnames(res$steiger)) {
    steiger_p  <- unname(res$steiger$steiger_pval[1])
    steiger_dir <- unname(res$steiger$correct_causal_direction[1])
  }

  ## MR-PRESSO 全局 P（结构：$`MR-PRESSO results`$`Global Test`$Pvalue）
  presso_p <- NA_real_
  if (!is.null(res$mrpresso) && is.list(res$mrpresso)) {
    presso_p <- tryCatch(
      unname(res$mrpresso[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]),
      error = function(e) NA_real_)
  }

  ## Egger 截距 P
  egger_int_p <- if (!is.null(res$pleiotropy) && nrow(res$pleiotropy) > 0)
    unname(res$pleiotropy$pval[1]) else NA_real_

  ## action=3 敏感性（剔回文）结果（若主循环已跑）
  sens <- res$sens_action3
  sens_b <- if (!is.null(sens)) unname(sens$b) else NA_real_
  sens_p <- if (!is.null(sens)) unname(sens$pval) else NA_real_
  sens_nsnp <- if (!is.null(sens)) unname(sens$nsnp) else NA_integer_

  data.frame(
    gene = res$gene, track = res$track,
    direction_OA = res$direction_OA, direction_muscle = res$direction_muscle,
    eqtl_source = res$eqtl_source,
    outcome = res$outcome, outcome_id = res$outcome_id,
    method = mr$method, nsnp = res$nsnp,
    b = mr$b, se = mr$se, pval = mr$pval,
    or = exp(mr$b),
    or_lci = exp(mr$b - 1.96 * mr$se), or_uci = exp(mr$b + 1.96 * mr$se),
    is_primary = mr$method == primary_method,
    f_min = res$fstat$min, f_mean = res$fstat$mean, f_pass = res$fstat$pass,
    iv_tier = res$iv_tier, iv_params = res$iv_params,
    steiger_p = steiger_p, steiger_correct_dir = steiger_dir,
    mrpresso_global_p = presso_p,
    het_q_pval_ivw = q_p_ivw, het_q_pval_egger = q_p_eg,
    egger_intercept_p = egger_int_p,
    harmonise_action = res$harmonise_action,
    sens_action3_b = sens_b, sens_action3_p = sens_p, sens_action3_nsnp = sens_nsnp,
    stringsAsFactors = FALSE)
}

## 多重校正（STROBE-MR 口径：同时报告原始 P、校正方法、校正后阈值）
## 检验族 = 基因 × 结局 × eQTL 源 → 恰好等于主方法行集合（每 key 一行）；
## Bonferroni 为主（保守），BH(FDR) 为敏感性，两者并列报告
apply_multiple_correction <- function(summary_df) {
  primary <- summary_df[summary_df$is_primary %in% TRUE, , drop = FALSE]
  n_tests <- nrow(primary)
  primary$n_tests_family <- n_tests
  primary$p_bonferroni  <- pmin(primary$pval * n_tests, 1)
  primary$p_fdr         <- p.adjust(primary$pval, method = "BH")
  primary$sig_bonferroni <- primary$p_bonferroni < CONFIG$alpha_bonferroni
  primary$sig_fdr        <- primary$p_fdr < CONFIG$alpha_fdr
  log_msg(sprintf("多重校正：检验族 = 基因×结局×eQTL源，共 %d 个主方法检验；",
                  n_tests),
          sprintf("Bonferroni 显著 %d 个，BH-FDR 显著 %d 个",
                  sum(primary$sig_bonferroni, na.rm = TRUE),
                  sum(primary$sig_fdr, na.rm = TRUE)))
  primary
}

# =============================================================================
# 5. 贝叶斯共定位（coloc-SuSiE 混合策略，PPH4 双阈值判定）
# =============================================================================
# 方法依据：Wallace 2021 PLoS Genet —— SuSiE 能识别可信集时优先 coloc-SuSiE
# （处理多因果变异），否则回退单因果变异 coloc.abf；coloc v6 已内置接口。
# 数据要求：共定位需区域内【全部】SNP 汇总统计（非仅显著 SNP）——
#   eQTL 侧：eQTLGen 全量 cis 文件（核查 2 已确认存在；缺失则整列跳过）；
#   结局侧：v2.8 起优先读本地镜像区域文件（CONFIG$coloc_region_dir，远程 VCF
#           服务已关停）；镜像缺失才回退 gwasvcf API 历史路径或 OUTCOMES$file
#           本地全量汇总分支（outcome_region_mode="local"）。
#   GTEx v8 无全配对文件（隐私），gtex_muscle 源不做 coloc。

## ---- OpenGWAS VCF → 区域汇总统计 data.frame --------------------------------
## OpenGWAS VCF 基因型字段：ES(beta), SE, LP(-log10 P), AF(效应等位基因频率), SS 等
vcf_to_region_df <- function(vcf) {
  if (is.null(vcf) || nrow(vcf) == 0) return(NULL)
  rr <- SummarizedExperiment::rowRanges(vcf)
  g <- VariantAnnotation::geno(vcf)
  get_geno <- function(nm) {
    if (nm %in% names(g)) as.numeric(g[[nm]][, 1]) else rep(NA_real_, length(rr))
  }
  lp <- get_geno("LP")
  pval <- 10^(-lp)
  alt <- VariantAnnotation::alt(vcf)
  alt_chr <- if (length(alt) > 0) as.character(unlist(alt)) else rep(NA_character_, length(rr))
  data.frame(
    SNP = names(rr),
    chr = sub("^chr", "", as.character(GenomicRanges::seqnames(rr))),
    pos = GenomicRanges::start(rr),
    effect_allele = alt_chr,
    other_allele = as.character(VariantAnnotation::ref(vcf)),
    beta = get_geno("ES"), se = get_geno("SE"), eaf = get_geno("AF"),
    pval = pval, stringsAsFactors = FALSE)
}

## ---- 结局侧区域数据：API 分支（gwasvcf 区域查询，分块 ≤10k SNP/次，
##      限速指数退避重试 ≤3 次，429/503/超时均重试）---------------------------
fetch_outcome_region_api <- function(oc_id, chr, start, end) {
  ## OpenGWAS 公开 VCF 下载地址（无需 API 鉴权；若远程 tabix 因网络/鉴权失败，
  ## with_retry 重试后返回 NULL，status 留痕并提示改用 local 模式）
  url <- sprintf("https://gwas.mrcieu.ac.uk/files/%s/%s.vcf.gz", oc_id, oc_id)
  fetch_block <- function(b_start, b_end, depth = 0L) {
    chrompos <- sprintf("%s:%d-%d", chr, b_start, b_end)
    vcf <- with_retry(
      function() gwasvcf::query_gwas(vcf = url, chrompos = chrompos),
      desc = sprintf("区域查询 %s %s", oc_id, chrompos))
    if (is.null(vcf)) return(NULL)
    df <- vcf_to_region_df(vcf)
    ## 单块 SNP 数触顶 → 二分递归（最深 4 层），保证单次 ≤ region_chunk_snp_cap
    if (!is.null(df) && nrow(df) >= CONFIG$region_chunk_snp_cap &&
        depth < 4L && (b_end - b_start) > 1) {
      mid <- floor((b_start + b_end) / 2)
      log_msg(sprintf("  块 %s SNP 数 %d 触顶 %d，二分拆分",
                      chrompos, nrow(df), CONFIG$region_chunk_snp_cap))
      return(rbind(fetch_block(b_start, mid, depth + 1L),
                   fetch_block(mid + 1, b_end, depth + 1L)))
    }
    df
  }
  bounds <- seq(start, end, by = CONFIG$region_block_bp)
  pieces <- list()
  for (bs in bounds) {
    be <- min(bs + CONFIG$region_block_bp - 1, end)
    pieces[[length(pieces) + 1]] <- fetch_block(bs, be)
    Sys.sleep(CONFIG$api_base_sleep)   # 限速：块间 sleep，避免触发 429
  }
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  if (length(pieces) == 0) return(NULL)
  region <- do.call(rbind, pieces)
  region <- region[!duplicated(region$SNP), , drop = FALSE]
  log_msg(sprintf("  结局区域数据 %s：%s:%d-%d，共 %d SNP",
                  oc_id, chr, start, end, nrow(region)))
  region
}

## ---- 结局侧区域数据：本地全量汇总文件分支（CONFIG$outcome_region_mode = "local"）----
## 【v2.8 改名】原 fetch_outcome_region_local → fetch_outcome_region_sumstats
## （OUTCOMES$file 旧分支，供 GO Consortium 等全量汇总文件使用；函数名让位给
## v2.8 的区域镜像分支），函数体未动
## 本地文件列要求：chr, pos, SNP, beta, se, maf（hg19/GRCh37 坐标）
fetch_outcome_region_sumstats <- function(oc, chr, start, end) {
  if (is.null(oc$file)) {
    log_msg("  [local分支] 结局 ", oc$id, " 未配置本地文件（OUTCOMES$file），跳过")
    return(NULL)
  }
  if (!file.exists(oc$file)) {
    log_msg("  [local分支] 结局本地文件缺失：", oc$file, "，跳过")
    return(NULL)
  }
  region <- cached_fread(oc$file)
  region <- region[region$chr == chr & region$pos >= start & region$pos <= end, ]
  if (nrow(region) == 0) return(NULL)
  data.frame(SNP = region$SNP, chr = region$chr, pos = region$pos,
             effect_allele = NA_character_, other_allele = NA_character_,
             beta = region$beta, se = region$se, eaf = region$maf,
             pval = NA_real_, stringsAsFactors = FALSE)
}

## ---- 【v2.8】结局侧区域数据：本地镜像分支（远程 VCF 关停后的主路径）--------
## 读取沙箱预提取的 region 文件（用户落位 CONFIG$coloc_region_dir），
## 统一 schema（tab 分隔，gzip，九列）：SNP, chr, pos, effect_allele,
## other_allele, eaf, beta, se, p——与 build_outcome_dataset 期望列名一致
fetch_outcome_region_local <- function(oc, gc) {
  if (is.null(oc$region_file)) return(NULL)
  f <- file.path(CONFIG$coloc_region_dir, oc$region_file)
  if (!file.exists(f)) return(NULL)   # 缺失由调度层回退 API（此处不重复告警）
  region <- cached_fread(f)
  region <- as.data.frame(region)     # data.table 变量子集禁用（v2.8 自查口径）
  need <- c("SNP", "chr", "pos", "effect_allele", "other_allele",
            "eaf", "beta", "se", "p")
  miss <- setdiff(need, colnames(region))
  if (length(miss) > 0) {
    stop("coloc 区域文件 ", f, " 缺列：", paste(miss, collapse = ", "),
         "（实际列名：", paste(colnames(region), collapse = ", "),
         "）——期望九列 schema：", paste(need, collapse = ", "))
  }
  region$effect_allele <- toupper(region$effect_allele)
  region$other_allele  <- toupper(region$other_allele)
  w <- CONFIG$cis_window_kb * 1000
  ## chr 列为数值型、gc$chr 为字符 → as.character 对齐比较（避免类型不一致静默全 FALSE）
  region <- region[as.character(region$chr) == as.character(gc$chr) &
                   region$pos >= gc$start - w & region$pos <= gc$end + w, , drop = FALSE]
  if (nrow(region) == 0) return(NULL)
  log_msg(sprintf("  结局区域数据（本地镜像）%s @ %s：chr%s:%d-%d，共 %d SNP",
                  oc$id, basename(f), gc$chr, gc$start - w, gc$end + w, nrow(region)))
  region
}

## 结局区域数据统一入口（v2.8 调度：本地镜像存在 → 优先 local；
## 缺失 → 回退 API 历史路径；outcome_region_mode="local" → OUTCOMES$file 旧分支）
fetch_outcome_region <- function(oc, gc) {
  w <- CONFIG$cis_window_kb * 1000
  if (!is.null(oc$region_file) &&
      file.exists(file.path(CONFIG$coloc_region_dir, oc$region_file))) {
    log_msg("  [coloc] 结局区域走本地镜像：", oc$region_file)
    return(fetch_outcome_region_local(oc, gc))
  }
  if (identical(CONFIG$outcome_region_mode, "local")) {
    fetch_outcome_region_sumstats(oc, gc$chr, gc$start - w, gc$end + w)
  } else {
    ## 【v2.8】远程 VCF 服务 gwas.mrcieu.ac.uk/files/ 已关停（301→opengwas.io
    ## →404），本路径预期失败，代码留作历史；失败时 run_coloc_for_pair 记
    ## cacheable=FALSE，镜像文件落位后重跑即可自动启用
    fetch_outcome_region_api(oc$id, gc$chr, gc$start - w, gc$end + w)
  }
}

## ---- 结局区域 data.frame → coloc dataset -----------------------------------
build_outcome_dataset <- function(region_df, oc) {
  ## 【v2.8 自查硬化】eaf 必须为开区间 (0,1)：本地镜像中 eaf=0/1 系频率四舍五入
  ## 或单态位点，coloc 的 MAF=pmin(eaf,1-eaf) 会得到 0，check_dataset/方差项不稳；
  ## 直接剔除（数量极少，见头部交集复核口径）。
  d <- region_df[!is.na(region_df$beta) & !is.na(region_df$se) &
                 !is.na(region_df$eaf) & region_df$se > 0 &
                 region_df$eaf > 0 & region_df$eaf < 1, , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  ## 与 API 分支一致的 rsID 唯一化：重复 rsID（多等位/源文件重复行）保留首个
  ## 有效行，避免 align_coloc_pair() 的 match() 在重复键中任意取行。
  if (anyDuplicated(d$SNP)) {
    n_dup <- sum(duplicated(d$SNP))
    d <- d[!duplicated(d$SNP), , drop = FALSE]
    log_msg("  [coloc] ", oc$id, " 结局区域剔除重复 rsID ", n_dup,
            " 行（保留首个有效行，与 API 分支去重口径一致）")
  }
  ## 【v2.8】替代源结局（coloc_n 非空，如步行速度 Neale v3）coloc 样本量用
  ## coloc_n 而非 oc$n（后者是 MR 主线 ukb-b-4711 的 n），并 log 透明声明
  n_use <- if (!is.null(oc$coloc_n)) oc$coloc_n else oc$n
  if (!is.null(oc$coloc_n))
    log_msg("  [coloc] ", oc$id, " 使用替代源样本量 coloc_n=", oc$coloc_n,
            "（≠ MR 主线 n=", oc$n, "；替代源声明见 OUTCOMES 注释与头部 v2.8 块）")
  ds <- list(beta = d$beta, varbeta = d$se^2,
             MAF = pmin(d$eaf, 1 - d$eaf),   # eaf 为效应等位基因频率，coloc 要 MAF
             snp = d$SNP, position = d$pos,
             ## 【N1 修复】结局侧等位基因列必须进入 ds，否则 align_coloc_pair()
             ## 的 has_alleles 恒 FALSE、等位基因协调分支永远不可达；
             ## local 分支 effect_allele/other_allele 为 NA → 走既有降级逻辑
             ea = d$effect_allele, oa = d$other_allele,
             type = ifelse(oc$type == "binary", "cc", "quant"), N = n_use)
  if (oc$type == "binary") ds$s <- oc$s   # 病例比例（OUTCOMES 配置）
  ds
}

## ---- eQTL 侧区域数据：eQTLGen 全量 cis 文件（核查 2 结论落实点）------------
## 全量文件含该基因全部受检 SNP（eQTLGen cis 定义为基因中点 ±1 Mb，与本模块
## cis_window_kb=1000 一致）；MAF 缺失的 SNP 直接剔除（coloc 对 MAF 敏感，
## 不允许用 0.5 兜底）。文件缺失则返回 status，coloc 整列跳过并留痕。
fetch_eqtlgen_full_region <- function(gene_symbol, gc) {
  if (!file.exists(EQTL_CONFIG$eqtlgen_full_file)) {
    msg <- paste0("eQTLGen 全量 cis 文件缺失（", EQTL_CONFIG$eqtlgen_full_file,
                  "）——coloc 跳过；核查 2 已确认官网提供该文件，",
                  "严禁用 FDR0.05 显著配对文件冒充全区域数据")
    log_msg("  [coloc跳过] ", msg)
    ## 【M3 修复】文件缺失属"环境未就绪"而非"生物学阴性"，不缓存（cacheable=FALSE），
    ## 与 API 失败同口径：文件到位后重跑即自动启用，不会被陈旧 NULL 缓存挡回
    return(list(ds = NULL, status = msg, cacheable = FALSE))
  }
  need_cols <- c("Pvalue", "SNP", "SNPChr", "SNPPos", "AssessedAllele",
                 "OtherAllele", "Zscore", "GeneSymbol")
  dat <- cached_fread(EQTL_CONFIG$eqtlgen_full_file, select = need_cols)
  dat <- dat[dat$GeneSymbol == gene_symbol, ]
  if (nrow(dat) == 0) {
    ## 【N3】基因在全量文件无记录属结构性阴性（非环境问题），缓存避免重复扫描 4GB 文件
    return(list(ds = NULL, status = paste0(gene_symbol, " 在 eQTLGen 全量文件中无记录"),
                cacheable = TRUE))
  }
  maf_tbl <- cached_fread(EQTL_CONFIG$eqtlgen_maf_file)
  ## 【fix4】MAF 表列名适配 + eaf 条件推导（共用 adapt_maf_table()；
  ## coloc 路径 na_policy="drop"，eaf 缺失/等位基因不一致者剔除并计数，
  ## 不允许 0.5 兜底）
  dat <- adapt_maf_table(maf_tbl, dat, gene_symbol,
                         context = "fetch_eqtlgen_full_region", na_policy = "drop")
  conv <- z_to_beta_se(dat$Zscore, dat$eaf, EQTL_CONFIG$eqtlgen_n)
  w <- CONFIG$cis_window_kb * 1000
  keep <- dat$SNPChr == gc$chr & dat$SNPPos >= gc$start - w & dat$SNPPos <= gc$end + w
  dat <- dat[keep, ]; conv <- list(beta = conv$beta[keep], se = conv$se[keep])
  if (nrow(dat) < CONFIG$coloc_min_snps) {
    return(list(ds = NULL, status = sprintf("区域内 SNP 数 %d < 阈值 %d",
                                            nrow(dat), CONFIG$coloc_min_snps),
                cacheable = TRUE))   # 结构性阴性，可缓存
  }
  ## 【fix4 要点 4】ds$MAF 改为 pmin(eaf, 1-eaf)（minor AF）——coloc 只需次等位
  ## 基因频率；且 align_coloc_pair 的 pal_mid 与 |ΔMAF|>0.3 比较对象
  ## gwas_ds$MAF 也是 minor（build_outcome_dataset 中 pmin(eaf,1-eaf)），
  ## 两侧均为 minor AF 口径，否则回文位点会被系统性误剔
  ds <- list(beta = conv$beta, varbeta = conv$se^2, MAF = pmin(dat$eaf, 1 - dat$eaf),
             snp = dat$SNP, position = dat$SNPPos,
             ## 【B2 修复】携带等位基因列供 align_coloc_pair() 协调
             ea = dat$AssessedAllele, oa = dat$OtherAllele,
             type = "quant", N = EQTL_CONFIG$eqtlgen_n)
  list(ds = ds, status = "ok")
}

## ---- 两侧 dataset 按 rsID 交集对齐 + 等位基因协调（B2 修复）-----------------
## coloc 要求相同 SNP 集合、顺序一致，且两侧 beta 相对同一效应等位基因。
## 协调规则（逐 SNP，按 eQTL 侧等位基因为基准；N2 修订版，保守优先）：
##   非回文 SNP：
##   1) 同向匹配（ea1==ea2 & oa1==oa2，含 GWAS 侧链互补后同向）：直接用；
##   2) 交换匹配（ea1==oa2 & oa1==ea2，含链互补后交换）：翻转 GWAS 侧 beta
##      符号（MAF=pmin(eaf,1-eaf) 对翻转不变；ds 未存 eaf，无需同步翻转）；
##   回文 SNP（AT/TA/CG/GC）单独分支（互补链即交换链，旧逻辑会误翻同向回文）：
##   3) pal & direct：不翻转直接用；两侧 MAF 一致性 sanity check
##     （|ΔMAF|>0.3 剔除并计数）；
##   4) pal & !direct：一律剔除并计数（无外部参考无法可靠判链，
##      宁可丢 SNP 不引入符号噪声）；
##   5) pal & MAF∈(0.42,0.58) 中间带：无论匹配形态均剔除（优先级最高）；
##   6) 其余无法协调（含等位基因 NA）：剔除并 log 留痕计数。
## local 结局分支等位基因为 NA → 降级：警告日志 + 仅按 rsID 对齐 + 该 key 的
## coloc 结果标注 allele_harmonised="no(local无等位基因列)"，解读时注意。
align_coloc_pair <- function(eqtl_ds, gwas_ds) {
  common <- intersect(eqtl_ds$snp, gwas_ds$snp)
  if (length(common) < CONFIG$coloc_min_snps) {
    return(list(ok = FALSE,
                status = sprintf("两侧 rsID 交集 SNP 数 %d < 阈值 %d",
                                 length(common), CONFIG$coloc_min_snps)))
  }
  idx1 <- match(common, eqtl_ds$snp)
  idx2 <- match(common, gwas_ds$snp)

  has_alleles <- !is.null(eqtl_ds$ea) && !is.null(gwas_ds$ea) &&
                 !all(is.na(eqtl_ds$ea)) && !all(is.na(gwas_ds$ea))
  harmonised_tag <- "yes"
  beta2 <- gwas_ds$beta[idx2]
  keep <- rep(TRUE, length(common))

  if (!has_alleles) {
    log_msg("  ★警告：一侧等位基因列缺失（local 结局分支常见）——本 key 仅按 rsID",
            "对齐、未做等位基因协调，PPH4 解读需谨慎（建议改用 api 模式或补齐等位基因列）")
    harmonised_tag <- "no(无等位基因列)"
  } else {
    comp <- c(A = "T", T = "A", C = "G", G = "C")
    e1 <- toupper(eqtl_ds$ea[idx1]); o1 <- toupper(eqtl_ds$oa[idx1])
    e2 <- toupper(gwas_ds$ea[idx2]); o2 <- toupper(gwas_ds$oa[idx2])
    maf1 <- eqtl_ds$MAF[idx1]; maf2 <- gwas_ds$MAF[idx2]
    na_allele <- is.na(e1) | is.na(o1) | is.na(e2) | is.na(o2)

    same <- !na_allele & e1 == e2 & o1 == o2
    swap <- !na_allele & e1 == o2 & o1 == e2
    e2c <- unname(comp[e2]); o2c <- unname(comp[o2])
    same_c <- !na_allele & !is.na(e2c) & !is.na(o2c) & e1 == e2c & o1 == o2c
    swap_c <- !na_allele & !is.na(e2c) & !is.na(o2c) & e1 == o2c & o1 == e2c
    direct <- same | same_c                   # 同向（含链互补后同向）
    flip <- (swap | swap_c) & !direct         # 交换（含链互补后交换）
    pal <- paste0(e1, o1) %in% c("AT", "TA", "CG", "GC")

    ## 【N2 修复】回文 SNP（AT/TA/CG/GC）单独分支，保守优先：
    ## 回文位点 direct 成立时 swap_c 恒成立（互补链即交换链），旧逻辑
    ## flip=swap|swap_c 会把同向回文误翻转；且带外回文无外部参考无法可靠
    ## 判链。规则：
    ##   ① pal & same（原始同向；链互补后同向 same_c 对回文与交换等价、
    ##     不予承认）：不翻转直接用；两侧 MAF 一致性 sanity check
    ##     （|ΔMAF|>0.3 视为频率异常/非同一变异，剔除并计数）；
    ##   ② pal & !same（仅互补/交换匹配或完全不符）：一律剔除并计数
    ##     （宁可丢 SNP 不引入符号噪声）；
    ##   ③ pal & MAF∈(0.42,0.58) 中间带：无论匹配形态均剔除（优先级最高）。
    pal_mid <- pal & !is.na(maf1) & maf1 > 0.42 & maf1 < 0.58
    ## 关键细化（数值模拟验证）：回文位点仅承认【原始同向】same（ea1==ea2 &
    ## oa1==oa2）；same_c（链互补后同向）对回文与 swap 数学等价、不可判别，
    ## 若承认 same_c 会把"仅交换/仅互补"的回文误判为可保留——故回文的
    ## direct 判定一律用 same，不用 same_c。
    pal_maf_bad <- pal & same & !pal_mid & !is.na(maf1) & !is.na(maf2) &
                   abs(maf1 - maf2) > 0.3
    ## 非回文：维持同向直接用、交换/互补翻转逻辑（含 NA 等位基因剔除）
    np_keep <- !pal & (direct | flip) & !na_allele
    ## 回文：仅 same（原始同向）且非中间带且 MAF 一致才保留，永不翻转
    pal_keep <- pal & same & !pal_mid & !pal_maf_bad & !na_allele
    flip_final <- flip & !pal & !na_allele    # 回文永不翻转
    keep <- np_keep | pal_keep
    beta2[flip_final & keep] <- -beta2[flip_final & keep]
    ## 注：ds 仅存 MAF（=pmin(eaf,1-eaf)，对符号翻转不变），无需同步翻转；
    ## 若未来 ds 携带 eaf，翻转处必须同步 1-eaf（当前已确认无此字段）。
    fail_drop <- !pal & !(direct | flip)      # 非回文无法协调（含 NA）
    pal_nondirect_drop <- pal & !pal_mid & !same  # 回文非原始同向（含 gwas 侧 NA）
    log_msg(sprintf(paste0("  等位基因协调：交集 %d → 保留 %d",
                           "（非回文同向 %d，非回文翻转 %d，回文同向保留 %d；",
                           "剔除：回文中间带 %d + 回文非同向 %d + 回文MAF不一致 %d",
                           " + 无法协调 %d）"),
                    length(common), sum(keep),
                    sum(direct & !pal & keep), sum(flip_final & keep),
                    sum(pal & direct & keep),
                    sum(pal_mid), sum(pal_nondirect_drop),
                    sum(pal_maf_bad), sum(fail_drop)))
    harmonised_tag <- "yes"
  }

  if (sum(keep) < CONFIG$coloc_min_snps) {
    return(list(ok = FALSE,
                status = sprintf("等位基因协调后 SNP 数 %d < 阈值 %d",
                                 sum(keep), CONFIG$coloc_min_snps)))
  }
  sub_ds <- function(ds, idx, keep, beta_override = NULL) {
    b <- if (is.null(beta_override)) ds$beta[idx][keep] else beta_override[keep]
    out <- list(beta = b, varbeta = ds$varbeta[idx][keep], MAF = ds$MAF[idx][keep],
                snp = ds$snp[idx][keep], position = ds$position[idx][keep],
                type = ds$type, N = ds$N)
    ## 【v2.9 修复】`$` 部分匹配陷阱：ds 无 "s" 字段但有 "snp" 时，ds$s
    ## 静默部分匹配返回 ds$snp（SNP 名字符向量）→ out$s 被污染 →
    ## check_dataset 报 "s must be between 0 and 1"（实跑全灭根因）。
    ## 改为 names() 精确判断 + [[ 精确取值（[[ 不做部分匹配）。
    if ("s" %in% names(ds) && !is.null(ds[["s"]])) out$s <- ds[["s"]]
    out
  }
  list(ok = TRUE,
       d1 = sub_ds(eqtl_ds, idx1, keep),
       d2 = sub_ds(gwas_ds, idx2, keep, beta_override = beta2),
       nsnps = sum(keep), allele_harmonised = harmonised_tag, status = "ok")
}

## ---- 混合策略主函数：SuSiE 优先，回退 coloc.abf ------------------------------
run_coloc_hybrid <- function(eqtl_dataset, gwas_dataset, gene_symbol, outcome_label) {
  chk <- tryCatch({
    coloc::check_dataset(eqtl_dataset); coloc::check_dataset(gwas_dataset); TRUE
  }, error = function(e) {
    log_msg("  check_dataset 未通过（", gene_symbol, " × ", outcome_label,
            "）：", conditionMessage(e))
    ## 【v2.9】失败时留 dataset 指纹（names 与 s 的 class/长度/值），便于诊断
    .fp <- function(d, tag) {
      log_msg("    [指纹] ", tag, " names: ", paste(names(d), collapse = ","))
      if ("s" %in% names(d))
        log_msg("    [指纹] ", tag, "$s class=", class(d[["s"]])[1],
                " len=", length(d[["s"]]),
                if (is.numeric(d[["s"]]) && length(d[["s"]]) == 1)
                  paste0(" val=", signif(d[["s"]], 4)) else "")
    }
    .fp(eqtl_dataset, "d1"); .fp(gwas_dataset, "d2")
    FALSE
  })
  if (!chk) return(NULL)

  ## 5.1 先尝试 SuSiE 精细定位（两个性状分别拟合；需两侧均有可信集）
  susie_q <- tryCatch(coloc::runsusie(eqtl_dataset), error = function(e) NULL)
  susie_g <- tryCatch(coloc::runsusie(gwas_dataset), error = function(e) NULL)
  if (!is.null(susie_q) && !is.null(susie_g) &&
      length(susie_q$sets$cs) > 0 && length(susie_g$sets$cs) > 0) {
    ## 【M1 修复】coloc.susie 显式传 p12，与 abf 分支口径一致
    res <- coloc::coloc.susie(susie_q, susie_g, p12 = CONFIG$coloc_prior_p12)
    attr(res, "branch") <- "coloc-SuSiE"
    return(res)
  }
  ## 5.2 回退：单因果变异 coloc.abf
  res <- coloc::coloc.abf(dataset1 = eqtl_dataset, dataset2 = gwas_dataset,
                          p12 = CONFIG$coloc_prior_p12)
  attr(res, "branch") <- "coloc.abf(单因果变异假设)"
  res
}

## ---- 结果判读：PPH4 双阈值（宽 0.8 判定口径主阈值；严 0.9 并列报告）---------
## 注意 coloc.abf 的 $summary 是命名向量，coloc.susie 的是 data.frame
## （每个可信集配对一行），必须分别处理；提取命名向量前一律 unname()
interpret_coloc <- function(coloc_res) {
  branch <- attr(coloc_res, "branch")
  s <- coloc_res$summary
  if (is.data.frame(s)) {
    if (nrow(s) == 0) {
      return(list(PPH4 = NA_real_, PPH3 = NA_real_, PPH4_by_cs_pair = numeric(0),
                  pass_strict = FALSE, pass_relaxed = FALSE,
                  n_cs_pairs = 0L, nsnps = NA_integer_, branch = branch))
    }
    pph4_all <- s$PP.H4.abf
    list(
      PPH4 = max(pph4_all, na.rm = TRUE),
      PPH3 = max(s$PP.H3.abf, na.rm = TRUE),
      PPH4_by_cs_pair = pph4_all,             # 留痕：全部 CS 配对的 PPH4
      pass_strict  = any(pph4_all > CONFIG$coloc_pph4_strict, na.rm = TRUE),
      pass_relaxed = any(pph4_all > CONFIG$coloc_pph4_relaxed, na.rm = TRUE),
      n_cs_pairs = nrow(s),
      nsnps = if ("nsnps" %in% colnames(s)) unname(s$nsnps[1]) else NA_integer_,
      branch = branch)
  } else {
    pph4 <- unname(s["PP.H4.abf"])
    list(
      PPH4 = pph4,
      PPH3 = unname(s["PP.H3.abf"]),
      pass_strict  = pph4 > CONFIG$coloc_pph4_strict,   # >0.9，对齐 Yin et al.
      pass_relaxed = pph4 > CONFIG$coloc_pph4_relaxed,  # >0.8，判定口径主阈值
      n_cs_pairs = 1L,
      nsnps = if ("nsnps" %in% names(s)) unname(s["nsnps"]) else NA_integer_,
      branch = branch)
  }
}

# =============================================================================
# 6. 已知因果基因对照（备忘录附录 D 清单做成 CONFIG 配置块）
# =============================================================================

## 跑完自动比对 10 标志物在双侧清单的命中情况（双侧同命中单独高亮），
## 并汇总每基因的判定口径达成状态（Bonferroni + Steiger + PPH4）
null_na <- function(v) if (is.null(v) || length(v) == 0) NA else unname(v)

check_known_causal_overlap <- function(mr_primary, coloc_summary) {
  kc <- CONFIG$known_causal
  rows <- lapply(CONFIG$genes, function(g) {
    in_sarc <- g %in% kc$sarcopenia
    in_oa   <- g %in% kc$oa
    mp <- mr_primary[mr_primary$gene == g, , drop = FALSE]
    sig_oc <- if (nrow(mp) > 0) mp$outcome[mp$sig_bonferroni %in% TRUE] else character(0)
    steiger_ok <- if (nrow(mp) > 0)
      mp$outcome[!is.na(mp$steiger_p) & mp$steiger_p < 0.05] else character(0)
    cs <- coloc_summary[coloc_summary$gene == g, , drop = FALSE]
    pph4_08 <- if (nrow(cs) > 0)
      unique(cs$outcome[!is.na(cs$PPH4) & cs$PPH4 > CONFIG$coloc_pph4_relaxed]) else character(0)
    pph4_09 <- if (nrow(cs) > 0)
      unique(cs$outcome[!is.na(cs$PPH4) & cs$PPH4 > CONFIG$coloc_pph4_strict]) else character(0)
    causal_ok <- Reduce(intersect, list(sig_oc, steiger_ok, pph4_08))
    data.frame(
      gene = g,
      in_sarcopenia_known_list = in_sarc,            # 肌少症侧（Yin et al. 2024）
      in_oa_known_list = in_oa,                      # OA 侧清单
      dual_hit_both_lists = in_sarc & in_oa,         # 双侧同命中 → 单独高亮
      highlight_priority = if (in_sarc & in_oa) "★双侧同命中：同一基因驱动两病最强证据，最高优先级呈现"
                           else if (in_sarc || in_oa) "单侧命中已发表清单"
                           else "",
      mr_bonferroni_sig_outcomes = paste(sig_oc, collapse = ";"),
      steiger_pass_outcomes = paste(steiger_ok, collapse = ";"),
      coloc_pph4_gt_0.8_outcomes = paste(pph4_08, collapse = ";"),
      coloc_pph4_gt_0.9_outcomes = paste(pph4_09, collapse = ";"),
      causal_established_outcomes = paste(causal_ok, collapse = ";"),
      causal_established_any = length(causal_ok) > 0,
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

## ---- 森林图（主方法效应；二分类为 OR，连续结局 exp(b) 仅供方向可视化）-------
plot_forest <- function(mr_primary) {
  if (is.null(mr_primary) || nrow(mr_primary) == 0) return(invisible(NULL))
  p <- ggplot(mr_primary,
              aes(x = or, y = gene, color = eqtl_source,
                  shape = sig_bonferroni %in% TRUE)) +
    geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
    geom_pointrange(aes(xmin = or_lci, xmax = or_uci),
                    position = position_dodge(width = 0.5)) +
    scale_x_log10() +
    facet_wrap(~outcome, scales = "free_x") +
    labs(x = "效应估计（log 轴；二分类=OR，连续结局为 exp(b) 仅供参考）", y = NULL,
         shape = "Bonferroni 显著", color = "eQTL 来源",
         title = "模块 H：10 锁定标志物主效应森林图（主方法 IVW/Wald）") +
    theme_bw(base_size = 11)
  ggplot2::ggsave(file.path(CONFIG$outdir, "moduleH_MR_forest_primary.pdf"),
                  p, width = 11, height = 7)
  log_msg("森林图已写出：moduleH_MR_forest_primary.pdf")
  invisible(p)
}


# =============================================================================
# 【fix11｜2026-08-29】模块 H SuSiE 敏感性分析（方案二：6 名义显著组合中
#  eQTLGen 侧 5 组合）——1000G EUR LD 注入 + 真 runsusie/coloc.susie
# =============================================================================
# 背景（备忘录 §43）：fix10 诊断石锤 24/24 组合走 coloc.abf 的根因是
#   dataset 缺 LD 矩阵、runsusie 报错被静默吞。本脚本对 5 个名义显著
#   MR 组合（eQTLGen 侧）补做 LD+SuSiE 精细定位敏感性分析：
#     FBLN5×ALM / CBR3×knee_oa / HEMK1×ALM / GDE1×low_grip / RNF14×walk_pace
#   （HEMK1×ALM-GTEx 侧：GTEx v8 起全配对文件不公开，无共定位数据，不在列）
#
# 两阶段运行（FIX11_MODE 切换）：
#   stage1（LD 文件未到位）：跑 5 组合完整构造链，导出对齐后 SNP 列表
#     （rsID/pos/ea/oa）到 results/moduleH/fix11_susie/snplist_*.csv；
#     交给沙箱侧算 1000G EUR LD 矩阵。
#   stage2（LD 文件 ld_<gene>_<outcome>.rds 落位 FIX11_DIR 后）：
#     读 LD + REF/ALT 元数据 → 等位基因符号对齐（四类匹配，不匹配剔除）
#     → 注入 d1/d2$LD → runsusie×2 → coloc.susie → 与 abf 结果对照汇总
#     fix11_susie_summary.csv。
#
# LD 符号对齐口径：LD 方阵符号基准=1000G 正链 ALT 等位基因计数；
#   d1$beta 基准=eQTLGen AssessedAllele（align_coloc_pair 未改 eQTL 侧符号，
#   d2 已翻转到与 d1 ea 一致）。对齐规则：ea==ALT → 不翻；ea==REF → 翻；
#   ea==comp(ALT) → 翻（链反）；ea==comp(REF) → 不翻；其余剔除。
# 第 1-1472 行与 fix9 逐字节一致（同 fix10），主流程替换为本段。
# =============================================================================

FIX11_MODE <- "stage2"          # fix11b：LD 已交付，直接 stage2 真跑
FIX11_DIR  <- file.path("results", "moduleH", "fix11_susie")
dir.create(FIX11_DIR, recursive = TRUE, showWarnings = FALSE)

FIX11_COMBOS <- list(
  list(gene = "FBLN5", outcome = "alm"),
  list(gene = "CBR3",  outcome = "knee_oa"),
  list(gene = "HEMK1", outcome = "alm"),
  list(gene = "GDE1",  outcome = "low_grip"),
  list(gene = "RNF14", outcome = "walk_pace")
)

## ---- 组合构造链（与 fix10 同函数同参数，返回对齐对 + SNP 元数据）-----------
fix11_build_pair <- function(gene0, oc_key) {
  oc0 <- OUTCOMES[[oc_key]]
  gc0 <- gene_coords(gene0)
  eq0 <- fetch_eqtlgen_full_region(gene0, gc0)
  if (!identical(eq0$status, "ok")) return(list(ok = FALSE, status = paste0("eQTL: ", eq0$status)))
  reg0 <- fetch_outcome_region(oc0, gc0)
  gw0  <- build_outcome_dataset(reg0, oc0)
  if (is.null(gw0)) return(list(ok = FALSE, status = "结局 dataset 为 NULL"))
  al0  <- align_coloc_pair(eq0$ds, gw0)
  if (!isTRUE(al0$ok)) return(list(ok = FALSE, status = paste0("align: ", al0$status)))
  ## d1$beta 符号基准 = eq0$ds 的 ea（AssessedAllele）；match 回原始区域取等位基因
  m <- match(al0$d1$snp, eq0$ds$snp)
  meta <- data.frame(snp = al0$d1$snp, pos = al0$d1$position,
                     ea = eq0$ds$ea[m], oa = eq0$ds$oa[m], stringsAsFactors = FALSE)
  list(ok = TRUE, d1 = al0$d1, d2 = al0$d2, meta = meta,
       nsnps = al0$nsnps, allele_harmonised = al0$allele_harmonised)
}

if (FIX11_MODE == "stage1") {
  log_msg("==== fix11 阶段 1：导出 5 组合对齐 SNP 列表 ====")
  for (cb in FIX11_COMBOS) {
    log_msg("---- ", cb$gene, " x ", cb$outcome, " ----")
    bp <- fix11_build_pair(cb$gene, cb$outcome)
    if (!bp$ok) { log_msg("  [跳过] ", bp$status); next }
    fn <- file.path(FIX11_DIR, sprintf("snplist_%s_%s.csv", cb$gene, cb$outcome))
    data.table::fwrite(bp$meta, fn)
    log_msg("  导出 ", fn, "：", nrow(bp$meta), " SNP（等位基因协调=",
            bp$allele_harmonised, "）")
  }
  log_msg("==== fix11 阶段 1 完成：把 snplist_*.csv 5 个文件发回 ====")
}

## =============================================================================
## 阶段 2：LD 注入 + 符号对齐 + runsusie/coloc.susie + abf 对照汇总
## =============================================================================
comp_base <- function(a) chartr("ATCG", "TAGC", a)   # 单碱基互补（仅 A/T/C/G）

if (FIX11_MODE == "stage2") {
  log_msg("==== fix11 阶段 2：LD 注入 + SuSiE 精细定位 ====")
  if (!requireNamespace("Matrix", quietly = TRUE))
    stop("阶段 2 需要 Matrix 包（nearPD 正定化）：install.packages('Matrix') 后重跑")
  ## abf 既有结果（对照列）
  abf_file <- file.path(CONFIG$outdir, "moduleH_coloc_summary.csv")
  abf <- if (file.exists(abf_file)) data.table::fread(abf_file) else NULL

  res_rows <- list()
  for (cb in FIX11_COMBOS) {
    tag <- sprintf("%s x %s", cb$gene, cb$outcome)
    log_msg("---- ", tag, " ----")
    ld_file <- file.path(FIX11_DIR, sprintf("ld_%s_%s.rds", cb$gene, cb$outcome))
    if (!file.exists(ld_file)) {
      log_msg("  [跳过] LD 文件缺失：", ld_file)
      res_rows[[tag]] <- data.frame(gene = cb$gene, outcome = cb$outcome,
                                    status = "ld_file_missing", stringsAsFactors = FALSE)
      next
    }
    ld_obj <- readRDS(ld_file)   # list(snp=, ref=, alt=, LD=<方阵，行序同 snp>)
    bp <- fix11_build_pair(cb$gene, cb$outcome)
    if (!bp$ok) { log_msg("  [跳过] ", bp$status)
      res_rows[[tag]] <- data.frame(gene = cb$gene, outcome = cb$outcome,
                                    status = paste0("build_fail: ", bp$status),
                                    stringsAsFactors = FALSE); next }

    ## ---- LD 子集化到对齐 SNP（顺序以 d1$snp 为准）----
    mm <- match(bp$d1$snp, ld_obj$snp)
    if (any(is.na(mm))) {
      log_msg("  [警告] ", sum(is.na(mm)), "/", length(mm),
              " 对齐 SNP 在 1000G EUR 面板无匹配（多为 indel/低频），剔除")
    }
    keep_ld <- !is.na(mm)
    d1s <- lapply(bp$d1[c("beta","varbeta","MAF","snp","position")], function(x) x[keep_ld])
    d1s$type <- bp$d1$type; d1s$N <- bp$d1$N
    d2s <- lapply(bp$d2[c("beta","varbeta","MAF","snp","position")], function(x) x[keep_ld])
    d2s$type <- bp$d2$type; d2s$N <- bp$d2$N
    if ("s" %in% names(bp$d2)) d2s$s <- bp$d2[["s"]]
    ea_keep <- bp$meta$ea[keep_ld]
    LD <- ld_obj$LD[mm[keep_ld], mm[keep_ld], drop = FALSE]
    ref_k <- ld_obj$ref[mm[keep_ld]]; alt_k <- ld_obj$alt[mm[keep_ld]]

    ## ---- 符号对齐（四类匹配；不匹配剔除）----
    s_sign <- ifelse(ea_keep == alt_k, 1,
              ifelse(ea_keep == ref_k, -1,
              ifelse(ea_keep == comp_base(alt_k), -1,
              ifelse(ea_keep == comp_base(ref_k), 1, NA))))
    n_bad <- sum(is.na(s_sign))
    if (n_bad > 0) log_msg("  [符号] ", n_bad, " SNP 等位基因无法与 1000G 匹配，剔除")
    ok2 <- !is.na(s_sign)
    d1s <- lapply(d1s[c("beta","varbeta","MAF","snp","position")], function(x) x[ok2])
    d1s$type <- "quant"; d1s$N <- EQTL_CONFIG$eqtlgen_n
    d2s <- lapply(d2s[c("beta","varbeta","MAF","snp","position")], function(x) x[ok2])
    d2s$type <- bp$d2$type; d2s$N <- bp$d2$N
    if ("s" %in% names(bp$d2)) d2s$s <- bp$d2[["s"]]
    sv <- s_sign[ok2]
    LD <- LD[ok2, ok2, drop = FALSE]
    ## 符号注入：beta 翻转 + LD 行列翻转
    d1s$beta <- d1s$beta * sv; d2s$beta <- d2s$beta * sv
    LD <- LD * (sv %o% sv)
    LD <- (LD + t(LD)) / 2                      # 强制对称（防脱手）
    diag(LD) <- 1
    log_msg("  [对齐] 最终 ", length(d1s$snp), " SNP；符号翻转 ", sum(sv < 0), " 个")

    ## ---- LD 正定兜底（503 样本面板秩亏常见）----
    if (any(!is.finite(LD))) stop("LD 含非有限值，中止（检查矩阵文件）")
    ## 【fix11c 修复 1】runsusie 硬性要求 LD dimnames 与 dataset$snp 对应——
    ## Python 端矩阵无名，实跑报错 "colnames in LD do not contain all SNPs"
    dimnames(LD) <- list(d1s$snp, d1s$snp)
    ## 【fix11c 修复 2】正定检查由 eigen（5057² 约 12 分钟）改为 chol 试探（秒级）
    pd_ok <- tryCatch({ chol(LD); TRUE }, error = function(e) FALSE)
    if (!pd_ok) {
      log_msg("  [LD] chol 失败（面板秩亏）→ Matrix::nearPD 正定化")
      LD <- as.matrix(Matrix::nearPD(LD, corr = TRUE)$mat)
      dimnames(LD) <- list(d1s$snp, d1s$snp)   # nearPD 后 dimnames 丢失需重设
    }

    d1s$LD <- LD; d2s$LD <- LD

    ## ---- 真 SuSiE（fix11d：直调 susie_rss + 外部 LD 分级加固）----------------
    ## 背景：coloc::runsusie 实跑报 "estimated prior variance is unreasonably
    ##   large"（强 eQTL z 信号 × 503 人面板 LD 噪声 → z-R 不一致，SuSiE RSS
    ##   外部 LD 场景的已知问题；HEMK1 d1 实跑实证，d2 弱信号侧成功）。
    ##   runsusie 不暴露 susie_rss 关键参数，故直调并分级降级：
    ##   L1: 原始 R + estimate_residual_variance=TRUE（官方外部 LD 推荐）；
    ##   L2: R shrinkage（0.9R+0.1I，吸收参考面板抽样噪声）+ est.resid；
    ##   L3: shrink + 固定残差与先验（最保守，必出结果）。
    run_susie_ext <- function(ds, tag) {
      z <- ds$beta / sqrt(ds$varbeta)
      R <- ds$LD
      mk <- function(fit, Rm) {
        fit$sets <- susieR::susie_get_cs(fit, Xcorr = Rm)
        fit
      }
      lev <- NULL
      fit <- tryCatch(
        { lev <- "L1(rawR+est.resid)"
          mk(susieR::susie_rss(z, R, n = ds$N, L = 10, max_iter = 1000,
                               estimate_residual_variance = TRUE), R) },
        error = function(e) e)
      if (inherits(fit, "error")) {
        log_msg("  [SuSiE-", tag, "] L1 失败：", conditionMessage(fit), " → L2 shrink")
        R2 <- 0.9 * R + 0.1 * diag(nrow(R)); dimnames(R2) <- dimnames(R)
        fit <- tryCatch(
          { lev <- "L2(shrink0.1+est.resid)"
            mk(susieR::susie_rss(z, R2, n = ds$N, L = 10, max_iter = 1000,
                                 estimate_residual_variance = TRUE), R2) },
          error = function(e) e)
        if (inherits(fit, "error")) {
          log_msg("  [SuSiE-", tag, "] L2 失败：", conditionMessage(fit), " → L3 固定先验")
          fit <- tryCatch(
            { lev <- "L3(shrink0.1+fixed)"
              mk(susieR::susie_rss(z, R2, n = ds$N, L = 10, max_iter = 1000,
                                   estimate_residual_variance = FALSE,
                                   estimate_prior_variance = FALSE), R2) },
            error = function(e) e)
        }
      }
      list(fit = fit,
           level = if (inherits(fit, "error")) paste0("FAIL: ", conditionMessage(fit)) else lev)
    }

    r1 <- run_susie_ext(d1s, "d1")
    r2 <- run_susie_ext(d2s, "d2")
    sq <- r1$fit; sg <- r2$fit

    ## ---- fix11e：收敛审计 ----
    ## IBSS 默认 max_iter=100；fix11d 警告 9 提示至少一个拟合未收敛。
    ## 本版 max_iter=1000 并显式记录 converged/niter，拟合对象存盘供复核。
    conv1 <- if (inherits(sq, "error")) NA else isTRUE(sq$converged)
    conv2 <- if (inherits(sg, "error")) NA else isTRUE(sg$converged)
    nit1  <- if (inherits(sq, "error")) NA_integer_ else as.integer(sq$niter)
    nit2  <- if (inherits(sg, "error")) NA_integer_ else as.integer(sg$niter)
    log_msg("  [收敛审计] d1 converged=", conv1, " niter=", nit1,
            " | d2 converged=", conv2, " niter=", nit2)
    saveRDS(list(d1_fit = sq, d2_fit = sg,
                 d1_level = r1$level, d2_level = r2$level,
                 d1_converged = conv1, d2_converged = conv2,
                 d1_niter = nit1, d2_niter = nit2),
            file.path(FIX11_DIR, paste0("fix11e_fits_", tag, ".rds")))
    d1_err <- if (inherits(sq, "error")) r1$level else NA_character_
    d2_err <- if (inherits(sg, "error")) r2$level else NA_character_
    ncs1 <- if (inherits(sq, "error")) NA_integer_ else length(sq$sets$cs)
    ncs2 <- if (inherits(sg, "error")) NA_integer_ else length(sg$sets$cs)
    log_msg("  [SuSiE] d1 可信集=", ifelse(is.na(ncs1), paste0("ERROR: ", d1_err), ncs1),
            "（", ifelse(is.na(ncs1), "-", r1$level), "）",
            " | d2 可信集=", ifelse(is.na(ncs2), paste0("ERROR: ", d2_err), ncs2),
            "（", ifelse(is.na(ncs2), "-", r2$level), "）")

    susie_pph4 <- NA_real_; susie_pph3 <- NA_real_; n_cs_pairs <- NA_integer_
    if (!inherits(sq, "error") && !inherits(sg, "error") && ncs1 > 0 && ncs2 > 0) {
      cs <- coloc::coloc.susie(sq, sg, p12 = CONFIG$coloc_prior_p12)
      sdf <- cs$summary
      if (is.data.frame(sdf) && nrow(sdf) > 0) {
        susie_pph4 <- max(sdf$PP.H4.abf, na.rm = TRUE)
        susie_pph3 <- max(sdf$PP.H3.abf, na.rm = TRUE)
        n_cs_pairs <- nrow(sdf)
      }
      log_msg("  [coloc.susie] CS 对数=", ifelse(is.na(n_cs_pairs), 0, n_cs_pairs),
              " | max PPH4=", signif(susie_pph4, 4), " | max PPH3=", signif(susie_pph3, 4))
    }

    ## ---- abf 对照列 ----
    abf_pph4 <- abf_pph3 <- NA_real_
    if (!is.null(abf)) {
      ## 三键匹配：outcome 列可能是 key（"alm"）、OpenGWAS id（"ebi-a-..."）或中文 label
      oc0 <- OUTCOMES[[cb$outcome]]
      ocol <- as.character(abf$outcome)
      hit <- abf[abf$gene == cb$gene &
                 (ocol == cb$outcome | ocol == oc0$id | ocol == oc0$label |
                  grepl(cb$outcome, oc0$id, ignore.case = TRUE)), ]
      if (nrow(hit) >= 1) { abf_pph4 <- hit$PPH4[1]; abf_pph3 <- hit$PPH3[1] }
      else log_msg("  [对照] abf 汇总表未匹配到本组合行（不影响 SuSiE 结果）")
    }

    res_rows[[tag]] <- data.frame(
      gene = cb$gene, outcome = cb$outcome, status = "ok",
      n_snp_final = length(d1s$snp),
      d1_cs = ncs1, d2_cs = ncs2, n_cs_pairs = n_cs_pairs,
      susie_PPH4 = susie_pph4, susie_PPH3 = susie_pph3,
      abf_PPH4 = abf_pph4, abf_PPH3 = abf_pph3,
      d1_level = if (inherits(sq, "error")) NA_character_ else r1$level,
      d2_level = if (inherits(sg, "error")) NA_character_ else r2$level,
      d1_converged = conv1, d1_niter = nit1,
      d2_converged = conv2, d2_niter = nit2,
      d1_error = d1_err, d2_error = d2_err, stringsAsFactors = FALSE)
  }

  res_tab <- do.call(rbind, res_rows)
  outf <- file.path(FIX11_DIR, "fix11e_susie_summary_convergence.csv")
  data.table::fwrite(res_tab, outf)
  log_msg("==== fix11e 收敛审计完成：", outf, " ====")
  print(res_tab)
}
