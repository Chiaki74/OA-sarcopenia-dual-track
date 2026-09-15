################################################################################
## OA × 肌少症整合生信研究 — 第 1 周任务①：R 环境自检 + 一键安装脚本
## 文件名：setup_R_environment.R
## 版本：v1.0（2026-08-18）  配套：实验计划书 v1.4a 附录 C / 部署操作指引 Step 1
##
## 功能：
##   1) 检查 R 版本（要求 >= 4.3.0，不满足则中止并给出升级指引）
##   2) 分三通道安装全部依赖：CRAN / Bioconductor / GitHub
##   3) 版本敏感项校验：Seurat >= 5.0、coloc >= 6.0（含 SuSiE 接口）
##   4) 安装后自检：逐包加载 + 关键函数存在性检查 + OpenGWAS JWT 预检
##   5) 输出环境报告：setup_report_<日期>.csv + sessionInfo 存档
##
## 用法（Windows / RStudio）：
##   方式一（推荐）：在 RStudio 中打开本文件，点击 "Source" 按钮；
##   方式二：命令行  Rscript setup_R_environment.R
##
## 注意：
##   - 首次完整安装约需 40–120 分钟（取决于网络与需编译的包数量）；
##   - 脚本可重复运行（幂等）：已安装且版本合格的包会自动跳过；
##   - CIBERSORT 不在任何公共仓库（需斯坦福官网注册获取源码），本脚本不安装，
##     请在第 3 周免疫浸润模块启动前按《部署操作指引》附录手动配置；
##   - 需要编译源码的包依赖 Rtools44（Windows），请先完成指引 Step 0。
################################################################################

## ============================================================================
## 0. 全局配置
## ============================================================================
MIN_R_VERSION  <- "4.3.0"
SEURAT_MIN     <- "5.0.0"     # 模块 G 单细胞分析要求 Seurat v5 架构
COLOC_MIN      <- "6.0.0"     # 模块 H 共定位要求 coloc v6+（coloc.susie / runsusie）
REPORT_FILE    <- paste0("setup_report_", format(Sys.Date(), "%Y%m%d"), ".csv")
SESSIONINFO_FILE <- paste0("sessionInfo_", format(Sys.Date(), "%Y%m%d"), ".txt")

options(timeout = 600)        # 大文件下载放宽到 10 分钟
## 先装 BiocManager 并把 CRAN+Bioconductor 仓库合并入 options(repos)，
## 使 CRAN 通道安装（如 WGCNA）能自动解析其 Bioconductor 依赖（impute 等）。
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
## 国内网络如镜像慢，可在下方加 BiocManager 镜像参数；默认使用官方仓库集合：
options(repos = BiocManager::repositories())

msg <- function(...) cat(sprintf(...), "\n")

## ============================================================================
## 1. R 版本检查（硬性门槛）
## ============================================================================
msg("============================================")
msg("[1/6] 检查 R 版本 ...")
msg("当前 R 版本：%s", R.version.string)

if (getRversion() < MIN_R_VERSION) {
  msg("【中止】当前 R 版本低于 %s。", MIN_R_VERSION)
  msg("请先到 https://cloud.r-project.org/bin/windows/base/ 下载安装最新版 R（4.4.x），")
  msg("安装后重启 RStudio，再重新运行本脚本。详见《部署操作指引》Step 0。")
  stop(sprintf("R version < %s", MIN_R_VERSION))
}
msg("✔ R 版本满足要求（>= %s）", MIN_R_VERSION)

## ============================================================================
## 2. 包清单（与实验计划书附录 C 对应）
## ============================================================================
## 2.1 CRAN 包
cran_pkgs <- c(
  ## 基础工具
  "data.table", "dplyr", "tidyr", "stringr", "readr", "openxlsx",
  "ggplot2", "pheatmap", "RColorBrewer", "ggrepel", "patchwork", "cowplot",
  "remotes", "usethis", "pkgbuild",
  ## 模块 B：共享基因识别（WGCNA 的 Bioc 依赖 impute/preprocessCore 见 bioc_pkgs）
  "RobustRankAggreg", "WGCNA",
  ## 模块 D：机器学习诊断标志物
  "glmnet", "e1071", "randomForest", "xgboost", "caret",
  ## 模块 D：模型评估（ROC + 校准 + DCA）
  "pROC", "rms", "ResourceSelection", "dcurves",
  ## 模块 F/G：免疫浸润与单细胞
  "xCell", "Seurat", "SeuratObject", "harmony",
  ## 模块 H：MR（CRAN 部分；coloc v6 仅在 GitHub，见 github_pkgs）
  "ieugwasr", "susieR",
  ## 待决定去留项（低成本先装，用不用在第 3–5 周决策点定）
  "SHAPforxgboost", "shapviz"
)

## 2.2 Bioconductor 包
bioc_pkgs <- c(
  ## 模块 A：数据获取与预处理
  "GEOquery", "limma", "DESeq2", "sva", "affy", "oligo",
  ## 模块 B：WGCNA 的 Bioconductor 依赖（务必经 Bioc 通道，勿走 CRAN）
  "impute", "preprocessCore",
  ## 模块 C：功能富集
  "clusterProfiler", "org.Hs.eg.db", "DOSE", "enrichplot", "ReactomePA",
  ## 模块 E/F：验证与免疫浸润
  "GSVA", "ComplexHeatmap", "SummarizedExperiment",
  ## 待决定去留项
  "UCell"
)

## 2.3 GitHub 包（命名向量：名字 = library 包名，值 = GitHub 仓库）
## 注意：rondolab/MR-PRESSO 安装后的包名是 MRPRESSO（无连字符）；
## coloc v6+ 未发布到 CRAN（CRAN 停留在 5.x），必须走 GitHub 才能满足 >= 6.0。
github_pkgs <- c(
  "TwoSampleMR" = "MRCIEU/TwoSampleMR",   # 模块 H 主 MR 框架
  "MRPRESSO"    = "rondolab/MR-PRESSO",   # 模块 H 多效性离群检测
  "coloc"       = "chr1swallace/coloc",   # 模块 H 共定位 v6+（coloc.susie/runsusie）
  "CellChat"    = "jinworks/CellChat"     # 模块 G 细胞通讯（v2.x，Seurat v5 兼容分支）
)

## ============================================================================
## 3. 安装函数（幂等：已装且版本合格则跳过）
## ============================================================================
need_install <- function(pkg, min_version = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(TRUE)
  if (!is.null(min_version)) {
    v <- tryCatch(as.character(packageVersion(pkg)), error = function(e) "0")
    if (utils::compareVersion(v, min_version) < 0) return(TRUE)
  }
  FALSE
}

msg("============================================")
msg("[2/6] 安装 Bioconductor 包（%d 个，先于 CRAN 通道以便解析依赖）...", length(bioc_pkgs))
for (p in bioc_pkgs) {
  if (need_install(p)) {
    msg("  → 安装 %s ...", p)
    tryCatch(BiocManager::install(p, update = FALSE, ask = FALSE, quiet = TRUE),
             error = function(e) msg("  【失败】%s：%s", p, conditionMessage(e)))
  } else {
    msg("  ✔ %s（%s）已就位", p, as.character(packageVersion(p)))
  }
}

msg("============================================")
msg("[3/6] 安装 CRAN 包（%d 个）...", length(cran_pkgs))
for (p in cran_pkgs) {
  min_v <- switch(p, "Seurat" = SEURAT_MIN, NULL)
  if (need_install(p, min_v)) {
    msg("  → 安装/升级 %s ...", p)
    tryCatch(install.packages(p, dependencies = TRUE, quiet = TRUE),
             error = function(e) msg("  【失败】%s：%s", p, conditionMessage(e)))
  } else {
    msg("  ✔ %s（%s）已就位", p, as.character(packageVersion(p)))
  }
}

msg("============================================")
msg("[4/6] 安装 GitHub 包（%d 个）...", length(github_pkgs))
for (p in names(github_pkgs)) {
  g <- github_pkgs[[p]]
  min_v <- switch(p, "coloc" = COLOC_MIN, NULL)   # coloc 必须 >= 6.0（v6 仅在 GitHub）
  if (need_install(p, min_v)) {
    msg("  → 安装 %s（GitHub: %s）...", p, g)
    tryCatch(remotes::install_github(g, upgrade = "never", quiet = TRUE),
             error = function(e) msg("  【失败】%s：%s（网络受限时可稍后在 RStudio 手动重试）", p, conditionMessage(e)))
  } else {
    msg("  ✔ %s（%s）已就位", p, as.character(packageVersion(p)))
  }
}

## ============================================================================
## 4. 版本敏感项强制校验
## ============================================================================
msg("============================================")
msg("[5/6] 版本敏感项校验 ...")
version_check <- function(pkg, min_v, why) {
  ok <- requireNamespace(pkg, quietly = TRUE) &&
        utils::compareVersion(as.character(packageVersion(pkg)), min_v) >= 0
  msg("  %s %s（要求 >= %s，当前 %s）—— %s",
      if (ok) "✔" else "✘", pkg, min_v,
      tryCatch(as.character(packageVersion(pkg)), error = function(e) "未安装"), why)
  ok
}
v_seurat <- version_check("Seurat", SEURAT_MIN, "模块 G 单细胞/单核分析框架")
v_coloc  <- version_check("coloc",  COLOC_MIN,  "模块 H 共定位 coloc-SuSiE 混合策略")

## ============================================================================
## 5. 全量自检：逐包加载 + 关键函数存在性 + OpenGWAS JWT 预检
## ============================================================================
msg("============================================")
msg("[6/6] 全量自检（加载测试 + 关键函数检查）...")

all_pkgs <- c(cran_pkgs, bioc_pkgs, names(github_pkgs))

## 关键函数存在性检查（函数 → 所属包 → 用途）
key_functions <- list(
  list(fun = "coloc.susie",         pkg = "coloc",       use = "共定位 SuSiE 主线（多因果变异）"),
  list(fun = "runsusie",            pkg = "coloc",       use = "SuSiE 精细定位封装"),
  list(fun = "coloc.abf",           pkg = "coloc",       use = "共定位回退方案（单因果变异）"),
  list(fun = "directionality_test", pkg = "TwoSampleMR", use = "Steiger 方向性检验"),
  list(fun = "mr_presso",           pkg = "MRPRESSO",    use = "多效性离群值检验"),
  list(fun = "api_status",          pkg = "ieugwasr",    use = "OpenGWAS 连接状态"),
  list(fun = "susie",               pkg = "susieR",      use = "SuSiE 底层引擎"),
  list(fun = "AggregateExpression", pkg = "Seurat",      use = "Seurat v5 伪 bulk 接口（v5 标志函数）")
)
fun_check <- function(fun, pkg) {
  requireNamespace(pkg, quietly = TRUE) && exists(fun, envir = asNamespace(pkg))
}

## OpenGWAS JWT 预检（模块 H 的咽喉要道；未配置不阻塞本脚本，但给出醒目提示）
jwt <- Sys.getenv("OPENGWAS_JWT")
jwt_ok <- nzchar(jwt)
if (jwt_ok) {
  msg("  ✔ 检测到 OPENGWAS_JWT 环境变量（长度 %d 字符）", nchar(jwt))
  msg("    正式连通性验证请按《部署操作指引》Step 2.4 执行（需联网调用 API）。")
} else {
  msg("  ⚠ 未检测到 OPENGWAS_JWT —— MR 模块（模块 H）将无法访问 OpenGWAS。")
  msg("    请按《部署操作指引》Step 2 注册账号、获取 JWT 并写入 .Renviron 后重启 R。")
}

## 汇总报告
report <- data.frame(
  package   = all_pkgs,
  channel   = c(rep("CRAN", length(cran_pkgs)),
                rep("Bioconductor", length(bioc_pkgs)),
                rep("GitHub", length(github_pkgs))),
  installed = vapply(all_pkgs, function(p) requireNamespace(p, quietly = TRUE), logical(1)),
  version   = vapply(all_pkgs, function(p)
    tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_), character(1)),
  load_ok   = NA,
  stringsAsFactors = FALSE
)
report$load_ok <- vapply(all_pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) return(FALSE)
  tryCatch({ suppressPackageStartupMessages(library(p, character.only = TRUE)); TRUE },
           error = function(e) FALSE)
}, logical(1))

fun_report <- data.frame(
  function_name = vapply(key_functions, `[[`, character(1), "fun"),
  package       = vapply(key_functions, `[[`, character(1), "pkg"),
  purpose       = vapply(key_functions, `[[`, character(1), "use"),
  available     = vapply(key_functions, function(x) fun_check(x$fun, x$pkg), logical(1))
)

write.csv(report, REPORT_FILE, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(fun_report, sub("\\.csv$", "_keyfun.csv", REPORT_FILE), row.names = FALSE, fileEncoding = "UTF-8")
writeLines(capture.output(sessionInfo()), SESSIONINFO_FILE)

## ============================================================================
## 6. 总结输出
## ============================================================================
n_fail  <- sum(!report$installed | !report$load_ok)
n_ffail <- sum(!fun_report$available)

msg("")
msg("==================== 环境搭建总结 ====================")
msg("包总数 %d：成功 %d，失败 %d", nrow(report), nrow(report) - n_fail, n_fail)
if (n_fail > 0) {
  msg("失败清单：")
  for (p in report$package[!report$installed | !report$load_ok]) msg("  ✘ %s", p)
  msg("→ 请先按《部署操作指引》Step 1.4 常见报错表处理后重跑本脚本（已装包会自动跳过）。")
}
msg("关键函数 %d 项：可用 %d，缺失 %d", nrow(fun_report), nrow(fun_report) - n_ffail, n_ffail)
if (n_ffail > 0) {
  for (i in which(!fun_report$available))
    msg("  ✘ %s::%s（%s）", fun_report$package[i], fun_report$function_name[i], fun_report$purpose[i])
}
msg("版本敏感项：Seurat %s / coloc %s", if (v_seurat) "✔" else "✘", if (v_coloc) "✔" else "✘")
msg("OpenGWAS JWT：%s", if (jwt_ok) "✔ 已配置" else "⚠ 未配置（Step 2 完成后再验证）")
msg("报告文件：%s", REPORT_FILE)
msg("=====================================================")

if (n_fail == 0 && n_ffail == 0 && v_seurat && v_coloc) {
  msg("★ 环境搭建【通过】，可进入下一步：Step 2（OpenGWAS JWT 配置）。")
} else {
  msg("★ 环境搭建【未完全通过】，请按上述失败项处理后重跑。")
}
