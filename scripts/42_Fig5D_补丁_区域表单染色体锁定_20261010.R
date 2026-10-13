# ============================================================
# 42_Fig5D_补丁_区域表单染色体锁定_20261010.R
# 用途：修复 Fig 5D 把 7 条染色体的位点压到同一坐标轴、轴标签错标 chr16 的问题
# 根因：region_sumstats_FBLN5_alm.csv 实为 7 个基因位点（chr1/3/5/10/14/16/21，
#       各约 4 Mb 窗口）的合并文件；脚本 21 panel_d 取 chr 列首值作轴标签。
# 修法（双保险，建议两个都做）：
#   A. 输入文件替换为 region_sumstats_FBLN5_alm_chr14only.csv（已交付，28,217 行，
#      chr14: 90.38–94.38 Mb，含 lead rs28529055 [92,428,216 bp, hg19]）
#   B. 在 21_fig5_MR_共定位_ClinRheum风格_20260915_via5.R 的 panel_d() 中插入
#      下面的防御代码块（插入点：cch <- pick_col(...) 之后、d <- data.frame(...) 之前）
# ============================================================

message("== 脚本版本：42_Fig5D_补丁_20261010（插入 21 号脚本 panel_d 函数内） ==")

## ---- 插入开始（panel_d 内，pick_col 列名探测之后）----
  ## --- [补丁 2026-10-10] 多染色体防御：仅保留 lead SNP 所在染色体 ---
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
## ---- 插入结束 ----

# 自检（插入后运行脚本 21 时控制台应出现）：
#   [Fig5][d] 区域表含 7 条染色体 -> 仅保留 lead 所在 chr14（28217/177822 行）
# 出图验收：x 轴 = "chr14 (Mb)"；lead rs28529055 菱形位于 ~92.4 Mb；
#   ~93.1 Mb 处独立高峰（-log10P 最高 61）保留且为灰色（LD NA）——
#   该形态直观支持正文 "linked but independent variants" 结论，无需回避。
