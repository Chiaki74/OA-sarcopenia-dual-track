## =============================================================================
## 00_setwd_根目录定位.R
## 用法：在每个脚本最开头 source 本文件，或直接把本段粘到脚本开头。
## 逻辑：① 若手工指定 ROOT_MANUAL 且存在 → 用之；
##       ② 否则从"当前脚本所在目录"（Rscript 与 RStudio 均兼容）向上最多 8 层，
##          找含标志目录（data / results / 01_raw / GEO / moduleI_data）的那一层；
##       ③ 都找不到则回退当前工作目录并警告。
## =============================================================================

## ---- ① 手工指定（最稳妥，推荐改成你的项目根目录绝对路径）-------------------
ROOT_MANUAL <- ""   ## 例：ROOT_MANUAL <- "D:/OA_sarcopenia/project"

## ---- ② 标志目录（你的项目根目录下存在的任意一个即可）-----------------------
ROOT_MARKERS <- c("data", "results", "01_raw", "GEO", "moduleI_data", "figures")

## 当前脚本所在目录的探测（Rscript --file= 与 RStudio Source 两种场景）
.script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- grep("--file=", a, value = TRUE)
  if (length(f) >= 1) return(dirname(normalizePath(sub("--file=", "", f[1]))))
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable() && rstudioapi::getActiveDocumentContext()$path != "")
    return(dirname(normalizePath(rstudioapi::getActiveDocumentContext()$path)))
  getwd()
}

locate_root <- function(markers = ROOT_MARKERS, max_up = 8) {
  if (nzchar(ROOT_MANUAL) && dir.exists(ROOT_MANUAL)) return(normalizePath(ROOT_MANUAL))
  d <- .script_dir()
  for (i in seq_len(max_up)) {
    if (any(dir.exists(file.path(d, markers)))) return(normalizePath(d))
    parent <- dirname(d)
    if (parent == d) break
    d <- parent
  }
  warning("未自动定位根目录（标志目录：", paste(markers, collapse = "/"),
          "），使用当前工作目录：", getwd(), "；建议改用 ROOT_MANUAL 手工指定。")
  normalizePath(getwd())
}

PROJ_ROOT <- locate_root()
setwd(PROJ_ROOT)
message("[setwd] 项目根目录：", PROJ_ROOT)
