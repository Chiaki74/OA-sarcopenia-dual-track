# Archive manifest — 2026-10-12 push

本清单记录 2026-10-12 归档推送时 scripts/ 目录的选取口径：**每个分析步骤仅收录最终封版（终版）脚本**，被取代的中间迭代留存于分析工作站本地归档，不入库。

## 终版选取要点（被取代版本 → 入库终版）

| 步骤 | 入库终版 | 被取代版本（留存本地） |
|---|---|---|
| 06 moduleE 外部验证 | `06_moduleE_validation_fix5_20260823.R` | 原版 + fix1–fix4 |
| 07 moduleF 免疫浸润 | `07_moduleF_immune_infiltration_fix3_20260827.R` | 原版 + fix1–fix2 |
| 08/08a/09 moduleG 单细胞 | fix3 / fix2 / fix3 各终版 | 各自前序 fix |
| 10–13e moduleH MR/共定位/SuSiE | fix8、fix9、fix10、fix11 + stage2/11c/11d/11e | 各自前序 fix |
| 11 moduleI 可成药 | `11_moduleI_druggable_fix1_20260826.R` | 原版 |
| 12 moduleI 列线图/分类器 | `12_*.R` + 12b + 12d + 12e | 12c（被 12d 带进度版取代） |
| 14 moduleK CellChat | `14c_moduleK_CellChat_分组标签与presto修正_20260914.R` | 14、14b（修正链中间版） |
| 16–21 主图绘制 | 16_via3、17_via9、18_via23、19_via12、20_via9、21_via14 | 全部前序 via 版 |
| 27/29/30/31 出图数据与拼图 | 27_via3、29_via8、30_via3、31_via3 | 各自前序 via 版 |
| 43 GSE168505 预测重建 | `43_..._via3.R`（复制脚本 23 管线 + 三重闸门校验） | 43 原版、via2（口径证伪，弃用） |
| fix14 CIBERSORTx 混合物 | `fix14v4_CIBERSORTx混合物制备_20260903.R` | fix14、v2、v3 |

## 整体排除项（有意不入库）

- `plot_fig2–fig10_*_20260903` 系列：2026-09-03 初版主图脚本，已被 16–21 号 ClinRheum 风格终版及 32–38 号重绘脚本整体取代。
- `moduleH_TwoSampleMR_coloc_pipeline.R`：moduleH 早期底稿，已由 10_moduleH 起的有日期编号链取代。
- `G0_aria2_fast_download_*.bat`：本机下载工具，与分析管线无关。
- 原始数据与 `results/` 产物：GEO 公共数据按再分发条款不入库；中间结果可由管线复算（见 .gitignore）。
- CIBERSORT 源代码与 LM22 矩阵：Stanford 许可限制，不入库（正文 Methods 已声明获取途径）。

## 封版表（tables/）

表 S1 源文件（moduleA 状态表）、S11、S12（2026-10-08 封版）、S13 附、S16、S19、S20 及 S20 附两件、fix15 CIBERSORT 双版对照表（S9 源文件）、GSM 分组对照表。表 S2–S8、S10、S14、S15、S17、S18 的载体为 `results/` 下逐队列/逐模块全量表，按 .gitignore 口径不入库；其来源文件路径见 Online Resource 1 第二节各表注。
