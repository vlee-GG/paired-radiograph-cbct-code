# 坐标重复读数与差异复核方案

## 1. 目的

评估同一观察者两次探测器中心坐标读数的重复性。两次读数间隔约一年，第二次读数时不显示第一次的结果。对二维 centre difference 大于 2 mm 的记录，由原观察者本人复核，并完整保存两次原始读数、后续复核读数、决定及理由。

## 2. 坐标定义

每个接收器中心记录两个距离，单位为 cm：第一列为中心到网格下边界的距离，第二列为中心到网格右边界的距离。边界指 20 × 20 cm 网格刻度线的边缘，而非透明板边缘。

## 3. 触发条件

对匹配图片计算二维 centre difference：

`sqrt((Δx)^2 + (Δy)^2)`，其中 cm 差值乘以 10 转换为 mm。

主要分析使用 2,244 条完整配对记录的冻结复核前版本。对 63 条结构化复核记录，从审计记录恢复其进入复核时的既有坐标；其余 2,181 条采用未发生变动的既有坐标。由此重建的复核前数据中恰有 63 条 centre difference 大于 2 mm，全部进入同一观察者的后续复核。

## 4. 同一观察者的后续差异复核

原观察者重新查看原始图片，对63组差异进行复核。结构化记录保留后续读数及其与前两次读数的二维距离。后续读数允许修改，每次修改均进入审计日志；若修改发生在决定后，旧决定自动失效并需重新确认。第二次读数对第一次结果设盲；后续差异复核不作为独立观察者或第三人仲裁验证。`third`、`third_reader_id`及相关文件名保留为历史程序字段，表示后续读数和复核阶段，不代表另一名观察者。

## 5. 最终裁决

- `retain_existing`：第三次读数支持既有读数，最终采用既有坐标。
- `correct_to_new`：第三次读数支持新读数，最终采用新坐标。
- `coordinate_uncertain`：三次读数仍有明显分歧，保留不确定标记。本次最终坐标版本采用同一观察者的后续复核读数，但不得在论文中将其描述为已消除不确定性。

原因使用标准英文枚举：`existing_reading_incorrect` 或 `new_reading_incorrect`。

## 6. 更正后的下游验证

对 `correct_to_new` 记录依次完成：重新生成 crop、文件 QC、crop 重提取验证和自动检索验证。四项结果分别记录，不以坐标差异小于 2 mm 代替下游验证。

## 7. 留痕与可追溯性

公开两次原始读数、第三次读数、决定、原因、操作者标识、时间戳、最终坐标和下游验证状态。审计日志采用 JSONL 追加式记录；发布包使用 SHA-256 清单验证文件完整性。

## 8. 论文报告

在 2,244 条冻结复核前完整配对记录中报告二维 centre difference 的均值与标准差、中位数与四分位距、横纵坐标平均绝对差、Bland–Altman 95% 一致性界限、双向混合效应、单次测量、绝对一致性ICC（ICC(A,1)），以及 ≤1 mm、>1–2 mm、>2 mm、>2–5 mm 和 >5 mm 的数量与比例。另报告由 >2 mm 阈值触发同一观察者复核的 63 条记录、触发比例、裁决结果、理由及 `coordinate_uncertain` 数量，并发布裁决后的最终 target。复核后的 target 与参与裁决的新读数之间的一致性不作为主要结果，因为该比较会受到纳入偏倚影响；如需计算，仅限内部 QC 并清楚标注。完成全部下游验证前，应明确说明相关验证仍在进行中。


## Released evidence

The two 335-file coordinate-reading snapshots are in `validation/coordinate_review/01_coordinates/`. The repeatability analysis uses `04_analysis/record_comparison_pre_review.csv`; the 63-case initial decisions are in `04_analysis/review_records.csv`, and the subsequent workflow export is in `02_review_records/third_review_records.json`. Current published coordinates remain in `source_cbct/`; combined release checks are in `validation/release_qc/qc_summary.json` (release-root paths). Historical website working copies and event logs are not distributed.


## 观察者与统计定义更正（2026-10-19）

作者确认两次原始读数及后续差异复核均由同一观察者完成。两次读数间隔约一年，第二次看不到第一次结果。此前“两名观察者”“第三人仲裁”的描述据此更正。原始数据和历史记录不改写；归档中的旧称谓属于历史标签，应结合当前方案解释。绝对一致性ICC点估计与原双向绝对一致性公式相同，当前设计应解释为观察者内重复性，不能改用仅评估一致趋势的ICC(3,1)。复核前差异及2 mm内比例均指探测器中心坐标读数差异，不代表解剖配准误差。
