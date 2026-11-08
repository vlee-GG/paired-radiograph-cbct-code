# Repeated coordinate reading and discrepancy-review data dictionary

## Reading design

The same observer performed the two reading sessions approximately one year apart. First-session results were concealed during the second session. The same observer subsequently reviewed the 63 pairs with detector-centre differences greater than 2 mm. Historical third-reader field names identify this review stage, not a separate observer. The frozen two-session readings, before any review-driven amendments, are the inputs to the within-observer repeatability analysis.

## Source and preservation

- `04_analysis/record_comparison_pre_review.csv` contains the frozen paired readings.
- `04_analysis/review_records.csv` contains the initial 63-case review.
- `02_review_records/third_review_records.json` contains the subsequent structured review export.
- Website working copies and the append-only event log are historical internal materials, outside the publication package.

## Record identity and readings

- `id`: stable `root/folder/image` identifier.
- `existing`: frozen first-session reading in centimetres.
- `new`: frozen second-session reading from the same observer, in centimetres.
- `third`: subsequent review measurement by the same observer, in centimetres. The historical field name denotes the additional reading, not a third person.
- `pre_review_centre_difference_mm`: Euclidean difference between existing and new readings.
- `third_reader_id`, `third_read_at`: provenance of the subsequent review; historical schema names are preserved.
- A submitted third reading may be revised. Each revision is appended to the audit log with the previous and new values. If a final decision already exists, revising the third reading invalidates that decision, restores the working target to the frozen existing coordinate and requires a new adjudication.

## Final decisions

- `retain_existing`: third reading supports the existing coordinates; working target is restored to the frozen existing value.
- `correct_to_new`: third reading supports the new coordinates; working target is changed to the frozen new value.
- `coordinate_uncertain`: the three readings remain materially discordant; no coordinate is forced and `final_coordinate` remains null.
- `decision_reason`, `decided_by`, `decided_at`: mandatory adjudication provenance. `decision_reason` is constrained to `existing_reading_incorrect` (Existing reading is incorrect) or `new_reading_incorrect` (New reading is incorrect) so that reasons can be counted consistently.

## Downstream validation after coordinate correction

For `correct_to_new`, all four Boolean fields must be completed separately:

- `crop_regenerated`
- `file_qc_passed`
- `crop_reextraction_validated`
- `automatic_retrieval_validated`

## Reporting rules

- Pre-review analysis uses the frozen paired-reading CSV and includes all 335 matched directories and 2244 matched records.
- Post-review analysis uses `working` after adjudication.
- `coordinate_uncertain` records must be reported explicitly and must not be altered solely to make centre difference ≤2 mm.
- The public export must retain both original readings, the third reading, distances from the third reading to each original reading, final decision, reason, final coordinate and downstream validation fields.

## Version scope and preserved workflow flags (2026-10-14)

`02_review_records/third_review_records.json` currently contains 46 retain-existing,
14 correct-to-new and 3 uncertain decisions after two reversions. Historical
`04_analysis/review_records.csv`, `04_analysis/analysis_summary.json` and
`01_coordinates/target_final/` reflect the initial 44/16/3 stage. Published coordinates
are the source-coordinate TXT files under the root `source_cbct/`, and published crop
origins are in `metadata/metadata.csv`. Current coverage is in `validation/release_qc/qc_summary.json`.

`final_coordinate=null` in the three uncertain adjudication records means that no
confirmed final coordinate is asserted by that table. The historical target export used
same-observer review values as stated in its generation manifest; this does not resolve uncertainty.
The original downstream flags retain their workflow-record meanings. A false flag is not
retroactively replaced by the later release-wide crop-origin check.

`_release_context` is documentation metadata added on 2026-10-14, not a new observation,
reading, adjudication or analysis result. It records the scope and links to current evidence.

## Auxiliary analysis CSV dictionaries

### record_comparison_pre_review.csv

The table has 2,244 frozen paired readings in a 13-column schema. The 63 triggered rows can be obtained using `centre_difference_mm > 2` with the numerical threshold tolerance described in the analysis. Numeric readings and
differences are populated for every row; no numeric missing-value code is used.

| Field | Type / unit | Definition |
|---|---|---|
| `id` | string | Historical photograph/record identity, not a released sample_id; links the reading and adjudication evidence. |
| `directory` | string | Historical acquisition-series directory identifying the source reading set. |
| `image` | string | Acquisition photograph filename used for the coordinate reading; the photograph is not included in the release. |
| `index` | integer / ordinal | 1-based row index within that coordinate reading set; not automatically the retained crop ordinal. |
| `existing_bottom_cm` | float / cm | Frozen existing distance to the grid bottom edge. |
| `existing_right_cm` | float / cm | Frozen existing distance to the grid right edge. |
| `new_bottom_cm` | float / cm | Second-session reading by the same observer of distance to the grid bottom edge. |
| `new_right_cm` | float / cm | Second-session reading by the same observer of distance to the grid right edge. |
| `delta_x_mm` | float / mm | 10 × (new_right_cm − existing_right_cm); sign follows the right-edge distance convention. |
| `delta_y_mm` | float / mm | 10 × (new_bottom_cm − existing_bottom_cm); sign follows the bottom-edge distance convention. |
| `abs_delta_x_mm` | nonnegative float / mm | Absolute value of delta_x_mm. |
| `abs_delta_y_mm` | nonnegative float / mm | Absolute value of delta_y_mm. |
| `centre_difference_mm` | nonnegative float / mm | sqrt(delta_x_mm² + delta_y_mm²); selects the triggered table using a strict >2 mm threshold. |

## Released evidence

The two 335-file coordinate-reading snapshots are in `validation/coordinate_review/01_coordinates/`. The repeatability analysis uses `04_analysis/record_comparison_pre_review.csv`; the 63-case initial decisions are in `04_analysis/review_records.csv`, and the subsequent workflow export is in `02_review_records/third_review_records.json`. Current published coordinates remain in `source_cbct/`; combined release checks are in `validation/release_qc/qc_summary.json` (release-root paths). Historical website working copies and event logs are not distributed.


## Current ICC interpretation

ICC(A,1) is reported using a two-way mixed-effects model, single measurements and absolute agreement. For n targets and k=2 sessions, the estimate is `(MS_target - MS_error) / (MS_target + (k-1)*MS_error + k*(MS_session - MS_error)/n)`. It is algebraically the same point-estimate formula as the previously labelled ICC(2,1); the study-design interpretation is within-observer repeatability. It is not the consistency-only ICC(3,1). Frozen historical outputs retain their original labels and numerical values.
