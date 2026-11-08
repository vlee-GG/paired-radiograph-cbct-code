# Technical validation

Four directories provide the protocols, record-level evidence and summaries used in the Data Descriptor.

| Step | Directory | Main evidence |
|---|---|---|
| 1. Detector-centre coordinate readings | `coordinate_review/` | 2,244 paired readings, agreement statistics and adjudication |
| 2. Full-release quality and reproducibility | `release_qc/` | Five-view series review, crop/coordinate reproduction, exclusions, file integrity and partition checks |
| 3. Blinded local anatomical correspondence | `anatomical_review/` | 250 original pairs and 50 substituted controls; two clinicians' first submissions, display choices and statistics |
| 4. Five-view independent retrieval | `multiview_retrieval/` | Fixed candidates, 222,900 scores, ranks, summaries and Figure 6 source data |

The anatomical-review directory now supports the two-clinician 250-original/50-control study described in the manuscript. First-submission scores and the completed case-level author follow-up of all 91 cases with at least one label-discordant or uncertain answer (`anatomical_review/author_review_records.csv`) are reported separately; no AI-assisted opinions are included.

The two original coordinate TXT snapshots (335 files each) are directly available under `coordinate_review/01_coordinates/new_readings/` and `coordinate_review/01_coordinates/target_final/`. The 90-degree-only compatibility CSV is replaced by filtering `multiview_retrieval/per_record.csv` on `angle = 90`; the separate interrater summary is included in `anatomical_review/review_stats.json`. Current published coordinates are the files in `source_cbct/`, not the archived `target_final` snapshot.

See `release_qc/qc_summary.json` for the combined quality-check summary and `MANIFEST_SHA256.txt` (paths relative to the repository root) for current checksums. The frozen supporting tools and curated inputs are retained in `release_qc/release_tools_snapshot.zip`; current published-layout programs are in `../scripts/`.
