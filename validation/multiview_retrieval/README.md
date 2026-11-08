# Five-view independent retrieval

Current analysis: **16 October 2026**. All 2,229 records were scored afresh at 60°, 75°, 90°, 105° and 120°, against the same 20 candidates per record (recorded crop plus 19 within-subset distractors). All 222,900 scores come from the current images, including 90°. No view fusion is used.

- `summary.csv`: counts, percentages, mean rank and score-discrimination AUC for all records and each anatomical subset. `angle` is the nominal angle in degrees. Top-1 intervals use 4,000 partition-unit cluster bootstrap resamples, seed 20261016.
- `per_record.csv`: 11,145 view-specific scores and ranks; rank is one plus the number of strictly higher-scoring distractors.
- `candidate_scores.npz`: `ids`, `angles`, `candidate_ids`, and `scores` with shape (2229, 5, 20). Candidate slot zero is the recorded assignment.
- `candidate_sets.csv`: fixed candidate assignments; `protocol.json`: processing and search settings. `development.json` preserves all three original parameter-development records, including their original UTF-8 text and SHA-256 values.
- `Figure_6_source_data.csv`: cumulative retrieval and subset Top-1 source data.
- `input_file_hashes.csv`, `all_views_pixel_qc.csv` and `../release_qc/qc_summary.json` document input integrity and numerical checks. The complete numerical verification is under `supporting_checks.retrieval_verification.report` in the QC summary. The frozen manifest fingerprint refers to the input snapshot before documentation synchronization. The original scoring snapshot is preserved. A later documentation update assigns distinct SOP instance UIDs and provenance headers to 115 already-reoriented DICOM files, with exact PixelData and decoded-array equality. `supporting_checks.scoring_input_bridge.report` in the QC summary connects the unchanged score tensor to current files; `metadata/radiograph_provenance.csv` (relative to the release root) supplies final whole-file hashes. Frozen input-file hashes and timestamps describe the scoring snapshot, not the updated headers.

Only this full-release analysis is distributed. The configuration retains the original 255-record 90° development and the 80-query angle-extension selection. Full-release statistics include development records. Scripts: `cbct_cut_1.py` and `five_view_retrieval.py` in the accompanying code repository. The detector reads the release images and fixed candidate sets; it does not import stored reference scores.

Original input-file fingerprints and scientific results are unchanged; the complete bridge report is in the QC summary. The separate 90-degree compatibility export is replaced by filtering `per_record.csv` on `angle = 90`.


## Data locations

Images and current consolidated metadata are in the separate PhysioNet dataset. Paths beginning `metadata/` or `source_cbct/` refer to that dataset root. Validation paths refer to this code repository. Pre-consolidation metadata are preserved in `../release_qc/metadata_snapshot/` for interpreting historical input fingerprints. Current file QC, UID normalization and radiograph provenance fields are consolidated into dataset `metadata/files.csv` (`qc_`, `uid_`, `provenance_`); group fields and complete group-to-animal links are retained in dataset `metadata/metadata.csv`.
