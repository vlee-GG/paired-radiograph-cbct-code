# Full-release quality and reproducibility

| File | Evidence |
|---|---|
| `manual_qc.csv` | Author review of all five views in 335 source-volume series; all accepted |
| `release_validation.csv` | Source-volume retention and the reconciliation of 2,244 targets to 2,229 released records |
| `exclusions.csv` | Final net exclusions, distinct from the historical action log |
| `qc_summary.json` | Crop/coordinate reproduction, file checks and the four validation steps; `supporting_checks` preserves eight merged reports, including the complete retrieval-verification JSON and all ten fold-pair comparisons with their original SHA-256 values |
| `figure_provenance.csv` | Figure record IDs, display settings and the stated verification status |
| `orientation_events.csv` | Per-view transformation provenance |
| `PROVENANCE.md` | UID/header normalization, final orientation and the relation between scored inputs and current files |

All 14,044 imaging and coordinate-file rows are indexed in the separate dataset's `metadata/files.csv`, including consolidated `qc_` fields. The original tables are preserved in `metadata_snapshot/`. Imaging and source-coordinate bytes were not modified by this reorganization. Crop checks evaluate extraction from arrays; coordinate-reading differences evaluate detector-centre localization.

Paths embedded in unchanged DICOM DerivationDescription fields describe the earlier release layout. In particular, the historical `validation/orientation_events.csv` now resolves to `validation/release_qc/orientation_events.csv`; the files themselves are not reserialized to edit that historical text.

In `qc_summary.json`, `supporting_checks.retrieval_verification.report` preserves all fields of the retrieval check. `supporting_checks.partition_overlap_check.rows` contains the ten fold-pair comparisons; `original_csv_utf8` retains the exact CSV text. `original_path` identifies the source before consolidation, not a separate downloadable file.

`supporting_checks.scoring_input_bridge.report` in `qc_summary.json` preserves all scoring-input bridge fields and its original SHA-256. The separate bridge file is no longer distributed.


## Data locations

Images and current consolidated metadata are in the separate PhysioNet dataset. Paths beginning `metadata/` or `source_cbct/` refer to that dataset root. Validation paths refer to this code repository. Pre-consolidation metadata are preserved in `../release_qc/metadata_snapshot/` for interpreting historical input fingerprints. Current file QC, UID normalization and radiograph provenance fields are consolidated into dataset `metadata/files.csv` (`qc_`, `uid_`, `provenance_`); group fields and complete group-to-animal links are retained in dataset `metadata/metadata.csv`.
