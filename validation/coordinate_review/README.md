# Detector-centre coordinate readings

The same observer read the coordinates twice, approximately one year apart, with the first-session results concealed during the second session. Each paired reading records the detector centre relative to the board lower and right boundaries (cm). Differences are two-dimensional Euclidean distances converted to mm. The paired-reading audit contains 2,244 paired readings across 335 source volumes: median difference 0.50 mm; 2,181 pairs (97.2%) within 2 mm; 63 pairs subsequently re-examined by the same observer.

| File | Purpose |
|---|---|
| `04_analysis/record_comparison_pre_review.csv` | Complete paired readings used for agreement analysis |
| `04_analysis/analysis_summary.json` | Agreement statistics and the initial adjudication analysis |
| `04_analysis/review_records.csv` | Initial 63-case discrepancy review: 44 retain-existing, 16 proposed corrections, 3 uncertain |
| `02_review_records/third_review_records.json` | Subsequent workflow export: 46 retain-existing, 14 correct-to-new, 3 uncertain |
| `03_protocol/REVIEW_PROTOCOL.md` | Repeated-reading and discrepancy-review procedure |
| `03_protocol/DATA_DICTIONARY.md` | Fields, units and interpretation |

The initial and subsequent decision tables refer to different stages. Of 16 initial proposed corrections, 15 concern released records: 13 applied and 2 reverted; one is unreleased. Uncertain final coordinates remain unconfirmed. Published coordinates are in `source_cbct/`; crop positions are in `metadata/metadata.csv` (both release-root paths). See `../release_qc/qc_summary.json` for current crop verification.

Two sets of original TXT snapshots are available directly in this directory:

| Directory | Files | Role |
|---|---:|---|
| `01_coordinates/new_readings/` | 335 | Preserved new-reading TXT snapshot |
| `01_coordinates/target_final/` | 335 | Initial adjudication snapshot (44/16/3 stage) |

File names, units, line order and comments are unchanged. The agreement analysis uses the frozen paired readings in `04_analysis/record_comparison_pre_review.csv`; subsequent edits to some working readings are documented in the preserved analysis. `target_final` is the initial snapshot; current released coordinates remain in `source_cbct/`.

Checksums for this code and validation repository are in `../MANIFEST_SHA256.txt`, with paths relative to the repository root. The separate imaging dataset has its own root `SHA256SUMS.txt`. The repeated-reading design assesses within-observer repeatability. The subsequent review was performed by the same observer. Historical names such as `third_reader_id` and `third_review_records.json` are retained for schema compatibility and identify the later review stage, not a third person. The current protocol states the within-observer study design.

Absolute-agreement ICCs use a two-way mixed-effects, single-measurement model (ICC(A,1)): 0.9896 for the right-edge coordinate and 0.9997 for the lower-edge coordinate. The point-estimate formula is unchanged from the earlier two-way absolute-agreement calculation. All repeatability statistics use the two frozen sessions before review; the 97.2% is a proportion of paired detector-centre readings, not anatomical registration accuracy.


## Data locations

Images and current consolidated metadata are in the separate PhysioNet dataset. Paths beginning `metadata/` or `source_cbct/` refer to that dataset root. Validation paths refer to this code repository. Pre-consolidation metadata are preserved in `../release_qc/metadata_snapshot/` for interpreting historical input fingerprints. Current file QC, UID normalization and radiograph provenance fields are consolidated into dataset `metadata/files.csv` (`qc_`, `uid_`, `provenance_`); group fields and complete group-to-animal links are retained in dataset `metadata/metadata.csv`.
