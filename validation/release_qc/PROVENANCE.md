# Image-processing provenance

## UID and header normalization

The 2026-10-13 normalization creates release-specific identifiers for all 11,145
orientation-standardised radiographs. It replaces 9,425 nonconformant instance/series
identifiers and also gives the remaining derived exports distinct release identifiers.
UIDs use the DICOM UUID-derived 2.25 root. A mapping from each old instance and original
series to its new identifier preserves grouping; StudyInstanceUID is unchanged.
The exports contain one distinct SeriesInstanceUID per radiograph. That existing grouping
is preserved. No series were merged, and no record, view or source-volume link changed.
The pre-normalization files contain no referenced-UID sequences.

The historical table `metadata_snapshot/dicom_uid_mapping.csv` contains one row per radiograph, identified by its path,
sample, source volume and nominal view. Old UID strings are historical text and may not
be valid DICOM UIDs. Old/new whole-file hashes and sizes support provenance. Stored-pixel
and decoded-array hashes are equal before and after normalization. The per-row JSON
header repair list records every corrected or omitted optional malformed field.

Known exporter VR errors were corrected while retaining valid values: VOILUTFunction
LINEAR is encoded as CS, empty PatientTelephoneNumbers as SH, and the empty scheduled
protocol-code field as SQ. Malformed optional ContextIdentifier, ContextUID and
PerformedProtocolType values are omitted where they cannot be interpreted as legal
values; their exact original VR and textual content remain in the mapping and backups.
New file-meta implementation identifiers describe the pydicom 3.0.2 writer, and
DerivationDescription records the normalization and points to the mapping.

All image pixels, window settings, matrix dimensions, recorded spacing, acquisition
settings, Study UIDs, transfer syntax and other non-repaired attributes are unchanged.
The original pre-normalization published files were backed up outside this release;
these backups are distinct from raw acquisition exports. See qc_summary.json#/supporting_checks/dicom_uid_normalization
for validation totals and scope. These checks do not constitute a complete IOD or PACS
certification and do not calibrate image geometry or resolve PixelSpacing limitations.

### Historical snapshot and final instances (2026-10-16)

This report and dicom_uid_mapping.csv describe the 2026-10-13 event. Their original hashes and counts are preserved. Later reorientation affects 115 instances; final file hashes and SOP UIDs are retained in the separate dataset's metadata/files.csv (sha256 and provenance_ fields); the original table is in metadata_snapshot/radiograph_provenance.csv. See the final-orientation section below for the cumulative 133-record / 665-view / 47-series status and exact snapshot recovery.

### Retained acquisition and specimen fields

The released headers retain the acquisition-system StudyDate, SeriesDate and StudyTime, and the internal specimen name in PatientName and PatientID. These are non-human identifiers that encode the same acquisition day already published as acquisition_day_id.

## Final orientation

The 2026-10-13 normalization mapping is an immutable historical event snapshot. It is not a table of final hashes after subsequent reorientation. `metadata_snapshot/radiograph_provenance.csv` records all 11,145 final instances and links them to that snapshot. Its fields are also consolidated with the `provenance_` prefix in the separate dataset's `metadata/files.csv`.

The original manual QC records 110 paired records (550 views) with an extra 180-degree rotation relative to the author-reported base transform. Comparison with the normalized pixel hashes identifies 23 additional mouse records (115 views): applying the inverse 180-degree rotation reproduces the historical stored-pixel SHA-256 exactly. The cumulative final status is **133 records (6.0%), 665 views, in 47 source series**. The original operation time and actor are unknown and remain blank in `orientation_events.csv`; 2026-10-16 is the documentation date.

Original `orientation_rot180_rounds` uses the coordinate-target ordinal, which can differ from retained-crop order after exclusions. The final fields in `manual_qc.csv` explicitly use coordinate-target ordinals. Original columns, review dates and review judgements are preserved. These records describe applied transformations and their evidence; they do not independently certify anatomical laterality or direction.

The 115 later-reoriented images receive distinct SOP instance UIDs with matching file-meta UIDs, a SourceImageSequence referencing the historical normalized instance, and a DerivationDescription. This header update preserves all PixelData bytes and decoded arrays. Study and series grouping are unchanged. No geometric axes were invented. Current file sizes/hashes are indexed in the separate dataset's metadata/files.csv and SHA256SUMS.txt. The pre-consolidation tables are retained in metadata_snapshot/ for historical interpretation.

## Scoring inputs and final released files

The five-view scoring snapshot and candidate_scores.npz are preserved. `qc_summary.json#/supporting_checks/scoring_input_bridge/report` connects the scored inputs to the final released headers and records pixel-array equality.

Portable processing commands, curated input annotations and verification records are available in the [versioned code archive](https://github.com/vlee-GG/paired-radiograph-cbct-code/blob/5a87fd6ab1e8ee4e65fd88b1ed94b7111175570d/archive/release_tools_20261016.zip). The original acquisition exports are not included; the base transformation remains author-reported.

### Current blinded study

The current 250-original/50-control study uses the final released orientations and no physical scale annotations. Its inputs, section choices and renderer are in `../anatomical_review/`.

## Historical document names

The former `DICOM_UID_NORMALIZATION.md` and `ORIENTATION_PROVENANCE.md` documents are consolidated here. Historical paths embedded in unchanged DICOM derivation fields identify the processing event; the released image files are not rewritten to update document names.
