# Methods and reproducibility

[Back to the repository](../README.md)

## Five-view retrieval

Run commands from the repository root. The main program is `five_view_retrieval.py`, which includes its numerical helpers.

```bash
python five_view_retrieval.py --dataset <release-root> --out <new-output> --workers 2
```

The output directory must be outside the dataset. Output includes `candidate_scores.npz` and `summary.csv`. Use `--workers` and `--fft-workers` to control CPU concurrency. Use a fresh output directory for each image snapshot; interrupted computation can resume only against the same frozen inputs.

All 2,229 records are evaluated separately at nominal angles 60°, 75°, 90°, 105° and 120°. Each uses its recorded CT assignment and 19 distractors from the same anatomical subset, with seed 20260920. All views use the identical candidate sets supplied in `validation/multiview_retrieval/candidate_sets.csv` in the data release. Candidate slot zero is the recorded assignment. No view fusion enters the reported results.

The CT descriptor averages the crop along array axis 2, retains `[10:110,20:140]`, transposes the image and applies a 1st-to-99.5th percentile window. Radiographs use the same percentile window. Both are resized to 139 rows by 114 columns and converted to normalized gradient magnitudes. Gradient NCC is maximized over eight rotation–magnification variants: 0°, ±3°, ±6° at unit magnification, and 0°, ±6° at 1.08 magnification. Horizontal displacement limits are ±20 working pixels; vertical limits are ±62, ±50, ±27, ±50 and ±62 pixels for the five angles.

The original 90° settings were developed on 255 records from a 511-record sample. The angle extension compared 13 configurations on 80 queries, 20 per anatomical subset, drawn from the 445 records in fold 0. The selected configuration is `angle_shift_axis1_factor1`. Candidates are drawn from the complete release, and full-release results include development records. Parameter comparisons and query identifiers accompany the data in `validation/multiview_retrieval/development.json`, which retains all original records and hashes.

Rank is one plus the number of distractors with a strictly higher score. Top-k is the fraction of queries with rank at most k. [Reference results](../reference/five_view_summary.csv) contain full-release counts and percentages. Complete scores have shape `(2229, 5, 20)`, with aligned `ids`, `angles` and `candidate_ids`; the dataset's `per_record.csv` contains 11,145 view-specific rows. These scores describe retrieval among the specified candidates. The displacement search is empirical and does not provide calibrated registration error.

`reference/baseline_90.npz` contains the 90° slice of the current 16 October 2026 analysis, for reference only. The detector recomputes every angle directly from the images; it does not load these reference scores. Input SHA-256 and sizes are checked against both metadata and the frozen release manifest. Output directories must be outside the dataset. The current run scored 222,900 comparisons; an independent angle-first loop checked 800 scores with a maximum difference of 2.23e-16.

## CBCT cropping

```bash
python cbct_cut_1.py --nii-root ../ct_nii --txt-root ../set --out ./crops
python cbct_cut_1.py --help
```

This program expects the acquisition-folder layout. Coordinate TXT records contain the reference voxel index, reference board coordinates in centimetres and photograph-derived target centres. `--match prefix3` (default) pairs directories by their first three digits; `--match dirname` pairs identical names recursively. `--prefix` restricts groups; `--force` recomputes existing outputs; `--fix-overflow` re-cuts oversized windows. Positioning images accompany the output crops.

At 0.25 mm isotropic spacing, the target crop is 120 × 160 × 140 voxels. In-plane centres are displaced from the reference by 40 voxels per centimetre; the second board-coordinate axis is reversed. Lower bounds are `int(i_center - 60)` and `int(j_center - 80)`; upper bounds add 120 and 160. Depth is `[160:300]` for both species. Out-of-bounds windows are clipped without padding.

In the released dataset, `crop_index_in_source` numbers retained crops and can differ from the original coordinate-target ordinal after exclusions. Use the crop origin and retained dimensions in `metadata/metadata.csv` to identify the corresponding target; see the data dictionary for the mapping.

## Data and review materials

Imaging data are supplied separately. PhysioNet publication and the dataset DOI are pending; the intended access procedure requires identification and contributor review, subject to platform approval. The current version 1.0.0 package contains 14,768 files (46.44 GiB), including 335 source CBCT volumes, 2,229 local crops and 11,145 radiographs. The supporting [release-tools archive](../archive/release_tools_20261016.zip) contains portable metadata/QC, coordinate-linkage, full crop-integrity, partition, checksum and historical normalization/review-replay code, with explicit curated inputs and tests. It rebuilds partition membership and replays the published fold allocation; it does not infer animal identities or an undocumented original randomization seed.

The coordinate-repeatability audit used 2,244 paired detector-centre readings from the same observer in two sessions approximately one year apart. First-session results were concealed during the second session. The same observer subsequently reviewed the 63 pairs exceeding 2 mm; 2,181 pairs (97.2%) were within 2 mm before review. Single-measurement, absolute-agreement ICCs use a two-way mixed-effects model. Current protocols and data dictionaries are in `validation/coordinate_review/03_protocol/`; historical field names such as `third_reader_id` denote the later review stage, not a separate observer. This audit is distinct from the two-practitioner anatomical review.

The current blinded study included 300 distinct pairs (250 original pairs and 50 substituted-CT controls), independently assessed by two clinicians. First-submission label agreement was 95.7% and 72.3%. One author subsequently reviewed all 91 cases with at least one label-discordant or uncertain answer, comprising 89 between-reader disagreements and two shared deviations. The author judged 82 original pairs approximately corresponding and nine controls non-corresponding, with reference labels and both readers' first-submission answers visible. The original blinded scores are retained.

The current display generator, per-case section indices, display windows and case-level assessments accompany the dataset in `validation/anatomical_review/`. Sections maximize the fraction of voxels above 300 rescaled units; display windows are 3–59,745 for radiographs and approximately −1,000–1,780 for CBCT. The current four-panel displays contain no physical scale bars or anatomical direction markers. The generator produces 320 distinct training/formal composites, including 20 training cases and 300 formal cases. The 15 hidden repeats per reader reuse formal-case images. CBCT values are not calibrated Hounsfield units.

Python scripts use the [MIT licence](../LICENSE). The dataset, `reference/baseline_90.npz` and `reference/five_view_summary.csv` use CC BY 4.0. Cite the accompanying data paper and code repository; formal DOI citations will be added when assigned.

## Final orientation provenance

The historical normalization mapping is preserved; final instances are listed in the dataset's `metadata/radiograph_provenance.csv`. Cumulative extra 180-degree rotations cover 133 records, 665 views and 47 source series, as documented in `validation/release_qc/PROVENANCE.md`. Round indices are original coordinate-target ordinals, which can differ from retained-crop ordinals. Distinct SOP instance UIDs and source-instance references were assigned to 115 later-reoriented views, with exact pixel-byte and decoded-array equality during this header update. The original scoring-input snapshot and score tensor are preserved; the scoring-input bridge under `supporting_checks.scoring_input_bridge.report` in `validation/release_qc/qc_summary.json` records exact equality of the scored pixels after header changes. Original review decisions and dates remain unchanged. The current `validation/anatomical_review/build_review_composites.py` generates the two-clinician study displays. Historical recovery/replay commands require additional inputs outside the current data package.

## Four validation evidence directories

The dataset organizes evidence under `validation/coordinate_review/`, `validation/release_qc/`, `validation/anatomical_review/` and `validation/multiview_retrieval/`. The 707 validation files include the two 335-file coordinate-reading snapshots. The release-root `SHA256SUMS.txt` is the single checksum list. The combined check summary is `validation/release_qc/qc_summary.json`; manual QC, crop/release reconciliation and exclusions are in the same directory. All ten fold-pair comparisons and the complete retrieval-verification report are embedded under `supporting_checks.partition_overlap_check` and `supporting_checks.retrieval_verification` in that summary. Only the partition-unit example split is distributed. The supporting release-tools archive reads the relocated manual QC table. Image/coordinate inputs, candidate scores and scientific parameters are unchanged. Two coordinate TXT snapshots are directly available in `validation/coordinate_review/01_coordinates/new_readings/` and `validation/coordinate_review/01_coordinates/target_final/` (335 files each). The latter preserves the initial review stage; current released coordinates remain under `source_cbct/`. Filter `validation/multiview_retrieval/per_record.csv` on `angle = 90` instead of reading the cancelled compatibility CSV; interrater statistics are in `validation/anatomical_review/review_stats.json`.
