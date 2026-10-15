# Methods and reproducibility

Detailed documentation for the accompanying cropping and retrieval scripts.
For installation and the shortest usage examples, see [README.md](README.md).

## What is not in this package

The scripts that construct the metadata, file-QC, partition and checksum tables are not
included; they read raw acquisition directories that are not part of the release. The data
record's Code Availability section states this. Everything those tables describe can be
checked against the release itself. `metadata/data_dictionary.csv` contains 228 definitions
for the core record, file and review tables; some auxiliary historical CSVs are not covered.

## Blinded-review materials in the data release

The data package contains 14,783 files, with 14,782 entries in its root checksum manifest,
11 metadata CSV tables and 722 validation files. Imaging and coordinate files number
14,044; documentation, metadata and validation account for the other 739 files.
These are dataset counts, separate from this code repository's `MANIFEST_SHA256.txt`.

The following files are included in the dataset, under `validation/anatomical_review/`:

| File | Purpose |
|---|---|
| `build_review_composites.py` | reconstruct the 300 composites used in the historical blinded review |
| `anat_choices.csv` | selected per-case section indices and full-precision display windows |

```bash
python ../physionet-paired-radiograph-cbct/validation/anatomical_review/build_review_composites.py --dataset ../physionet-paired-radiograph-cbct --out ./review_cases
```

Sections maximise the fraction of voxels above 300 rescaled units along the corresponding
array axis. The fixed windows are 3 to 59,745 for radiographs and -1000.0001220703 to
1780.0089111328 rescaled units for CBCT. These CBCT values are not calibrated Hounsfield units.

The generator preserves the images seen by the reviewers, including historical annotations
labelled 10 mm whose rendered lengths were calculated incorrectly from image spacing and
display scaling. Those annotations are not valid physical scale bars and must not be used
for distance measurements. The reported review assesses qualitative compatibility of local
anatomy. Any subsequent images with corrected scale bars must be identified as a revised
display version rather than presented as the original reviewed material. The published
review decisions and statistics are unchanged.

## `cbct_cut_1.py` — crop a local volume

Reads a reconstructed source CBCT volume and its coordinate record, works out the window
the radiographs cover, and writes the local crop plus a positioning-check image.

```bash
python cbct_cut_1.py --help

python cbct_cut_1.py --nii-root ../ct_nii --txt-root ../set --out ../ct_cut
python cbct_cut_1.py --nii-root ../ct_nii --txt-root ../set --out ../ct_cut --prefix 812,621
python cbct_cut_1.py --nii-root ../nii --txt-root ../set_all --out ../ct_cut --match dirname
python cbct_cut_1.py --out ../ct_cut --fix-overflow
```

| Option | Default | Meaning |
|---|---|---|
| `--nii-root` / `--txt-root` | `ct_nii` / `set` | roots holding the source volumes and the coordinate records |
| `--out` | `ct_cut` | output root |
| `--prefix` | all | comma-separated group prefixes to process |
| `--match` | `prefix3` | `prefix3` pairs directories whose names start with the same three digits; `dirname` pairs identical names, searched recursively |
| `--force` | off | recompute, ignoring existing outputs |
| `--fix-overflow` | off | re-cut only boxes that came out oversized; add `--no-vis` to skip redrawing |

Output goes to `<out>/<txt directory name>/<txt file stem>/<n>.nii.gz`, with check images
under `visualizations/`. Existing outputs are skipped unless `--force` is given.

### Geometry

The window reproduces the detector footprint: 30 x 40 mm on a 640 x 640 matrix spanning
16 x 16 cm, so the in-plane scale is 40 voxels per cm.

```
px_per_cm = 640 / 16 = 40
cube_w, cube_h = 3 cm * 40 = 120 ,  4 cm * 40 = 160
dx = x_target_cm - x_reference_cm
dy = y_reference_cm - y_target_cm          # the second axis is reversed
i_c = i_reference + dx * 40
j_c = j_reference + dy * 40
i0 = int(i_c - 60),  i1 = i0 + 120         # round the width, then derive the upper
j0 = int(j_c - 80),  j1 = j0 + 160         # bound, so the box cannot exceed 120 x 160
k0, k1 = 160, 300                          # 140 slices, both species
crop = volume[clip(i0,0,640):clip(i1,0,640),
              clip(j0,0,640):clip(j1,0,640), k0:k1]
```

Rounding both bounds independently produced 161-voxel boxes when a centre landed on a
floating-point boundary such as 243.99999999999997. Windows that run past the array edge
are clipped and are correspondingly smaller.

Each coordinate record holds, in order, the reference voxel index, the reference board
coordinate in centimetres, and the photograph-derived target centres. Note that
`crop_index_in_source` in `metadata.csv` is the ordinal of the **retained crop**, which
equals the ordinal of the target in the coordinate file only where no acquisition round was
dropped; `metadata/data_dictionary.csv` records the 11 volumes where the two differ and how
to recover the target ordinal.

## `correspondence_check.py` — automated retrieval check

Every released record is compared with parallel-beam projections of its own recorded crop
window and of 19 distractor crops drawn from the same subset, with every candidate scored
identically.

```bash
python correspondence_check.py --help

python correspondence_check.py --dataset ../physionet-paired-radiograph-cbct --out .
```

| Option | Default | Meaning |
|---|---|---|
| `--dataset` | `physionet-paired-radiograph-cbct` | released dataset directory |
| `--out` | `.` | output directory |
| `--distractors` | `19` | distractors per record |
| `--seed` | `20260920` | distractor sampling seed |

Writes `correspondence_check.csv`, `.npz` and `.txt`. Runtime is about 23 minutes over
2,229 records on one core; most of it is reading the crops and the radiographs.

### Method

| Stage | Detail |
|---|---|
| Radiograph | DICOM -> 1-99.5 percentile window -> 1/6 scale (139 x 114) -> gradient magnitude |
| Crop | NIfTI -> mean along the depth axis -> a 25 x 30 mm window, chosen so its field matches the 24.1 x 29.2 mm radiograph -> transpose -> eight (rotation, magnification) variants -> 1/6 scale -> gradient magnitude |
| Similarity | normalised cross-correlation by FFT, zero-padded to twice the image size, maximised over integer shifts inside a +/-27 x +/-20 frequency mask |
| Translation reach | +/-27 and +/-20 at 1/6 scale = +/-162 x +/-120 full-resolution pixels = +/-5.67 x +/-4.20 mm |
| Variants | `(rotation in degrees, magnification)`: `(0, 1.00)`, `(±3, 1.00)`, `(±6, 1.00)`, `(0, 1.08)`, `(±6, 1.08)` |
| Rank | 1 + the number of distractors scoring above the record's own crop; chance is 5% first, mean rank 10.5 |

Gradient magnitude is what makes the two modalities comparable: a radiograph records
transmission and the crop projection averages voxels, so their absolute intensities have no
common scale, while their edge structure does.

Note the rotation search spans only six degrees. The check is therefore **not** indifferent
to the 180-degree orientation corrections applied to 41 series: those records score lower
against their own crop than the rest, and scoring the released radiograph against its
pre-correction copy shows the released orientation is the one the check prefers, for 108 of
the 110 corrected records (median 0.34 against 0.11).

### Result on the released data

| | |
|---|---|
| Records | 2,229, covering all 335 source volumes |
| Ranked first | **85.4%** (1,904), 95% bootstrap CI 83.9-86.9 |
| Within the top three | **95.6%** |
| Mean rank | **1.52** against a chance mean rank of 10.5 |
| AUC | **0.958**, 95% CI 0.953-0.963, over 42,351 mismatched pairs |
| By subset | rabbit limb 89.4%, maxilla 78.9%, mandible 79.8%, mouse 60.3% |

This describes the release as it stands rather than estimating performance on unseen data:
the search variants were developed on 255 records from a 511-record sample and then retained for the release-level analysis.

## Scope

Both programs were written for the authors' working tree and are published as they were run;
retarget them with the options above.

The original correspondence_check.py covers the 90-degree view. The five-view extension below evaluates every nominal angle separately. It is an automated consistency measure:
it involves no reviewer decisions, provides no registration-error bound, and does not
confirm every pairing. It is not a substitute for the blinded anatomical-correspondence
review reported in the data record. The measure answers how *distinguishable* a crop is,
whereas anatomical correspondence requires *identity*. The cropping program yields
approximate local anatomical correspondence; it computes no projection matrix and does not
calibrate CBCT values to Hounsfield units.


## Five-view independent retrieval

Run `python five_view_retrieval.py --dataset <release-root> --out <new-output> --workers 2`. The output directory must be outside the dataset. `five_view_core.py` provides the original gradient NCC and display preprocessing. `baseline_90.npz` preserves the original full-precision 90° reference scores, with distractors mapped to the frozen candidate order. This analysis file follows the dataset's CC BY 4.0 licence; Python scripts follow MIT.

Each of 2,229 query records is evaluated independently at nominal angles 60°, 75°, 90°, 105° and 120°. Twenty candidates contain its recorded assignment and 19 distractors from the same anatomical subset, sampled with seed 20260920. All views use exactly the same candidates, supplied as `validation/multiview_retrieval/candidate_sets.csv`. Every candidate receives the same transformations and displacement search for a given view. No fusion is used in the manuscript results.

The CT reference is the mean of `crop[10:110,20:140,:]` along array axis 2, followed by transpose and a 1st-to-99.5th percentile window. Radiographs use the same percentile window. Images are resized to 139 rows by 114 columns and converted to normalized gradient magnitudes. Candidate transforms are rotations 0, ±3°, ±6° at unit magnification and 0, ±6° at 1.08 magnification. Gradient NCC is maximized over these variants and translations. Horizontal limits are ±20 working pixels; vertical limits are ±62, ±50, ±27, ±50 and ±62 pixels at the respective five angles. These limits implement an empirical angle-dependent displacement search, rather than a calibrated source trajectory or object-space scale.

The original 90° settings were developed on 255 records from a 511-record sample. For the four additional views, 13 reference/projection/search configurations were compared on 80 queries, 20 per anatomical subset, drawn from the 445 records in fold 0. The selected configuration is `angle_shift_axis1_factor1`. Candidate pools draw from the complete release. Full-release results include the development records. See `development_queries.json`, `development_model_protocol.json` and `development_model_selection.json` in the data release for all model comparisons and query identifiers. The separate finite-source experiment is not used for the manuscript results.

Rank equals one plus the number of distractors with strictly higher score. Top-k is the fraction of queries with rank at most k. Full-release Top-1 percentages are 34.1, 61.1, 85.4, 51.1 and 26.6; Top-3 percentages are 60.0, 85.2, 95.6, 78.2 and 51.5. `candidate_scores.npz` has shape (2229,5,20), aligned `ids`, `angles` and `candidate_ids`; slot zero is the recorded assignment. `summary.csv` includes anatomical subsets and `per_record.csv` contains all 11,145 view-specific ranks. These are descriptive retrieval outcomes against recorded assignments.

Reproduce Figure 6 with `python plot_five_view_retrieval.py --dataset <release-root> --out <figure-output>`. It reads the score tensor and metadata, producing vector PDF/SVG and PNG preview files. Panel a shows cumulative retrieval by rank; panel b shows Top-1 by anatomical subset.
