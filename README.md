# Code for *Paired physically acquired five-view radiographs and cone-beam computed tomography of rabbits and mice*

Two programs: the local cropping program that produced the released crops, and the automated
retrieval check that produced `validation/correspondence_check.csv`.

> Dataset: `physionet-paired-radiograph-cbct` version 1.0.0, CC BY 4.0.
> 14,762 files, 46.41 GiB, 2,229 paired records over 335 source CBCT volumes.
> DOI: *to be inserted once the deposition is issued.*

Code licence: MIT, see [LICENSE](LICENSE). The dataset is released separately under CC BY 4.0.

## Contents

| File | |
|---|---|
| `cbct_cut_1.py` | crop a local volume from a source CBCT volume and its coordinate record |
| `correspondence_check.py` | automated retrieval check over all released records |
| `requirements.txt` | dependencies |
| `LICENSE` | MIT |
| `MANIFEST_SHA256.txt` | checksums for the files above |

Both programs run against the released dataset. Neither needs anything that is not in the
release, and neither writes into it unless you point it there.

## What is not in this package

The scripts that construct the metadata, file-QC, partition and checksum tables are not
included; they read raw acquisition directories that are not part of the release. The data
record's Code Availability section states this. Everything those tables describe can be
checked against the release itself, and `metadata/data_dictionary.csv` defines every field.

## Install

Python 3.9 or later.

```bash
pip install -r requirements.txt
python -c "import numpy, nibabel, pydicom, PIL, matplotlib; print('ok')"
```

`matplotlib` is used only by `cbct_cut_1.py`, for its positioning-check images.

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
the search variants and the shift window are fixed constants, so nothing is fitted to these
records.

## Scope

Both programs were written for the authors' working tree and are published as they were run;
retarget them with the options above.

The retrieval check covers the 90-degree view only. It is an automated consistency measure:
it involves no reviewer decisions, provides no registration-error bound, and does not
confirm every pairing. It is not a substitute for the blinded anatomical-correspondence
review reported in the data record. The measure answers how *distinguishable* a crop is,
whereas anatomical correspondence requires *identity*. The cropping program yields
approximate local anatomical correspondence; it computes no projection matrix and does not
calibrate CBCT values to Hounsfield units.

## Licence

Code: **MIT**, see [LICENSE](LICENSE). Copyright (c) 2026  Hong Teng.

Data: **CC BY 4.0**, as stated in the data record. The two are licensed separately; the MIT
licence covers the scripts in this repository, not the dataset.

## Citing

Please cite the data record for the dataset and this repository for the code. A
`CITATION.cff` and a Zenodo DOI for the archived release will be added once the dataset
deposition is issued.
