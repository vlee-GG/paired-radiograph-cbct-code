# Paired radiograph–CBCT code

Code for *Paired physically acquired five-view radiographs and cone-beam computed tomography of rabbits and mice*.

## Scripts

- [`cbct_cut_1.py`](cbct_cut_1.py): generate local CBCT crops from source volumes and coordinate records.
- [`correspondence_check.py`](correspondence_check.py): reproduce the automated radiograph–CBCT retrieval check.

## Quick start

Requires **Python 3.10 or later**.

```bash
git clone https://github.com/vlee-GG/paired-radiograph-cbct-code.git
cd paired-radiograph-cbct-code
pip install -r requirements.txt
```

Run the retrieval check using the dataset directory:

```bash
python correspondence_check.py --dataset ../physionet-paired-radiograph-cbct --out ./results
```

Generate crops using source NIfTI volumes and coordinate TXT files arranged in the acquisition-folder layout:

```bash
python cbct_cut_1.py --nii-root ../ct_nii --txt-root ../set --out ./crops
```

Use `python cbct_cut_1.py --help` or `python correspondence_check.py --help` for options. See [METHODS.md](METHODS.md) for input layout, crop-index mapping and algorithm details.

## Dataset and validation

The dataset contains **2,229 paired records**, five radiographic views per record and **335 source CBCT volumes**. Imaging data are distributed separately from this code repository.

PhysioNet publication and the dataset DOI are pending. The intended access procedure requires user identification and contributor review, subject to PhysioNet approval.

With the default settings, the retrieval check uses 19 distractors per record and yields **Top-1 85.4%**, **Top-3 95.6%** and **AUC 0.958**.

Image-quality review covers all five views; this automated retrieval check uses the **90° view**. It measures consistency of approximate anatomical correspondence, without calibrated pixel-to-voxel registration. CBCT intensities are not calibrated Hounsfield units. Historical review images contain invalid scale annotations; see [review display limitations](METHODS.md#blinded-review-materials-in-the-data-release).

## Licence and citation

Code: [MIT](LICENSE). Dataset: CC BY 4.0, licensed separately.

Please cite the accompanying data paper and this repository. Formal dataset and code-release citations will be added when their DOIs are available.
