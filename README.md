# Paired radiograph–CBCT code

Code for *Paired physically acquired five-view radiographs and cone-beam computed tomography of rabbits and mice*.

## Scripts

- [`cbct_cut_1.py`](cbct_cut_1.py): generate local CBCT crops from source volumes and coordinate records.
- [`five_view_retrieval.py`](five_view_retrieval.py): reproduce full-release retrieval separately for all five views.
- [`plot_five_view_retrieval.py`](plot_five_view_retrieval.py): reproduce Figure 6.
- [`correspondence_check.py`](correspondence_check.py): original 90° reference analysis.

## Quick start

Requires **Python 3.10 or later**.

```bash
git clone https://github.com/vlee-GG/paired-radiograph-cbct-code.git
cd paired-radiograph-cbct-code
pip install -r requirements.txt
```

Run the retrieval check using the dataset directory:

```bash
python five_view_retrieval.py --dataset ../physionet-paired-radiograph-cbct --out ./results
```

Generate crops using source NIfTI volumes and coordinate TXT files arranged in the acquisition-folder layout:

```bash
python cbct_cut_1.py --nii-root ../ct_nii --txt-root ../set --out ./crops
```

Use `python cbct_cut_1.py --help` or `python correspondence_check.py --help` for options. See [METHODS.md](METHODS.md) for input layout, crop-index mapping and algorithm details.

## Dataset and validation

The dataset contains **2,229 paired records**, five radiographic views per record and **335 source CBCT volumes**. Imaging data are distributed separately from this code repository.

PhysioNet publication and the dataset DOI are pending. The intended access procedure requires user identification and contributor review, subject to PhysioNet approval.

All 2,229 records are evaluated separately at 60°, 75°, 90°, 105° and 120°, using the same 20 CT candidates per record. Full-release Top-1 rates are **34.1%, 61.1%, 85.4%, 51.1% and 26.6%**, respectively. See [five_view_summary.csv](five_view_summary.csv) and [METHODS.md](METHODS.md).

Full scores, candidate assignments and record-level ranks are distributed with the dataset under `validation/multiview_retrieval/`. Parameter selection and the approximate-correspondence scope are documented in METHODS.md.

## Licence and citation

Code: [MIT](LICENSE). Dataset: CC BY 4.0, licensed separately.

Please cite the accompanying data paper and this repository. Formal dataset and code-release citations will be added when their DOIs are available.
