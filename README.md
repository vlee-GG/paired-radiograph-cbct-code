# Paired radiograph–CBCT code

Two programs for the rabbit and mouse dataset: **CBCT cropping** and **five-view retrieval**.

## Install

Python **3.10+**:

```bash
pip install -r requirements.txt
```

## Run

**Five-view retrieval** — all 2,229 records, evaluated separately at five angles:

```bash
python five_view_retrieval.py --dataset ../physionet-paired-radiograph-cbct --out ./results
```

**CBCT cropping** — source volumes and coordinate TXT files in acquisition-folder layout:

```bash
python cbct_cut_1.py --nii-root ../ct_nii --txt-root ../set --out ./crops
```

Use `--help` for options. [Methods and input details](docs/METHODS.md) · [Reference results](reference/five_view_summary.csv).

Imaging data are distributed separately; the PhysioNet deposit is pending. Full scores and ranks accompany the dataset under `validation/multiview_retrieval/`.

Code: [MIT](LICENSE). Dataset and reference analysis files: CC BY 4.0. Please cite the accompanying data paper and this repository. Checksums: [reference/MANIFEST_SHA256.txt](reference/MANIFEST_SHA256.txt), with paths relative to the repository root.
