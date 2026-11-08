# Paired radiograph–CBCT: code and validation

Code and technical-validation evidence for 2,229 paired rabbit/mouse records. Imaging data are deposited separately on PhysioNet (project under preparation).

```text
scripts/       CBCT cropping and five-view retrieval programs
validation/
  release_qc/            Release checks and provenance
  coordinate_review/    Within-observer detector-centre repeatability
  multiview_retrieval/   Independent retrieval at five angles
  anatomical_review/    Blinded clinician review and author re-examination
```

## Run

Python 3.10+:

```bash
pip install -r scripts/requirements.txt
python scripts/five_view_retrieval.py --dataset /path/to/dataset --out /path/to/results
python scripts/cbct_cut_1.py --dataset /path/to/dataset --out /path/to/crops
```

The retrieval program uses the reorganized published dataset layout and fixed candidate sets in this repository. It recomputes image features and scores, using output outside the dataset. The cropping program accepts the published dataset through `--dataset`: it re-extracts the recorded voxel boxes as float32 NIfTI files with translated source affines, without overwriting files. These outputs preserve scaled source intensities and are not byte-identical copies of the quantized published crops. Without `--dataset`, the original acquisition-folder interface is retained. `scripts/METHODS.md` describes its coordinate and axis conventions.

The four validation directories contain protocols, inputs, per-record results and summaries. They replace validation evidence previously bundled in the imaging deposit. `release_qc/metadata_snapshot/` preserves the pre-consolidation tables for interpretation of frozen validation inputs. Current dataset metadata consolidate these fields into `metadata.csv` and `files.csv`. Scientific results are unchanged by relocation. Historical input fingerprints describe their recorded snapshots; current repository checksums are in `validation/MANIFEST_SHA256.txt`.

Code license: MIT (`LICENSE`). Dataset-derived validation records: CC BY 4.0 (`validation/DATA_LICENSE.txt`). Cite the dataset and the fixed code/validation release used. No archival DOI has yet been assigned.
