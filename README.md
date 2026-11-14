# Paired radiograph–CBCT: code and validation

Code and technical-validation evidence for 2,229 paired rabbit/mouse records: five physically acquired radiographs and one local CBCT crop per record. The imaging data and consolidated reuse metadata are plan to deposite separately in Zenodo. The pairs describe approximate local anatomical correspondence.

```text
scripts/       CBCT cropping and five-view retrieval programs
validation/
  release_qc/            Release checks and provenance
  coordinate_review/    Within-observer detector-centre repeatability
  multiview_retrieval/   Independent retrieval at five angles
  anatomical_review/    Blinded clinician review and author re-examination
```

## Imaging-data access

The version 1.0.0 imaging-data package is plan to deposite in zenodo. Its access policy is open download under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

The Zenodo deposit contains one archive, `paired-radiograph-cbct.zip` (49,125,707,075 bytes; approximately 49.13 GB). Extract it before running the programs and use the extracted `paired-radiograph-cbct/` directory as the dataset root. The archive contains 14,052 files: 11,145 DICOM radiographs, 2,229 local NIfTI crops, 335 source NIfTI volumes, 335 coordinate TXT files, seven metadata/documentation files and the root README.

Verify the downloaded ZIP against the MD5 checksum displayed on the Zenodo record. Verify individual imaging and coordinate files against the paths, byte sizes and SHA-256 values in `metadata/files.csv`. The current English dictionary is `metadata/data_dictionary.csv`. Keep linked observations together using `metadata/partition_units.csv`; `metadata/splits/` provides an example five-fold allocation.

Code is licensed under [MIT](LICENSE). Imaging data and [dataset-derived validation records](validation/DATA_LICENSE.txt) are licensed under CC BY 4.0.

## Run

Python 3.10+:

```bash
pip install -r scripts/requirements.txt
python scripts/five_view_retrieval.py --dataset /path/to/dataset --out /path/to/results
python scripts/cbct_cut_1.py --dataset /path/to/dataset --out /path/to/crops
```

The retrieval program uses the reorganized published dataset layout and fixed candidate sets in this repository. It recomputes image features and scores, using output outside the dataset. The cropping program accepts the published dataset through `--dataset`: it re-extracts the recorded voxel boxes as float32 NIfTI files with translated source affines, without overwriting files. These outputs preserve scaled source intensities and are not byte-identical copies of the quantized published crops. Without `--dataset`, the original acquisition-folder interface is retained. `scripts/METHODS.md` describes its coordinate and axis conventions.

The four validation directories contain protocols, inputs, per-record results and summaries. They replace validation evidence previously bundled in the imaging deposit. `release_qc/metadata_snapshot/` preserves the pre-consolidation tables for interpretation of frozen validation inputs. Current dataset metadata consolidate these fields into `metadata.csv` and `files.csv`. Scientific results are unchanged by relocation. Historical input fingerprints describe their recorded snapshots; current repository checksums are in `validation/MANIFEST_SHA256.txt`.

## Citation and fixed code release

Teng, H., Huang, Y., Gao, W. & Liu, W. *Paired physically acquired five-view radiographs and cone-beam computed tomography of rabbits and mice*. Version 1.0.0. Zenodo. https://doi.org/10.5281/zenodo.24156428

Cite the dataset alongside the code/validation version used. The fixed code and validation snapshot is [v1.1.1](https://github.com/vlee-GG/paired-radiograph-cbct-code/releases/tag/v1.1.1), commit `90a0b26ead05107e20ab1e0a487149d4af48bed5`. The dataset DOI identifies the imaging deposit; it is not a DOI for the code release. Current documentation updates on `main` do not change that fixed snapshot.

