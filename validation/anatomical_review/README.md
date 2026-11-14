# Blinded assessment of local anatomical correspondence

The current manuscript analysis uses two clinicians' first confirmed submissions for **300 distinct cases: 250 original pairs and 50 substituted-CT controls**. Reader 1 is an orthopaedic clinician with three years of clinical experience at a county-level hospital; reader 2 practises in a dental clinic.

| Result | Reader 1 | Reader 2 |
|---|---:|---:|
| Agreement with constructed reference labels | 287/300 (95.7%) | 217/300 (72.3%) |
| Original pairs rated yes | 246/250 (98.4%) | 168/250 (67.2%) |
| Substituted controls rated no | 41/50 (82.0%) | 49/50 (98.0%) |
| Uncertain answers | 10/300 | 58/300 |
| Agreement on 15 hidden repeat pairs | 15/15 | 11/15 |

Between-reader agreement was 211/300 (70.3%); 89 cases had different reader answers. Subsequent author review covered all 91 cases with at least one label-discordant or uncertain answer: these 89 cases plus two shared deviations. All 82 original pairs were judged approximately corresponding, and all nine substituted-CT controls were judged non-corresponding. The author assessments agreed with the recorded reference labels in all 91 cases and are documented separately from first-submission scores. Similar bony contours or trajectories were recorded for 75 original pairs and compatible joint morphology or relative anatomical relationships for 34, with both selected in 27; these criteria are non-exclusive. Of the two shared deviations, c0065 was a substituted control rated yes by both readers and judged non-corresponding at follow-up; c0280 was an original pair rated no by both readers and judged approximately corresponding.

## Files

- `PROTOCOL.md`: sampling, display, reader task and analysis definitions.
- `review_records.csv`: 300 original/control mappings, both first-submission answers and the 91-case author-review scope.
- `presentation_records.csv`: frozen presentation order and first-submission answers, including training and repeats.
- `review_stats.json`: counts, agreement and cluster-bootstrap intervals.
- `display_choices.csv`, `input_hashes.csv`, `build_review_composites.py`: fixed display inputs and rendering.
- `analyze_review.py`: recompute counts, agreement and intervals from released tables.
- `author_review_records.csv`: all 91 case-level author conclusions and selected structural criteria.
- Release-root `SHA256SUMS.txt`: checksums for these materials and the rest of the release. Follow-up methods and record interpretation are in `PROTOCOL.md`.

No AI-assisted review opinions are included in the current manuscript or these current assessment materials.


## Data locations

Images and current consolidated metadata are in the separate [Zenodo dataset](https://doi.org/10.5281/zenodo.24156428) (version 1.0.0; CC BY 4.0; open-download access policy). Extract `paired-radiograph-cbct.zip` and use its `paired-radiograph-cbct/` directory as the dataset root. Paths beginning `metadata/` or `source_cbct/` refer to that dataset root. Validation paths refer to this code repository. Pre-consolidation metadata are preserved in `../release_qc/metadata_snapshot/` for interpreting historical input fingerprints. Current file QC, UID normalization and radiograph provenance fields are consolidated into dataset `metadata/files.csv` (`qc_`, `uid_`, `provenance_`); group fields and complete group-to-animal links are retained in dataset `metadata/metadata.csv`.
