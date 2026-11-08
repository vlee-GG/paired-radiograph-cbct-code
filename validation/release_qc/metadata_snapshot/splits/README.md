# Example partitions

One example five-fold partition is provided in `example_5fold_by_partition_unit.csv`.

## Partition units and independence

Folds are assigned to **partition units**, not to individual data groups. A partition unit is a
connected component of records linked by any of: shared `data_group_id`, shared
`source_volume_id`, shared `specimen_id`, or shared **verified** `animal_id`. The 108 data groups
collapse into **100 partition units**.

The ten fold-pair checks are in `../../validation/release_qc/qc_summary.json`, under `supporting_checks.partition_overlap_check.rows`:

| Fold pair | shared groups | shared source volumes | shared verified animals |
|---|---:|---:|---:|
| all 10 pairs | 0 | 0 | 0 |

Fold sizes: 19 / 19 / 21 / 20 / 21 units, i.e. 445 / 444 / 449 / 444 / 447 records.

**Residual limitation.** 29 units (408 records) belong to the rabbit-head subset whose animal
identity is not verifiable. Overlap between those units cannot be excluded. Use this split for
animal-independent evaluation of the rabbit-limb and mouse subsets; treat the rabbit-head subset
as a separately reported, exploratory subset.
