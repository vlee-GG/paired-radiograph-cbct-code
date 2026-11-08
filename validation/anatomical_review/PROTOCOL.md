# Two-clinician blinded anatomical assessment

## Sample and reference classes

The frozen study was prepared on 16 October 2026 (sampling seed 2026101603). Source-balanced sampling selected 250 original pairs (170 rabbit limb, 30 maxilla, 20 mandible and 30 mouse) and 50 substituted controls (30, 8, 6 and 6, respectively), covering 251 source volumes and 95 query partition units. The CT crop in each control came from the same source volume and anatomical subset, with zero three-dimensional crop-box overlap and a box-centre distance of at least 30 mm in array space. These are constructed record-level reference classes; a non-overlapping box can contain similar repeated anatomy.

Twenty training cases (14 original, six substituted) came from source volumes excluded from the formal sample. The formal orders contained 315 presentations per reader: 300 distinct cases and 15 hidden repeat presentations (12 original, three substituted), separated by at least 40 positions. Previously reviewed records were excluded during sampling. Case membership and presentation orders are frozen in the released CSV tables.

## Readers and task

Two clinicians independently reviewed the same cases in separately randomised orders. Reader 1 had three years of clinical experience in orthopaedics at a county-level hospital. Reader 2 practised in a dental clinic. These background descriptions were supplied by the authors.

Each case displayed the released 90-degree radiograph (A), a whole-crop CBCT maximum-intensity projection (B), and two orthogonal CBCT sections (C, D). Readers assessed approximate correspondence of local bony structures, joint morphology and anatomical relationships, using yes, no or uncertain. Pixel-level overlap or identical size and projection appearance was not required. Constructed class labels, record identities and the other reader's answers were concealed. Reasons could be entered for uncertain assessments; confidence scores were not collected.

## Display generation

Radiographs were displayed in the released orientation. The CT panels used `cube.max(axis=2).T`, `cube[:,:,k].T` and `cube[i,:,:]`. Section indices maximised the fraction of voxels above 300 rescaled units, with the centre index used when none exceeded the threshold. Radiograph and CBCT windows were [3, 59745] and [-1000.0001220703, 1780.0089111328], respectively. Every composite used the same 1560 by 1560 layout without physical scale bars or anatomical direction markers. Display choices and original PNG hashes are in `display_choices.csv`; input file hashes are in `input_hashes.csv`. The browser used WebP derivatives of these composites; captions were subsequently changed from 平片 to X光. Panel image content and class assignments were preserved.

## First-submission analysis

For each composite, the first actual exposure in each reader's formal order defines the primary observation. Answers are those present at that reader's first confirmed final submission, reconstructed from saved events. Edits made before confirmation are included; later rechecks do not replace the primary scores. Training and extra repeat exposures are excluded from the denominator of 300.

Label agreement is `(original rated yes + substituted rated no) / 300`. Uncertain responses remain in the denominator and are also reported separately. Counts by reference class, three-category between-reader agreement, Cohen's kappa and Gwet's AC1 are reported. Hidden repeats assess within-reader response agreement. Confidence intervals use 20,000 percentile bootstrap resamples of query partition units (seed 20261020), with source-volume clustering as a sensitivity analysis. They are conditional on the sampled cases and constructed labels.

## Subsequent author review

One author subsequently reviewed all 91 cases in which at least one reader's answer differed from the reference label, including uncertain responses. These comprised 89 between-reader disagreements and two shared deviations. Reference labels and both readers' first-submission answers were available during follow-up. The principal displays were the 90-degree radiograph and the whole-crop CBCT maximum-intensity projection. Case-level conclusions, structural criteria and actual review dates are recorded in `author_review_records.csv`. The original 89 cases were reviewed on 24 October 2026. The two shared deviations, c0065 and c0280, were added and reviewed on 25 October 2026. The criterion for c0053 was corrected on author instruction on 25 October 2026. The follow-up was unblinded; original labels and first-submission answers were retained.

## Follow-up record interpretation

`author_review_records.csv` contains all 91 case mappings, the actual reviewer and review date, final conclusions, selected structural criteria, displayed/confirmed materials and save revision. `review_records.csv` identifies the follow-up scope using `author_review_in_scope`, independently of `readers_disagree`. The original reference labels and clinician answers are unchanged; corresponding counts are under `author_followup` in `review_stats.json`.

Materials marked `shown_in_page` record the panels displayed and do not imply that an external CT volume was viewed. Author follow-up was unblinded and remains separate from the first-submission scores.
