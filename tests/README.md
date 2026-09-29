# Review regression tests

From the repository root, activate an analysis environment containing the
packages in `workflow/envs/te_analysis.yaml` and put Snakemake on `PATH`, then run:

```bash
python -m unittest discover -s tests -v
```

If Snakemake is installed in a separate environment, set `SNAKEMAKE` to its
executable path. The active `python` and `bedtools` must be available on `PATH`.
Missing Snakemake or bedtools causes the corresponding tests to skip; inspect
the test summary before treating this as an integration pass.

The workflow tests copy source into temporary directories and generate two
synthetic transcripts. They exercise real GTF parsing, mode-specific sequence
preparation and length normalization, chunk splitting, output merging and
feature extraction. A **deterministic RepeatMasker stub** generates known hits
and explicit no-hit outputs; no repeat library is downloaded or searched.
Three chunks for two transcripts additionally exercise the empty-chunk path.
These tests check workflow behavior and feature accounting, not RepeatMasker
sensitivity, biological validity, cluster execution or environment resolution.

No existing paper inputs, outputs or Snakemake metadata are read or modified.
