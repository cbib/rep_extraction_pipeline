# Paper revision release fixes

Review branch: `review/paper-release-blockers`, based on `0148e59` (`dev`).
Parent review branch: `review/te-paper-release`, based on `80f5d9a`. These changes
were committed only on the review branches; the original `dev` and parent `tests`
branches were not modified by this work. No remote push is included.

## Changes

- Commit the three unspliced review configs. Preserve the v49 exclusive cohort,
  replacing its machine-specific paths with relative sibling-project paths.
- Track sequence and classification files as workflow inputs; fix the feature-only
  target to produce CSVs; validate modes and chunk counts.
- Namespace subset-ID intermediates by annotation and cohort paths; keep stderr
  out of annotation and length output files.
- Normalize and merge RepeatMasker output without dropping rows by fixed header
  position; handle empty/no-hit chunks, validate output/GFF agreement, and preserve
  GFF versions. Do not reuse stale raw outputs after a failed/incomplete invocation.
- Preserve zero-valued global feature columns for entirely no-hit feature exports.
- Remove the hardcoded recovery rule; document ordinary chunk reruns with
  `--rerun-incomplete`. Preserve existing recovery directories and ignore them.
- Declare statsmodels and GNU awk, keep pandas in the tested 2.x series, and supply
  separate local and SLURM driver environments/profiles.
- Correct startup commands, rule/target descriptions, normalization terminology,
  reference requirements, and provenance limitations.

## Validation

- Eleven regression tests passed with no skips: six RepeatMasker output boundary
  tests and five workflow tests. Real GTF parsing, sequence preparation, splitting,
  merging and feature extraction run in temporary directories; RepeatMasker is a
  deterministic stub. Both spliced and unspliced lengths, strand-suffixed IDs,
  all-no-hit exports, empty chunks, FASTA cohort subsetting, and same-path input
  content-change invalidation are covered.
- Test environment: Python 3.12, pandas 2.3.3, numpy 2.3.5, scipy 1.16.3,
  scikit-learn 1.7.2, statsmodels 0.14.5; workflow driver Snakemake 9.13.2.
- The updated paper v47 unspliced workflow builds in a full forced dry-run with
  the portable default profile. No paper outputs were rerun or rewritten.
- Parsed all 3,021,590 rows of an existing RepeatMasker chunk and verified that
  its GFF contains the same number of features.
- All three final environment specifications (analysis, local driver and SLURM
  driver) resolved successfully in online Conda dry-run solves on Linux x86-64.
  No environment was installed.
- Repository pre-commit checks and whitespace checks passed for code changes.

## Limits and release handoff

These checks do not constitute a new full-cohort RepeatMasker run, a fresh
installed-environment execution, or SLURM submission. Environment specifications
are not immutable lockfiles. Exact historical repeat-library provenance cannot
be recovered by changing code; preserve and publish the available logs, reference
checksums, cohort files and environment/library metadata with the paper artifacts.

The separate `fix/features_mem` branch and three existing stashes are deliberately
not applied. They require a separate output-equivalence review; no work was dropped.
The cbib remote advertised only preprint `main` at `7f3602a` when checked. Publish
this submodule review branch before publishing a parent revision that references
its new commit. Local commits alone are not available to remote fresh clones.
