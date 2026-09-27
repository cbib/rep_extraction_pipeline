# TE pipeline

This Snakemake workflow annotates transcript sequences with RepeatMasker, extracts transcript-level TE features, and compares protein-coding with lncRNA transcripts. Its default `all` target produces the full feature table, statistical reports, and plots. Class-specific and combined feature tables require explicit targets.

## Inputs

Run from `te_pipeline/`. The published templates are:
- `config/config.default.spliced.yaml` to scan spliced transcripts.
- `config/config.default.unspliced.yaml` to scan unspliced transcripts (or full gene bodies).

Reference files are absent from a fresh checkout and should be obtained from GENCODE. Place them at the configured paths or override those paths. GTF and genome contig names must match.

| Input | Format and naming rule |
| --- | --- |
| `source_gtf` | Full GENCODE GTF with `transcript_id` attributes. When set, the workflow selects classification IDs into `gencode_gtf`; otherwise supply `gencode_gtf` directly. Use a distinct subset output path for each cohort. |
| `gencode_fasta` | Mature transcript FASTA for `sequence_mode: spliced`; headers must identify transcripts in the GTF. |
| `genome_fasta` | Genome FASTA for `sequence_mode: unspliced`; bedtools extracts strand-aware genomic transcript spans, including introns. |
| `pc_transcripts_fasta`, `lncrna_transcripts_fasta` | Separate class FASTAs for `classification_mode: fasta`. The first header token before `|` is used as the ID. |
| `pc_transcript_ids_file`, `lncrna_transcript_ids_file` | Alternative plain-text files, one transcript ID per line, for `classification_mode: id_file`. |

## Requirements

Conda or Mamba creates two kinds of environments: `te_driver` runs Snakemake 9 and starts the workflow; Snakemake creates separate environments for individual rules. Those rule environments supply RepeatMasker, bedtools, samtools, GNU awk, and Python analysis packages. RepeatMasker needs a configured search engine and compatible repeat library. An initial Conda solve may need network access; reference downloads are separate. No GPU is required. Full-cohort runs need substantial CPU, memory, storage, and time. SLURM is optional.

## Run

From `te_pipeline/`, create and activate the Conda environment that runs Snakemake, then copy both published templates. Edit each local config for its matching references and cohort. The spliced run sends the transcript FASTA directly to RepeatMasker; the unspliced run extracts genomic transcript spans using the subset GTF and genome FASTA, avoiding a RepeatMasker run on the full genome sequence.

```bash
conda env create -f workflow/envs/te_driver.yaml
conda activate te_driver
cp config/config.default.spliced.yaml config/config.spliced.local.yaml
cp config/config.default.unspliced.yaml config/config.unspliced.local.yaml
# Spliced
snakemake --snakefile Snakefile --configfile config/config.spliced.local.yaml \
  --profile profiles/default --cores 8

# Unspliced
snakemake --snakefile Snakefile --configfile config/config.unspliced.local.yaml \
  --profile profiles/default --cores 8
```

Add `--dry-run` to either command to check its DAG before execution. For SLURM, create and activate the Snakemake environment from `workflow/envs/te_driver_slurm.yaml` instead, configure your site's account, partition, and resource limits, then replace `--profile profiles/default --cores 8` with `--profile profiles/slurm`.

For a smaller check, activate an environment with the packages in `workflow/envs/te_analysis.yaml` and bedtools on `PATH`, then run the synthetic suite from this directory:

```bash
export SNAKEMAKE="$(conda run -n te_driver which snakemake)"
python -m unittest discover -s tests -v
```

The tests stub RepeatMasker. `te_driver` supplies Snakemake, while the active `python` must supply rule dependencies such as `tqdm`. Set `XDG_CACHE_HOME` to a writable directory if your home cache is read-only.

## Configuration

Always pass `--configfile`: the Snakefile has no default dataset config. Values below are from the published unspliced template; the spliced template differs in sequence source, mode, chunk count, and output directory. Numbered configs are local run records: `config/config_47.yaml` and `config/config_49.yaml` (spliced), `config/config_47_unspliced.yaml` (one chunk), `config/config_47_unspliced_chunks.yaml` (eight-chunk paper run), and `config/config_49_unspliced.yaml` (exclusive review cohort). Edit copied templates or use `--config` for overrides.

### Path parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `source_gtf` | `resources/annotation/gencode.v47.annotation.gtf` | Full GTF for cohort subsetting; optional if `gencode_gtf` already exists. |
| `gencode_gtf` | `resources/annotation/gencode.v47.annotation.unspliced.subset.gtf` | Selected GTF output; the spliced template uses `.spliced.subset.gtf`. Only the unspliced run uses its coordinates to create RepeatMasker input. |
| `genome_fasta` | `resources/genome/GRCh38.primary_assembly.genome.fa` | Unspliced sequence source; absent from the spliced template. |
| `gencode_fasta` | Not set | Spliced template uses `resources/transcripts/gencode.v47.transcripts.fa` directly, without GTF filtering. |
| `pc_transcripts_fasta`, `lncrna_transcripts_fasta` | `resources/transcripts/gencode.v47.pc_transcripts.fa`, `resources/transcripts/gencode.v47.lncRNA_transcripts.fa` | Class membership inputs; only their IDs are needed for unspliced extraction. |
| `pc_transcript_ids_file`, `lncrna_transcript_ids_file` | Not set | Supply both paths when switching to `id_file` mode. |

### Functioning parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `datasets` | `[gencode.v47]` | Output dataset name; unspliced mode appends `_unspliced`. |
| `sequence_mode` | `unspliced` | `unspliced` uses genomic transcript spans; `spliced` uses mature transcript FASTA. |
| `classification_mode` | `fasta` | Choose `fasta` or `id_file` for coding/lncRNA membership. |
| `n_chunks` | `8` | Positive integer number of RepeatMasker FASTA chunks. |
| `repeatmasker_species` | `human` | Species passed to RepeatMasker. |
| Analysis thresholds and model settings | Not set | Fixed in the analysis scripts; legacy fields in numbered configs are unused. |

For example, add `--config n_chunks=1 genome_fasta=/path/to/GRCh38.fa` to an unspliced Snakemake command to override those two values.

### Visualization customization

| Parameter | Default | Meaning |
| --- | --- | --- |
| Plot output directory | `results/gencode.v47_unspliced/plots/` | Derived from `datasets` and `sequence_mode`; spliced output is `results/gencode.v47/plots/`. |
| Plot selection | Presence, volcano, PCA | Fixed by `generate_visualizations`; no config selection. |
| Figure format and colors | PNG and script defaults | Fixed by the plotting script; legacy config fields are unused. |

### Environment

| Parameter | Default | Meaning |
| --- | --- | --- |
| `threads` | `20` | CPU threads requested per RepeatMasker chunk; Snakemake scales to available cores. |
| Local profile | `profiles/default`: local executor, `cores: 1` | Example command overrides cores with `--cores 8`; rule Conda environments are enabled. |
| SLURM profile | `profiles/slurm`: `jobs: 30` | Cluster job limit and resource overrides; requires the SLURM executor plugin and site settings. |
| `latency-wait` | `45` seconds in both profiles | Wait for output files on shared storage. |

## Outputs

For the published templates, `{dataset}` is `gencode.v47` (spliced) or `gencode.v47_unspliced` (unspliced). Paths are relative to `te_pipeline/`.

| Main result | Meaning |
| --- | --- |
| `results/{dataset}/repeatmasker/all_transcripts.out` and `.out.gff` | Merged RepeatMasker annotations. A no-hit chunk contributes no annotations. |
| `results/{dataset}/features/all_transcripts_te_features.csv` | Full TE feature table, including transcripts without RepeatMasker hits (zero counts); use this for downstream analysis. |
| `results/{dataset}/analysis/univariate_tests.csv`, `categorical_tests.csv`, `pca_scores.csv`, `summary_report.txt` | Comparisons and summary for coding and lncRNA transcripts; unclassified transcripts are excluded. |
| `results/{dataset}/plots/hit_presence_comparison.png`, `volcano_plot.png`, `pca_plot.png` | Tracked comparison plots. |

The default target names result files, not a completion marker. Target `pipeline_all_features` to stop at the full feature table. Class-specific CSVs under `features/` and `combined/classified_te_features.csv` are optional explicit targets. Scripts also write auxiliary files not all declared as Snakemake outputs; see `te_pipeline_manifest.yaml`.

## Reproducibility and limits

Environment YAMLs constrain major versions but are not solved lockfiles. Saved v47 unspliced logs report RepeatMasker 4.2.2, NCBI/RMBlast 2.14.1+, and Dfam 3.9; original reference checksums and the solved environment are unknown. The historical feature table was not regenerated after current workflow fixes. In the paper workspace, Snakemake 9.22.0 dry-runs resolved the v47 spliced, v47 eight-chunk unspliced, and v49 unspliced DAGs; all 11 unit and synthetic workflow tests passed.
