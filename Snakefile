# Flexible Snakemake workflow for TE feature extraction
# Supports: GTF annotation + flexible classification (FASTA or ID file)
# MODIFIED: Full GTF support, flexible classification, biotype tracking
# Author: Generated for bioinformatics pipeline
# Date: 2025-11-06

import os
import hashlib
import json
from pathlib import Path


# ============================================================================
# Configuration
# ============================================================================



# ============================================================================
# Wildcards
# ============================================================================
wildcard_constraints:
    dataset="[a-zA-Z0-9._-]+"


# ============================================================================
# Global variables
# ============================================================================
THREADS = config.get("threads", 4)
N_CHUNKS = config.get("n_chunks", 1)
if isinstance(N_CHUNKS, bool) or not isinstance(N_CHUNKS, int) or N_CHUNKS < 1:
    raise ValueError("n_chunks must be a positive integer")
CHUNKS = list(range(N_CHUNKS))

# Input files
GENCODE_GTF = config["gencode_gtf"]  # GENCODE v47 annotation - REQUIRED
SOURCE_GTF = config.get("source_gtf", "")  # Full GTF to subset from (optional)
GENCODE_FASTA = config.get("gencode_fasta", "")  # Full transcript FASTA (optional)

# Classification approach
CLASSIFICATION_MODE = config.get("classification_mode", "fasta")  # "fasta" or "id_file"

if CLASSIFICATION_MODE not in {"fasta", "id_file"}:
    raise ValueError("classification_mode must be fasta or id_file")

# Conditional inputs based on classification mode
if CLASSIFICATION_MODE == "fasta":
    PC_TRANSCRIPTS_FA = config.get("pc_transcripts_fasta", "")
    LNCRNA_TRANSCRIPTS_FA = config.get("lncrna_transcripts_fasta", "")
elif CLASSIFICATION_MODE == "id_file":
    PC_TRANSCRIPT_IDS = config.get("pc_transcript_ids_file", "")
    LNCRNA_TRANSCRIPT_IDS = config.get("lncrna_transcript_ids_file", "")

# RepeatMasker parameters
RM_SPECIES = config.get("repeatmasker_species", "human")

# Sequence mode: "spliced" (exonic FASTA) or "unspliced" (full genomic span via bedtools getfasta)
SEQUENCE_MODE = config.get("sequence_mode", "spliced")
GENOME_FASTA = config.get("genome_fasta", "")
if SEQUENCE_MODE not in {"spliced", "unspliced"}:
    raise ValueError("sequence_mode must be spliced or unspliced")

# ============================================================================
# Validation and setup
# ============================================================================
if SOURCE_GTF:
    if not os.path.exists(SOURCE_GTF):
        raise ValueError(f"source_gtf not found: {SOURCE_GTF}")
    # GENCODE_GTF will be created by the subset_gtf rule
elif not os.path.exists(GENCODE_GTF):
    raise ValueError(f"GENCODE GTF not found: {GENCODE_GTF}")

if SEQUENCE_MODE == "unspliced":
    if not GENOME_FASTA:
        raise ValueError("sequence_mode 'unspliced' requires genome_fasta in config")
    if not os.path.exists(GENOME_FASTA):
        raise ValueError(f"Genome FASTA not found: {GENOME_FASTA}")
elif SEQUENCE_MODE == "spliced":
    if not GENCODE_FASTA or not os.path.exists(GENCODE_FASTA):
        raise ValueError("sequence_mode 'spliced' requires gencode_fasta in config")

if CLASSIFICATION_MODE == "fasta":
    if not PC_TRANSCRIPTS_FA or not LNCRNA_TRANSCRIPTS_FA:
        raise ValueError(
            "fasta mode requires both pc_transcripts_fasta and lncrna_transcripts_fasta"
        )
    if not os.path.exists(PC_TRANSCRIPTS_FA) or not os.path.exists(LNCRNA_TRANSCRIPTS_FA):
        raise ValueError("Transcript FASTA files not found")

elif CLASSIFICATION_MODE == "id_file":
    if not PC_TRANSCRIPT_IDS or not LNCRNA_TRANSCRIPT_IDS:
        raise ValueError(
            "id_file mode requires both pc_transcript_ids_file and lncrna_transcript_ids_file"
        )
    if not os.path.exists(PC_TRANSCRIPT_IDS) or not os.path.exists(LNCRNA_TRANSCRIPT_IDS):
        raise ValueError("Transcript ID files not found")


# ============================================================================
# Target rules  — rule all MUST be first so it is the default Snakemake target
# ============================================================================
# Get list of datasets from config; append sequence mode suffix for non-default modes
_datasets_raw = config.get("datasets", ["default"])
_mode_suffix = "" if SEQUENCE_MODE == "spliced" else f"_{SEQUENCE_MODE}"
DATASETS = [f"{d}{_mode_suffix}" for d in _datasets_raw]

rule all:
    input:
        expand(
            [
                # GTF parsing
                "results/{dataset}/annotation/transcripts_from_gtf.bed",
                "results/{dataset}/annotation/transcript_biotypes.txt",
                "results/{dataset}/annotation/transcript_lengths.txt",
                # Classification (depends on mode)
                "results/{dataset}/annotation/pc_transcript_ids.txt",
                "results/{dataset}/annotation/lncrna_transcript_ids.txt",
                # RepeatMasker on full database
                "results/{dataset}/repeatmasker/all_transcripts.out.gff",
                "results/{dataset}/repeatmasker/all_transcripts.out",
                # Feature extraction on full database
                "results/{dataset}/features/all_transcripts_te_features.csv",
                "results/{dataset}/analysis/univariate_tests.csv",
                # Statistical analysis and visualization
                "results/{dataset}/analysis/summary_report.txt",
                "results/{dataset}/plots/hit_presence_comparison.png",
            ],
            dataset=DATASETS,
        )


# ============================================================================
# GTF subsetting subdag (only active when source_gtf is set in config)
# ============================================================================
if SOURCE_GTF:
    _pc_input   = PC_TRANSCRIPTS_FA   if CLASSIFICATION_MODE == "fasta" else PC_TRANSCRIPT_IDS
    _lnc_input  = LNCRNA_TRANSCRIPTS_FA if CLASSIFICATION_MODE == "fasta" else LNCRNA_TRANSCRIPT_IDS

    _subset_key = hashlib.sha256(json.dumps([os.path.abspath(p) for p in
        (GENCODE_GTF, SOURCE_GTF, _pc_input, _lnc_input)] + [CLASSIFICATION_MODE]).encode()).hexdigest()[:16]
    SUBSET_IDS = f"resources/annotation/subset_transcript_ids.{_subset_key}.txt"

    rule extract_fasta_ids_for_subset:
        """Collect transcript IDs from classification inputs for GTF subsetting."""
        input:
            pc=_pc_input,
            lncrna=_lnc_input,
        output:
            ids=SUBSET_IDS,
        log:
            f"logs/subset_gtf/{_subset_key}.extract_ids.log",
        run:
            ids = set()
            for source in (input.pc, input.lncrna):
                with open(source) as handle:
                    for line in handle:
                        if CLASSIFICATION_MODE == "fasta":
                            if line.startswith(">"):
                                ids.add(line[1:].split()[0].split("|")[0])
                        elif line.strip():
                            ids.add(line.strip())
            with open(output.ids, "w") as handle:
                handle.writelines(identifier + "\n" for identifier in sorted(ids))
            Path(log[0]).write_text(f"Collected {len(ids)} transcript IDs\n")

    rule subset_gtf:
        """Filter full GTF to only transcripts present in the classification inputs."""
        input:
            gtf=SOURCE_GTF,
            ids=SUBSET_IDS,
        output:
            gtf=GENCODE_GTF,
        conda:
            "workflow/envs/te_analysis.yaml"
        log:
            f"logs/subset_gtf/{_subset_key}.subset_gtf.log",
        shell:
            r"""
            awk -F'\t' 'BEGIN {{ while ((getline id < "{input.ids}") > 0) ids[id]=1 }}
                 /^#/ {{ print; next }}
                 match($9, /transcript_id "([^"]+)"/, a) && a[1] in ids {{ print }}' \
              {input.gtf:q} > {output.gtf:q} 2> {log:q}
            """

rule create_subset_gtf:
    """Standalone target: create the subset GTF (run before the main pipeline)."""
    input:
        GENCODE_GTF,


# ============================================================================
# Step 1: Parse GENCODE GTF - Extract transcripts and biotypes
# ============================================================================
rule parse_gencode_gtf:
    """Parse GENCODE GTF to extract transcript information."""
    input:
        gtf=GENCODE_GTF,
    output:
        bed="results/{dataset}/annotation/transcripts_from_gtf.bed",
        biotypes="results/{dataset}/annotation/transcript_biotypes.txt",
    conda:
        "workflow/envs/te_analysis.yaml"
    log:
        "logs/{dataset}/parse_gtf.log",
    resources:
        mem_mb=16000,
    shell:
        """
        python workflow/scripts/parse_gencode_gtf.py \
            --gtf {input.gtf} \
            --output-bed {output.bed} \
            --output-biotypes {output.biotypes} \
            2>&1 | tee {log}
        """

# ============================================================================
# Step 1.5: Extract transcript lengths from GTF (for later use in feature extraction)
# ============================================================================
rule extract_transcript_lengths:
    """Extract transcript lengths: exon sum (spliced) or genomic span (unspliced)."""
    input:
        gtf=GENCODE_GTF,
        bed="results/{dataset}/annotation/transcripts_from_gtf.bed",
    output:
        lengths="results/{dataset}/annotation/transcript_lengths.txt",
    conda:
        "workflow/envs/te_analysis.yaml"
    log:
        "logs/{dataset}/extract_lengths.log",
    resources:
        mem_mb=16000,
    run:
        if SEQUENCE_MODE == "unspliced":
            # ponytail: col4=transcript_id, end-start gives genomic span (BED is 0-based)
            shell("awk '{{print $4 \"\\t\" $3-$2}}' {input.bed:q} > {output.lengths:q} 2> {log:q}")
        else:
            shell(
                """
                awk -F"\\t" '$3 == "exon" {{
                    match($9, /transcript_id "([^"]+)"/, arr);
                    L[arr[1]] += $5 - $4 + 1
                }} END {{
                    for (t in L) print t "\\t" L[t]
                }}' {input.gtf:q} > {output.lengths:q} 2> {log:q}
                """
            )


# ============================================================================
# Step 2: Handle Classification - FASTA or ID file mode
# ============================================================================
# MODE A: Extract IDs from FASTA files
rule extract_ids_from_fasta:
    """Extract transcript IDs from FASTA files (FASTA mode)."""
    input:
        pc_fa=PC_TRANSCRIPTS_FA if CLASSIFICATION_MODE == "fasta" else PC_TRANSCRIPT_IDS,
        lncrna_fa=LNCRNA_TRANSCRIPTS_FA if CLASSIFICATION_MODE == "fasta" else LNCRNA_TRANSCRIPT_IDS,
    output:
        pc_ids="results/{dataset}/annotation/pc_transcript_ids.txt",
        lncrna_ids="results/{dataset}/annotation/lncrna_transcript_ids.txt",
    conda:
        "workflow/envs/te_analysis.yaml"
    log:
        "logs/{dataset}/extract_ids.log",
    run:
        if CLASSIFICATION_MODE == "fasta":
            shell(
                """
            # Extract IDs from pc_transcripts.fa
            grep "^>" {input.pc_fa} | \
                sed 's/^>//g' | \
                cut -d'|' -f1 > {output.pc_ids}

            # Extract IDs from lncrna_transcripts.fa
            grep "^>" {input.lncrna_fa} | \
                sed 's/^>//g' | \
                cut -d'|' -f1 > {output.lncrna_ids}
            """
            )
        else:
            # If using ID file mode, just copy the files
            shell("cp {input.pc_fa:q} {output.pc_ids:q}")
            shell("cp {input.lncrna_fa:q} {output.lncrna_ids:q}")


# ============================================================================
# Step 2.5: Prepare transcript FASTA (spliced or unspliced)
# ============================================================================
rule prepare_transcript_fasta:
    """Create all_transcripts.fa: symlink spliced FASTA or extract unspliced sequences."""
    input:
        bed="results/{dataset}/annotation/transcripts_from_gtf.bed",
        sequence=GENOME_FASTA if SEQUENCE_MODE == "unspliced" else GENCODE_FASTA,
    output:
        fa="results/{dataset}/annotation/all_transcripts.fa",
    conda:
        "workflow/envs/te_analysis.yaml"
    log:
        "logs/{dataset}/prepare_transcript_fasta.log",
    run:
        if SEQUENCE_MODE == "unspliced":
            shell(
                "bedtools getfasta -fi {input.sequence:q} -bed {input.bed}"
                " -nameOnly -s -fo {output.fa} 2>&1 | tee {log}"
            )
        else:
            shell("ln -sf $(realpath {input.sequence:q}) {output.fa} 2>&1 | tee {log}")

# ============================================================================
# Step 3: Index and prepare for RepeatMasker
# ============================================================================
rule index_transcripts:
    """Index transcript FASTA for RepeatMasker."""
    input:
        fa="results/{dataset}/annotation/all_transcripts.fa",
    output:
        fai="results/{dataset}/annotation/all_transcripts.fa.fai",
    conda:
        "workflow/envs/te_analysis.yaml"
    log:
        "logs/{dataset}/index_transcripts.log",
    shell:
        """
        samtools faidx {input.fa} 2>&1 | tee {log}
        """


# ============================================================================
# Step 4: Run RepeatMasker on full database
# ============================================================================
rule check_fasta_headers:
    """Check and clean FASTA headers for RepeatMasker compatibility."""
    input:
        fa="results/{dataset}/annotation/all_transcripts.fa",
    output:
        cleaned_fa="results/{dataset}/annotation/all_transcripts_headers_checked.fa",
    conda:
        "workflow/envs/te_analysis.yaml"
    log:
        "logs/{dataset}/check_fasta_headers.log",
    shell:
        """
        python workflow/scripts/check_fasta_headers.py \
            --input {input.fa} \
            --output {output.cleaned_fa} \
            2>&1 | tee {log}
        """


rule split_fasta_for_repeatmasker:
    """Split transcript FASTA into N_CHUNKS pieces for parallel RepeatMasker."""
    input:
        fa="results/{dataset}/annotation/all_transcripts_headers_checked.fa",
    output:
        expand("results/{{dataset}}/repeatmasker/chunks/chunk_{chunk}.fa", chunk=CHUNKS),
    log:
        "logs/{dataset}/split_fasta.log",
    run:
        total = sum(1 for line in open(input.fa) if line.startswith('>'))
        if total == 0:
            raise ValueError("Transcript FASTA contains no sequences")
        chunk_size = (total + N_CHUNKS - 1) // N_CHUNKS
        handles = [open(p, 'w') for p in output]
        idx, seq_count = 0, 0
        with open(input.fa) as f:
            for line in f:
                if line.startswith('>'):
                    if seq_count > 0 and seq_count % chunk_size == 0:
                        idx = min(idx + 1, N_CHUNKS - 1)
                    seq_count += 1
                handles[idx].write(line)
        for h in handles:
            h.close()
        with open(log[0], 'w') as lf:
            lf.write(f"Split {total} sequences into {N_CHUNKS} chunks (chunk_size={chunk_size})\n")


rule run_repeatmasker_chunk:
    """Run RepeatMasker on one FASTA chunk."""
    input:
        fa="results/{dataset}/repeatmasker/chunks/chunk_{chunk}.fa",
        helper="workflow/scripts/repeatmasker_outputs.py",
    output:
        gff="results/{dataset}/repeatmasker/chunks/chunk_{chunk}.out.gff",
        out="results/{dataset}/repeatmasker/chunks/chunk_{chunk}.out",
    conda:
        "workflow/envs/te_analysis.yaml"
    params:
        species=RM_SPECIES,
        outdir=lambda wc, output: str(Path(output.gff).parent),
    threads: THREADS
    resources:
        mem_mb=46000,
        runtime="3d",
    log:
        "logs/{dataset}/repeatmasker_chunk_{chunk}.log",
    benchmark:
        "benchmarks/{dataset}/repeatmasker_chunk_{chunk}.txt",
    shell:
        """
        if [ -s {input.fa:q} ]; then
            # Do not mistake stale raw files for successful output on a rerun.
            rm -f -- {input.fa:q}.out {input.fa:q}.out.gff
            RepeatMasker -species {params.species:q} -pa {threads} -gff \
                -dir {params.outdir:q} -s {input.fa:q} 2>&1 | tee {log:q}
            python {input.helper:q} normalize \
                --raw-out {input.fa:q}.out --raw-gff {input.fa:q}.out.gff \
                --out {output.out:q} --gff {output.gff:q}
        else
            echo "Empty chunk: no RepeatMasker execution needed" > {log:q}
            python {input.helper:q} normalize --empty \
                --out {output.out:q} --gff {output.gff:q}
        fi
        """


rule merge_repeatmasker_chunks:
    """Merge per-chunk RepeatMasker outputs into a single file."""
    input:
        gffs=expand("results/{{dataset}}/repeatmasker/chunks/chunk_{chunk}.out.gff", chunk=CHUNKS),
        outs=expand("results/{{dataset}}/repeatmasker/chunks/chunk_{chunk}.out", chunk=CHUNKS),
        helper="workflow/scripts/repeatmasker_outputs.py",
    output:
        gff="results/{dataset}/repeatmasker/all_transcripts.out.gff",
        out="results/{dataset}/repeatmasker/all_transcripts.out",
    params:
        n_chunks=N_CHUNKS,
    log:
        "logs/{dataset}/merge_repeatmasker.log",
    conda:
        "workflow/envs/te_analysis.yaml"
    shell:
        """
        python {input.helper:q} merge \
            --outs {input.outs:q} --gffs {input.gffs:q} \
            --out {output.out:q} --gff {output.gff:q} 2>&1 | tee {log:q}
        """


# ============================================================================
# Step 5: Extract TE features for full database
# ============================================================================
rule extract_all_features:
    """Extract TE features for all transcripts."""
    input:
        repeatmasker="results/{dataset}/repeatmasker/all_transcripts.out",
        bed="results/{dataset}/annotation/transcripts_from_gtf.bed",
        biotypes="results/{dataset}/annotation/transcript_biotypes.txt",
        lengths="results/{dataset}/annotation/transcript_lengths.txt",
        pc_ids="results/{dataset}/annotation/pc_transcript_ids.txt",
        lncrna_ids="results/{dataset}/annotation/lncrna_transcript_ids.txt",
    output:
        features="results/{dataset}/features/all_transcripts_te_features.csv",
        summary="results/{dataset}/features/all_transcripts_te_summary.txt",
    conda:
        "workflow/envs/te_analysis.yaml"
    log:
        "logs/{dataset}/extract_features_all.log",
    shell:
        """
        python workflow/scripts/te_feature_extractor.py \
            --repeatmasker {input.repeatmasker} \
            --transcripts {input.bed} \
            --biotypes {input.biotypes} \
            --lengths {input.lengths} \
            --pc-ids {input.pc_ids} \
            --lnc-ids {input.lncrna_ids} \
            --output-prefix results/{wildcards.dataset}/features/all_transcripts \
            2>&1 | tee {log}
        """


# ============================================================================
# Step 6: Filter features by classification
# ============================================================================
rule filter_features_by_class:
    """Filter TE features for pc vs lncRNA classes."""
    input:
        features="results/{dataset}/features/all_transcripts_te_features.csv",
        pc_ids="results/{dataset}/annotation/pc_transcript_ids.txt",
        lncrna_ids="results/{dataset}/annotation/lncrna_transcript_ids.txt",
    output:
        pc_features="results/{dataset}/features/pc_transcripts_te_features.csv",
        lncrna_features="results/{dataset}/features/lncrna_transcripts_te_features.csv",
    conda:
        "workflow/envs/te_analysis.yaml"
    log:
        "logs/{dataset}/filter_features.log",
    shell:
        """
        python workflow/scripts/filter_features_by_class.py \
            --features {input.features} \
            --pc-ids {input.pc_ids} \
            --lncrna-ids {input.lncrna_ids} \
            --output-pc {output.pc_features} \
            --output-lncrna {output.lncrna_features} \
            2>&1 | tee {log}
        """


# ============================================================================
# Step 7: Combine features and prepare for analysis
# ============================================================================
rule combine_classified_features:
    """Combine pc and lncRNA features with group labels."""
    input:
        pc_features="results/{dataset}/features/pc_transcripts_te_features.csv",
        lncrna_features="results/{dataset}/features/lncrna_transcripts_te_features.csv",
    output:
        combined="results/{dataset}/combined/classified_te_features.csv",
    conda:
        "workflow/envs/te_analysis.yaml"
    log:
        "logs/{dataset}/combine_features.log",
    run:
        import pandas as pd

        pc_df = pd.read_csv(input.pc_features)
        lncrna_df = pd.read_csv(input.lncrna_features)

        pc_df["group"] = "Coding"
        lncrna_df["group"] = "lncRNA"

        combined = pd.concat([pc_df, lncrna_df], ignore_index=True)
        combined.to_csv(output.combined, index=False)
        print("Combined features saved to {output.combined}")


# ============================================================================
# Step 8: Statistical analysis
# ============================================================================
rule statistical_analysis:
    """Perform statistical analysis on classified features."""
    input:
        features="results/{dataset}/features/all_transcripts_te_features.csv",
    output:
        univariate="results/{dataset}/analysis/univariate_tests.csv",
        categorical="results/{dataset}/analysis/categorical_tests.csv",
        pca_scores="results/{dataset}/analysis/pca_scores.csv",
        summary="results/{dataset}/analysis/summary_report.txt",
    conda:
        "workflow/envs/te_analysis.yaml"
    log:
        "logs/{dataset}/statistical_analysis.log",
    shell:
        """
        python workflow/scripts/te_statistical_analyzer.py \
            --features {input.features} \
            --output-prefix results/{wildcards.dataset}/analysis \
            2>&1 | tee {log}
        """


# ============================================================================
# Step 9: Visualizations
# ============================================================================
rule generate_visualizations:
    """Generate comprehensive visualizations."""
    input:
        features="results/{dataset}/features/all_transcripts_te_features.csv",
        tests="results/{dataset}/analysis/univariate_tests.csv",
        pca="results/{dataset}/analysis/pca_scores.csv",
    output:
        presence="results/{dataset}/plots/hit_presence_comparison.png",
        volcano="results/{dataset}/plots/volcano_plot.png",
        pca="results/{dataset}/plots/pca_plot.png",
    conda:
        "workflow/envs/te_analysis.yaml"
    log:
        "logs/{dataset}/visualization.log",
    shell:
        """
        python workflow/scripts/te_visualizer.py \
            --features {input.features} \
            --test-results {input.tests} \
            --pca-scores {input.pca} \
            --output-dir results/{wildcards.dataset}/plots \
            2>&1 | tee {log}
        """


# ============================================================================
# Clean-up rules
# ============================================================================
rule clean:
    """Remove all output files for a dataset."""
    conda:
        "workflow/envs/te_analysis.yaml"
    log:
        "logs/{dataset}/clean.log",
    shell:
        """
        rm -rf results/{wildcards.dataset}/*
        echo "All output files removed for dataset {wildcards.dataset}." 2>&1 | tee {log}
        """

# ===========================================================================
# Additional rule targets
# ===========================================================================

rule pipeline_all_features:
    input:
        expand(rules.extract_all_features.output.features, dataset=DATASETS),
