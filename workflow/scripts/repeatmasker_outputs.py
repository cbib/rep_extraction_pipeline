#!/usr/bin/env python3
"""Normalize and merge RepeatMasker outputs without assuming header length.

A missing GFF is allowed only for an empty chunk or a recognized .out table
with no hits (including RepeatMasker's no-repeat message). A missing .out for
a nonempty chunk remains an error even after a successful process exit.
"""
import argparse
from pathlib import Path

HEADER = (
    "   SW   perc perc perc  query position in query matching repeat position in repeat\n"
    "score   div. del. ins.  sequence begin end (left) repeat class/family begin end (left) ID\n"
    "\n"
)
NO_HITS = "There were no repetitive sequences detected"


def data_rows(path):
    """Yield data rows and reject unrecognized non-header text."""
    recognized = False
    with open(path) as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            fields = stripped.split()
            if fields[0].isdigit():
                if len(fields) not in (15, 16):
                    raise ValueError(
                        f"Malformed RepeatMasker row in {path}: {stripped}"
                    )
                for index in (0, 5, 6, 14):
                    int(fields[index])
                for index in (1, 2, 3):
                    float(fields[index])
                recognized = True
                yield line if line.endswith("\n") else line + "\n"
            elif NO_HITS.lower() in stripped.lower() or fields[0] in {"SW", "score"}:
                recognized = True
            else:
                raise ValueError(
                    f"Unrecognized RepeatMasker output in {path}: {stripped}"
                )
    if not recognized:
        raise ValueError(f"Empty or unrecognized RepeatMasker output: {path}")


def normalize_outputs(raw_out, raw_gff, output_out, output_gff, empty=False):
    output_out, output_gff = Path(output_out), Path(output_gff)
    output_out.parent.mkdir(parents=True, exist_ok=True)
    output_gff.parent.mkdir(parents=True, exist_ok=True)
    hits = 0
    with output_out.open("w") as handle:
        handle.write(HEADER)
        if not empty:
            for row in data_rows(raw_out):
                handle.write(row)
                hits += 1
    if not empty and hits:
        if raw_gff is None or not Path(raw_gff).is_file():
            raise ValueError("RepeatMasker produced hits but no GFF output")
        with open(raw_gff) as source, output_gff.open("w") as target:
            n_features = 0
            for line in source:
                target.write(line)
                if line.strip() and not line.startswith("#"):
                    n_features += 1
        if n_features != hits:
            raise ValueError(
                f"RepeatMasker GFF/.out hit counts differ: {n_features} vs {hits}"
            )
    else:
        output_gff.write_text("##gff-version 2\n")


def merge_outputs(outs, gffs, output_out, output_gff):
    if len(outs) != len(gffs) or not outs:
        raise ValueError("Expected matching nonempty lists of chunk outputs")
    # Parsing every row preserves all hits even for headerless .out files.
    with open(output_out, "w") as target:
        target.write(HEADER)
        for path in outs:
            for row in data_rows(path):
                target.write(row)
    versions = set()
    for path in gffs:
        version = None
        with open(path) as source:
            for line in source:
                if line.startswith("##gff-version"):
                    version = line.strip()
                elif line.strip() and not line.startswith("#"):
                    if version is None:
                        raise ValueError(
                            f"GFF features lack version declaration: {path}"
                        )
                    versions.add(version)
                    break
    if len(versions) > 1:
        raise ValueError("Cannot merge different GFF format versions")
    with open(output_gff, "w") as target:
        target.write(next(iter(versions), "##gff-version 2") + "\n")
        for path in gffs:
            with open(path) as source:
                for line in source:
                    if not line.startswith("##gff-version"):
                        target.write(line)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    normalize = sub.add_parser("normalize")
    normalize.add_argument("--raw-out")
    normalize.add_argument("--raw-gff")
    normalize.add_argument("--empty", action="store_true")
    merge = sub.add_parser("merge")
    merge.add_argument("--outs", nargs="+", required=True)
    merge.add_argument("--gffs", nargs="+", required=True)
    for command in (normalize, merge):
        command.add_argument("--out", required=True)
        command.add_argument("--gff", required=True)
    args = parser.parse_args()
    if args.command == "normalize":
        if not args.empty and not args.raw_out:
            parser.error("--raw-out is required unless --empty is used")
        normalize_outputs(args.raw_out, args.raw_gff, args.out, args.gff, args.empty)
    else:
        merge_outputs(args.outs, args.gffs, args.out, args.gff)


if __name__ == "__main__":
    main()
