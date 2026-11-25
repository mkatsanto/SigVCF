"""Command-line entry point for SigVCF (SBS96 only)."""

import argparse
import os
import re

import pandas as pd

from .sbs96 import build_sbs96_channels, compute_sbs96
from .io_formats import load_matrix_order, save_in_format


def _sample_name_from_path(vcf_path):
    base = os.path.basename(vcf_path)
    return re.sub(r"\.vcf(\.gz)?$", "", base)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="sigvcf",
        description="Extract an SBS96 mutation-count matrix directly from "
                    "VCFs, streaming -- no whole-genome or whole-VCF load "
                    "into memory. Runs fully offline.",
    )
    parser.add_argument("--matrix-type", choices=["sbs96"], default="sbs96",
                        help="Only sbs96 is supported; kept so existing callers keep working.")
    parser.add_argument("--vcf", nargs="+", required=True,
                        help="One or more .vcf/.vcf.gz files, one per sample.")
    parser.add_argument("--sample-names", nargs="+", default=None,
                        help="Column names, one per --vcf, in the same order. "
                             "Defaults to each VCF's filename stem.")
    parser.add_argument("--fasta", required=True,
                        help="Reference FASTA with a matching .fai index alongside it.")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--output-prefix", default="counts_matrix",
                        help="Output files are written as {output-prefix}_{format}.dat")
    parser.add_argument("--output-formats", nargs="+", default=["tab_with_index"],
                        help="Any name; only 'sigfit' and 'sigLASSO' get special "
                             "delimiter/quoting, everything else is tab-separated.")
    parser.add_argument("--matrix-order", default=None,
                        help="Optional file with one channel label per line; "
                             "reindexes the output matrix to this row order.")
    return parser


def _resolve_sample_names(vcf_paths, sample_names):
    if sample_names is None:
        return [_sample_name_from_path(p) for p in vcf_paths]
    if len(sample_names) != len(vcf_paths):
        raise ValueError(
            f"--sample-names has {len(sample_names)} entries but --vcf has {len(vcf_paths)}"
        )
    return sample_names


def main(argv=None):
    args = build_parser().parse_args(argv)

    channels = build_sbs96_channels()
    sample_names = _resolve_sample_names(args.vcf, args.sample_names)

    columns = {}
    for vcf_path, sample_name in zip(args.vcf, sample_names):
        counts = compute_sbs96(vcf_path, args.fasta)
        columns[sample_name] = [counts[c] for c in channels]

    matrix = pd.DataFrame(columns, index=pd.Index(channels, name="Type"))

    if args.matrix_order:
        matrix = matrix.reindex(load_matrix_order(args.matrix_order))

    for format_name in args.output_formats:
        output_path = os.path.join(args.output_dir, f"{args.output_prefix}_{format_name}.dat")
        save_in_format(matrix, output_path, format_name)
        print(f"Created {output_path}")


if __name__ == "__main__":
    main()
