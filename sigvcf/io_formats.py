"""Matrix I/O matching the conventions the rest of this ecosystem expects.

Output shape mirrors every other counts_*.dat producer in this project
(e.g. SigFitBench/scripts/prepare_real_samples.py): a "Type"-indexed,
tab- or comma-separated table, mutation channels as rows and samples
as columns, so a matrix built here plugs into the existing fitting
tool rules with no format-specific handling on their side.
"""

import os


def save_in_format(counts, output_path, format_name):
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    if format_name == "sigfit":
        counts.to_csv(output_path, sep=",", index=True, header=True)
    elif format_name == "sigLASSO":
        counts.to_csv(output_path, sep=",", index=True, header=True, quoting=1)
    else:
        counts.to_csv(output_path, sep="\t", index=True, header=True)


def load_matrix_order(path):
    """One channel label per line; tolerates the double-quoted style some
    existing matrix-order files use (e.g. '"A[C>A]A"')."""
    with open(path) as fh:
        return [line.strip().strip('"') for line in fh if line.strip()]
