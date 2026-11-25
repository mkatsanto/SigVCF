"""SBS96: single-base substitutions in trinucleotide context.

Row order and channel labels (e.g. "A[C>A]A") match the COSMIC/
SigProfilerMatrixGenerator/MutationalPatterns convention, so a matrix
built here lines up 1:1 with theirs.
"""

from collections import Counter

from .core import BASES, iter_sbs96_events

_PYR_SUBS = [("C", "A"), ("C", "G"), ("C", "T"), ("T", "A"), ("T", "C"), ("T", "G")]


def build_sbs96_channels():
    """The 96 canonical channel labels, e.g. "A[C>A]A", in COSMIC row order."""
    return [
        f"{five}[{ref}>{alt}]{three}"
        for ref, alt in _PYR_SUBS
        for five in BASES
        for three in BASES
    ]


def compute_sbs96(vcf_path, fasta_path):
    """Return a Counter of channel label -> mutation count for one VCF."""
    counts = Counter(dict.fromkeys(build_sbs96_channels(), 0))
    for event in iter_sbs96_events(vcf_path, fasta_path):
        if event.channel in counts:
            counts[event.channel] += 1
    return counts
