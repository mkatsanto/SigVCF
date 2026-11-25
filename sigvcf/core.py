"""VCF streaming and trinucleotide-context primitives for SBS96."""

from typing import NamedTuple

import pysam

BASES = "ACGT"
_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")
PYRIMIDINES = frozenset("CT")


def revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


class SNV(NamedTuple):
    chrom: str
    pos: int  # 1-based, VCF convention
    ref: str
    alt: str


def stream_snvs(vcf_path):
    """Yield single-base substitutions from a VCF, one record at a time.

    Uses ``pysam.VariantFile``'s htslib-backed iterator, so a
    multi-gigabyte VCF is never held in memory at once. Indels, MNVs,
    non-ACGT alleles, and non-variant ALT==REF records are skipped;
    multiallelic records yield one SNV per valid single-base ALT.
    """
    with pysam.VariantFile(vcf_path) as vf:
        for record in vf:
            ref = record.ref
            if not ref or len(ref) != 1 or ref not in BASES:
                continue
            for alt in record.alts or ():
                if len(alt) == 1 and alt in BASES and alt != ref:
                    yield SNV(record.chrom, record.pos, ref, alt)


def stream_isolated_snvs(vcf_path):
    """Yield only SNVs with no SNV at the adjacent position on either side.

    Two or more SNVs at consecutive positions are a dinucleotide (DBS) or
    multi-nucleotide (MNV) event, not independent single-base substitutions.
    MutationalPatterns and SigProfilerMatrixGenerator both keep them out of
    the SBS96 matrix, so counting them here would inflate SBS96 by one
    mutation per base of every DBS/MNV. Uses a small same-chromosome run
    buffer, so memory stays O(longest run). Assumes the VCF is position-
    sorted (standard for VCF).
    """
    run = []
    for snv in stream_snvs(vcf_path):
        if run and run[-1].chrom == snv.chrom and snv.pos == run[-1].pos + 1:
            run.append(snv)
        else:
            if len(run) == 1:
                yield run[0]
            run = [snv]
    if len(run) == 1:
        yield run[0]


def normalize_snv_to_pyrimidine(ref, alt, trinucleotide_context):
    """Fold a substitution + its 3bp context onto the pyrimidine (C/T) reference strand.

    Returns (pyr_ref, pyr_alt, pyr_context). If ``ref`` is already a
    pyrimidine, returns the inputs unchanged; otherwise reverse-
    complements all three.
    """
    if ref in PYRIMIDINES:
        return ref, alt, trinucleotide_context
    return revcomp(ref), revcomp(alt), revcomp(trinucleotide_context)


class SBS96Event(NamedTuple):
    snv: SNV
    channel: str  # e.g. "A[C>A]A"
    on_plus_strand: bool  # False -> the pyrimidine fold needed a revcomp to reach it


def iter_sbs96_events(vcf_path, fasta_path):
    """Stream every isolated SNV with its pyrimidine-folded trinucleotide channel."""
    with pysam.FastaFile(fasta_path) as fasta:
        for snv in stream_isolated_snvs(vcf_path):
            try:
                context = fasta.fetch(snv.chrom, snv.pos - 2, snv.pos + 1).upper()
            except KeyError:
                continue  # snv.chrom isn't a contig in this FASTA
            except ValueError:
                continue  # snv.pos is 1 or 2: the 3bp window starts before position 1
            if len(context) != 3 or context[1] != snv.ref:
                continue  # short read past a contig end, or REF doesn't match the FASTA
            pyr_ref, pyr_alt, pyr_context = normalize_snv_to_pyrimidine(snv.ref, snv.alt, context)
            channel = f"{pyr_context[0]}[{pyr_ref}>{pyr_alt}]{pyr_context[2]}"
            yield SBS96Event(snv, channel, pyr_ref == snv.ref)
