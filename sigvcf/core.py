"""VCF streaming and trinucleotide-context primitives for SBS96."""

from collections import Counter
from typing import NamedTuple

import pysam

BASES = "ACGT"
_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")
PYRIMIDINES = frozenset("CT")

# Why a mutation in a VCF was left out of the SBS96 counts.
# The first two are expected (they are not single-base substitutions).
# The last three point to a mismatch between the VCF and the genome FASTA.
NOT_SNV = "not_snv"                      # indel, MNV, non-ACGT allele, or ALT same as REF
NEXT_TO_SNV = "next_to_another_snv"      # part of a DBS/MNV: SNVs at back-to-back positions
CONTIG_NOT_IN_GENOME = "contig_not_in_genome"
AT_CONTIG_EDGE = "at_contig_edge"        # no base on one side, so no 3-base context
REF_MISMATCH = "ref_mismatch"            # REF in the VCF is not the base in the FASTA

EXPECTED_SKIPS = (NOT_SNV, NEXT_TO_SNV)
GENOME_MISMATCH_SKIPS = (CONTIG_NOT_IN_GENOME, AT_CONTIG_EDGE, REF_MISMATCH)
SKIP_REASONS = EXPECTED_SKIPS + GENOME_MISMATCH_SKIPS


class SigVCFInputError(Exception):
    """A VCF that SigVCF cannot turn into correct counts."""


class SkipReport:
    """Counts, for one VCF, the mutations used and the mutations left out (by reason)."""

    def __init__(self):
        self.counted = 0
        self.skipped = Counter(dict.fromkeys(SKIP_REASONS, 0))
        self.missing_contigs = Counter()  # contig name -> SNVs on it

    def genome_mismatches(self):
        return sum(self.skipped[r] for r in GENOME_MISMATCH_SKIPS)

    def genome_mismatch_fraction(self):
        """Share of the SNVs checked against the FASTA that did not fit it."""
        checked = self.counted + self.genome_mismatches()
        return self.genome_mismatches() / checked if checked else 0.0


def revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


class SNV(NamedTuple):
    chrom: str
    pos: int  # 1-based, VCF convention
    ref: str
    alt: str


def stream_snvs(vcf_path, report=None):
    """Yield single-base substitutions from a VCF, one record at a time.

    Uses ``pysam.VariantFile``'s htslib-backed iterator, so a
    multi-gigabyte VCF is never held in memory at once. Indels, MNVs,
    non-ACGT alleles, and non-variant ALT==REF records are skipped;
    multiallelic records yield one SNV per valid single-base ALT.

    SigVCF does not read genotypes: every record counts for the file's one
    sample. A VCF with several sample columns would mix their mutations
    into one column, so it is refused.
    """
    report = report if report is not None else SkipReport()
    with pysam.VariantFile(vcf_path) as vf:
        samples = list(vf.header.samples)
        if len(samples) > 1:
            raise SigVCFInputError(
                f"{vcf_path} has {len(samples)} sample columns ({', '.join(samples[:5])}"
                f"{', ...' if len(samples) > 5 else ''}). SigVCF counts every record "
                "for one sample, so give it one VCF per sample (for example with "
                "'bcftools view -s <sample> -c 1')."
            )
        for record in vf:
            ref = record.ref
            alts = record.alts or ()
            if not ref or len(ref) != 1 or ref not in BASES:
                report.skipped[NOT_SNV] += max(len(alts), 1)
                continue
            if not alts:
                report.skipped[NOT_SNV] += 1
                continue
            for alt in alts:
                if len(alt) == 1 and alt in BASES and alt != ref:
                    yield SNV(record.chrom, record.pos, ref, alt)
                else:
                    report.skipped[NOT_SNV] += 1


def stream_isolated_snvs(vcf_path, report=None):
    """Yield only SNVs with no SNV at the adjacent position on either side.

    Two or more SNVs at consecutive positions are a dinucleotide (DBS) or
    multi-nucleotide (MNV) event, not independent single-base substitutions.
    MutationalPatterns and SigProfilerMatrixGenerator both keep them out of
    the SBS96 matrix, so counting them here would inflate SBS96 by one
    mutation per base of every DBS/MNV. Uses a small same-chromosome run
    buffer, so memory stays O(longest run). Assumes the VCF is position-
    sorted (standard for VCF).
    """
    report = report if report is not None else SkipReport()
    run = []
    for snv in stream_snvs(vcf_path, report):
        if run and run[-1].chrom == snv.chrom and snv.pos == run[-1].pos + 1:
            run.append(snv)
        else:
            if len(run) == 1:
                yield run[0]
            else:
                report.skipped[NEXT_TO_SNV] += len(run)
            run = [snv]
    if len(run) == 1:
        yield run[0]
    else:
        report.skipped[NEXT_TO_SNV] += len(run)


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


def iter_sbs96_events(vcf_path, fasta_path, report=None):
    """Stream every isolated SNV with its pyrimidine-folded trinucleotide channel.

    If ``report`` (a SkipReport) is given, it is filled in as the VCF is read.
    """
    report = report if report is not None else SkipReport()
    with pysam.FastaFile(fasta_path) as fasta:
        for snv in stream_isolated_snvs(vcf_path, report):
            try:
                context = fasta.fetch(snv.chrom, snv.pos - 2, snv.pos + 1).upper()
            except KeyError:
                report.skipped[CONTIG_NOT_IN_GENOME] += 1
                report.missing_contigs[snv.chrom] += 1
                continue
            except ValueError:
                report.skipped[AT_CONTIG_EDGE] += 1  # the 3-base window starts before position 1
                continue
            if len(context) != 3:
                report.skipped[AT_CONTIG_EDGE] += 1  # the 3-base window runs past the contig end
                continue
            if context[1] != snv.ref:
                report.skipped[REF_MISMATCH] += 1
                continue
            pyr_ref, pyr_alt, pyr_context = normalize_snv_to_pyrimidine(snv.ref, snv.alt, context)
            channel = f"{pyr_context[0]}[{pyr_ref}>{pyr_alt}]{pyr_context[2]}"
            report.counted += 1
            yield SBS96Event(snv, channel, pyr_ref == snv.ref)
