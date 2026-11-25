# SigVCF

Extracts an SBS96 mutation-count matrix directly from VCFs -- no
upstream matrix-generation step needed. (SBS192 and DBS78 were removed;
only SBS96 is supported for now.)

Everything streams: SNVs are read one VCF record at a time via
`pysam.VariantFile`, and trinucleotide context is pulled from an
indexed FASTA (`.fa` + `.fai`) via `pysam.FastaFile`, so no genome is
held in RAM. No network access and no bundled reference data -- built
to run inside an offline Apptainer/Singularity container.

## Install

```
pip install -e .
```

## Usage

```
sigvcf --vcf sample1.vcf sample2.vcf.gz \
    --fasta /ref/GRCh38.fa \
    --output-dir out/ \
    --output-prefix counts_myrun \
    --output-formats tab_with_index sigfit
```

`--matrix-type sbs96` is accepted (and is the default) so existing
callers keep working.

Output is written as `{output-prefix}_{format}.dat` per requested
`--output-formats` entry: a tab- or comma-separated table, the 96
channels as rows (index name `Type`, labels like `A[C>A]A`, COSMIC
order) and one column per `--vcf` (named from `--sample-names` or each
file's stem). This is the shape [SigFitBench](../SigFitBench)'s
fitting-tool rules already consume (see
`SigFitBench/workflow/rules/vcf_input.smk`).

Pass `--matrix-order some_order.dat` (one channel label per line) to
reindex the output rows to match a specific reference-signature file.

## Behaviour

- Only single-base substitutions are counted; indels, MNVs and non-ACGT
  alleles are skipped. SNVs at consecutive positions (a DBS, or an MNV of
  3+) are excluded from SBS96, as in MutationalPatterns and
  SigProfilerMatrixGenerator. The VCF must be position-sorted. Multiallelic records count once per valid ALT.
- Substitutions are folded onto the pyrimidine (C/T) strand.
- A record is skipped if its contig is not in the FASTA, it sits within
  1 bp of a contig start/end, or its REF disagrees with the FASTA.

## Credits and references

SigVCF is a new, independent implementation. It does not include or
import code from the tools below. It follows their published method for
building an SBS96 matrix so that its output can be used in their place.
If you use SigVCF, please also cite the work it is based on.

**The SBS96 categories** (the 96 types, their labels such as `A[C>A]A`,
their order, and recording each substitution on the C/T strand) come
from the COSMIC mutational signatures work:

- Alexandrov LB, et al. Signatures of mutational processes in human
  cancer. *Nature* 500, 415–421 (2013).
  https://doi.org/10.1038/nature12477
- Alexandrov LB, et al. The repertoire of mutational signatures in human
  cancer. *Nature* 578, 94–101 (2020).
  https://doi.org/10.1038/s41586-020-1943-3

**The tools whose counting rules SigVCF follows** (for example, leaving
out SNVs at back-to-back positions), and whose output it was checked
against:

- SigProfilerMatrixGenerator: Bergstrom EN, et al.
  SigProfilerMatrixGenerator: a tool for visualizing and exploring
  patterns of small mutational events. *BMC Genomics* 20, 685 (2019).
  https://doi.org/10.1186/s12864-019-6041-2
- MutationalPatterns: Blokzijl F, Janssen R, van Boxtel R, Cuppen E.
  MutationalPatterns: comprehensive genome-wide analysis of mutational
  processes. *Genome Medicine* 10, 33 (2018).
  https://doi.org/10.1186/s13073-018-0539-0
- MutationalPatterns (current version): Manders F, et al.
  MutationalPatterns: the one stop shop for the analysis of mutational
  processes. *BMC Genomics* 23, 134 (2022).
  https://doi.org/10.1186/s12864-022-08357-3
