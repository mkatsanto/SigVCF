"""SigVCF: SBS96 mutation-count matrix extraction directly from VCFs.

Streams the VCF and an indexed FASTA (no genome or VCF held in memory)
and runs fully offline, so it can ship inside an Apptainer/Singularity
image on an air-gapped cluster.
"""

__version__ = "0.1.0"
