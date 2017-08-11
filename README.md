# PacbioAssemble

This code takes L1 casette sequencing data and generates final polshed sequences.

Current PacBio FL-L1 Analysis Pipeline:

1. cp Raw file to Raw Files DIR
2. Align to human with blasr using command like:
- blasr -sa /path/to/phaseI.fasta.sa -sam -nproc 8 -clipping soft -out /path/to/out.sam /path/to/bas.h5 /path/to/phaseI.fasta
- Note: Fasta file is a 1000 bp cutout around the expected insertion site (500 bp on each side)
3. Align to LINE1 with blasr using command like:
- blasr -sa /path/to/LINE1.fa.sa -sam -nproc 8 -clipping soft -out /path/to/out.sam /path/to/bas.h5 /path/to/LINE1.fa
4. Run Pacbio.jar like:
- Pacbio.jar reference.fa aligned_blasr.bam aligned_LINE1.bam working_dir/ fofn.txt
- where fofn is simply a file with the name of the bas.h5 file in it (weird convention of ConsenseTools)
5. Align consensus sequences to LINE1 ref using bowtie like:
- bowtie2 --local -f -x /home/eugene.gardner/MELT/me_refs/line_ref/LINE1 -U consensus_seqs.fa -S aligned.sam
6. Manually curate in IGV

7. Get ORFS/Sequences/etc:
- Run LIPAC/FinishingScript/L1_analysis like:
- perl L1_analysis.pl <working_dir>/

- Additionally look at lengths of reads and distribution
- Use ‘_’ in the names for the small contigs. DO NOT USE ANYTHING ELSE.
