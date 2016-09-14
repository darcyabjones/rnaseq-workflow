# rnaseq-workflow
A makefile that takes adapter trimmed fastq files from an RNA-seq experiment and aligns them to a genomes, performs basic QC, and counts transcript abundance.

Currently, the alignment tool implemented is HISAT2.
I would like to use STAR as well for comparison of results.

BAM QC is performed using RSeQC and MultiQC.
I might look at some tools in Picard later.

Fragment counts for later DGE is performed using HTSeq-count.
I will also implement StringTie and Cufflinks later.
