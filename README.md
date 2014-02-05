ICGC-TCGA-PCAP
==============

NGS reference implementations and helper code for the ICGC/TCGA Pan-Cancer Analysis Project

This repository contains code to run genomic alignments of paired end data
and subsequent calling algorithms under a fire and forget model with
automatic resume of jobs from the last successfully processed step.

The intention is to provide reference implementations but also provide
simple executable wrappers that are useful for the scientific community
who may have little IT support.

Some elements of this repository are specific to the ICGC/TCGA PanCancer project
data submissions process.

###Programs

####bwa_aln.pl
Perform paired end alignment of data using bwa backtrack algorithm - Requires BWA 0.6.2 (final stable release for this method)

####bwa_mem.pl
Perform paired end alignment of data using bwa mem algorithm - Requires BWA 0.7.6a+

####diff_bams.pl
Compare two bam files ignoring differences in header information which have no effect on mapping result. @SQ number and order are expected to match.

###Utilities
####monitor.py
Utility script which can be used to monitor CPU, memory etc for any program, e.g.

    monitor.py bwa_aln.pl ...

Please be aware that all commands under this are also prefixed with:

    numactl --interleave=all

###PanCancer specific
####bam\_to\_sra\_sub.pl
Generate SRA XML required for GNOS submissions.
