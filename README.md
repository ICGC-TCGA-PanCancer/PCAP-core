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

---

###Dependencies
Some of the code included in this package has dependencies on [biobambam](https://github.com/gt1/biobambam).  This also requires the library [libmaus](https://github.com/gt1/libmaus).  Please ensure that you compile libmaus with **snappy** and **io_lib** to enable the relevant features.  These are described on the libmaus page.

Please see the `INSTALL` file for full list of dependencies.

---

###Programs
####bwa_aln.pl
Perform paired end alignment of data using bwa backtrack algorithm - Requires [BWA 0.6.2](https://github.com/lh3/bwa/archive/0.6.2.tar.gz) (final stable release for this method).

####bwa_mem.pl
Perform paired end alignment of data using bwa mem algorithm - Requires BWA 0.7.6a+.
Please see the [BWA releases](https://github.com/lh3/bwa/releases) for the current stable release.

###PanCancer specific
####bam\_to\_sra\_sub.pl
Generate SRA XML required for GNOS submissions.

###Utilities
####diff_bams.pl
Compare two bam files ignoring differences in header information which have no effect on mapping result. @SQ number and order are expected to match.

####monitor.py
Utility script which can be used to monitor CPU, memory etc for any program, e.g.

    monitor.py bwa_aln.pl ...

Please be aware that all commands under this are also prefixed with:

    numactl --interleave=all

---

##Creating a release
####Preparation
* Commit all relevant changes.
* Pull a clean version of the repo and use this for the following steps.

####Cutting the release
1. Update `lib/PCAP.pm` to the correct version (adding rc/beta to end if applicable).
2. Run `./prerelease.sh`
3. Check all tests and coverage reports are acceptable.
4. Commit the updated docs tree and updated module/version.
5. Push commits.
6. Use the GitHub tools to draft a release.
