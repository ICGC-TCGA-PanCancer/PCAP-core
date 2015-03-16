ICGC-TCGA-PCAP
==============

NGS reference implementations and helper code for the ICGC/TCGA Pan-Cancer Analysis Project

This repository contains code to run genomic alignments of paired end data
and subsequent calling algorithms.

The intention is to provide reference implementations and simple to execute wrappers
that are useful for the scientific community who may have little IT/bioinformatic support.

Please see the [wiki](https://github.com/ICGC-TCGA-PanCancer/PCAP-core/wiki) for further details.

---

###Dependencies/Install
Some of the code included in this package has dependencies on external packages:

 * [biobambam](https://github.com/gt1/biobambam)
 * [bwa](https://github.com/lh3/bwa)
 * [samtools](https://github.com/samtools/samtools)
 * [kentUtils](https://github.com/ENCODE-DCC/kentUtils)

And various perl modules.

Please see the respective licence for each before use.

Please use `setup.sh` to install the dependencies.  Please be aware that this expects basic C
compilation libraries and tools to be available, most are listed in `INSTALL`.

---

###Programs

Please see the [wiki](https://github.com/ICGC-TCGA-PanCancer/PCAP-core/wiki) for details of programs.

---

##Creating a release
####Preparation
* Commit/push all relevant changes.
* Pull a clean version of the repo and use this for the following steps.

####Cutting the release
1. Update `lib/PCAP.pm` to the correct version.
2. Update `c/Makefile` to the correct version.
3. Ensure upgrade path for new version number is added to `lib/PCAP.pm`.
4. Update `CHANGES.md` to show major items.
5. Run `./prerelease.sh`
6. Check all tests and coverage reports are acceptable.
7. Commit the updated docs tree and updated module/version.
8. Push commits.
9. Use the GitHub tools to draft a release.
