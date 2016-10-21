### 3.3.0
* adds `map_threads|mt` option to `bwa_mem.pl` to allow more control of parallel processing in one shot submission.
* adds `bwa_pl|l` option to `bwa_mem.pl` to allow preload of different malloc libraries.

### 3.2.0
* Move from legacy kent bigwig manipulation code and to [`cgpBigWig`](https://github.com/cancerit/cgpBigWig)
  * Faster and handles the huge number of contigs in many new reference builds.
  * Resulting changes to underlying installed tools is `bwcat` now `bwjoin` to be more descriptive of actual function.

### 3.1.1
Handle recent changes to BioPerl structure

### 3.1.0
* Use BWA default for `-T` previously hard coded to `-T 0`.
  *  Can be passed through `bwa_mem.pl` other ags to bwa via the `-b` option.
*  Fix `bam2bedgraph` compilation since changes to underlying libraries
*  `bamToBw.pl` - expose read flag filters
*  Drop dependancy on `Bio::DB::HTS INSTALL.pl` as can't fix to known good version.
*  Added travis CI
*  Add support for output directly to CRAM
*  `bwa` version upgraded to [0.7.15](https://github.com/lh3/bwa/releases/tag/v0.7.15)

### 3.0.0
* Threading module now converts currently running step to bash script for following reasons:
  * Changes logging to use file redirects instead of Capture::Tiny - prevent log bleed into wrong files
  * Commands for failed jobs remain after shutdown for easy debug/testing
* Log and progress file names simplified so more portable.

### 2.4.1
* Modified reheadSQ to be more robust.

### 2.1.0
* Adds `xam_coverage_bins.pl` which calculates fraction of targets covered at various depths (BAM/CRAM), using BED/GFF3 as target bait file.

### 2.0.0
* bwa_mem.pl
    * allow user to specify BWA mapping parameters
    * now accepts CRAM as input
* bamToBw.pl - now accepts CRAM as input.
* bam_stats - Adds 2 new stats:
    * `#_mapped_pairs`
    * `#_inter_chr_pairs`
* Dependancy changes
    * __WARNING__: ensure all related tools handle these updates
    * [samtools](https://github.com/samtools/samtools/releases), now only uses htslib based versions (1.3+, handling deprecated use of sort)
    * [Bio::DB::HTS](http://search.cpan.org/~rishidev/Bio-DB-HTS/lib/Bio/DB/HTS.pm) htslib bindings replacing Bio::DB::Sam

### 1.14.0
* bwa_mem.pl - Option to disable duplicate marking

### 1.13.0
* bam_stats - Unit tests for C code
* bam_stats - Fix to median insert size calculation

### 1.11.0
* bam_stats - new rna switch to give more appropriate insert size stats
* bam_stats - more robust handling of optional RG header entries
* bam_stats - allows streaming IO (thanks to @jenniferliddle)
* bwa_mem.pl - Handle `'` in RG header line/IDs
* Generally improved version handling and updated versions of some tools.

### 1.9.1
* Changed final log folder to include sample name and analysis type, prevents clash when lots of data to same output loc.

### 1.9.0
* Fix bugs #52 and #53
* Modified bwa_mem.pl to accept multi-readgroup BAM as input

### 1.7.1
* Turns out BWA mem still requires fixmates to get proper isize distributions
* bumped biobambam to [0.0.191](https://github.com/gt1/biobambam/releases/tag/0.0.191-release-20150401083643)

### 1.7.0
* Switched to bam_stats C in ``bwa_mem.pl``.
* Updates to ``bam_to_sra.pl`` to prevent bad SM values in unaligned BAM uploads.

### 1.6.0
* Adding local file mode for sites that cannot download from GNOS when the xml_to_bas.pl script runs
* [gnos_pull.pl](https://github.com/ICGC-TCGA-PanCancer/PCAP-core/wiki/Scripts-PanCancer-Specific#gnos_pullpl) - see linked docs

### 1.5.4
* bam_stats C
    * Reference file parameter is now optional to replicate bam_stats.pl functionality.
    * Warnings in help, and when a cram file is given as reference from header may not be found, and bam stats will fail.

### 1.5.3
* bam_stats C - changed array for khash in insert size calculations in order to make code more robust.
* Header RG line reading now reads anything not a tab or newline as it should when determining what the values of tags are.

### 1.5.2
bamToBw.pl fixes
* Pull actual binaries from jkent_util not html page associated
* process name corrections in bamToBw.pm command line args

### 1.5.1
* bam_stats c now has CRAM support.
* Also dropped need for samtools v1.x api as can be handled by htslib on it's own.

### 1.5.0
* bamToBw.pl and new biobambam dep

No changes to old tools, just additions and prep for handling CRAM input.

### 1.4.0
bam_stats in C, less than 2 hours to generate stats on a sample level BAM file of ~120GB.

### 1.3.0
* bam_stats.pl is now multi-threaded, can get ~50% runtime reduction with 3-4 threads, memory still <500MB.
* Upgrades biobambam to 0.0.185 (and dependencies).

### 1.2.3
xml_to_bas.pl - detect readgroup id clashes and attempt to reconcile, #54

### 1.2.2
Fixed bug in bwa_mem.pl when using '-f' option on paired fastq.

### 1.2.1
Makes xml_to_bas.pl more robust on AWS. Retrieved XML was being truncated on some network configurations.

### 1.2.0
Modifications made to the bwa_mem.pl code to split a lane of data into fragments to reduce failure recovery time. Primarily added to handle X10 data better.

Also updated samtools to 0.1.20, last version that is currently compatible with Bio::DB::Sam.

### 1.1.1
Fix missing dependancy and build a relocatable version of biobambam suitable for use in artifactory.

### 1.1.0
* Minor enhancement to bwa_mem.pl to automatically generate the *.bas file.
* Added xml_to_bas.pl for pancancer users, see the wiki for details.
* Fixed a few minor issues, #36, #37, #39

### 1.0.4
* Install biobambam 0.0.142 to prevent over-counting of duplicates when multiple
libraries, also required libmaus 0.0.124.
* Improve install for those working with multiple perl installs.
* Improve version inheritance, less code

### 1.0.3
* Corrected issue from dynamic de-reference of hash, issue for pre 5.14 perl
and potentially unstable in future.
* Added missing project code to cv terms.
* Bug-fixed upgrade path, still needs better solution.
* Cleaned up messaging in Threaded module.

### 1.0.2
* Upgrade install to pull biobambam 0.0.138
    * fastqtobam option 'pairedfile' for where readnames don't have trailing '/1' or '/2'.
    * fastqtobam option to relax qscore validation without turning off... careful
* Upgrade install to pull BWA 0.7.8
    * performance improvements for short read alignment (100bp)

### 1.0.1
* Upgrade install to pull biobambam 0.0.135
    * fastqtobam supports Casava v1.8
    * bamsort supports NM/MD correction during sam->bam/merge process
* Minor enhancement to BAS reader module.
* Sample name from command line passed through to SM of RG header in bwa_mem.pl
* SRA.pm - check that rg id is unique within run of code (thanks to Junjun Zhang)
* Threads.pm - join interval is now configurable.

### 1.0.0
* bam_stats.pl actually installed now.
* Basic *.bas perl access module.
* Upgraded libmaus/biobambam to resolve patch and CentOS install issue.
* Reference implementations ensure unique RG:ID between files.

### 0.3.0
* Changes for the re-worked PanCancer submission SOP.
* Patch for libmaus issue as not going to be a release in time.
* Bug fix for *.info files (bam_to_sra_sub.pl).
* Added bam_stats.pl.
* Project is now defaulted when not provided (bam_to_sra_sub.pl).

### 0.2.0 14-Mar-2014
* Updated biobambam version
* Documented additional dependencies
* Improved install implementation

### 0.0.2-beta2 04-Feb-2014
* Updated module naming in preparation for publication to GitHub.
* Added license boiler plate.
* bam_to_sra_sub.pl generates valid XML for GNOS, some features disabled until modifications to GNOS can be made (warnings indicate this on execution)

### 0.0.2-beta 29-Jan-2014
* Pre release with basic SRA XML generation (GENOS)
* Updated requirements for biobambam of 0.0.120
* Tests update to reflect change in biobambam requirement

### 0.0.1 11-Sep-2013
* Initial codebase for PanCancer alignment with BWA 0.6.2.
