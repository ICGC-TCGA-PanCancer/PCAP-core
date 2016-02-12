#!/bin/bash

SOURCE_BWA="https://github.com/lh3/bwa/archive/0.7.12.tar.gz"
# for bio db sam
SOURCE_SAMTOOLS="https://github.com/samtools/samtools/archive/0.1.20.tar.gz"

# for bamstats
SOURCE_HTSLIB="https://github.com/samtools/htslib/archive/1.2.1.tar.gz"

# for bigwig
SOURCE_JKENT_BIN="https://github.com/ENCODE-DCC/kentUtils/raw/master/bin/linux.x86_64"

# for biobambam
SOURCE_BBB_BIN_DIST="https://github.com/gt1/biobambam2/releases/download/2.0.25-release-20151105154334/biobambam2-2.0.25-release-20151105154334-x86_64-etch-linux-gnu.tar.gz"

get_distro () {
  if hash curl 2>/dev/null; then
    curl -sS -o $1.tar.gz -L $2
  else
    wget -nv -O $1.tar.gz $2
  fi
  mkdir -p $1
  tar --strip-components 1 -C $1 -zxf $1.tar.gz
}

get_file () {
# output, source
  if hash curl 2>/dev/null; then
    curl -sS -o $1 -L $2
  else
    wget -nv -O $1 $2
  fi
}

if [ "$#" -ne "1" ] ; then
  echo "Please provide an installation path  such as /opt/ICGC"
  exit 0
fi

set -e

CPU=`cat /proc/cpuinfo | egrep "^processor" | wc -l`
echo "Max compilation CPUs set to $CPU"

INST_PATH=$1

# get current directory
INIT_DIR=`pwd`

# cleanup inst_path
mkdir -p $INST_PATH
cd $INST_PATH
INST_PATH=`pwd`
mkdir -p $INST_PATH/bin
cd $INIT_DIR

# make sure that build is self contained
unset PERL5LIB
ARCHNAME=`perl -e 'use Config; print $Config{archname};'`
PERLROOT=$INST_PATH/lib/perl5
export PERL5LIB="$PERLROOT"
export PATH="$INST_PATH/bin:$PATH"

#create a location to build dependencies
SETUP_DIR=$INIT_DIR/install_tmp
mkdir -p $SETUP_DIR

## grab cpanm and stick in final bin:
rm -f $INST_PATH/bin/cpanm
get_file $INST_PATH/bin/cpanm https://cpanmin.us/
chmod +x $INST_PATH/bin/cpanm
CPANM=`which cpanm`
echo $CPANM

if ! ( perl -MExtUtils::MakeMaker -e 1 >/dev/null 2>&1); then
    echo
    echo "WARNING: Your Perl installation does not seem to include a complete set of core modules.  Attempting to cope with this, but if installation fails please make sure that at least ExtUtils::MakeMaker is installed.  For most users, the best way to do this is to use your system's package manager: apt, yum, fink, homebrew, or similar."
fi

perlmods=( "File::ShareDir" "File::ShareDir::Install" "Const::Fast" )
for i in "${perlmods[@]}" ; do
  perl $CPANM --mirror http://cpan.metacpan.org -l $INST_PATH $i
done

# figure out the upgrade path
COMPILE=`echo 'nothing' | perl -I lib -MPCAP -ne 'print PCAP::upgrade_path($_);'`
if [ -e "$INST_PATH/lib/perl5/PCAP.pm" ]; then
  COMPILE=`perl -I $INST_PATH/lib/perl5 -MPCAP -e 'print PCAP->VERSION,"\n";' | perl -I lib -MPCAP -ne 'print PCAP::upgrade_path($_);'`
fi

#Need to add CaVEMan stuff here... will depend on samtools too (for now).

echo -n "Building jkentUtils ..."
if [ -e $SETUP_DIR/jkentUtils.success ]; then
  echo -n " previously installed ...";
else
  cd $SETUP_DIR
  if [[ `uname -m` == x86_64 ]] ; then
    get_file $INST_PATH/bin/wigToBigWig $SOURCE_JKENT_BIN/wigToBigWig
    chmod +x $INST_PATH/bin/wigToBigWig
    get_file $INST_PATH/bin/bigWigMerge $SOURCE_JKENT_BIN/bigWigMerge
    chmod +x $INST_PATH/bin/bigWigMerge
  else
    if [ ! -e $INST_DIR/bin/bigWigMerge ]; then
      echo "Binaries only available for x86_64, please install bigWigMerge from kentUtils: https://github.com/ENCODE-DCC/kentUtils"
      exit 1
    fi
    if [ ! -e $INST_DIR/bin/wigToBigWig ]; then
      echo "Binaries only available for x86_64, please install wigToBigWig from kentUtils: https://github.com/ENCODE-DCC/kentUtils"
      exit 1
    fi
  fi
fi
echo

echo -n "Building htslib ..."
if [ -e $SETUP_DIR/htslib.success ]; then
  echo -n " previously installed ...";
else
  cd $SETUP_DIR
  if [ ! -e htslib ]; then
    get_distro "htslib" $SOURCE_HTSLIB
  fi
  make -C htslib -j$CPU
  touch $SETUP_DIR/htslib.success
fi
echo

export HTSLIB="$SETUP_DIR/htslib"

echo -n "Building bam_stats ..."
if [ -e $SETUP_DIR/bam_stats.success ]; then
  echo -n " previously installed ...";
else
  cd $INIT_DIR
  make -C c clean
  make -C c -j$CPU
  cp bin/bam_stats $INST_PATH/bin/.
  touch $SETUP_DIR/bam_stats.success
  # need to clean up as will clash with other version
  rm -rf $SAMTOOLS
  make -C c clean
fi
echo

cd $SETUP_DIR
if [[ ",$COMPILE," == *,bwa,* ]] ; then
  echo -n "Building BWA ..."
  if [ -e $SETUP_DIR/bwa.success ]; then
    echo -n " previously installed ..."
  else
      get_distro "bwa" $SOURCE_BWA
      cd $SETUP_DIR/bwa
      make -j$CPU
      cp bwa $INST_PATH/bin/.
      touch $SETUP_DIR/bwa.success
  fi
  echo
else
  echo "BWA - No change between PCAP versions"
fi

if [[ ",$COMPILE," == *,biobambam,* ]] ; then
  echo -n "Building biobambam ..."
  if [ -e $SETUP_DIR/biobambam.success ]; then
    echo -n " previously installed ..."
  else
      cd $SETUP_DIR
      mkdir -p biobambam
      get_distro "biobambam" $SOURCE_BBB_BIN_DIST
      mkdir -p $INST_PATH/bin $INST_PATH/etc $INST_PATH/lib $INST_PATH/share
      rm -f biobambam/bin/curl # don't let this file in SSL doesn't work
      cp -r biobambam/bin/* $INST_PATH/bin/.
      cp -r biobambam/etc/* $INST_PATH/etc/.
      cp -r biobambam/lib/* $INST_PATH/lib/.
      cp -r biobambam/share/* $INST_PATH/share/.
      touch $SETUP_DIR/biobambam.success
  fi
  export PERL5LIB="$PERLROOT"
  echo
else
  echo "biobambam - No change between PCAP versions"
fi

cd $INIT_DIR

#perl $CPANM --mirror http://cpan.metacpan.org --notest -l $INST_PATH BioPerl

if [[ ",$COMPILE," == *,samtools,* ]] ; then
  echo -n "Building samtools v0.x ..."
  if [ -e $SETUP_DIR/samtools.success ]; then
    echo -n " previously installed ...";
  else
    cd $SETUP_DIR
    if [ ! -e samtools ]; then
      get_distro "samtools" $SOURCE_SAMTOOLS
      perl -i -pe 's/^CFLAGS=\s*/CFLAGS=-fPIC / unless /\b-fPIC\b/' samtools/Makefile
    fi
    make -C samtools -j$CPU
    cp samtools/samtools $INST_PATH/bin/.
    export SAMTOOLS="$SETUP_DIR/samtools"
    echo
    echo -n "Building Bio::DB::Sam..."
    cd $SETUP_DIR
    mkdir -p bioDbSam
    get_distro "bioDbSam" https://github.com/GMOD/GBrowse-Adaptors/archive/release-1_42.tar.gz
    cd $SETUP_DIR/bioDbSam/Bio-SamTools
    perl $CPANM --mirror http://cpan.metacpan.org --notest -l $INST_PATH .
    touch $SETUP_DIR/samtools.success
    echo
  fi
  echo
else
  echo "samtools - No change between PCAP versions"
fi

cd $INIT_DIR

echo -n "Installing Perl prerequisites ..."
perl $CPANM --mirror http://cpan.metacpan.org --notest -l $INST_PATH --installdeps .
echo

echo -n "Installing PCAP ..."
  cd $INIT_DIR
  perl $CPANM --mirror http://cpan.metacpan.org --notest -l $INST_PATH .
echo

# cleanup all junk
rm -rf $SETUP_DIR
rm -rf $INIT_DIR/bin/biobambam



echo
echo
echo "Please add the following to beginning of path:"
echo "  $INST_PATH/bin"
echo "Please add the following to beginning of PERL5LIB:"
echo "  $PERLROOT"
echo

exit 0
