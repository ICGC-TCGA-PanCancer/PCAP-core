#!/bin/bash

SOURCE_BWA="https://github.com/lh3/bwa/archive/0.7.12.tar.gz"

SOURCE_SAMTOOLS="https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2"

# for bamstats
SOURCE_HTSLIB="https://github.com/samtools/htslib/archive/1.3.1.tar.gz"

# for bigwig
SOURCE_JKENT_BIN="https://github.com/ENCODE-DCC/kentUtils/raw/master/bin/linux.x86_64"
# for Bio::DB::BigWig
SOURCE_KENTSRC="http://hgdownload.cse.ucsc.edu/admin/jksrc.zip"
# for fast merging of per-chr BW files
SOURCE_LIB_BW="https://github.com/dpryan79/libBigWig/archive/0.1.6.tar.gz"

# for biobambam
SOURCE_BBB_BIN_DIST="https://github.com/gt1/biobambam2/releases/download/2.0.50-release-20160705161609/biobambam2-2.0.50-release-20160705161609-x86_64-etch-linux-gnu.tar.gz"

BIODBHTS_INSTALL="https://raw.githubusercontent.com/Ensembl/Bio-HTS/master/INSTALL.pl"

get_distro () {
  EXT=""
  if [[ $2 == *.tar.bz2* ]] ; then
    EXT="tar.bz2"
  elif [[ $2 == *.zip* ]] ; then
    EXT="zip"
  elif [[ $2 == *.tar.gz* ]] ; then
    EXT="tar.gz"
  else
    echo "I don't understand the file type for $1"
    exit 1
  fi
  if hash curl 2>/dev/null; then
    curl -sS -o $1.$EXT -L $2
  else
    wget -nv -O $1.$EXT $2
  fi
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

CPU=`grep -c ^processor /proc/cpuinfo`
if [ $? -eq 0 ]; then
  if [ "$CPU" -gt "6" ]; then
    CPU=6
  fi
else
  CPU=1
fi
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

## grab cpanm and stick in workspace, then do a self upgrade into bin:
get_file $SETUP_DIR/cpanm https://cpanmin.us/
perl $SETUP_DIR/cpanm -l $INST_PATH App::cpanminus
CPANM=`which cpanm`
echo $CPANM

if ! ( perl -MExtUtils::MakeMaker -e 1 >/dev/null 2>&1); then
    echo
    echo "WARNING: Your Perl installation does not seem to include a complete set of core modules.  Attempting to cope with this, but if installation fails please make sure that at least ExtUtils::MakeMaker is installed.  For most users, the best way to do this is to use your system's package manager: apt, yum, fink, homebrew, or similar."
fi

perlmods=( "File::ShareDir" "File::ShareDir::Install" "Const::Fast" "File::Which" )
for i in "${perlmods[@]}" ; do
  $CPANM --mirror http://cpan.metacpan.org -l $INST_PATH $i
done

# figure out the upgrade path
COMPILE=`echo 'nothing' | perl -I lib -MPCAP -ne 'print PCAP::upgrade_path($_);'`
if [ -e "$INST_PATH/lib/perl5/PCAP.pm" ]; then
  COMPILE=`perl -I $INST_PATH/lib/perl5 -MPCAP -e 'print PCAP->VERSION,"\n";' | perl -I lib -MPCAP -ne 'print PCAP::upgrade_path($_);'`
fi
COMPILE=",$COMPILE,"

#Need to add CaVEMan stuff here... will depend on samtools too (for now).

echo -n "Building jkentUtils ..."
if [ -e $SETUP_DIR/jkentUtils.success ]; then
  echo -n " previously installed ...";
else
  cd $SETUP_DIR
  if [[ `uname -m` == x86_64 ]] ; then
    get_file $INST_PATH/bin/wigToBigWig $SOURCE_JKENT_BIN/wigToBigWig
    chmod +x $INST_PATH/bin/wigToBigWig
  else
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
  get_distro "htslib" $SOURCE_HTSLIB
  mkdir -p htslib
  tar --strip-components 1 -C htslib -zxf htslib.tar.gz
  make -C htslib -j$CPU
  touch $SETUP_DIR/htslib.success
fi
echo

echo -n "Building libBigWig ..."
if [ -e $SETUP_DIR/libBigWig.success ]; then
  echo -n " previously installed ...";
else
  cd $SETUP_DIR
  get_distro "libBigWig" $SOURCE_LIB_BW
  mkdir -p libBigWig
  tar --strip-components 1 -C libBigWig -zxf libBigWig.tar.gz
  make -C libBigWig -j$CPU install prefix=$INST_PATH
  rm -f $INST_PATH/lib/libBigWig.so
  touch $SETUP_DIR/libBigWig.success
fi
echo

export HTSLIB="$SETUP_DIR/htslib"

echo -n "Building bam_stats ..."
if [ -e $SETUP_DIR/bam_stats.success ]; then
  echo -n " previously installed ...";
else
  cd $INIT_DIR
  make -C c clean
  make -C c -j$CPU prefix=$INST_PATH
  cp bin/bam_stats $INST_PATH/bin/.
  cp bin/bwcat $INST_PATH/bin/.
  cp bin/reheadSQ $INST_PATH/bin/.
  cp bin/bam2bedgraph $INST_PATH/bin/.
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
      mkdir -p bwa
      tar --strip-components 1 -C bwa -zxf bwa.tar.gz
      make -C bwa -j$CPU
      cp bwa/bwa $INST_PATH/bin/.
      touch $SETUP_DIR/bwa.success
  fi
  echo
else
  echo "BWA - No change between PCAP versions"
fi

if [[ ",$COMPILE," == *,biobambam,* ]] ; then
  echo -n "Building biobambam2 ..."
  if [ -e $SETUP_DIR/biobambam2.success ]; then
    echo " previously installed2 ..."
  else
    cd $SETUP_DIR
    get_distro "biobambam2" $SOURCE_BBB_BIN_DIST
    mkdir -p biobambam2
    tar --strip-components 1 -C biobambam2 -zxf biobambam2.tar.gz
    mkdir -p $INST_PATH/bin $INST_PATH/etc $INST_PATH/lib $INST_PATH/share
    rm -f biobambam2/bin/curl # don't let this file in SSL doesn't work
    cp -r biobambam2/bin/* $INST_PATH/bin/.
    cp -r biobambam2/etc/* $INST_PATH/etc/.
    cp -r biobambam2/lib/* $INST_PATH/lib/.
    cp -r biobambam2/share/* $INST_PATH/share/.
    touch $SETUP_DIR/biobambam2.success
    echo
  fi
else
  echo "biobambam - No change between PCAP versions"
fi

cd $INIT_DIR

if [[ ",$COMPILE," == *,samtools,* ]] ; then
  echo -n "Building samtools ..."
  if [ -e $SETUP_DIR/samtools.success ]; then
    echo -n " previously installed ...";
  else
    cd $SETUP_DIR
    get_distro "samtools" $SOURCE_SAMTOOLS
    mkdir -p samtools
    tar --strip-components 1 -C samtools -xjf samtools.tar.bz2
    cd samtools
    ./configure --enable-plugins --enable-libcurl --prefix=$INST_PATH
    make all all-htslib
    make install install-htslib
    touch $SETUP_DIR/samtools.success
  fi
  echo
else
  echo "samtools - No change between PCAP versions"
fi

cd $INIT_DIR

if [[ ",$COMPILE," == *,samtools,* ]] ; then
  echo -n "Building Bio::DB::HTS ..."
  if [ -e $SETUP_DIR/biohts.success ]; then
    echo " previously installed ...";
  else
    cd $SETUP_DIR
    $CPANM --mirror http://cpan.metacpan.org --notest -l $INST_PATH Module::Build Bio::Root::Version
    # now Bio::DB::HTS
    get_file "INSTALL.pl" $BIODBHTS_INSTALL
    perl -I $PERL5LIB INSTALL.pl --prefix $INST_PATH --static
    rm -f INSTALL.pl
    touch $SETUP_DIR/biohts.success
  fi
  echo
else
  echo "Bio::DB::HTS - No change between PCAP versions" # based on samtools tag
fi

echo -n "Building kentsrc + Bio::DB::BigFile ..."
if [ -e $SETUP_DIR/kentsrc.success ]; then
  echo " previously installed ...";
else
  cd $SETUP_DIR
  get_distro "kentsrc" $SOURCE_KENTSRC
  unzip -q kentsrc.zip
  perl -pi -e 's/(\s+CFLAGS=)$/${1}-fPIC/' kent/src/inc/common.mk
  cd kent/src/lib
  export MACHTYPE=i686    # for a 64-bit system
  make
  cd ../
  export KENT_SRC=`pwd`
  cd $SETUP_DIR
  $CPANM --mirror http://cpan.metacpan.org -l $INST_PATH Bio::DB::BigFile
fi

cd $INIT_DIR

echo -n "Installing Perl prerequisites ..."
$CPANM --mirror http://cpan.metacpan.org --notest -l $INST_PATH --installdeps .
echo

echo -n "Installing PCAP ..."
$CPANM --mirror http://cpan.metacpan.org -l $INST_PATH .
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
