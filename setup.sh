#!/bin/bash

SOURCE_BWA="https://github.com/lh3/bwa/archive/v0.7.15.tar.gz"

SOURCE_SAMTOOLS="https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2"

# for bamstats and Bio::DB::HTS
SOURCE_HTSLIB="https://github.com/samtools/htslib/releases/download/1.3.1/htslib-1.3.1.tar.bz2"

# Bio::DB::HTS
SOURCE_BIOBDHTS="https://github.com/Ensembl/Bio-HTS/archive/2.3.tar.gz"

# for Bio::DB::BigWig
SOURCE_KENTSRC="ftp://ftp.sanger.ac.uk/pub/cancer/legacy-dependancies/jksrc.v334.zip"
SOURCE_BIGFILE="http://www.cpan.org/authors/id/L/LD/LDS/Bio-BigFile-1.07.tar.gz"

# for biobambam
SOURCE_BBB_BIN_DIST="https://github.com/gt1/biobambam2/releases/download/2.0.54-release-20160802163650/biobambam2-2.0.54-release-20160802163650-x86_64-etch-linux-gnu.tar.gz"

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
  rm -f $1.$EXT
  if hash curl 2>/dev/null; then
    curl --retry 10 -sS -o $1.$EXT -L $2
  else
    wget --tries=10 -nv -O $1.$EXT $2
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

# check for cgpBigWig
if [ -e "$INST_PATH/bin/bwjoin" ]; then
  echo -e "\tcgpBigWig installation found";
else
  echo -e "\tERROR: Please see README.md and install cgpBigWig";
  exit 1
fi

## grab cpanm and stick in workspace, then do a self upgrade into bin:
get_file $SETUP_DIR/cpanm https://cpanmin.us/
perl $SETUP_DIR/cpanm -l $INST_PATH App::cpanminus
CPANM=`which cpanm`
echo $CPANM

if ! ( perl -MExtUtils::MakeMaker -e 1 >/dev/null 2>&1); then
    echo
    echo "WARNING: Your Perl installation does not seem to include a complete set of core modules.  Attempting to cope with this, but if installation fails please make sure that at least ExtUtils::MakeMaker is installed.  For most users, the best way to do this is to use your system's package manager: apt, yum, fink, homebrew, or similar."
fi

if [ -e $SETUP_DIR/basePerlDeps.success ]; then
  echo "Previously installed base perl deps..."
else
  perlmods=( "ExtUtils::CBuilder" "Module::Build~0.42" "File::ShareDir" "File::ShareDir::Install" "Const::Fast" "File::Which" "LWP::UserAgent" "Bio::Root::Version@1.006924")
  for i in "${perlmods[@]}" ; do
    $CPANM --no-interactive --notest --mirror http://cpan.metacpan.org -l $INST_PATH $i
  done
  touch $SETUP_DIR/basePerlDeps.success
fi


# figure out the upgrade path
COMPILE=`echo 'nothing' | perl -I lib -MPCAP -ne 'print PCAP::upgrade_path($_);'`
if [ -e "$INST_PATH/lib/perl5/PCAP.pm" ]; then
  COMPILE=`perl -I $INST_PATH/lib/perl5 -MPCAP -e 'print PCAP->VERSION,"\n";' | perl -I lib -MPCAP -ne 'print PCAP::upgrade_path($_);'`
fi
COMPILE=",$COMPILE,"

echo -n "Get htslib ..."
if [ -e $SETUP_DIR/htslibGet.success ]; then
  echo " already staged ...";
else
  echo
  cd $SETUP_DIR
  get_distro "htslib" $SOURCE_HTSLIB
  touch $SETUP_DIR/htslibGet.success
fi

echo -n "Building Bio::DB::HTS ..."
if [ -e $SETUP_DIR/biohts.success ]; then
  echo " previously installed ...";
else
  echo
  cd $SETUP_DIR
  rm -rf bioDbHts
  get_distro "bioDbHts" $SOURCE_BIOBDHTS
  mkdir -p bioDbHts/htslib
  tar --strip-components 1 -C bioDbHts -zxf bioDbHts.tar.gz
  tar --strip-components 1 -C bioDbHts/htslib -jxf $SETUP_DIR/htslib.tar.bz2
  cd bioDbHts/htslib
  perl -pi -e 'if($_ =~ m/^CFLAGS/ && $_ !~ m/\-fPIC/i){chomp; s/#.+//; $_ .= " -fPIC -Wno-unused -Wno-unused-result\n"};' Makefile
  make -j$CPU
  rm -f libhts.so*
  cd ../
  env HTSLIB_DIR=$SETUP_DIR/bioDbHts/htslib perl Build.PL --install_base=$INST_PATH
  ./Build test
  ./Build install
  cd $SETUP_DIR
  rm -f bioDbHts.tar.gz
  touch $SETUP_DIR/biohts.success
fi

echo -n "Building htslib ..."
if [ -e $SETUP_DIR/htslib.success ]; then
  echo " previously installed ...";
else
  echo
  mkdir -p htslib
  tar --strip-components 1 -C htslib -jxf htslib.tar.bz2
  cd htslib
  ./configure --enable-plugins --enable-libcurl --prefix=$INST_PATH
  make -j$CPU
  make install
  cd $SETUP_DIR
  touch $SETUP_DIR/htslib.success
fi

export HTSLIB=$INST_PATH

cd $INIT_DIR

if [[ ",$COMPILE," == *,samtools,* ]] ; then
  echo -n "Building samtools ..."
  if [ -e $SETUP_DIR/samtools.success ]; then
    echo " previously installed ...";
  else
  echo
    cd $SETUP_DIR
    rm -rf samtools
    get_distro "samtools" $SOURCE_SAMTOOLS
    mkdir -p samtools
    tar --strip-components 1 -C samtools -xjf samtools.tar.bz2
    cd samtools
    ./configure --enable-plugins --enable-libcurl --prefix=$INST_PATH
    make -j$CPU all all-htslib
    make install all all-htslib
    cd $SETUP_DIR
    rm -f samtools.tar.bz2
    touch $SETUP_DIR/samtools.success
  fi
else
  echo "samtools - No change between PCAP versions"
fi

cd $SETUP_DIR
if [[ ",$COMPILE," == *,bwa,* ]] ; then
  echo -n "Building BWA ..."
  if [ -e $SETUP_DIR/bwa.success ]; then
    echo " previously installed ..."
  else
    echo
    get_distro "bwa" $SOURCE_BWA
    mkdir -p bwa
    tar --strip-components 1 -C bwa -zxf bwa.tar.gz
    make -C bwa -j$CPU
    cp bwa/bwa $INST_PATH/bin/.
    rm -f bwa.tar.gz
    touch $SETUP_DIR/bwa.success
  fi
else
  echo "BWA - No change between PCAP versions"
fi

if [[ ",$COMPILE," == *,biobambam,* ]] ; then
  echo -n "Building biobambam2 ..."
  if [ -e $SETUP_DIR/biobambam2.success ]; then
    echo " previously installed2 ..."
  else
  echo
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
    rm -f biobambam2.tar.gz
    touch $SETUP_DIR/biobambam2.success
  fi
else
  echo "biobambam - No change between PCAP versions"
fi

cd $INIT_DIR

echo -n "Building kentsrc + Bio::DB::BigFile ..."
if [ -e $SETUP_DIR/kentsrc.success ]; then
  echo " previously installed ...";
else
  echo
  cd $SETUP_DIR
  rm -rf kent kentsrc.zip
  get_distro "kentsrc" $SOURCE_KENTSRC
  unzip -q kentsrc.zip
  perl -pi -e 's/(\s+CFLAGS=)$/${1}-fPIC/' kent/src/inc/common.mk
  cd kent/src/lib
  export MACHTYPE=i686    # for a 64-bit system
  make -j$CPU
  cd ../
  export KENT_SRC=`pwd`
  cd $SETUP_DIR
  mkdir -p bigfile
  get_distro "bigfile" $SOURCE_BIGFILE
  tar --strip-components 1 -C bigfile -zxf bigfile.tar.gz
  cd bigfile
  chmod u+w Build.PL
  patch -p1 Build.PL < $INIT_DIR/dists/patch/Bio-BigFile_build.patch
  perl Build.PL --install_base=$INST_PATH
  ./Build
  ./Build test
  ./Build install
  rm -f kentsrc.zip
  touch $SETUP_DIR/kentsrc.success
fi

cd $INIT_DIR

echo -n "Building PCAP-c ..."
if [ -e $SETUP_DIR/bam_stats.success ]; then
  echo " previously installed ...";
else
  echo
  cd $INIT_DIR
  make -C c clean
  env HTSLIB=$SETUP_DIR/htslib make -C c -j$CPU prefix=$INST_PATH
  cp bin/bam_stats $INST_PATH/bin/.
  cp bin/reheadSQ $INST_PATH/bin/.
  touch $SETUP_DIR/bam_stats.success
  make -C c clean
fi

cd $INIT_DIR

echo -n "Building PCAP_perlPrereq ..."
if [ -e $SETUP_DIR/PCAP_perlPrereq.success ]; then
  echo "PCAP_perlPrereq previously installed ...";
else
  echo
  $CPANM --no-interactive --notest --mirror http://cpan.metacpan.org --notest -l $INST_PATH --installdeps .
  touch $SETUP_DIR/PCAP_perlPrereq.success
fi

echo -n "Installing PCAP ..."
$CPANM -v --no-interactive --mirror http://cpan.metacpan.org -l $INST_PATH .
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
