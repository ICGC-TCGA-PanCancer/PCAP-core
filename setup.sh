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

done_message () {
    if [ $? -eq 0 ]; then
        echo " done."
        if [ "x$1" != "x" ]; then
            echo $1
        fi
    else
        echo " failed.  See setup.log file for error messages." $2
        echo "    Please check INSTALL file for items that should be installed by a package manager"
        exit 1
    fi
}

get_distro () {
  if hash curl 2>/dev/null; then
    curl --insecure -sS -o $1.tar.gz -L $2
  else
    wget -nv -O $1.tar.gz $2
  fi
  mkdir -p $1
  tar --strip-components 1 -C $1 -zxf $1.tar.gz
}

get_file () {
# output, source
  if hash curl 2>/dev/null; then
    curl --insecure -sS -o $1 -L $2
  else
    wget -nv -O $1 $2
  fi
}

if [ "$#" -ne "1" ] ; then
  echo "Please provide an installation path  such as /opt/ICGC"
  exit 0
fi

CPU=`cat /proc/cpuinfo | egrep "^processor" | wc -l`
echo "Max compilation CPUs set to $CPU"

INST_PATH=$1

# get current directory
INIT_DIR=`pwd`

# cleanup inst_path
mkdir -p $INST_PATH/bin
cd $INST_PATH
INST_PATH=`pwd`
cd $INIT_DIR

# make sure that build is self contained
unset PERL5LIB
ARCHNAME=`perl -e 'use Config; print $Config{archname};'`
PERLROOT=$INST_PATH/lib/perl5
PERLARCH=$PERLROOT/$ARCHNAME
export PERL5LIB="$PERLROOT:$PERLARCH"

#create a location to build dependencies
SETUP_DIR=$INIT_DIR/install_tmp
mkdir -p $SETUP_DIR

# re-initialise log file
echo > $INIT_DIR/setup.log

# log information about this system
(
    echo '============== System information ===='
    set -x
    lsb_release -a
    uname -a
    sw_vers
    system_profiler
    grep MemTotal /proc/meminfo
    set +x
    echo; echo
) >>$INIT_DIR/setup.log 2>&1

## grab cpanm:
rm -f $SETUP_DIR/cpanm
get_file $SETUP_DIR/cpanm https://cpanmin.us/
chmod +x $SETUP_DIR/cpanm

perlmods=( "File::ShareDir" "File::ShareDir::Install" "Const::Fast" )

set -e
for i in "${perlmods[@]}" ; do
  echo -n "Installing build prerequisite $i..."
  (
    set -x
    $SETUP_DIR/cpanm -v --mirror http://cpan.metacpan.org -l $INST_PATH $i
    set +x
    echo; echo
  ) >>$INIT_DIR/setup.log 2>&1
  done_message "" "Failed during installation of $i."
done

# figure out the upgrade path
COMPILE=`echo 'nothing' | perl -I lib -MPCAP -ne 'print PCAP::upgrade_path($_);'`
if [ -e "$INST_PATH/lib/perl5/PCAP.pm" ]; then
  COMPILE=`perl -I $INST_PATH/lib/perl5 -MPCAP -e 'print PCAP->VERSION,"\n";' | perl -I lib -MPCAP -ne 'print PCAP::upgrade_path($_);'`
fi

#add bin path for PCAP install tests
export PATH="$INST_PATH/bin:$PATH"

#Need to add CaVEMan stuff here... will depend on samtools too (for now).

echo -n "Building jkentUtils ..."
if [ -e $SETUP_DIR/jkentUtils.success ]; then
  echo -n " previously installed ...";
else
  cd $SETUP_DIR
  if [[ `uname -m` == x86_64 ]] ; then
    (
    get_file $INST_PATH/bin/wigToBigWig $SOURCE_JKENT_BIN/wigToBigWig
    chmod +x $INST_PATH/bin/wigToBigWig
    get_file $INST_PATH/bin/bigWigMerge $SOURCE_JKENT_BIN/bigWigMerge
    chmod +x $INST_PATH/bin/bigWigMerge
    ) >>$INIT_DIR/setup.log 2>&1
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
done_message "" "Failed to build jkentUtils."

echo -n "Building htslib ..."
if [ -e $SETUP_DIR/htslib.success ]; then
  echo -n " previously installed ...";
else
  cd $SETUP_DIR
  (
  set -xe
  if [ ! -e htslib ]; then
    get_distro "htslib" $SOURCE_HTSLIB
  fi
  make -C htslib -j$CPU
  touch $SETUP_DIR/htslib.success
  )>>$INIT_DIR/setup.log 2>&1
fi
done_message "" "Failed to build htslib."

export HTSLIB="$SETUP_DIR/htslib"

echo -n "Building bam_stats ..."
if [ -e $SETUP_DIR/bam_stats.success ]; then
  echo -n " previously installed ...";
else
  cd $INIT_DIR
  (
  set -xe
  make -C c clean
  make -C c -j$CPU
  cp bin/bam_stats $INST_PATH/bin/.
  touch $SETUP_DIR/bam_stats.success
  # need to clean up as will clash with other version
  rm -rf $SAMTOOLS
  make -C c clean
  )>>$INIT_DIR/setup.log 2>&1
fi
done_message "" "Failed to build bam_stats."

cd $SETUP_DIR
if [[ ",$COMPILE," == *,bwa,* ]] ; then
  echo -n "Building BWA ..."
  if [ -e $SETUP_DIR/bwa.success ]; then
    echo -n " previously installed ..."
  else
    (
      get_distro "bwa" $SOURCE_BWA
      cd $SETUP_DIR/bwa
      make -j$CPU
      cp bwa $INST_PATH/bin/.
      touch $SETUP_DIR/bwa.success
    ) >>$INIT_DIR/setup.log 2>&1
  fi
  done_message "" "Failed to build bwa."
else
  echo "BWA - No change between PCAP versions"
fi

if [[ ",$COMPILE," == *,biobambam,* ]] ; then
  echo -n "Building biobambam ..."
  if [ -e $SETUP_DIR/biobambam.success ]; then
    echo -n " previously installed ..."
  else
    (
      cd $SETUP_DIR
      mkdir -p biobambam
      get_distro "biobambam" $SOURCE_BBB_BIN_DIST
      mkdir -p $INST_PATH/bin $INST_PATH/include $INST_PATH/lib $INST_PATH/share
      rm -f biobambam/bin/curl* # breaks OS installs
      cp -r biobambam/bin/* $INST_PATH/bin/.
      cp -r biobambam/include/* $INST_PATH/include/.
      cp -r biobambam/lib/* $INST_PATH/lib/.
      cp -r biobambam/share/* $INST_PATH/share/.
      touch $SETUP_DIR/biobambam.success
    ) >>$INIT_DIR/setup.log 2>&1
  fi
  export PERL5LIB="$PERLROOT:$PERLARCH"
  done_message "" "Failed to build biobambam."
else
  echo "biobambam - No change between PCAP versions"
fi

cd $INIT_DIR

if [[ ",$COMPILE," == *,samtools,* ]] ; then
  echo -n "Building samtools v0.x ..."
  if [ -e $SETUP_DIR/samtools.success ]; then
    echo -n " previously installed ...";
  else
    cd $SETUP_DIR
    (
    set -x
    if [ ! -e samtools ]; then
      get_distro "samtools" $SOURCE_SAMTOOLS
      perl -i -pe 's/^CFLAGS=\s*/CFLAGS=-fPIC / unless /\b-fPIC\b/' samtools/Makefile
    fi
    make -C samtools -j$CPU
    set +x
    cp samtools/samtools $INST_PATH/bin/.
    touch $SETUP_DIR/samtools.success
    )>>$INIT_DIR/setup.log 2>&1
  fi
  done_message "" "Failed to build samtools v0.x."
else
  echo "samtools - No change between PCAP versions"
fi

export SAMTOOLS="$SETUP_DIR/samtools"

cd $INIT_DIR

echo -n "Installing Perl prerequisites ..."
if ! ( perl -MExtUtils::MakeMaker -e 1 >/dev/null 2>&1); then
    echo
    echo "WARNING: Your Perl installation does not seem to include a complete set of core modules.  Attempting to cope with this, but if installation fails please make sure that at least ExtUtils::MakeMaker is installed.  For most users, the best way to do this is to use your system's package manager: apt, yum, fink, homebrew, or similar."
fi
(
  set -x
  $SETUP_DIR/cpanm -v --mirror http://cpan.metacpan.org --notest -l $INST_PATH/ --installdeps . < /dev/null
  set +x
) >>$INIT_DIR/setup.log 2>&1
done_message "" "Failed during installation of core dependencies."

echo -n "Installing PCAP ..."
(
  cd $INIT_DIR
  perl Makefile.PL INSTALL_BASE=$INST_PATH
  make
  make test
  make install
) >>$INIT_DIR/setup.log 2>&1
done_message "" "PCAP install failed."

# cleanup all junk
rm -rf $SETUP_DIR
rm -rf $INIT_DIR/bin/biobambam



echo
echo
echo "Please add the following to beginning of path:"
echo "  $INST_PATH/bin"
echo "Please add the following to beginning of PERL5LIB:"
echo "  $PERLROOT"
echo "  $PERLARCH"
echo

exit 0
