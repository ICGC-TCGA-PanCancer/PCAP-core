#!/bin/bash

SOURCE_BWA="https://github.com/lh3/bwa/archive/0.7.10.tar.gz"
SOURCE_SAMTOOLS="https://github.com/samtools/samtools/archive/0.1.19.tar.gz"

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
    curl -sS -o $1.tar.gz -L $2
  else
    wget -nv -O $1.tar.gz $2
  fi
  mkdir -p $1
  tar --strip-components 1 -C $1 -zxf $1.tar.gz
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

perlmods=( "File::ShareDir" "File::ShareDir::Install" "Const::Fast" )

set -e
for i in "${perlmods[@]}" ; do
  echo -n "Installing build prerequisite $i..."
  (
    set -x
    $INIT_DIR/bin/cpanm -v --mirror http://cpan.metacpan.org -l $INST_PATH $i
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
  unset PERL5LIB
  echo -n "Building biobambam ..."
  if [ -e $SETUP_DIR/biobambam.success ]; then
    echo -n " previously installed ..."
  else
    (
      cd $SETUP_DIR
      $INIT_DIR/bin/build_biobambam_relocatable.sh
      cp -r biobambam/* $INST_PATH/.
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
  echo -n "Building samtools ..."
  if [ -e $SETUP_DIR/samtools.success ]; then
    echo -n " previously installed ...";
  else
    cd $SETUP_DIR
    if( [ "x$SAMTOOLS" == "x" ] ); then
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
  fi
  done_message "" "Failed to build samtools."
else
  echo "samtools - No change between PCAP versions"
fi

export SAMTOOLS="$SETUP_DIR/samtools"

#add bin path for PCAP install tests
export PATH="$INST_PATH/bin:$PATH"

cd $INIT_DIR

echo -n "Installing Perl prerequisites ..."
if ! ( perl -MExtUtils::MakeMaker -e 1 >/dev/null 2>&1); then
    echo
    echo "WARNING: Your Perl installation does not seem to include a complete set of core modules.  Attempting to cope with this, but if installation fails please make sure that at least ExtUtils::MakeMaker is installed.  For most users, the best way to do this is to use your system's package manager: apt, yum, fink, homebrew, or similar."
fi
(
  set -x
  $INIT_DIR/bin/cpanm -v --mirror http://cpan.metacpan.org -notest -l $INST_PATH/ --installdeps . < /dev/null
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

echo
echo
echo "Please add the following to beginning of path:"
echo "  $INST_PATH/bin"
echo "Please add the following to beginning of PERL5LIB:"
echo "  $PERLROOT"
echo "  $PERLARCH"
echo

exit 0
