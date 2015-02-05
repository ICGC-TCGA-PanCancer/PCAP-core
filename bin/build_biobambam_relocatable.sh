#! /bin/bash
LIBMAUSVERSION=0.0.191-release-20150129231009
BIOBAMBAMVERSION=0.0.185-release-20150116094537
SNAPPYVERSION=1.1.2
IOLIBVERSION=1.13.9
CHRPATHVERSION=0.16
BUILDDIR=${PWD}
INSTALLDIR=${BUILDDIR}/install-dir
TOOLSDIR=${BUILDDIR}/tools-dir
PAR=`cat /proc/cpuinfo | egrep "^processor" | wc -l`

# get chrpath
if [ ! -f chrpath-${CHRPATHVERSION}.tar.gz ] ; then
	wget -O - "https://alioth.debian.org/frs/download.php/file/3979/chrpath-${CHRPATHVERSION}.tar.gz" > chrpath-${CHRPATHVERSION}.tar.gz
fi

# get iolib
if [ ! -f io_lib-${IOLIBVERSION}.tar.gz ] ; then
	wget -O - "http://downloads.sourceforge.net/project/staden/io_lib/${IOLIBVERSION}/io_lib-${IOLIBVERSION}.tar.gz?&use_mirror=kent" \
		> io_lib-${IOLIBVERSION}.tar.gz
fi

# get snappy
if [ ! -f snappy-${SNAPPYVERSION}.tar.gz ] ; then
  cp $BUILDDIR/../dists/snappy-${SNAPPYVERSION}.tar.gz snappy-${SNAPPYVERSION}.tar.gz
#	wget -O - https://snappy.googlecode.com/files/snappy-${SNAPPYVERSION}.tar.gz > snappy-${SNAPPYVERSION}.tar.gz
fi

# get libmaus
if [ ! -f libmaus-${LIBMAUSVERSION}.tar.gz ] ; then
	wget -O - https://github.com/gt1/libmaus/archive/${LIBMAUSVERSION}.tar.gz > libmaus-${LIBMAUSVERSION}.tar.gz
fi

# get biobambam
if [ ! -f biobambam-${BIOBAMBAMVERSION}.tar.gz ] ; then
	wget -O - https://github.com/gt1/biobambam/archive/${BIOBAMBAMVERSION}.tar.gz > biobambam-${BIOBAMBAMVERSION}.tar.gz
fi

tar xzvf chrpath-${CHRPATHVERSION}.tar.gz
mv chrpath-${CHRPATHVERSION} chrpath-${CHRPATHVERSION}-src
mkdir -p chrpath-${CHRPATHVERSION}-build
cd chrpath-${CHRPATHVERSION}-build
${BUILDDIR}/chrpath-${CHRPATHVERSION}-src/configure --prefix=${TOOLSDIR}
make -j${PAR}
make -j${PAR} install
cd ..
rm -fR chrpath-${CHRPATHVERSION}-src chrpath-${CHRPATHVERSION}-build

rm -fR ${INSTALLDIR}
mkdir -p ${INSTALLDIR}

# build iolib
tar xzvf io_lib-${IOLIBVERSION}.tar.gz
mv io_lib-${IOLIBVERSION} io_lib-${IOLIBVERSION}-src
mkdir -p io_lib-${IOLIBVERSION}-build
cd io_lib-${IOLIBVERSION}-build
LDFLAGS="-Wl,-rpath=XORIGIN/../lib" ${BUILDDIR}/io_lib-${IOLIBVERSION}-src/configure --prefix=${INSTALLDIR}
make -j${PAR}
make -j${PAR} install
cd ..
rm -fR io_lib-${IOLIBVERSION}-src io_lib-${IOLIBVERSION}-build

# build snappy
tar xzvf snappy-${SNAPPYVERSION}.tar.gz
mv snappy-${SNAPPYVERSION} snappy-${SNAPPYVERSION}-src
mkdir -p snappy-${SNAPPYVERSION}-build
cd snappy-${SNAPPYVERSION}-build
LDFLAGS="-Wl,-rpath=XORIGIN/../lib" ${BUILDDIR}/snappy-${SNAPPYVERSION}-src/configure --prefix=${INSTALLDIR}
make -j${PAR}
make -j${PAR} install
cd ..
rm -fR snappy-${SNAPPYVERSION}-src snappy-${SNAPPYVERSION}-build

# build libmaus
tar xzvf libmaus-${LIBMAUSVERSION}.tar.gz
mv libmaus-${LIBMAUSVERSION} libmaus-${LIBMAUSVERSION}-src
mkdir -p libmaus-${LIBMAUSVERSION}-build
cd libmaus-${LIBMAUSVERSION}-build
LDFLAGS="-Wl,-rpath=XORIGIN/../lib" ${BUILDDIR}/libmaus-${LIBMAUSVERSION}-src/configure --prefix=${INSTALLDIR} \
	--with-snappy=${INSTALLDIR} \
	--with-io_lib=${INSTALLDIR}
make -j${PAR}
make -j${PAR} install
cd ..
rm -fR libmaus-${LIBMAUSVERSION}-src libmaus-${LIBMAUSVERSION}-build

# build biobambam
tar xzvf biobambam-${BIOBAMBAMVERSION}.tar.gz
mv biobambam-${BIOBAMBAMVERSION} biobambam-${BIOBAMBAMVERSION}-src
mkdir -p biobambam-${BIOBAMBAMVERSION}-build
cd biobambam-${BIOBAMBAMVERSION}-build
LDFLAGS="-Wl,-rpath=XORIGIN/../lib" ${BUILDDIR}/biobambam-${BIOBAMBAMVERSION}-src/configure --prefix=${INSTALLDIR} \
	--with-libmaus=${INSTALLDIR}
make -j${PAR}
make -j${PAR} install
cd ..
rm -fR biobambam-${BIOBAMBAMVERSION}-src biobambam-${BIOBAMBAMVERSION}-build

for i in `find ${INSTALLDIR} -name \*.so\*` ; do
	ORIG=`objdump -x ${i} | grep RPATH | awk '{print $2}'`
	MOD=`echo "$ORIG" | sed "s/XORIGIN/\\$ORIGIN/"`
	${TOOLSDIR}/bin/chrpath -r "${MOD}" ${i}
done

for i in ${INSTALLDIR}/bin/* ; do
	if [ ! -z `LANG=C file ${i} | egrep "ELF.*executable" | awk '{print $1}' | perl -p -e "s/://"` ] ; then
		ORIG=`objdump -x ${i} | grep RPATH | awk '{print $2}'`
		MOD=`echo "$ORIG" | sed "s/XORIGIN/\\$ORIGIN/"`
		${TOOLSDIR}/bin/chrpath -r "${MOD}" ${i}
	fi
done

#mv ${INSTALLDIR} biobambam-${BIOBAMBAMVERSION}
#tar czvf biobambam-${BIOBAMBAMVERSION}-`uname -p`-`uname -s`.tar.gz biobambam-${BIOBAMBAMVERSION}
#rm -fR biobambam-${BIOBAMBAMVERSION}

rm -fR ${TOOLSDIR}

# my additions
mv ${INSTALLDIR} biobambam
rm -f chrpath-${CHRPATHVERSION}.tar.gz
rm -f io_lib-${IOLIBVERSION}.tar.gz
rm -f snappy-${SNAPPYVERSION}.tar.gz
rm -f libmaus-${LIBMAUSVERSION}.tar.gz
rm -f biobambam-${BIOBAMBAMVERSION}.tar.gz
