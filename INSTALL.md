#To install run:

`./setup.sh /install/to/here`

'/install/to/here' is where you want the bin/lib folders to be created.

#Perl:
  Minimum version: 5.10.1
  Tested version: 5.16.3

#OS:
  This distribution will only work on *NIX type systems at present.

#Other Software
  For installation to proceed you require the following packages:


##For Ubuntu (tested with 14.04)
```
apt-get
build-essential zlib1g-dev libncurses5-dev libcurl4-gnutls-dev libssl-dev libexpat1-dev libgd-dev nettle-dev
```

##For Amazon Linux AMI 2016.03.0 x86_64

```
yum install
make glibc-devel gcc patch ncurses-devel expat-devel gd-devel perl-core openssl-devel libcurl-devel gnutls-devel libtasn1-devel p11-kit-devel gmp-devel nettle-devel
```

**Should nettle-devel not exist**
Use the following

```
yum install autoconf
wget https://git.lysator.liu.se/nettle/nettle/repository/archive.tar.gz?ref=nettle_3.2_release_20160128 -O nettle.tar.gz
mkdir -p nettle
tar --strip-components 1 -C nettle -zxf nettle.tar.gz
cd nettle
./.bootstrap
./configure && \
sudo make && \
sudo make check && \
sudo make install
```

setup.sh will install
-biobambam
-bwa
-samtools

**NOTE:** 
bwa_aln.pl will only function when 0.6.x installed
(you will need to make this available on path manually)
bwa_mem.pl will only function when 0.7.x installed

