#!/bin/bash 
rm -rf temp.`date +%Y%m%d`
mkdir -p temp.`date +%Y%m%d`
cd temp.`date +%Y%m%d`
NWVER='7.0.0-release'
TOPDIR="nwchem-7.0.0"
git clone --depth 1 --shallow-submodules -b release-7-0-0 https://github.com/nwchemgit/nwchem $TOPDIR
cd $TOPDIR/src/tools
rm -f *tar*gz
./get-tools-github
cd ../util
./util_ga_version.bash
./util_nwchem_version.bash
cd ..
make nwchem_config NWCHEM_MODULES=all\ python
export EACCSD=1
export IPCCSD=1
export CCSDTQ=1
export MRCC_METHODS=1
make 64_to_32   USE_INTERNALBLAS=y
#rm `find . -name dependencies`
#rm `find . -name include_stamp`
#rm `find peigs -name peigs_stamp.*`
cd ..
rm -rf bin lib
REV700="$(git describe)"
cd ..
echo 'revision ' $REV700
pwd
rm -f *md5 *tar*
tar --exclude=".git" -czvf nwchem-"${NWVER}".revision-"${REV700}"-src.`date +%Y-%m-%d`.tar.gz $TOPDIR/*
md5sum nwchem-"${NWVER}".revision-"${REV700}"-src.`date +%Y-%m-%d`.tar.gz > nwchem-"${NWVER}".revision"${REV700}"-src.`date +%Y-%m-%d`.tar.gz.md5
tar --exclude=".git" -cjvf nwchem-"${NWVER}".revision-"${REV700}"-src.`date +%Y-%m-%d`.tar.bz2 $TOPDIR/*
md5sum nwchem-"${NWVER}".revision-"${REV700}"-src.`date +%Y-%m-%d`.tar.bz2 >  nwchem-"${NWVER}".revision"${REV700}"-src.`date +%Y-%m-%d`.tar.bz2.md5
tar --exclude=".git" -cjvf nwchem-"${NWVER}".revision-"${REV700}"-srconly.`date +%Y-%m-%d`.tar.bz2 $TOPDIR/src/*
md5sum nwchem-"${NWVER}".revision-"${REV700}"-srconly.`date +%Y-%m-%d`.tar.bz2>  nwchem-"${NWVER}".revision"${REV700}"-srconly.`date +%Y-%m-%d`.tar.bz2.md5
tar --exclude=".git" --exclude="src" -cjvf nwchem-"${NWVER}".revision-"${REV700}"-nonsrconly.`date +%Y-%m-%d`.tar.bz2 $TOPDIR/*
md5sum nwchem-"${NWVER}".revision-"${REV700}"-nonsrconly.`date +%Y-%m-%d`.tar.bz2 > nwchem-"${NWVER}".revision"${REV700}"-nonsrconly.`date +%Y-%m-%d`.tar.bz2.md5
echo 'upload to http://192.101.105.206/'
