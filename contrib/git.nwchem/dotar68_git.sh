#!/bin/bash 
rm -rf temp.`date +%Y%m%d`
mkdir -p temp.`date +%Y%m%d`
cd temp.`date +%Y%m%d`
NWVER='6.8-release'
git clone  -b release-6-8 https://github.com/nwchemgit/nwchem nwchem-6.8
cd nwchem-6.8/src/tools
######REV68="$(svn info -r HEAD |grep Revision: |cut -c11-)"
rm -f *tar*gz
./get-tools-github
cd ../util
touch util_module_avail.F
sh ./util_ga_version.bash
sh ./util_nwchem_version.bash
#sed -i 's/N\/A/v5.6.2/'  util_ga_version.F
cat util_ga_version.F
cd ../..
REV68="$(git describe)"
cd ..
echo 'revision ' $REV68 
pwd
rm -f *md5 *tar*
tar --exclude=".git" -czvf nwchem-"${NWVER}".revision-"${REV68}"-src.`date +%Y-%m-%d`.tar.gz nwchem-6.8/*
md5sum nwchem-"${NWVER}".revision-"${REV68}"-src.`date +%Y-%m-%d`.tar.gz > nwchem-"${NWVER}".revision"${REV68}"-src.`date +%Y-%m-%d`.tar.gz.md5
tar --exclude=".git" -cjvf nwchem-"${NWVER}".revision-"${REV68}"-src.`date +%Y-%m-%d`.tar.bz2 nwchem-6.8/*
md5sum nwchem-"${NWVER}".revision-"${REV68}"-src.`date +%Y-%m-%d`.tar.bz2 >  nwchem-"${NWVER}".revision"${REV68}"-src.`date +%Y-%m-%d`.tar.bz2.md5
tar --exclude=".git" -cjvf nwchem-"${NWVER}".revision-"${REV68}"-srconly.`date +%Y-%m-%d`.tar.bz2 nwchem-6.8/src/*
md5sum nwchem-"${NWVER}".revision-"${REV68}"-srconly.`date +%Y-%m-%d`.tar.bz2>  nwchem-"${NWVER}".revision"${REV68}"-srconly.`date +%Y-%m-%d`.tar.bz2.md5
tar --exclude=".git" --exclude="src" -cjvf nwchem-"${NWVER}".revision-"${REV68}"-nonsrconly.`date +%Y-%m-%d`.tar.bz2 nwchem-6.8/*
md5sum nwchem-"${NWVER}".revision-"${REV68}"-nonsrconly.`date +%Y-%m-%d`.tar.bz2 > nwchem-"${NWVER}".revision"${REV68}"-nonsrconly.`date +%Y-%m-%d`.tar.bz2.md5
echo 'upload to http://192.101.105.206/'
