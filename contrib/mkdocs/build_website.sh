#!/bin/bash
MYPWD=`pwd`
if [[ -z "${NWCHEM_TOP}" ]]; then
    DIRMKDOCS=`dirname "$0"`
    NWCHEM_TOP=`echo ${MYPWD}/${DIRMKDOCS} | sed -e 's/\/contrib\/mkdocs.*//' `
fi
echo NWCHEM_TOP is ${NWCHEM_TOP}
if [[ "${MYPWD}" != ${NWCHEM_TOP}/contrib/mkdocs ]]; then
    echo wrong dir
    echo please cd to ${NWCHEM_TOP}/contrib/mkdocs
    exit
fi
rsync -av archivedforum/Special_AWCforum docs/.
cd docs
while read fname; do
    ls "$fname"
    rm -f tmptmp.txt
    ../remove_svg.sh $fname > tmptmp.txt
    cp $fname "$fname".tmp
    mv tmptmp.txt $fname
done <../mathfiles.txt
cd ..
if [[ -z "${MKDOCS_SERVE}" ]]; then
rm -rf nwchemgit.github.io_temp
git clone --depth 1 https://github.com/nwchemgit/nwchemgit.github.io  nwchemgit.github.io_temp
cd nwchemgit.github.io_temp
mkdocs -v gh-deploy --config-file ../mkdocs.yml --remote-branch master
cd ..
rm -rf  nwchemgit.github.io_temp
else
    mkdocs serve
fi    
# restore svg bits
cd docs
while read fname; do
    ls "$fname"
    mv "$fname".tmp $fname
done <../mathfiles.txt
