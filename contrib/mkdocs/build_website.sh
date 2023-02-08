#!/usr/bin/env bash
#MYPWD=`pwd`
#if [[ -z "${NWCHEM_TOP}" ]]; then
#    DIRMKDOCS=`dirname "$0"`
#    NWCHEM_TOP=`echo ${MYPWD}/${DIRMKDOCS} | sed -e 's/\/contrib\/mkdocs.*//' `
#fi
#echo NWCHEM_TOP is ${NWCHEM_TOP}
#if [[ "${MYPWD}" != ${NWCHEM_TOP}/contrib/mkdocs ]]; then
#    echo wrong dir
#    echo please cd to ${NWCHEM_TOP}/contrib/mkdocs
#    exit
#fi
#fresh clone of wiki
if [ -d "docs" ]; then
    echo ' WARNING'
    echo ' since the docs directory is already present, '
    echo ' it will not be updated'
    echo ' '
    cd docs
else    
    rm -rf docs #archivedforum
#git clone --depth 1 git@github.com:nwchemgit/nwchem-wiki.git docs
wget -q https://github.com/nwchemgit/nwchem-wiki/tarball/master -O - | tar -xz
mv nwchemgit-nwchem-wiki-* docs
cd docs
wget -q  https://github.com/nwchemgit/archivedforum/tarball/master -O - | tar -xz --wildcards   nwchemgit-archivedforum-*/Special_AWCforum/*
mv nwchemgit-archivedforum-*/Special_AWCforum .
#git clone --depth 1 git@github.com:nwchemgit/archivedforum.git
#mv archivedforum/Special_AWCforum .
fi
cd ..
if [[ -z "${MKDOCS_SERVE}" ]]; then
    #git clone --depth 1 https://github.com/nwchemgit/nwchemgit.github.io  nwchemgit.github.io_temp
rm -rf nwchemgit.github.io_temp
git clone --depth 1 git@github.com:nwchemgit/nwchemgit.github.io.git nwchemgit.github.io_temp
cd nwchemgit.github.io_temp
mkdocs  gh-deploy --config-file ../mkdocs.yml --remote-branch master
#echo "********"
#echo remember to apply preload.patch to nwchemgit.github.io
#echo using the preload.sh script
#echo "********"
cd ..
rm -rf  nwchemgit.github.io_temp
elif [[ "${MKDOCS_SERVE}" == "B" ]]; then
    echo 'building'
    mkdocs build
else
    mkdocs serve
fi    
# restore svg bits
rm -rf docs
#cd docs
#while read fname; do
#    ls "$fname"
#    mv "$fname".tmp $fname
#done <../mathfiles.txt
