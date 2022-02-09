#!/usr/bin/env bash
if [[ $1 == "" ]]; then
    user=nwchemgit
else
    user=$1
fi
if [[ $2 == "" ]]; then
    repos=nwchem
else
    repos=$2
fi
echo downloads for $user/$repos
echo https://api.github.com/repos/$user/$repos/releases
curl -s https://api.github.com/repos/$user/$repos/releases|egrep -i down|sed 's/\"browser\_download\_url\":\ \"https:\/\/github.com\/$user\/$repos\/releases\/download\///g' | sed 's/\"//g'
if [[ `uname -s` != 'Darwin' ]]; then
    echo Total number of downloads
    curl -s https://api.github.com/repos/$user/$repos/releases | egrep 'download_count'  | cut '-d:' -f 2 | sed 's/,/+/' | xargs echo | xargs  -n 1 -I{} echo {} 0 | bc
fi
# curl -s https://api.github.com/repos/nwchemgit/nwchem/releases|egrep dow|sed 's/\"browser\_download\_url\":\ \"https:\/\/github.com\/nwchemgit\/nwchem\/releases\/download\//g'
