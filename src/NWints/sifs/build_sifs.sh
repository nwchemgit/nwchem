#!/usr/bin/bash


rm -rf columbu* *.F
curl -LJO https://gitlab.com/api/v4/projects/36816383/repository/archive?path=Columbus/source/gcfci/colib/sifs
curl -LJO https://gitlab.com/api/v4/projects/36816383/repository/archive?path=Columbus/source/gcfci/colib/pack
curl -LJO https://gitlab.com/api/v4/projects/36816383/repository/archive?path=Columbus/source/gcfci/colib/basicio
for f in columbus-master*gz; do tar xzf "$f"; done
mv $(find columbus-master*/ -name "*F") .
wget -O bummer.F https://gitlab.com/columbus-program-system/columbus/-/raw/master/Columbus/source/gcfci/colib/humanio/bummer.F?inline=false
wget -O izero.F https://gitlab.com/columbus-program-system/columbus/-/raw/master/Columbus/source/gcfci/colib/linear_algebra/izero.F?inline=false
