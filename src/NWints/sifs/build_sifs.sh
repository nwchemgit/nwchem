#!/usr/bin/env bash
rm -rf columbu* *.F
curl -LJO https://gitlab.com/api/v4/projects/36816383/repository/archive?path=Columbus/source/gcfci/colib/sifs
ls -l columbus-master*gz
tar xzf columbus-master*gz
mv $(find columbus-master*s -name "*F") .
wget -O bummer.F https://gitlab.com/columbus-program-system/columbus/-/raw/master/Columbus/source/gcfci/colib/humanio/bummer.F?inline=false
