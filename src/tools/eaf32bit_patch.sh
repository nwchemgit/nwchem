#!/usr/bin/env bash
rm -f eaf32bit.patch
wget -O - https://github.com/GlobalArrays/ga/commit/2554b4431087c43f411f3e59426b83f6e9c48bd1.patch  | sed -e 's/ a\/pario/ '$1'\/pario/' | sed -e 's/ b\/pario/ '$1'\/pario/' > eaf32bit.patch
patch -p0 -s -N < eaf32bit.patch
echo eaf32bit.patch applied
