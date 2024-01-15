#!/usr/bin/env bash
rm -f o.tgz
wget -O o.tgz https://gitlab.com/api/v4/projects/36816383/repository/archive?path=Columbus/source/colib 2> /dev/null
gzip -d -t o.tgz
if [ "$?" -eq 0 ]; then
    echo OK
fi
