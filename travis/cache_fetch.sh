#!/bin/bash
set -v
          ls -Rla ~/cache
          [ -d ~/cache/libext ] && rsync -av ~/cache/libext/* src/libext/.  
          [ -d ~/cache/simint/simint_install/lib ] && \
          mkdir -p src/NWints/simint/libsimint_source/simint_install && \
          rsync -av ~/cache/simint/simint_install/* src/NWints/simint/libsimint_source/simint_install/.   || true
          [ -f ~/cache/tarballs/dftd3.tgz ] && rsync -av ~/cache/tarballs/dftd3.tgz  src/nwpw/nwpwlib/nwpwxc/dftd3.tgz  || true
          [ -f ~/cache/libxc/install/lib/libxc.a ] && rsync -av ~/cache/libxc/install src/libext/libxc/. || true
          echo "cache fetched"
	  echo " mpich debug "
	  ls -lart /home/runner/work/nwchem/nwchem/src/libext/mpich/mpich/../../include || true
	  ls -lart /home/runner/work/nwchem/nwchem/src/libext/include || true
	  ls -lart /home/runner/work/nwchem/nwchem/src/libext/mpich || true
