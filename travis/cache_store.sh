#!/bin/bash
          mkdir -p ~/cache/libext/lib  ~/cache/libext/bin  ~/cache/libext/include ~/cache/libext/mpich || true
          mkdir -p ~/cache/libext/mpich/mpich || true
          mkdir -p ~/cache/libext/lib  ~/cache/libext/bin  ~/cache/libext/include || true
          mkdir -p ~/cache/simint/simint_install || true
          mkdir -p ~/cache/tarballs || true
          mkdir -p ~/cache/libxc || true
          ls -la src/libext ||true
          [ -d "src/libext/lib" ] && rsync -av src/libext/lib/*  ~/cache/libext/lib/. || true
          [ -d "src/libext/bin" ] && rsync -av src/libext/bin/*  ~/cache/libext/bin/. || true
          [ -d "src/libext/include" ] && rsync -av src/libext/include/*  ~/cache/libext/include/.  || true
          [ -d "src/NWints/simint/libsimint_source/simint_install/lib" ] && rsync -av src/NWints/simint/libsimint_source/simint_install/* ~/cache/simint/simint_install/.  || true
          [ -f "src/nwpw/nwpwlib/nwpwxc/dftd3.tgz" ] && rsync -av  src/nwpw/nwpwlib/nwpwxc/dftd3.tgz ~/cache/tarballs/dftd3.tgz || true
          [ -f "src/libext/libxc/install/lib/libxc.a" ] && rsync -av  src/libext/libxc/install ~/cache/libxc/. || true
          echo "cache stored"
          ls -Rla ~/cache
	  du -sh ~/cache
	  du -sh ~/cache/* || true
	  du -sh ~/cache/*/* || true
