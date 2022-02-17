#!/usr/bin/env bash
rm -f nompif.patch
cat > nompif.patch <<EOF
--- $1/global/src/scalapack.F
+++ $1/global/src/scalapack.F
@@ -3605,7 +3605,7 @@ c     broadcast evals
 #endif
 #endif
       implicit none
-#include "mpif.h"
+c#include "mpif.h"
 #include "mafdecls.fh"
 #include "global.fh"
 #include "scalapack.fh"
EOF
patch -p0 -s -N < nompif.patch
echo nompif.patch applied
