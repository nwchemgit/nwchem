#!/usr/bin/env bash
rm -f peigstubs.patch
cat > peigstubs.patch <<EOF
--- $1/global/src/peigstubs.c
+++ $1/global/src/peigstubs.c
@@ -11,9 +11,9 @@
 #       define gai_diag_       F77_FUNC_(gai_diag,GAI_DIAG)
 #       define gai_diag_std_   F77_FUNC_(gai_diag_std,GAI_DIAG_STD)
 #       define gai_diag_reuse_ F77_FUNC_(gai_diag_reuse,GAI_DIAG_REUSE)
-extern gai_diag_(Integer*,Integer*,Integer*,DoublePrecision*);
-extern gai_diag_std_(Integer*,Integer*,DoublePrecision*);
-extern gai_diag_reuse_(Integer*,Integer*,Integer*,Integer*,DoublePrecision*);
+extern void gai_diag_(Integer*,Integer*,Integer*,DoublePrecision*);
+extern void gai_diag_std_(Integer*,Integer*,DoublePrecision*);
+extern void gai_diag_reuse_(Integer*,Integer*,Integer*,Integer*,DoublePrecision*);
 #   else
 #   endif
 #else
EOF
patch -p0 -s -N < peigstubs.patch
echo strdup.patch applied

