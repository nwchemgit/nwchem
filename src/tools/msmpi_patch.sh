#!/bin/bash
rm -f comex_configure.patch
cat > comex_configure.patch <<EOF
--- $1/comex/configure.orig	2020-06-04 11:16:29.506196033 -0700
+++ $1/comex/configure	2020-06-04 11:19:01.610898545 -0700
@@ -7501,11 +7501,19 @@
 #ifdef __cplusplus
 extern "C"
 #endif
+#ifdef WIN32
+char __stdcall MPI_Finalize ();
+#else
 char MPI_Init ();
+#endif
 int
 main ()
 {
+#ifdef WIN32
+return MPI_Finalize ();
+#else
 return MPI_Init ();
+#endif
   ;
   return 0;
 }
EOF
patch -p0 < comex_configure.patch

