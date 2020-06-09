#!/bin/bash
rm -f win64_ma_nxtval.patch
cat > win64_ma_nxtval.patch <<EOF
--- $1/ma/ma.c.org	2020-06-08 16:55:52.821749074 -0700
+++ $1/ma/ma.c	2020-06-08 16:56:01.765785136 -0700
@@ -114,7 +114,11 @@
  **/
 
 typedef unsigned int Guard;   /* for detection of memory trashing */
+#ifdef WIN64
+typedef unsigned long long ulongi; /* for brevity */
+#else
 typedef unsigned long ulongi; /* for brevity */
+#endif
 
 /* allocation request for a block */
 typedef struct _AR
--- $1/tcgmsg/tcgmsg-mpi/nxtval-armci.c.org	2020-06-08 16:56:58.610007026 -0700
+++ $1/tcgmsg/tcgmsg-mpi/nxtval-armci.c	2020-06-08 16:57:05.114032467 -0700
@@ -66,12 +66,17 @@
         if (*mproc > 0) {
 #if   SIZEOF_F77_INTEGER == SIZEOF_INT
             int op = ARMCI_FETCH_AND_ADD;
+            rc = ARMCI_Rmw(op,(void*)&local,(void*)pnxtval_counter,1,server);
 #elif SIZEOF_F77_INTEGER == SIZEOF_LONG
+            rc = ARMCI_Rmw(op,(void*)&local,(void*)pnxtval_counter,1,server);
             int op = ARMCI_FETCH_AND_ADD_LONG;
 #else
+#ifdef WIN64
+#else
+                Error("nxtval: not implemented",0);
 #   error
 #endif
-            rc = ARMCI_Rmw(op,(void*)&local,(void*)pnxtval_counter,1,server);
+#endif
         }
     } else {
         /* Not running in parallel ... just do a simulation */
EOF
patch -p0 < win64_ma_nxtval.patch

