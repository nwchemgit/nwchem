#!/usr/bin/env bash
rm -f ga_ma_align.patch
cat > ga_ma_align.patch <<EOF
--- $1/global/src/base.c.org	2023-10-19 15:49:17
+++ $1/global/src/base.c	2023-10-19 15:49:04
@@ -76,7 +76,7 @@
 /*#define CHECK_MA yes */
 
 /*uncomment line below to verify if MA base address is alligned wrt datatype*/
-#if !(defined(LINUX) || defined(CRAY) || defined(CYGWIN))
+#if !(defined(LINUX) || defined(CRAY) || defined(CYGWIN) || defined(MACX))
 #define CHECK_MA_ALGN 1
 #else
 #endif
EOF
patch -p0 -s -N < ga_ma_align.patch
echo ga_ma_align.patch applied
