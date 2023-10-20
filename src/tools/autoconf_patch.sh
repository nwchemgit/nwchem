#!/usr/bin/env bash
rm -f autoconf.patch
cat > autoconf.patch <<EOF
--- $1/configure.org	2022-11-03 14:47:40.000000000 -0700
+++ $1/configure	2023-10-20 15:36:29.516676155 -0700
@@ -21699 +21698,0 @@
-  s/\/.*\/lib.*\.a//p ;
@@ -21811 +21809,0 @@
-  s/\/.*\/lib.*\.a//p ;
@@ -22145 +22142,0 @@
-  s/\/.*\/lib.*\.a//p ;
@@ -22257 +22253,0 @@
-  s/\/.*\/lib.*\.a//p ;
EOF
patch -p0 -s -N < ./autoconf.patch
echo autoconf.patch applied
