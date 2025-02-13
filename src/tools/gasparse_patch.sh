#!/usr/bin/env bash
rm -f gasparse.patch
cat > gasparse.patch <<EOF
--- $1/global/src/sparse.array.c
+++ $1/global/src/sparse.array.c
@@ -2738,7 +2738,7 @@
   /* Set up global arrays to hold distributed indices and non-zero values */
   {
     int64_t isize = (rowdim+1)*nblocks;
-    int64_t totalsize = 0;
+    Integer totalsize = 0;
     Integer ndim = 1;
     Integer *offset = (Integer*)malloc(nprocs*sizeof(Integer));
     Integer *tmp = (Integer*)malloc(nprocs*sizeof(Integer));
EOF
patch -p0 -s -N < gasparse.patch
echo gasparse.patch applied
