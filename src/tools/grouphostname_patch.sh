#!/usr/bin/env bash
rm -f grouphostname.patch
cat > grouphostname.patch <<EOF
--- $1/comex/src-mpi-pr/groups.c
+++ $1/comex/src-mpi-pr/groups.c
@@ -453,7 +453,7 @@ void comex_group_init(MPI_Comm comm)
     /* create node comm */
     /* MPI_Comm_split requires a non-negative color,
      * so sort and sanitize */
-    sorted = (long*)malloc(sizeof(host_name_t) * g_state.size);
+    sorted = (host_name_t*)malloc(sizeof(host_name_t) * g_state.size);
     (void)memcpy(sorted, g_state.host, sizeof(host_name_t)*g_state.size);
     qsort(sorted, g_state.size, sizeof(host_name_t), cmpname);
     /* count is number of distinct host IDs that are lower than
EOF
patch -p0 -s -N < grouphostname.patch
echo grouphostname.patch applied
