#!/usr/bin/env bash                                                                                                                              
rm -f pt.patch
cat > pt.patch <<EOF
--- $1/comex/src-mpi-pt/comex.c
+++ $1/comex/src-mpi-pt/comex.c
@@ -457,8 +457,8 @@ int comex_finalize()
     free(fence_array);
 
     free(nb_state);
-#ifdef DEBUG
-    printf(" %d freed nb_state ptr %p \n", g_state.rank, nb_state);
+#if DEBUG
+    printf(" %d freed nb_state ptr %p \n", g_state.rank, nb_state);
 #endif
 
     MPI_Barrier(g_state.comm);
EOF
patch -p0 -s -N < pt.patch
echo pt.patch applied
