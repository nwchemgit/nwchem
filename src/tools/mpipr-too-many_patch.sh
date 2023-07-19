#!/usr/bin/env bash
rm -f mpipr-too-many.patch
cat > mpipr-too-many.patch <<EOF
--- $1/comex/src-mpi-pr/comex.c.org	2023-07-12 19:10:15.711084258 -0700
+++ $1/comex/src-mpi-pr/comex.c	2023-07-12 19:10:21.851117110 -0700
@@ -358,6 +358,7 @@
 static int devshm_initialized = 0;
 static long devshm_fs_left = 0;
 static long devshm_fs_initial = 0;
+STATIC void count_open_fds(void);
 
 int comex_init()
 {
@@ -7215,6 +7216,7 @@
 	    g_state.rank, g_state.node_size, devshm_fs_initial/CONVERT_TO_M, (long) ufs_statfs.f_bsize, (long)  g_state.node_size);
 #endif
   }
+  count_open_fds();
   //  if (size > 0) {
     newspace = (long) ( size*(g_state.node_size -1));
     //  }else{
@@ -7248,3 +7250,24 @@
 #endif
 #endif
 }
+
+STATIC void count_open_fds(void) {
+  FILE *f = fopen("/proc/sys/fs/file-nr", "r");
+
+  long nfiles, unused, maxfiles;
+  fscanf(f, "%ld %ld %ld", &nfiles, &unused, &maxfiles);
+#ifdef DEBUGSHM
+  if(nfiles % 1000 == 0) fprintf(stderr," %d: no. open files = %ld maxfiles = %ld\n", g_state.rank, nfiles, maxfiles);
+#endif
+  long mylimit = (maxfiles/100)*90;
+    if(nfiles > (maxfiles/100)*90) {
+      printf(" %d: running out of files; files = %ld  maxfiles = %ld\n", g_state.rank, nfiles, maxfiles);
+#if PAUSE_ON_ERROR
+    fprintf(stderr,"%d(%d): too many open files\n",
+            g_state.rank,  getpid());
+    pause();
+#endif
+    comex_error("count_open_fds: too many open files", -1);
+  }
+  fclose(f);
+}
EOF
patch -p0 -s -N < mpipr-too-many.patch
echo mpipr-too-many.patch applied
