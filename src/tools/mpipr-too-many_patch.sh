#!/usr/bin/env bash
rm -f mpipr-too-many.patch
cat > mpipr-too-many.patch <<EOF
--- $1/comex/src-mpi-pr/comex.c.org	2023-07-12 19:10:15.711084258 -0700
+++ $1/comex/src-mpi-pr/comex.c	2023-07-12 19:10:21.851117110 -0700
@@ -51,7 +51,7 @@
 
 #define XSTR(x) #x
 #define STR(x) XSTR(x)
-
+#define MIN(a, b) (((b) < (a)) ? (b) : (a))
 /* data structures */
 
 typedef enum {
@@ -358,6 +358,8 @@
 static int devshm_initialized = 0;
 static long devshm_fs_left = 0;
 static long devshm_fs_initial = 0;
+static long counter_open_fds = 0;
+STATIC void count_open_fds(void);
 
 int comex_init()
 {
@@ -4746,6 +4748,7 @@
 
     /* set the size of my shared memory object */
     check_devshm(fd, size);
+    count_open_fds();
     retval = ftruncate(fd, size);
     if (-1 == retval) {
         perror("_shm_create: ftruncate");
@@ -7248,3 +7251,29 @@
 #endif
 #endif
 }
+
+STATIC void count_open_fds(void) {
+#ifdef __linux__
+  /* check only every 100 ops && rank == 1 */
+  counter_open_fds += 1;
+  if (counter_open_fds % 100 == 0 && g_state.rank == MIN(1,g_state.node_size)) {
+    FILE *f = fopen("/proc/sys/fs/file-nr", "r");
+
+    long nfiles, unused, maxfiles;
+    fscanf(f, "%ld %ld %ld", &nfiles, &unused, &maxfiles);
+#ifdef DEBUGSHM
+    if(nfiles % 1000 == 0) fprintf(stderr," %d: no. open files = %ld maxfiles = %ld\n", g_state.rank, nfiles, maxfiles);
+#endif
+    if(nfiles > (maxfiles/100)*80) {
+      printf(" %d: running out of files; files = %ld  maxfiles = %ld \n", g_state.rank, nfiles, maxfiles);
+#if PAUSE_ON_ERROR
+      fprintf(stderr,"%d(%d): too many open files\n",
+	      g_state.rank,  getpid());
+      pause();
+#endif
+      comex_error("count_open_fds: too many open files", -1);
+  }
+    fclose(f);
+  }
+#endif
+}
EOF
patch -p0 -s -N < ./mpipr-too-many.patch
echo mpipr-too-many.patch applied
