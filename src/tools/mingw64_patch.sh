#!/bin/bash
rm -f win64_ma_nxtval.patch
cat > win64_ma_nxtval.patch <<EOF
--- $1/ma/ma.c.org	2020-02-28 11:39:16.000000000 -0800
+++ $1/ma/ma.c	2020-06-26 17:03:05.009072842 -0700
@@ -98,11 +98,11 @@
 #if defined(BGQ)
 #define ALIGNMENT	32
 #else
-#define ALIGNMENT	sizeof(long)
+#define ALIGNMENT	sizeof(size_t)
 #endif
 
 /* min size of block split and placed on free list */
-#define MINBLOCKSIZE mai_round((long)(ALIGNMENT + BLOCK_OVERHEAD_FIXED), \
+#define MINBLOCKSIZE mai_round((size_t)(ALIGNMENT + BLOCK_OVERHEAD_FIXED), \
         (ulongi)ALIGNMENT)
 
 /* signatures for guard words */
@@ -114,7 +114,7 @@
  **/
 
 typedef unsigned int Guard;   /* for detection of memory trashing */
-typedef unsigned long ulongi; /* for brevity */
+typedef size_t ulongi; /* for brevity */
 
 /* allocation request for a block */
 typedef struct _AR
@@ -180,7 +180,7 @@
 private void ma_preinitialize(char *caller);
 private Boolean mh2ad(Integer memhandle, AD **adout, BlockLocation location, char *caller);
 private void mh_free(AD *ad);
-private long mai_round(long value, ulongi unit);
+private size_t mai_round(size_t value, ulongi unit);
 private void str_ncopy(char *to, char *from, int maxchars);
 
 /* foreign routines */
@@ -364,7 +364,7 @@
 
 /* compute 0-based index for client_space from AD */
 #define client_space_index(ad) \
-    ((MA_AccessIndex)((long)((ad)->client_space - ma_base[(ad)->datatype]) / \
+    ((MA_AccessIndex)((size_t)((ad)->client_space - ma_base[(ad)->datatype]) / \
                  ma_sizeof[(ad)->datatype]))
 
 /* compute address of guard from AD */
@@ -601,9 +601,9 @@
         (void)printf("unknown");
     else
         (void)printf("%ld",
-            (long)memhandle);
+            (size_t)memhandle);
     (void)printf(", address 0x%lx",
-        (long)ad);
+        (size_t)ad);
 }
 
 /* ------------------------------------------------------------------------- */
@@ -646,11 +646,11 @@
 
     L_client_space = ar->nelem * ma_sizeof[datatype];
 
-    L_gap1 = ((long)B_base
-        - (long)B_address
-        - (long)sizeof(AD)
-        - (long)sizeof(Guard))
-        % (long)ma_sizeof[datatype];
+    L_gap1 = ((size_t)B_base
+        - (size_t)B_address
+        - (size_t)sizeof(AD)
+        - (size_t)sizeof(Guard))
+        % (size_t)ma_sizeof[datatype];
 
     if (L_gap1 < 0)
         L_gap1 += ma_sizeof[datatype];
@@ -668,8 +668,8 @@
      */
 
     if (ma_numalign > 0) {
-      unsigned long mask = (1<<ma_numalign)-1;
-      int diff = ((unsigned long) B_client_space) & mask;
+      size_t mask = (1<<ma_numalign)-1;
+      int diff = ((size_t) B_client_space) & mask;
       
       /* Check that the difference is a multiple of the type size.
        * If so, then we can shift the client space which is already
@@ -701,11 +701,11 @@
      *        - address
      */
 
-    L_gap2 = ((long)B_address
-        - (long)B_client_space
-        - (long)L_client_space
-        - (long)sizeof(Guard))
-        % (long)ALIGNMENT;
+    L_gap2 = ((size_t)B_address
+        - (size_t)B_client_space
+        - (size_t)L_client_space
+        - (size_t)sizeof(Guard))
+        % (size_t)ALIGNMENT;
 
     if (L_gap2 < 0)
         L_gap2 += ALIGNMENT;
@@ -784,8 +784,8 @@
      */
 
     if (ma_numalign > 0) {
-      unsigned long mask = (1<<ma_numalign)-1;
-      int diff = ((unsigned long) B_client_space) & mask;
+      size_t mask = (1<<ma_numalign)-1;
+      int diff = ((size_t) B_client_space) & mask;
       
       /* Check that the difference is a multiple of the type size.
        * If so, then we can shift the client space which is already
@@ -970,9 +970,9 @@
 #define NUMADFIELDS 7
 
     char    *fn[NUMADFIELDS];    /* field names */
-    long    fa[NUMADFIELDS];    /* field addresses */
+    size_t    fa[NUMADFIELDS];    /* field addresses */
     int        i;            /* loop index */
-    long    address;        /* other addresses */
+    size_t    address;        /* other addresses */
 
     /* set field names */
     fn[0] = "datatype";
@@ -984,13 +984,13 @@
     fn[6] = "checksum";
 
     /* set field addresses */
-    fa[0] = (long)(&(ad->datatype));
-    fa[1] = (long)(&(ad->nelem));
-    fa[2] = (long)(&(ad->name));
-    fa[3] = (long)(&(ad->client_space));
-    fa[4] = (long)(&(ad->nbytes));
-    fa[5] = (long)(&(ad->next));
-    fa[6] = (long)(&(ad->checksum));
+    fa[0] = (size_t)(&(ad->datatype));
+    fa[1] = (size_t)(&(ad->nelem));
+    fa[2] = (size_t)(&(ad->name));
+    fa[3] = (size_t)(&(ad->client_space));
+    fa[4] = (size_t)(&(ad->nbytes));
+    fa[5] = (size_t)(&(ad->next));
+    fa[6] = (size_t)(&(ad->checksum));
 
     /* print AD fields to stderr */
     (void)fprintf(stderr, "debug_ad_print:\n");
@@ -1003,19 +1003,19 @@
             fn[i]);
 
     /* print other addresses to stderr */
-    address = (long)guard1(ad);
+    address = (size_t)guard1(ad);
     (void)fprintf(stderr, "\t0x%lx  mod4,8,16=%d,%d,%-2d  guard1\n",
         address,
         address % 4,
         address % 8,
         address % 16);
-    address = (long)ad->client_space;
+    address = (size_t)ad->client_space;
     (void)fprintf(stderr, "\t0x%lx  mod4,8,16=%d,%d,%-2d  client_space\n",
         address,
         address % 4,
         address % 8,
         address % 16);
-    address = (long)guard2(ad);
+    address = (size_t)guard2(ad);
     (void)fprintf(stderr, "\t0x%lx  mod4,8,16=%d,%d,%-2d  guard2\n",
         address,
         address % 4,
@@ -1328,11 +1328,11 @@
         (void)printf("\ttype of elements:\t\t%s\n",
             ma_datatype[ad->datatype]);
         (void)printf("\tnumber of elements:\t\t%ld\n",
-            (long)ad->nelem);
+            (size_t)ad->nelem);
         (void)printf("\taddress of client space:\t0x%lx\n",
-            (long)ad->client_space);
+            (size_t)ad->client_space);
         (void)printf("\tindex for client space:\t\t%ld\n",
-            (long)(client_space_index(ad) + index_base));
+            (size_t)(client_space_index(ad) + index_base));
         (void)printf("\ttotal number of bytes:\t\t%lu\n",
             ad->nbytes);
     }
@@ -1641,7 +1641,7 @@
     {
         (void)sprintf(ma_ebuf,
             "invalid block address (0x%lx) for memhandle %ld",
-            (long)ad, (long)memhandle);
+            (size_t)ad, (size_t)memhandle);
         ma_error(EL_Nonfatal, ET_External, caller, ma_ebuf);
         return MA_FALSE;
     }
@@ -1652,7 +1652,7 @@
         {
             (void)sprintf(ma_ebuf,
                 "invalid checksum for memhandle %ld (name: '%s')",
-                (long)memhandle, ad->name);
+                (size_t)memhandle, ad->name);
             ma_error(EL_Nonfatal, ET_External, caller, ma_ebuf);
             return MA_FALSE;
         }
@@ -1664,7 +1664,7 @@
         {
             (void)sprintf(ma_ebuf,
                 "invalid guard(s) for memhandle %ld (name: '%s')",
-                (long)memhandle, ad->name);
+                (size_t)memhandle, ad->name);
             ma_error(EL_Nonfatal, ET_External, caller, ma_ebuf);
             return MA_FALSE;
         }
@@ -1676,7 +1676,7 @@
         {
             (void)sprintf(ma_ebuf,
                 "memhandle %ld (name: '%s') not in heap",
-                (long)memhandle, ad->name);
+                (size_t)memhandle, ad->name);
             ma_error(EL_Nonfatal, ET_External, caller, ma_ebuf);
             return MA_FALSE;
         }
@@ -1687,7 +1687,7 @@
         {
             (void)sprintf(ma_ebuf,
                 "memhandle %ld (name: '%s') not in stack",
-                (long)memhandle, ad->name);
+                (size_t)memhandle, ad->name);
             ma_error(EL_Nonfatal, ET_External, caller, ma_ebuf);
             return MA_FALSE;
         }
@@ -1699,7 +1699,7 @@
         {
             (void)sprintf(ma_ebuf,
                 "memhandle %ld (name: '%s') not in stack",
-                (long)memhandle, ad->name);
+                (size_t)memhandle, ad->name);
             ma_error(EL_Nonfatal, ET_External, caller, ma_ebuf);
             return MA_FALSE;
         }
@@ -1709,7 +1709,7 @@
         {
             (void)sprintf(ma_ebuf,
                 "memhandle %ld (name: '%s') not top of stack",
-                (long)memhandle, ad->name);
+                (size_t)memhandle, ad->name);
             ma_error(EL_Nonfatal, ET_External, caller, ma_ebuf);
             return MA_FALSE;
         }
@@ -1720,7 +1720,7 @@
         {
             (void)sprintf(ma_ebuf,
                 "memhandle %ld (name: '%s') not in heap or stack",
-                (long)memhandle, ad->name);
+                (size_t)memhandle, ad->name);
             ma_error(EL_Nonfatal, ET_External, caller, ma_ebuf);
             return MA_FALSE;
         }
@@ -1747,7 +1747,7 @@
     {
         (void)sprintf(ma_ebuf,
             "cannot find memhandle for block address 0x%lx",
-            (long)ad);
+            (size_t)ad);
         ma_error(EL_Nonfatal, ET_Internal, "mh_free", ma_ebuf);
     }
     else
@@ -1761,14 +1761,14 @@
  */
 /* ------------------------------------------------------------------------- */
 
-private long mai_round(value, unit)
-    long    value;        /* to round */
+private size_t mai_round(value, unit)
+    size_t    value;        /* to round */
     ulongi    unit;        /* to round to */
 {
     /* voodoo ... */
     unit--;
     value += unit;
-    value &= ~(long)unit;
+    value &= ~(size_t)unit;
     return value;
 }
 
@@ -1825,7 +1825,7 @@
     {
         (void)sprintf(ma_ebuf,
             "invalid datatype: %ld",
-            (long)datatype);
+            (size_t)datatype);
         ma_error(EL_Nonfatal, ET_Internal, "MAi_inform_base", ma_ebuf);
         return MA_FALSE;
     }
@@ -1982,7 +1982,7 @@
     {
         (void)sprintf(ma_ebuf,
             "block '%s', invalid datatype: %ld",
-            name, (long)datatype);
+            name, (size_t)datatype);
         ma_error(EL_Nonfatal, ET_External, "MA_allocate_heap", ma_ebuf);
         return MA_FALSE;
     }
@@ -1992,7 +1992,7 @@
     {
         (void)sprintf(ma_ebuf,
             "block '%s', invalid nelem: %ld",
-            name, (long)nelem);
+            name, (size_t)nelem);
         ma_error(EL_Nonfatal, ET_External, "MA_allocate_heap", ma_ebuf);
         return MA_FALSE;
     }
@@ -2156,7 +2156,7 @@
     {
         (void)sprintf(ma_ebuf,
             "memhandle %ld (name: '%s') not on heap used list",
-            (long)memhandle, ad->name);
+            (size_t)memhandle, ad->name);
         ma_error(EL_Nonfatal, ET_Internal, "MA_free_heap", ma_ebuf);
         return MA_FALSE;
     }
@@ -2216,7 +2216,7 @@
     {
         (void)sprintf(ma_ebuf,
             "block '%s', invalid nelem: %ld",
-            ad->name, (long)nelem);
+            ad->name, (size_t)nelem);
         ma_error(EL_Nonfatal, ET_External, "MA_free_heap_piece", ma_ebuf);
         return MA_FALSE;
     }
@@ -2228,7 +2228,7 @@
 
     if (ma_trace) 
     (void)printf("MA: freeing %ld elements of '%s'\n",
-            (long)nelem, ad->name);
+            (size_t)nelem, ad->name);
 
     /* set AR for data to keep */
     ar.datatype = ad->datatype;
@@ -2327,7 +2327,7 @@
     {
         (void)sprintf(ma_ebuf,
             "invalid datatype: %ld",
-            (long)datatype);
+            (size_t)datatype);
         ma_error(EL_Fatal, ET_External, "MA_get_mbase", ma_ebuf);
         return NULL;
     }
@@ -2472,7 +2472,7 @@
     {
         (void)sprintf(ma_ebuf,
             "invalid datatype: %ld",
-            (long)datatype);
+            (size_t)datatype);
         ma_error(EL_Nonfatal, ET_External, "MA_init", ma_ebuf);
         return MA_FALSE;
     }
@@ -2490,7 +2490,7 @@
         heap_bytes = (nominal_heap * ma_sizeof[datatype]) +
             (DEFAULT_REQUESTS_HEAP * max_block_overhead(datatype));
     }
-    heap_bytes = (unsigned long)mai_round((long)heap_bytes, (ulongi)ALIGNMENT);
+    heap_bytes = (size_t)mai_round((size_t)heap_bytes, (ulongi)ALIGNMENT);
 
     /* compute # of bytes in stack */
     if (nominal_stack < 0)
@@ -2502,7 +2502,7 @@
         stack_bytes = (nominal_stack * ma_sizeof[datatype]) +
             (DEFAULT_REQUESTS_STACK * max_block_overhead(datatype));
     }
-    stack_bytes = (unsigned long)mai_round((long)stack_bytes, (ulongi)ALIGNMENT);
+    stack_bytes = (size_t)mai_round((size_t)stack_bytes, (ulongi)ALIGNMENT);
 
     /* segment consists of heap and stack */
     total_bytes = heap_bytes + stack_bytes;
@@ -2609,7 +2609,7 @@
 
 public Integer MA_inquire_avail(Integer datatype)
 {
-    long    gap_length;    /* # of bytes between heap and stack */
+    size_t    gap_length;    /* # of bytes between heap and stack */
     Integer    nelem_gap;    /* max elements containable in gap */
 
 #ifdef STATS
@@ -2635,7 +2635,7 @@
     {
         (void)sprintf(ma_ebuf,
             "invalid datatype: %ld",
-            (long)datatype);
+            (size_t)datatype);
         ma_error(EL_Fatal, ET_External, "MA_inquire_avail", ma_ebuf);
         return DONTCARE;
     }
@@ -2648,7 +2648,7 @@
      */
 
     /* try space between heap and stack */
-    gap_length = (long)(ma_sp - ma_hp);
+    gap_length = (size_t)(ma_sp - ma_hp);
     if (gap_length > 0)
         nelem_gap = ma_nelem(ma_hp, (ulongi)gap_length, datatype);
     else
@@ -2673,7 +2673,7 @@
 
 public Integer MA_inquire_heap(Integer datatype)
 {
-    long    gap_length;    /* # of bytes between heap and partition */
+    size_t    gap_length;    /* # of bytes between heap and partition */
     Integer    nelem_gap;    /* max elements containable in gap */
     Integer    nelem_frag;    /* max elements containable in any frag */
 
@@ -2700,7 +2700,7 @@
     {
         (void)sprintf(ma_ebuf,
             "invalid datatype: %ld",
-            (long)datatype);
+            (size_t)datatype);
         ma_error(EL_Fatal, ET_External, "MA_inquire_heap", ma_ebuf);
         return DONTCARE;
     }
@@ -2713,7 +2713,7 @@
      */
 
     /* try space between heap and partition */
-    gap_length = (long)(ma_partition - ma_hp);
+    gap_length = (size_t)(ma_partition - ma_hp);
     if (gap_length > 0)
         nelem_gap = ma_nelem(ma_hp, (ulongi)gap_length, datatype);
     else
@@ -2744,7 +2744,7 @@
 
 public Integer MA_inquire_heap_check_stack(Integer datatype)
 {
-    long    gap_length;    /* # of bytes between heap and partition */
+    size_t    gap_length;    /* # of bytes between heap and partition */
     Integer    nelem_gap;    /* max elements containable in gap */
     Integer    nelem_frag;    /* max elements containable in any frag */
 
@@ -2771,7 +2771,7 @@
     {
         (void)sprintf(ma_ebuf,
             "invalid datatype: %ld",
-            (long)datatype);
+            (size_t)datatype);
         ma_error(EL_Fatal, ET_External, "MA_inquire_heap_check_stack", ma_ebuf);
         return DONTCARE;
     }
@@ -2784,7 +2784,7 @@
      */
 
     /* try space between heap and partition or heap and stack */
-    gap_length = min((long)(ma_partition - ma_hp), (long)(ma_sp - ma_hp));
+    gap_length = min((size_t)(ma_partition - ma_hp), (size_t)(ma_sp - ma_hp));
     if (gap_length > 0)
         nelem_gap = ma_nelem(ma_hp, (ulongi)gap_length, datatype);
     else
@@ -2814,7 +2814,7 @@
 
 public Integer MA_inquire_heap_no_partition(Integer datatype)
 {
-    long    gap_length;    /* # of bytes between heap and partition */
+    size_t    gap_length;    /* # of bytes between heap and partition */
     Integer    nelem_gap;    /* max elements containable in gap */
     Integer    nelem_frag;    /* max elements containable in any frag */
 
@@ -2841,7 +2841,7 @@
     {
         (void)sprintf(ma_ebuf,
             "invalid datatype: %ld",
-            (long)datatype);
+            (size_t)datatype);
         ma_error(EL_Fatal, ET_External, "MA_inquire_heap_no_partition", ma_ebuf);
         return DONTCARE;
     }
@@ -2854,7 +2854,7 @@
      */
 
     /* try space between heap and stack */
-    gap_length = (long)(ma_sp - ma_hp);
+    gap_length = (size_t)(ma_sp - ma_hp);
     if (gap_length > 0)
         nelem_gap = ma_nelem(ma_hp, (ulongi)gap_length, datatype);
     else
@@ -2882,7 +2882,7 @@
 
 public Integer MA_inquire_stack(Integer datatype)
 {
-    long    gap_length;    /* # of bytes between partition and stack */
+    size_t    gap_length;    /* # of bytes between partition and stack */
     Integer    nelem_gap;    /* max elements containable in gap */
 
 #ifdef STATS
@@ -2908,7 +2908,7 @@
     {
         (void)sprintf(ma_ebuf,
             "invalid datatype: %ld",
-            (long)datatype);
+            (size_t)datatype);
         ma_error(EL_Fatal, ET_External, "MA_inquire_stack", ma_ebuf);
         return DONTCARE;
     }
@@ -2921,7 +2921,7 @@
      */
 
     /* try space between partition and stack */
-    gap_length = (long)(ma_sp - ma_partition);
+    gap_length = (size_t)(ma_sp - ma_partition);
     if (gap_length > 0)
         nelem_gap = ma_nelem(ma_partition, (ulongi)gap_length, datatype);
     else
@@ -2949,7 +2949,7 @@
 
 public Integer MA_inquire_stack_check_heap(Integer datatype)
 {
-    long    gap_length;    /* # of bytes between partition and stack */
+    size_t    gap_length;    /* # of bytes between partition and stack */
     Integer    nelem_gap;    /* max elements containable in gap */
 
 #ifdef STATS
@@ -2975,7 +2975,7 @@
     {
         (void)sprintf(ma_ebuf,
             "invalid datatype: %ld",
-            (long)datatype);
+            (size_t)datatype);
         ma_error(EL_Fatal, ET_External, "MA_inquire_stack_check_heap", ma_ebuf);
         return DONTCARE;
     }
@@ -2988,7 +2988,7 @@
      */
 
     /* try space between partition and stack or heap and stack */
-    gap_length = min((long)(ma_sp - ma_partition), (long)(ma_sp - ma_hp));
+    gap_length = min((size_t)(ma_sp - ma_partition), (size_t)(ma_sp - ma_hp));
     if (gap_length > 0)
         nelem_gap = ma_nelem(ma_partition, (ulongi)gap_length, datatype);
     else
@@ -3018,7 +3018,7 @@
 
 public Integer MA_inquire_stack_no_partition(Integer datatype)
 {
-    long    gap_length;    /* # of bytes between heap and partition */
+    size_t    gap_length;    /* # of bytes between heap and partition */
     Integer    nelem_gap;    /* max elements containable in gap */
 
 #ifdef STATS
@@ -3044,7 +3044,7 @@
     {
         (void)sprintf(ma_ebuf,
             "invalid datatype: %ld",
-            (long)datatype);
+            (size_t)datatype);
         ma_error(EL_Fatal, ET_External, "MA_inquire_stack_no_partition", ma_ebuf);
         return DONTCARE;
     }
@@ -3057,7 +3057,7 @@
      */
 
     /* try space between heap and stack */
-    gap_length = (long)(ma_sp - ma_hp);
+    gap_length = (size_t)(ma_sp - ma_hp);
     if (gap_length > 0)
         nelem_gap = ma_nelem(ma_hp, (ulongi)gap_length, datatype);
     else
@@ -3101,7 +3101,7 @@
     {
         (void)sprintf(ma_ebuf,
             "memhandle %ld (name: '%s') not on stack used list",
-            (long)memhandle, ad->name);
+            (size_t)memhandle, ad->name);
         ma_error(EL_Nonfatal, ET_Internal, "MA_pop_stack", ma_ebuf);
         return MA_FALSE;
     }
@@ -3255,7 +3255,7 @@
     {
         (void)sprintf(ma_ebuf,
             "block '%s', invalid datatype: %ld",
-            name, (long)datatype);
+            name, (size_t)datatype);
         ma_error(EL_Nonfatal, ET_External, "MA_push_stack", ma_ebuf);
         return MA_FALSE;
     }
@@ -3265,7 +3265,7 @@
     {
         (void)sprintf(ma_ebuf,
             "block '%s', invalid nelem: %ld",
-            name, (long)nelem);
+            name, (size_t)nelem);
         ma_error(EL_Nonfatal, ET_External, "MA_push_stack", ma_ebuf);
         return MA_FALSE;
     }
@@ -3411,7 +3411,7 @@
     {
         (void)sprintf(ma_ebuf,
             "invalid alignment: %ld",
-            (long)value);
+            (size_t)value);
         ma_error(EL_Nonfatal, ET_External, "MA_set_numalign", ma_ebuf);
         return MA_FALSE;
     }
@@ -3451,7 +3451,7 @@
     {
         (void)sprintf(ma_ebuf,
             "invalid datatype: %ld",
-            (long)datatype1);
+            (size_t)datatype1);
         ma_error(EL_Fatal, ET_External, "MA_sizeof", ma_ebuf);
         return DONTCARE;
     }
@@ -3461,7 +3461,7 @@
     {
         (void)sprintf(ma_ebuf,
             "invalid nelem: %ld",
-            (long)nelem1);
+            (size_t)nelem1);
         ma_error(EL_Fatal, ET_External, "MA_sizeof", ma_ebuf);
         return DONTCARE;
     }
@@ -3471,7 +3471,7 @@
     {
         (void)sprintf(ma_ebuf,
             "invalid datatype: %ld",
-            (long)datatype2);
+            (size_t)datatype2);
         ma_error(EL_Fatal, ET_External, "MA_sizeof", ma_ebuf);
         return DONTCARE;
     }
@@ -3520,7 +3520,7 @@
     {
         (void)sprintf(ma_ebuf,
             "invalid datatype: %ld",
-            (long)datatype);
+            (size_t)datatype);
         ma_error(EL_Fatal, ET_External, "MA_sizeof_overhead", ma_ebuf);
         return DONTCARE;
     }
--- ga-5.7.2/ma/memcpy.h.org	2020-02-28 11:39:16.000000000 -0800
+++ ga-5.7.2/ma/memcpy.h	2020-06-26 17:02:36.464960502 -0700
@@ -22,7 +22,7 @@
  **/
 
 /* allocate bytes */
-#define bytealloc(nbytes)    malloc((unsigned long)(nbytes))
+#define bytealloc(nbytes)    malloc((size_t)(nbytes))
 
 /* deallocate bytes */
 #define bytefree(pointer)    (void)free((char *)(pointer))

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

