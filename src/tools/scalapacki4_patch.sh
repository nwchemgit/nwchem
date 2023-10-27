#!/usr/bin/env bash
rm -f scalapacki4.patch
cat > scalapacki4.patch <<EOF
--- $1/global/src/scalapack.F.org	2023-10-26 16:35:47.951995538 -0700
+++ $1/global/src/scalapack.F	2023-10-26 16:36:44.840284758 -0700
@@ -1076,13 +1076,12 @@
 c****
       character*1 uplo         ! (input) 'U' or 'L'
       integer     g_A          ! (input/output)
-      logical status
       integer dimA1, dimA2, typeA
       integer me, nproc
       integer n
       integer i, j, hBUF
+      integer j0,j1,i0,i1
       MA_ACCESS_INDEX_TYPE adrBUF
-
 c**** Check Environment
       nproc = ga_nnodes()
       me    = ga_nodeid()
@@ -1101,27 +1100,38 @@
       n = dimA1
       
 c**** Allocate BUF
-      status = ma_push_get(MT_DBL, n, 'BUF', hBUF, adrBUF)
-      if (.not.status)
+      if(.not.ma_push_get(MT_DBL, n, 'BUF', hBUF, adrBUF))
      &     call ga_error(' ga_zeroUL: mem alloc failed BUF ', -1)
+      do i=0,n-1
+         dbl_mb(adrBUF+i)=0d0
+      enddo
       
       call ga_sync()
 
-      do i = me+1, n, nproc
-         call ga_get(g_A, 1, n, i, i, dbl_mb(adrBUF), n)
+      i0=me+1
+      i1=n
+      if (uplo.eq.'L') then
+         i0=me+2
+      elseif (uplo.eq.'U') then 
+         i1=n-1   
+      else
+         call ga_error('ga_symUL: uplo must be L or U ', 1)
+      endif
+      do i = i0, i1, nproc
          if (uplo.eq.'L') then
-c****       case L: make zero the upper triangle            
-            call dcopy(i-1,0.0d0,0, dbl_mb(adrBUF),1)
+c**** case L: make zero the upper triangle            
+            j0=1
+            j1=i-1
          elseif (uplo.eq.'U') then
-c****       case U: make zero the lower triangle            
-            call dcopy(n-i,0.0d0,0, dbl_mb(adrBUF+i),1)
-         else
-            call ga_error('ga_symUL: uplo must be L or U ', 1)
+c**** case U: make zero the lower triangle            
+            j0=i+1
+            j1=n
          endif
-         call ga_put(g_A, 1, n, i, i, dbl_mb(adrBUF), n)
-      end do    !i
+         call ga_put(g_a, j0, j1, i, i, dbl_mb(adrBUF), n)
+      end do                    !i
 c
-      status = ma_pop_stack(hBUF)
+      if(.not.ma_pop_stack(hBUF)) call
+     c     ga_error(' ga_zeroUL: pop_stack failed  ',-1)
       call ga_sync()
       end
 
EOF
patch -p0 -s -N < ./scalapacki4.patch
echo scalapacki4.patch applied
