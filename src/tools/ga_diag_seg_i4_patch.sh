#!/usr/bin/env bash
rm -f ga_diag_seg_i4.patch
cat > ga_diag_seg_i4.patch <<EOF
--- $1/global/src/ga_diag_seq.F.org	2023-10-30 11:07:22.963036506 -0700
+++ $1/global/src/ga_diag_seq.F	        2023-10-28 20:21:34.000000000 -0700
@@ -28,11 +28,12 @@
       integer l_wrk
       MA_ACCESS_INDEX_TYPE k_wrk
       integer n2
+      INTGR4 n2_4
 #endif
       integer l_a, l_s
       MA_ACCESS_INDEX_TYPE k_a, k_s
       integer dim1, dim2, type, me
-      INTGR4 n4,ierr
+      INTGR4 n4,ierr,one4
       logical status
 c
 c
@@ -94,8 +95,10 @@
          call rsg(n, n, dbl_mb(k_a), dbl_mb(k_s), evals, 1,
      $        dbl_mb(k_v), dbl_mb(k_fv1), dbl_mb(k_fv2), ierr)
 #else
-         call dsygv(1,'V','U',n4,dbl_mb(k_a),n,dbl_mb(k_s),n4,
-     $              evals,dbl_mb(k_wrk),n2, ierr)
+         one4=1
+         n2_4=n2
+         call dsygv(one4,'V','U',n4,dbl_mb(k_a),n4,dbl_mb(k_s),n4,
+     $              evals,dbl_mb(k_wrk),n2_4, ierr)
          if (ierr.ne.0)
      $       call ga_error('ga_diag_seq: dsygv failed',ierr)
 c We used to copy to preserve code symmetry with EISPACK
EOF
patch -p0 -s -N < ./ga_diag_seg_i4.patch
echo ga_diag_seg_i4.patch applied
