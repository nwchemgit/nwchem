      subroutine movecs_lock(geom, basis, g_old, g_new, evals, occ)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
      integer geom, basis
      integer g_old             ! [input] Old movecs unchanged
      integer g_new             ! [input/output] New movecs reordered
      double precision evals(*) ! [input/output] Evals reordered
      double precision occ(*)   ! [input/output] Occupation reordered
c
c     Reorder the new MOs and evals so that the MOs correspond to
c     the old MOs according to a maximum overlap test
c
c     Form Sold-new = Cold*S*Cnew and find max element in each column.
c
      integer type, nbf, nmo
      integer g_s, g_tmp, g_sno
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c
      integer k_e, l_e, k_v, l_v, k_o, l_o,l_t,k_t
      integer i, j
      integer idamax
      external idamax
c
c     Get the overlap matrix
c
      g_s = ga_create_atom_blocked(geom, basis, 'movecs_lock:s')
      call ga_zero(g_s)
      call int_1e_ga(basis, basis, g_s, 'overlap',.false.)
c
c     Allocate workspace and form Snew-old
c
      call ga_inquire(g_old, type, nbf, nmo)
*ga:1:0
      if (.not. ga_create(mt_dbl, nbf, nmo, 'movecs_lock:tmp',
     $     nbf, 0, g_tmp)) call errquit
     $     ('movecs_lock: ga_create of tmp failed',nbf*nmo, GA_ERR)
*ga:1:0
      if (.not. ga_create(mt_dbl, nmo, nmo, 'movecs_lock:sno',
     $     nmo, 0, g_sno)) call errquit
     $     ('movecs_lock: ga_create of Sno failed',nmo*nmo, GA_ERR)
      if (.not. ma_push_get(mt_dbl, nmo, 'movecs_lock:e',l_e, k_e))
     $     call errquit('movecs_lock: insufficient memory',nmo, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nmo, 'movecs_lock:o',l_o, k_o))
     $     call errquit('movecs_lock: insufficient memory',nmo, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf, 'movecs_lock:v',l_v, k_v))
     $     call errquit('movecs_lock: insufficient memory',nbf, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf, 'movecs_lock:t',l_t, k_t))
     $     call errquit('movecs_lock: insufficient memory',nbf, MA_ERR)
c
      call two_index_transf(g_s, g_new, g_old, g_tmp, g_sno)
c
c     Loop thru columns of sno, determine max overlapping vector and copy
c
      call ga_copy(g_new, g_tmp)
      call dcopy(nmo, evals, 1, dbl_mb(k_e), 1)
      call dcopy(nmo, occ, 1, dbl_mb(k_o), 1)
      call ga_sync
c
      if(ga_nodeid().eq.0) then
      do j = 1, nmo 
         call ga_get(g_sno, 1, nmo, j, j, dbl_mb(k_v), 1)
         i = idamax(nmo, dbl_mb(k_v), 1)
c
c     new vector i corresponds to old vector j
c
         if (i. ne. j) then
C
C         invalidate column that you have taken
C
          call dcopy(nmo, 0.0d0,0, dbl_mb(k_t), 1)
         call ga_put(g_sno, i, i, 1, nmo, dbl_mb(k_t), 1)

            if (util_print('movecs_lock', print_default)) then
               write(6,*) ' movecs_lock', i, j
               call util_flush(6)
            endif
         endif
         call ga_get(g_tmp, 1, nbf, i, i, dbl_mb(k_v), 1)
         call ga_put(g_new, 1, nbf, j, j, dbl_mb(k_v), 1)
         evals(j) = dbl_mb(k_e+i-1)
         occ(j)   = dbl_mb(k_o+i-1)
         dbl_mb(k_e+i-1) = 1.0d99 ! Larger than any valid eval
      enddo
c
c     Verify that all vectors have been used
c
      do i = 1, nmo
         if (dbl_mb(k_e+i-1) .ne. 1.0d99) call errquit
     $        ('movecs_lock: maximum overlap test neglected vectors',i,
     &       UNKNOWN_ERR)
      enddo
            if (util_print('movecs_lock', print_default)) then
               call util_flush(6)
            endif
      endif
c
c     Clean up
c
      if (.not. ma_chop_stack(l_e)) call errquit('movecs_lock: ma?',0,
     &       MA_ERR)
      if (.not. ga_destroy(g_sno)) call errquit('movecs_lock: ga?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_tmp)) call errquit('movecs_lock: ga?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_s))   call errquit('movecs_lock: ga?',0,
     &       GA_ERR)
c
      call ga_sync
c
      end
