      subroutine movecs_lock(geom, basis, g_old, g_new, evals, occ)
*
* $Id: movecs_lock.F,v 1.6 1999-06-16 00:10:13 d3e129 Exp $
*
      implicit none
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
      integer k_e, l_e, k_v, l_v, k_o, l_o
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
     $     ('movecs_lock: ga_create of tmp failed',nbf*nmo)
*ga:1:0
      if (.not. ga_create(mt_dbl, nmo, nmo, 'movecs_lock:sno',
     $     nmo, 0, g_sno)) call errquit
     $     ('movecs_lock: ga_create of Sno failed',nmo*nmo)
      if (.not. ma_push_get(mt_dbl, nmo, 'movecs_lock:e',l_e, k_e))
     $     call errquit('movecs_lock: insufficient memory',nmo)
      if (.not. ma_push_get(mt_dbl, nmo, 'movecs_lock:o',l_o, k_o))
     $     call errquit('movecs_lock: insufficient memory',nmo)
      if (.not. ma_push_get(mt_dbl, nbf, 'movecs_lock:v',l_v, k_v))
     $     call errquit('movecs_lock: insufficient memory',nbf)
c
      call two_index_transf(g_s, g_new, g_old, g_tmp, g_sno)
c
c     Loop thru columns of sno, determine max overlapping vector and copy
c
      call ga_copy(g_new, g_tmp)
      call dcopy(nmo, evals, 1, dbl_mb(k_e), 1)
      call dcopy(nmo, occ, 1, dbl_mb(k_o), 1)
c
      do j = 1+ga_nodeid(), nmo, ga_nnodes() ! Should really use locality
         call ga_get(g_sno, 1, nmo, j, j, dbl_mb(k_v), 1)
         i = idamax(nmo, dbl_mb(k_v), 1)
c
c     new vector i corresponds to old vector j
c
         if (i. ne. j) then
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
      call ga_sync
      call ga_dgop(1, dbl_mb(k_e), nmo, 'max')
      do i = 1, nmo
         if (dbl_mb(k_e+i-1) .ne. 1.0d99) call errquit
     $        ('movecs_lock: maximum overlap test neglected vectors',i)
      enddo
c
c     Clean up
c
      if (.not. ma_pop_stack(l_v)) call errquit('movecs_lock: ma?',0)
      if (.not. ma_pop_stack(l_o)) call errquit('movecs_lock: ma?',0)
      if (.not. ma_pop_stack(l_e)) call errquit('movecs_lock: ma?',0)
      if (.not. ga_destroy(g_sno)) call errquit('movecs_lock: ga?',0)
      if (.not. ga_destroy(g_tmp)) call errquit('movecs_lock: ga?',0)
      if (.not. ga_destroy(g_s))   call errquit('movecs_lock: ga?',0)
c
      call ga_sync
c
      end
