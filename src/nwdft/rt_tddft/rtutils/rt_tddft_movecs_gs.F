C
C     rt_tddft_movecs_gs.F
C
C     Compute the *ground state* MO vectors for the supplied Fock matrix
C     (in MO basis).  Note, you can either send the closed shell total
C     Fock mat, or the individual alpha or beta spin Fock matrices.
C
C     Note we throw away the eigenvalues (can easily be returned if
C     desired).
C
C
      subroutine rt_tddft_movecs_gs (params, g_zfock_mo, g_movecs_gs)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "geom.fh"
#include "util.fh"
#include "rt_tddft.fh"
#include "matutils.fh"

C     == Inputs ==
      type (rt_params_t), intent(in) :: params
      integer, intent(in)            :: g_zfock_mo      !complex MO Fock
      integer, intent(in)            :: g_movecs_gs     !real-valued MO eigenvectors


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_movecs_gs: "

      
C     == Variables ==
      integer levals, ievals
      integer g_fockre_mo, g_fockim_mo
      integer dtype, dim1, dim2
C      integer i, g_q, g_l, g_tmp1, g_tmp2

      
C
C     Allocations.
C
      if (.not. ga_create (mt_dbl, params%ns_mo, params%ns_mo, "Re[F]",
     $     0, 0, g_fockre_mo))
     $     call errquit (pname//"alloc failed", 0, GA_ERR)

      if (.not. ga_create (mt_dbl, params%ns_mo, params%ns_mo, "Im[F]",
     $     0, 0, g_fockim_mo))
     $     call errquit (pname//"alloc failed", 0, GA_ERR)

c$$$      if (.not. ga_create (mt_dbl, params%ns_mo, params%ns_mo, "tmp",
c$$$     $     0, 0, g_tmp1)) !!!REMOV
c$$$     $     call errquit (pname//"alloc failed", 0, GA_ERR)
c$$$
c$$$      if (.not. ga_create (mt_dbl, params%ns_mo, params%ns_mo, "tmp",
c$$$     $     0, 0, g_tmp2)) !!REMOVE
c$$$     $     call errquit (pname//"alloc failed", 0, GA_ERR)
c$$$
      if (.not.ma_push_get (mt_dbl,params%ns_mo,"evals",levals,ievals))
     $     call errquit (pname//"failed to push evals", 0, GA_ERR)


C
C     Diagonalize.
C     
      call ga_zero (g_fockre_mo)
      call ga_zero (g_fockim_mo)
      call convert_z2d (g_zfock_mo, g_fockre_mo, g_fockim_mo)
      call ga_diag_std (g_fockre_mo, g_movecs_gs, dbl_mb(ievals))


C     XXX
C     XXX TESTING INVERSION
C     XXX
c$$$      call ga_zero (g_tmp1)  !will have evecs of C
c$$$      call ga_diag_std (g_movecs_gs, g_tmp1, dbl_mb(ievals))  !evals overwritten


C     (check that eigenvector matrix of the eigenvector matrix is unitary)
c$$$      call ga_zero (g_zfock_mo) !TMP
c$$$      call convert_d2z (1d0, g_tmp1, 0d0, g_tmp1, g_zfock_mo)
c$$$      if (.not. mat_is_unitary (g_zfock_mo, 1d-5))
c$$$     $     call errquit (pname//"Q is not unitary",0,0)


      
C     (diagonal mat of eigenvalues)
c$$$      call ga_zero (g_tmp2)
c$$$      do i = 1, params%ns_mo
c$$$         call ga_put (g_tmp2, i, i, i, i, -1d0*dbl_mb(ievals+i-1), 1)
c$$$      enddo
c$$$
c$$$      call ga_zero (g_fockre_mo) !recycle !XXX DELETE
c$$$      call ga_zero (g_fockim_mo) !recycle !XXX DELETE
c$$$
c$$$      call ga_dgemm ("N", "T", params%ns_mo, params%ns_mo, params%ns_mo,
c$$$     $     1d0, g_tmp2, g_tmp1, 0d0, g_fockre_mo)  ! L Q^T
c$$$
c$$$      call ga_dgemm ("N", "N", params%ns_mo, params%ns_mo, params%ns_mo,
c$$$     $     1d0, g_tmp1, g_fockre_mo, 0d0, g_fockim_mo) ! Q L Q^T
c$$$
c$$$      if (.not. mat_is_same (g_movecs_gs, g_fockim_mo, 1d-4))
c$$$     $     call errquit ("NOT SAME",0,0)
c$$$      call halt ()


      if (.not. ga_destroy (g_fockre_mo))
     $     call errquit (pname//"destroy failed", 0, GA_ERR)

      if (.not. ga_destroy (g_fockim_mo))
     $     call errquit (pname//"destroy failed", 0, GA_ERR)

      if (.not.ma_pop_stack (levals))
     $     call errquit (pname//"failed to pop evals", 0, GA_ERR)

      end subroutine
c $Id$
