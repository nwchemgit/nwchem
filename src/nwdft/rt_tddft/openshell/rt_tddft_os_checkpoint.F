C====================================================================
C
C     rt_tddft_os_checkpoint.F
C
C     Open-shell checking
C
C     Does a series of checks including:
C     - Hermicity of F and P
C     - Tr[P] = N_elec
C     - P^2 = P
C
      subroutine rt_tddft_os_checkpoint (params, g_zfock_ao, g_zdens_ao,
     $     g_zfock_mo, g_zdens_mo)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "matutils.fh"
#include "rt_tddft.fh"
      
      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_zdens_ao(2), g_zfock_ao(2) !ns_ao x ns_ao; alpha,beta
      integer, intent(in)           :: g_zdens_mo(2), g_zfock_mo(2) !ns_mo x ns_mo; alpha,beta

 
C     == Parameters ==
      character(*), parameter   :: pname = "rt_tddft_os_checkpoint: "


C     == Variables ==
      integer g_zscr            !complex scratch array, ns_mo x ns_mo
      

      call rt_tddft_os_confirm (params,'rt_tddft_os_checkpoint.F')


      if (.not. ga_duplicate(g_zfock_mo, g_zscr, "zscr"))
     $     call errquit(pname//"GA allocation failed", 0,GA_ERR)
      
      
C
C     Check symmetries of F and P in both AO and MO basis.
C
      if (.not. mat_is_hermitian(g_zdens_ao(1), params%tol_zero))
     $     call errquit (pname//"P alpha not Hermitian", 0, 0)
      if (.not. mat_is_hermitian(g_zdens_ao(2), params%tol_zero))
     $     call errquit (pname//"P beta not Hermitian", 0, 0)

      if (.not. mat_is_hermitian(g_zfock_ao(1), params%tol_zero))
     $     call errquit (pname//"F alpha not Hermitian", 0, 0)
      if (.not. mat_is_hermitian(g_zfock_ao(2), params%tol_zero))
     $     call errquit (pname//"F beta not Hermitian", 0, 0)

      if (.not. mat_is_hermitian(g_zdens_mo(1), params%tol_zero))
     $     call errquit (pname//"P' alpha not Hermitian", 0, 0)
      if (.not. mat_is_hermitian(g_zdens_mo(2), params%tol_zero))
     $     call errquit (pname//"P' beta not Hermitian", 0, 0)

      if (.not. mat_is_hermitian(g_zfock_mo(1), params%tol_zero))
     $     call errquit (pname//"F' alpha not Hermitian", 0, 0)
      if (.not. mat_is_hermitian(g_zfock_mo(2), params%tol_zero))
     $     call errquit (pname//"F' beta not Hermitian", 0, 0)

      
C
C     Check that Tr[P' alpha] = Nelec alpha, same for beta.
C
      if (.not. mat_is_tracex(g_zdens_mo(1), dble(params%nalpha),
     $     params%tol_zero))
     $     call errquit (pname//"Tr[P'alpha] /= Nalpha", 0, 0)

      if (.not. mat_is_tracex(g_zdens_mo(2), dble(params%nbeta),
     $     params%tol_zero))
     $     call errquit (pname//"Tr[P'beta] /= Nbeta", 0, 0)


C      
C     Check that P' alpha,beta are idempotent.
C     
      call ga_zero (g_zscr)
      call ga_copy (g_zdens_mo(1), g_zscr)
      if (.not. mat_is_idempotent(g_zscr, params%tol_zero))
     $     call errquit (pname//"P' alpha not idempotent", 0, 0)

      call ga_zero (g_zscr)
      call ga_copy (g_zdens_mo(2), g_zscr)
      if (.not. mat_is_idempotent(g_zscr, params%tol_zero))
     $     call errquit (pname//"P' beta not idempotent", 0, 0)

      if (.not. ga_destroy (g_zscr))
     $     call errquit (pname//"GA destroy failed", 0, GA_ERR)
     
      end subroutine




c $Id$
