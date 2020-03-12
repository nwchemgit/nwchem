C====================================================================
C
C     rt_tddft_cs_checkpoint.F
C
C     Closed-shell checking
C
C     Does a series of checks including:
C     - Hermicity of F and P
C     - Tr[P] = N_elec
C     - P^2 = P
C
      subroutine rt_tddft_cs_checkpoint (params, g_zfock_ao, g_zdens_ao,
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
      integer, intent(in)           :: g_zdens_ao, g_zfock_ao !ns_ao x ns_ao
      integer, intent(in)           :: g_zdens_mo, g_zfock_mo !ns_mo x ns_mo

 
C     == Parameters ==
      character(*), parameter   :: pname = "rt_tddft_cs_checkpoint: "


C     == Variables ==
      integer g_zscr            !complex scratch array, ns_mo x ns_mo
      

      call rt_tddft_cs_confirm (params,'rt_tddft_cs_checkpoint.F')


      if (.not. ga_duplicate(g_zfock_mo, g_zscr, "zscr"))
     $     call errquit(pname//"GA allocation failed", 0,GA_ERR)
      
      
C
C     Check symmetries of F and P in both AO and MO basis.
C
      if (.not. mat_is_hermitian(g_zdens_ao, params%tol_zero))
     $     call errquit (pname//"P not Hermitian", 0, 0)

      if (.not. mat_is_hermitian(g_zfock_ao, params%tol_zero))
     $     call errquit (pname//"F not Hermitian", 0, 0)

      if (.not. mat_is_hermitian(g_zdens_mo, params%tol_zero))
     $     call errquit (pname//"P' not Hermitian", 0, 0)

      if (.not. mat_is_hermitian(g_zfock_mo, params%tol_zero))
     $     call errquit (pname//"F' not Hermitian", 0, 0)

      
C
C     Check that Tr[P'] = Nelec .
C
      if (.not. mat_is_tracex(g_zdens_mo, dble(params%ntotel),
     $     params%tol_zero))
     $     call errquit (pname//"Tr[P'] /= Nelec", 0, 0)


C      
C     Check that 1/2 P' is idempotent, see Szabo and Ostlund pg 139.
C     Note this is only true for closed shell case.
C     
      call ga_zero (g_zscr)
      call ga_copy (g_zdens_mo, g_zscr)
      call ga_scale (g_zscr, zhalf)
      
      if (.not. mat_is_idempotent(g_zscr, params%tol_zero))
     $     call errquit (pname//"P' not idempotent", 0, 0)

      if (.not. ga_destroy (g_zscr))
     $     call errquit (pname//"GA destroy failed", 0, GA_ERR)
      
     
      end subroutine




c $Id$
