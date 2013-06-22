C
C     canorg_check.F
C
C     Checks that canonical transformations are correct:
C
C     (1) X^+ S X = 1_m    (X orthogonalizes S)
C     (2) Y^+ X = 1_m      (Y^+ is left inverse of X)
C     (3) X^+ Y = 1_m      (Y is the right inverse of X^+)
C
C     True but not checked:
C         Y X^+ = 1_n
C         X Y^+ = 1_n
C
C     !!! DISABLED (4) Tr[P'] = Nelec   (P->P' ok)  (only for alpha spin part)
C     !!! DISABLED (5) P' -> P works
C
C     The m subscript means the identity matrix is of that size.
C
      subroutine canorg_check (params, g_s)
      implicit none


#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "matutils.fh"
#include "rt_tddft.fh"


C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_s
C      integer, intent(in)           :: g_ddens_ao  !dens mat in AO bas; alpha, beta, or tot
C      integer, intent(in)           :: numelec     !number of elec stored in g_ddens_so 

      
C     == Parameters ==
      character(*), parameter :: pname = "canorg_check: "


C     == Variables ==
      integer g_zdens_ao        !n x n
      integer g_zsr_nn          !n x n scratch
      integer g_zsr_mm          !m x m scratch
      integer g_zsr_nm          !n x m scratch



C
C     Allocate scratch space.
C
      if (.not. ga_create (mt_dcpl, params%ns_mo, params%ns_mo,
     $     "scr mxm", 0, 0, g_zsr_mm))
     $     call errquit (pname//"scr_mm alloc failed", 0, GA_ERR)

      if (.not. ga_create (mt_dcpl, params%ns_ao, params%ns_mo,
     $     "scr nxm", 0, 0, g_zsr_nm))
     $     call errquit (pname//"scr_nm alloc failed", 0, GA_ERR)

      if (.not. ga_create (mt_dcpl, params%ns_ao, params%ns_ao,
     $     "scr nxn", 0, 0, g_zsr_nn))
     $     call errquit (pname//"scr_nn alloc failed", 0, GA_ERR)

      if (.not. ga_create (mt_dcpl, params%ns_ao, params%ns_ao,
     $     "zdens_ao", 0, 0, g_zdens_ao))
     $     call errquit (pname//"zdens_ao alloc failed", 0, GA_ERR)


C
C     (1) X^+ S X = 1_m
C
      call ga_zero (g_zsr_mm)
      call convert_d2z (1d0, g_s, 0d0, g_s, g_zsr_nn)
      call canorg_trans (params, "S", "AO->MO", g_zsr_nn, g_zsr_mm)
      if (.not. mat_is_ident (g_zsr_mm, params%tol_zero))
     $     call errquit (pname // "X failed to orthz S", 0, 0)


C
C     (2) Y^+ X = 1_m
C     
      call ga_zero (g_zsr_mm)
      call ga_zgemm ("T", "N", params%ns_mo, params%ns_mo, params%ns_ao,
     $     z1, params%g_zY, params%g_zX, z0, g_zsr_mm)
      if (.not. mat_is_ident (g_zsr_mm, params%tol_zero))
     $     call errquit (pname // "Y^+ X /= 1_m", 0, 0)


C
C     (3) X^+ Y = 1_m
C
      call ga_zero (g_zsr_mm)
      call ga_zgemm ("T", "N", params%ns_mo, params%ns_mo, params%ns_ao,
     $     z1, params%g_zX, params%g_zY, z0, g_zsr_mm)
      if (.not. mat_is_ident (g_zsr_mm, params%tol_zero))
     $     call errquit (pname // "X^+ Y /= 1_m", 0, 0)


C
C     (4) Tr[P'(alpha)] = Nelec(alpha), note that P'(alpha) is full dens
C     mat if doing closed shell (beta part is zero in this case).
C     
c$$$      call ga_zero (g_zsr_mm)
c$$$      call ga_zero (g_zsr_nn)
c$$$
c$$$      call convert_d2z (1d0, g_ddens_ao, 0d0, g_ddens_ao, g_zdens_ao)
c$$$      call canorg_trans (params, "P", "AO->MO", g_zdens_ao, g_zsr_mm)
c$$$
c$$$      if (.not.mat_is_tracex(g_zsr_mm, dble (numelec),params%tol_zero))
c$$$     $     call errquit (pname//
c$$$     $     "Tr[P'] /= Nelec; bad canorg transforms", 0, 0)


C     
C     (5) Check that we can transform P'->P to get back ground state
C     dens mat.  At this point, g_zsr_mm still holds P' from previous
C     transform.
C
c$$$       call ga_zero (g_zsr_nn)
c$$$       call canorg_trans (params, "P", "MO->AO", g_zsr_mm, g_zsr_nn)
c$$$
c$$$       if (.not. mat_is_same (g_zdens_ao, g_zsr_nn, params%tol_zero))
c$$$     $      call errquit (pname//"failed to convert P'->P", 0, 0)

       
C
C     Clean up.
C
      if (.not. ga_destroy (g_zsr_nn))
     $     call errquit (pname//"destroy zscr_nn failed", 0, GA_ERR)

      if (.not. ga_destroy (g_zsr_mm))
     $     call errquit (pname//"destroy zscr_mm failed", 0, GA_ERR)

      if (.not. ga_destroy (g_zsr_nm))
     $     call errquit (pname//"destroy zscr_nm failed", 0, GA_ERR)

      if (.not. ga_destroy (g_zdens_ao))
     $     call errquit (pname//"destroy zdens_ao failed", 0, GA_ERR)

      end subroutine
c $Id$
