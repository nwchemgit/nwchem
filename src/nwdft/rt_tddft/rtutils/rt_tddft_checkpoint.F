C====================================================================
C
C     rt_tddft_cs_checkpoint.F
C
C     Does a series of checks on the fock and density matrices,
C     including:
C      
C     - Hermicity of F and P
C     - Tr[P] = N_elec
C     - P^2 = P
C
C     Note, if this is for an open shell calculation, you can call this
C     independently for the alpha and beta matrices.
C
C
      logical function rt_tddft_checkpoint (params, nmats, 
     $     g_zfock_ao, g_zdens_ao, g_zfock_mo, g_zdens_mo)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "matutils.fh"
#include "rt_tddft.fh"
      
      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: nmats
      integer, intent(in)           :: g_zdens_ao(nmats)
      integer, intent(in)           :: g_zfock_ao(nmats)         !ns_ao x ns_ao
      integer, intent(in)           :: g_zdens_mo(nmats)
      integer, intent(in)           :: g_zfock_mo(nmats)         !ns_mo x ns_mo
C      integer, intent(in)           :: nelec                   !# e- in this dens mat

 
C     == Parameters ==
      character(*), parameter   :: pname = "rt_tddft_checkpoint: "
      character(len=*), parameter :: fmt = "(a,es8.2,a)"

C     == Variables ==
      integer g_zscr            !complex scratch array, ns_mo x ns_mo
      logical status
      integer imat
      integer me
      

      status = .true.
      me = ga_nodeid ()

c$$$      if (.not. ga_duplicate(g_zfock_mo, g_zscr, "zscr"))
c$$$     $     call errquit(pname//"GA allocation failed", 0,GA_ERR)
      
      
C
C     Check symmetries of F and P in both AO and MO basis.
C
c$$$      if (.not. mat_is_hermitian(g_zdens_ao, params%tol_zero))
c$$$     $     call errquit (pname//"P not Hermitian", 0, 0)
c$$$
c$$$      if (.not. mat_is_hermitian(g_zfock_ao, params%tol_zero))
c$$$     $     call errquit (pname//"F not Hermitian", 0, 0)
c$$$
c$$$      if (.not. mat_is_hermitian(g_zdens_mo, params%tol_zero))
c$$$     $     call errquit (pname//"P' not Hermitian", 0, 0)
c$$$
c$$$      if (.not. mat_is_hermitian(g_zfock_mo, params%tol_zero))
c$$$     $     call errquit (pname//"F' not Hermitian", 0, 0)

      do imat = 1, nmats
         if (.not. mat_is_hermitian(g_zdens_ao(imat),
     $        params%tol_zero)) then
            status = .false.
            if (me.eq.0) write (luout, fmt) pname//
     $           "P not Hermitian vs tolerance ",
     $           params%tol_zero, ""
         endif

         if (.not. mat_is_hermitian(g_zfock_ao(imat),
     $        params%tol_zero)) then
            status = .false.
            if (me.eq.0) write (luout, fmt) pname//
     $           "F not Hermitian vs tolerance ",
     $           params%tol_zero, ""
         endif
         
         if (.not. mat_is_hermitian(g_zdens_mo(imat),
     $        params%tol_zero)) then
            status = .false.
            if (me.eq. 0) write (luout, fmt) pname//
     $           "P' not Hermitian vs tolerance ",
     $           params%tol_zero, ""
         endif
         
         if (.not. mat_is_hermitian(g_zfock_mo(imat),
     $        params%tol_zero)) then
            status = .false.
            if (me.eq.0) write (luout, fmt) pname//
     $           "F' not Hermitian vs tolerance ",
     $           params%tol_zero, ""
         endif
      enddo
      
      
      rt_tddft_checkpoint = status


CXXX  [KAL]: DISABLE THESE CHECKS FOR NOW
C
C     Check that Tr[P'] = Nelec .
C
c$$$      if (.not. mat_is_tracex(g_zdens_mo, dble(nelec),
c$$$     $     params%tol_zero))
c$$$     $     call errquit (pname//"Tr[P'] /= Nelec", 0, 0)


C      
C     Check idempotency of the density matrix in MO basis.
C
C     For a closed shell case, 1/2 P' is idempotent (see Szabo and
C     Ostlund pg 139).
C
C     For an open shell case, P' is idempotent.
C     
c$$$      call ga_zero (g_zscr)
c$$$      call ga_copy (g_zdens_mo, g_zscr)
c$$$      call ga_scale (g_zscr, zhalf)
c$$$      
c$$$      if (.not. mat_is_idempotent(g_zscr, params%tol_zero))
c$$$     $     call errquit (pname//"P' not idempotent", 0, 0)

c$$$      if (.not. ga_destroy (g_zscr))
c$$$     $     call errquit (pname//"GA destroy failed", 0, GA_ERR)
      
     
      end function




c $Id$
