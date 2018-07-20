C
C     zfock_cs_exchim.F
C
C     Computes the imaginary part of the exchange for a complex Fock
C     matrix.
C
      subroutine zfock_cs_exchim (params, g_densim, Excim, g_fxim)
      implicit none

#include "bas.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "util.fh"
#include "cdft.fh"
#include "matutils.fh"
#include "rt_tddft.fh"


C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_densim

      
C     == Outputs ==
      double precision, intent(out) :: Excim(2)
      integer, intent(in)           :: g_fxim            !imag part exchange potential

      
C     == Parameters ==
      character(*), parameter     :: pname = "zfock_cs_exchim: "
      

C     == Variables ==
      double precision jfac, kfac
C      double precision tol2e
      double precision elapsed


      call rt_tddft_cs_confirm (params,'zfock_cs_exchim.F')

      
      if (params%prof) call prof_start (elapsed)


C     == Initializations ==
C      tol2e = 10.d0**(-itol2e)
      Excim(1) = 0d0

      call ga_zero (g_fxim)

C      
C     Note that for imaginary part of exchange we pass antisymm (.true.)
C     to the Fock builder, and we explicitly antisymmetrize the result.
C

      jfac = 0.0d0
      kfac = -0.5d0*xfac(1)

      call fock_2e(geom, ao_bas_han, 1, jfac, kfac,
     $     params%tol2e_im, oskel, g_densim, g_fxim, .true.)

C      call asym_fock2e (g_fxim)  !xxx 
C      call ga_scale (g_fxim, -1d0)  !FIXED: Proper canorg transform now, no need for this incorrect fudge

      call ga_antisymmetrize (g_fxim)
      call ga_scale (g_fxim, -1d0)

      
C
C     Check symmetries (disable for speed).
C
      if (params%checklvl .ge. 2) then
         if (.not. mat_is_symmetric (g_densim, "A", params%tol_zero))
     $        call errquit (pname//"Im[P] not antisymm", 0, 0)

         if (.not. mat_is_symmetric (g_fxim, "A", params%tol_zero))
     $        call errquit (pname//"Im[Vx] not antisymm", 0, 0)
      endif


C
C     Calculate exchange energy.
C
      Excim(1)= 0.5d0 * ga_ddot (g_densim, g_fxim)

      if (params%prof) call prof_end (elapsed, "Fock CS imag exch")
      
      end subroutine


c $Id$
