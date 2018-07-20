C
C     zfock_cs_exch.F
C
C     Compute both real and imaginary parts of exact exchange for a
C     closed shell complex Fock matrix.
C
C     We must call fock_2e() twice--once for the real (symmetric) part
C     of the exchange and once for the imaginary (antisymmetric) part.
C     If doing range separated exchange, we can't combine the real part
C     of exchange with Coulomb fock_2e() call, so that means in total to
C     build a CAM Fock matrix we need three separate fock_2e() calls.
C
C     Note that we only use this routine if doing CAM exchange or if
C     doing CD fitting for the Coulomb part.  Otherwise, for speed use
C     zfock_cs_coul_exchre() and zfock_cs_exchim().
C
C
      subroutine zfock_cs_exch (params, g_densre, g_densim, lcam, 
     $     Exc, g_fxre, g_fxim)
      implicit none

#include "bas.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "util.fh"
#include "cdft.fh"
#include "matutils.fh"
#include "case.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_densre, g_densim
      logical, intent(in)           :: lcam                !false=normal exch, true=CAM exch

      
C     == Outputs ==
      double precision, intent(out) :: Exc
      integer, intent(in)           :: g_fxre, g_fxim  !real and imag part of Vx

      
C     == Parameters ==
      character(*), parameter     :: pname = "zfock_cs_exch: "


C     == Variables ==
      double precision jfac, kfac
C      double precision tol2e
      double precision elapsed
      character*100 desc

      
      call rt_tddft_cs_confirm (params,'zfock_cs_exch.F')

      
      if (params%prof) call prof_start (elapsed)

      
C      tol2e=10.d0**(-itol2e)
      Exc = 0.d0

      call ga_zero (g_fxre)
      call ga_zero (g_fxim)
      
      jfac = 0.0d0
      kfac = -0.5d0*xfac(1)

C
C     Enable attenuation if doing CAM exchange.
C
      if (lcam) call case_setflags (.true.)
      
C
C     Compute real and imag parts of exchange.  Note ".false." means
C     *not* antisymmetric and ".true." means symmetric.  We explicitly
C     antisymmetrize the imag part afterwards.
C
      call fock_2e (geom, ao_bas_han, 1, jfac, kfac,
     &     params%tol2e_re, oskel, g_densre, g_fxre, .false.)

      call fock_2e (geom, ao_bas_han, 1, jfac, kfac,
     &     params%tol2e_im, oskel, g_densim, g_fxim, .true.)

C      call asym_fock2e (g_fxim)
C      call ga_scale (g_fxim, -1d0) !XXX NOTE NEGATIVE SIGN FUDGE
      
      call ga_antisymmetrize (g_fxim)
      call ga_scale (g_fxim, -1d0)



C
C     Calculate HF exchange energy.
C
      Exc = 0.5d0 * ga_ddot (g_densre, g_fxre) +
     $     0.5d0 * ga_ddot (g_densim, g_fxim)


C
C     Disable attenuation if we did CAM.
C      
      if (lcam) call case_setflags (.false.)


C
C     Profiling.
C
      if (lcam) then 
         desc = "Fock CS CAM exch"
      else
         desc = "Fock CS standard exch"
      endif
      
      if (params%prof) call prof_end (elapsed, trim(desc))
      end subroutine


c $Id$
