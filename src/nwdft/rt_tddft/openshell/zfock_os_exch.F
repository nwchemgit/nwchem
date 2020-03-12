C
C     zfock_os_exch.F
C
C     Compute both real and imaginary parts of exact exchange for a
C     open shell complex Fock matrix.
C
C     We must call fock_2e() twice--once for the real (symmetric) part
C     of the exchange and once for the imaginary (antisymmetric) part.
C     If doing range separated exchange, we can't combine the real part
C     of exchange with Coulomb fock_2e() call, so that means in total to
C     build a CAM Fock matrix we need three separate fock_2e() calls.
C
C     Note that for closed shell, the Falpha depends only on Palpha,
C     etc.
C
C
      subroutine zfock_os_exch (params, g_densre, g_densim, lcam, 
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
      type(rt_params_t), intent(in)     :: params
      integer, intent(in)               :: g_densre(2), g_densim(2)
      logical, intent(in)               :: lcam   !false=normal exchange, true=CAM exchange


      
C     == Outputs ==
      double precision, intent(out)     :: Exc
      integer, intent(in)               :: g_fxre(2), g_fxim(2)  !re and im Vx; alpha,beta

      
C     == Parameters ==
      character(*), parameter     :: pname = "zfock_os_exch: "


C     == Variables ==
      integer me
      double precision jfac(2), kfac(2)    !alpha, beta
C      double precision tol2e
      double precision elapsed
      character*100 desc


C      call rt_tddft_os_confirm (params,'junk')

      if (params%prof) call prof_start (elapsed)

      me = ga_nodeid()
      
C      tol2e=10.d0**(-itol2e)
      Exc = 0.d0

      call ga_zero (g_fxre(1))
      call ga_zero (g_fxre(2))
      call ga_zero (g_fxim(1))
      call ga_zero (g_fxim(2))
      
C
C     Note 1 in xfac(1) refers to xc factor *not* alpha,beta.  Also note
C     -1d0 prefactor *not* -0.5d0 like in closed shell.
C
      jfac(1) = 0.0d0
      jfac(2) = 0.0d0
      kfac(1) = -1d0*xfac(1)  
      kfac(2) = -1d0*xfac(1)

C
C     Enable attenuation if doing CAM exchange.
C
      if (lcam) call case_setflags (.true.)
      
C
C     Compute real and imag parts of exchange.  Note ".false." means
C     *not* antisymmetric and ".true." means symmetric.  We explicitly
C     antisymmetrize the imag part afterwards.
C
      call fock_2e (geom, ao_bas_han, 2, jfac, kfac,
     &     params%tol2e_re, oskel, g_densre, g_fxre, .false.)

      call fock_2e (geom, ao_bas_han, 2, jfac, kfac,
     &     params%tol2e_im, oskel, g_densim, g_fxim, .true.)

      
C
C     Antisymmetrize alpha, beta parts of Fock mat
C      
c$$$      call asym_fock2e (g_fxim(1))
c$$$C      call ga_scale (g_fxim(1), -1d0) !XXX NOTE NEGATIVE SIGN FUDGE
c$$$
c$$$      call asym_fock2e (g_fxim(2))
c$$$C      call ga_scale (g_fxim(2), -1d0) !XXX NOTE NEGATIVE SIGN FUDGE


      call ga_antisymmetrize (g_fxim(1))
      call ga_scale (g_fxim(1), -1d0)
         
      call ga_antisymmetrize (g_fxim(2))
      call ga_scale (g_fxim(2), -1d0)
      

C
C     Calculate HF exchange energy.
C
      Exc = 0.5d0 * (ga_ddot (g_densre(1), g_fxre(1)) +
     $     ga_ddot (g_densre(2), g_fxre(2)) +
     $     ga_ddot (g_densim(1), g_fxim(1)) +
     $     ga_ddot (g_densim(2), g_fxim(2)))   !note half


C
C     Disable attenuation if we did CAM.
C      
      if (lcam) call case_setflags (.false.)


C
C     Profiling.
C
      if (lcam) then 
         desc = "Fock OS CAM exch"
      else
         desc = "Fock OS standard exch"
      endif
      
      if (params%prof) call prof_end (elapsed, trim(desc))
      end subroutine

c $Id$
