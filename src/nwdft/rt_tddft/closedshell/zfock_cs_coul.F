C
C     zfock_cs_coul.F
C
C     Calculates 2 electron part of a closed shell complex Fock matrix.
C     This routine is used when you *only* need to compute the Coulomb
C     term (i.e., pure DFT, with no Hartree-Fock exchange).  Note that
C     input is the real part of the density matrix, and output is real
C     part of the Fock matrix.
C
      subroutine zfock_cs_coul (params, g_densre, Ecoul, g_fc)
      implicit none

#include "bas.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "util.fh"
#include "cdft.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_densre

      
C     == Outputs ==
      double precision, intent(out) :: Ecoul
      integer, intent(in)           :: g_fc

      
C     == Parameters ==
      character(*), parameter     :: pname = "zfock_cs_coul:"
      

C     == Variables ==
      double precision jfac, kfac
C      double precision tol2e
      double precision elapsed


      call rt_tddft_cs_confirm (params,'zfock_cs_coul.F')

      if (params%prof) call prof_start (elapsed)

      
C     == Initializations ==
C      tol2e=10.d0**(-itol2e)
      Ecoul = 0.d0

      call ga_zero (g_fc)

      jfac = 1d0
      kfac = 0d0
      

C     == Calculate Coulomb part of F and energy ==
      call fock_2e (geom, ao_bas_han, 1, jfac, kfac,
     $     params%tol2e_re, oskel, g_densre, g_fc, .false.)

      Ecoul = 0.5d0 * ga_ddot (g_densre, g_fc)


      if (params%prof) call prof_end (elapsed,"Fock CS standard coul")
      end subroutine
c $Id$
