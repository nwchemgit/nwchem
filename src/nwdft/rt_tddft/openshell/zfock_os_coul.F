C
C     zfock_os_coul.F
C
C     Calculates the Coulomb term of the open shell Fock matrix.
C     exchange).  Note that input is the real part of the density matrix
C     (with alpha, beta spin part), and output is real part of the Fock
C     matrix.
C
      subroutine zfock_os_coul (params, g_densre, Ecoul, g_fc)
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
      type(rt_params_t), intent(in)     :: params
      integer, intent(in)               :: g_densre(2)       !alpha, beta

      
C     == Outputs ==
      double precision, intent(out)     :: Ecoul
      integer, intent(in)               :: g_fc(2)           !alpha, beta

      
C     == Parameters ==
      character(*), parameter           :: pname = "zfock_os_coul:"
      

C     == Variables ==
      double precision jfac, kfac
C      double precision tol2e
      double precision elapsed
      integer g_densre_tot      !alpha+beta


C      call rt_tddft_os_confirm (params,'junk')


      if (params%prof) call prof_start (elapsed)

      
      if (.not.ga_duplicate (g_densre(1), g_densre_tot, "Ptot"))
     $     call errquit (pname//"cannot allocate Ptot", 1, GA_ERR)


C
C     Compute total dens mat.
C
      call ga_add (1d0, g_densre(1), 1d0, g_densre(2), g_densre_tot)


C      tol2e=10.d0**(-itol2e)
      Ecoul = 0.d0

      call ga_zero (g_fc(1))
      call ga_zero (g_fc(2))

      jfac = 1d0
      kfac = 0d0
      

C     Calculate Coulomb part of F and energy--both alpha, beta part of
C     Fock mat have same value (it depends on full density matrix).  The
C     result will be stored in g_fc(1) (the alpha part) and then copied
C     to the beta part.
      
      call fock_2e (geom, ao_bas_han, 1, jfac, kfac,
     $     params%tol2e_re, oskel, g_densre_tot, g_fc(1), .false.)

      
C
C     Copy to beta part.
C
      call ga_copy (g_fc(1), g_fc(2))


      Ecoul = 0.5d0 * ga_ddot (g_densre_tot, g_fc(1))   !note half

C      Ecoul = 0.5d0 * (ga_ddot (g_densre(1), g_fc(1)) +
C     $     ga_ddot (g_densre(2), g_fc(2)))
      

      if (.not.ga_destroy (g_densre_tot))
     $     call errquit (pname//"cannot destroy Ptot", 1, GA_ERR)


      if (params%prof) call prof_end (elapsed,"Fock OS standard coul")
      end subroutine
c $Id$
