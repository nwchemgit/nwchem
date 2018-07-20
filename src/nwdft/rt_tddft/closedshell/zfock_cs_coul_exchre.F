C     
C     zfock_cs_coul_exchre.F
C
C     Computes the Coulomb (which is pure real) and real parts of the
C     exact exchange for a closed shell complex Fock matrix.
C
      subroutine zfock_cs_coul_exchre (params, g_densre,
     $     Ecoul, Excre, g_fcxre)
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
      integer, intent(in)           :: g_densre

      
C     == Outputs ==
      double precision, intent(out) :: Ecoul, Excre(2)
      integer, intent(in)           :: g_fcxre         !contains Coul and re part of exch

      
C     == Parameters ==
      character(*), parameter     :: pname = "zfock_cs_coul_exchre: "
      

C     == Variables ==
      double precision jfac(2), kfac(2)
C      double precision tol2e
      integer g_dens_hans(2)
      integer g_tmps(2)         ! work arrays (nbf x nbf)
      double precision elapsed

      double precision JUNK
      

      call rt_tddft_cs_confirm (params,'zfock_cs_coul_exchre.F')

      if (params%prof) call prof_start (elapsed)

      
C     == Initializations ==

      call ga_zero (g_fcxre)

C      tol2e=10.d0**(-itol2e)
      Ecoul = 0d0
      Excre(1) = 0d0

      
      if (.not.ga_duplicate (g_fcxre, g_tmps(1), "work mat 1"))
     $     call errquit (pname//"cannot allocate g_tmps1", 1, GA_ERR)

      if (.not.ga_duplicate (g_fcxre, g_tmps(2), "work mat 2"))
     $     call errquit (pname//"cannot allocate g_tmps2", 1, GA_ERR)

      call ga_zero (g_tmps(1))
      call ga_zero (g_tmps(2))


C
C     The first facs refer to the coulomb part, and the second facs
C     refer to the real part of the exact exchange.

      g_dens_hans(1) = g_densre
      g_dens_hans(2) = g_densre

      jfac(1) = 1.0d0
      kfac(1) = 0.0d0

      jfac(2) = 0.0d0
      kfac(2) = -0.5d0*xfac(1)

      
C
C     Compute potentials in AO basis.
C     
      call fock_2e (geom, ao_bas_han, 2, jfac, kfac,
     $     params%tol2e_re, oskel, g_dens_hans, g_tmps, .false.)


C      
C     Calculate Coulomb energy and accumulate in Fock matrix.
C
      Ecoul = 0.5d0 * ga_ddot (g_densre, g_tmps(1))
      call ga_dadd (1d0, g_tmps(1), 1d0, g_fcxre, g_fcxre)


C
C     Calculate HF exchange energy and accumulate in Fock matrix.
C
      Excre(1) = 0.5d0 * ga_ddot (g_densre, g_tmps(2))
      call ga_dadd (1d0, g_tmps(2), 1d0, g_fcxre, g_fcxre)



C     == Clean up ==
      if (.not.ga_destroy (g_tmps(1)))
     $     call errquit (pname//"cannot destroy g_tmps(1)", 1, GA_ERR)
      if (.not.ga_destroy (g_tmps(2)))
     $     call errquit (pname//"cannot destroy g_tmps(2)", 1, GA_ERR)

      if (params%prof)
     $     call prof_end (elapsed,"Fock CS standard coul + real exch")

      end subroutine

c $Id$
