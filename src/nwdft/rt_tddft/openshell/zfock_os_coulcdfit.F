C
C     zfock_os_coulcdfit.F
C
C     Computes the Coulomb part of the complex open shell Fock Matrix
C     via charge density fitting.
C
      subroutine zfock_os_coulcdfit (params, g_densre, Ecoul, g_vc)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "cdft.fh"
#include "matutils.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_densre(2) !alpha, beta

      
C     == Ouputs ==
      double precision, intent(out) :: Ecoul
      integer, intent(in)           :: g_vc(2)

      
C     == Parameters ==
      character(*), parameter :: pname = "zfock_os_coulcdfit: "
      integer, parameter      :: iVcoul_opt = 1 !XXX HARDCODED


C     == Variables ==
      integer lcd_coef, icd_coef
      double precision elapsed
      integer g_densre_tot

      
C      call rt_tddft_os_confirm (params,'junk')

      
      if (params%prof) call prof_start (elapsed)

      
C
C     Allocations.
C
      if (.not.ma_Push_Get (mt_dbl,nbf_cd,"cd_coef",lcd_coef, icd_coef))
     $     call errquit (pname//"cannot allocate cd_coef", 0, MA_ERR)

      if (.not.ga_duplicate (g_densre(1), g_densre_tot, "Ptot"))
     $     call errquit (pname//"cannot allocate Ptot", 1, GA_ERR)


C
C     Compute total dens mat.
C
      call ga_add (1d0, g_densre(1), 1d0, g_densre(2), g_densre_tot)


C
C     Fit charge density with CD basis.
C
      call dft_fitcd (1, dbl_mb(icd_coef), dbl_mb(params%k_3ceri),
     $     Ecoul, g_densre_tot, params%ntotel, params%n_batch,
     $     params%n3c_int, int_mb(params%k_3cwhat), params%n3c_dbl,
cold     $     params%iwhat_max, params%n_semi_bufs, params%fd,
     $     params%iolgc, params%natoms, .false., 0d0, .false.)

      
C
C     Compute Vcoul using least squares fitting; store in alpha part.
C
      call dft_getvc(dbl_mb(icd_coef), dbl_mb(params%k_3ceri),
     $     Ecoul, g_vc(1), iVcoul_opt, params%n_batch,
     $     params%n3c_int, int_mb(params%k_3cwhat), params%n3c_dbl,
cold     $     params%iwhat_max, params%n_semi_bufs, params%fd,
     $     params%iolgc, .false., 1)


C
C     Copy to coulomb part of Fock mat to beta part.
C
      call ga_copy (g_vc(1), g_vc(2))

      
C
C     Clean up.
C      	
      if (.not. ma_chop_stack (lcd_coef))
     $     call errquit (pname//"failed to chop stack", 0, MA_ERR)

      if (.not.ga_destroy (g_densre_tot))
     $     call errquit (pname//"cannot destroy Ptot", 1, GA_ERR)

      
      if (params%prof) call prof_end(elapsed, "Fock OS CD fitting coul")

      end subroutine


c $Id$
