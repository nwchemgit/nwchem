C
C     rt_tddft_os_tdfock.F
C
C     Wrapper routine for building time-dependent open shell Fock matrix
C     including building Fock matrix from density matrix in AO basis,
C     and calculating and adding dipole interation with external uniform
C     E-field.
C
      logical function rt_tddft_os_tdfock (params, tt, g_zdens_ao,
     $     energies, g_zfock_ao)

      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "matutils.fh"
#include "rt_tddft.fh"


C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      double precision, intent(in)  :: tt            !current time
      integer, intent(in)           :: g_zdens_ao(2) !dens mat, ns_ao x ns_ao; alpha,beta

      
C     == Outputs ==
      type(rt_energies_t), intent(out) :: energies
      integer, intent(in)              :: g_zfock_ao(2) !fock mat, ns_ao x ns_ao


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_os_tdfock: "


C     == Variables ==
      type(rt_vector_t) field
      integer g_zscr(2)            !scratch, ns_ao x ns_ao; alpha, beta
      integer n


      call rt_tddft_os_confirm (params,'rt_tddft_os_tdfock.F')


      n = params%ns_ao            ! alias


      if (.not.ga_create(mt_dcpl,n ,n ,"zscr1", 0, 0, g_zscr(1)))
     $     call errquit ("failed to create zscr1", 0, GA_ERR)
      if (.not.ga_create(mt_dcpl,n ,n ,"zscr2", 0, 0, g_zscr(2)))
     $     call errquit ("failed to create zscr2", 0, GA_ERR)

      
C        
C     Build new complex fock mat from complex dens mat; this also
C     calculates energies.  Note, the input g_zdens_ao is in AO basis,
C     and the output g_zfock_ao is also in AO basis.
C     
      call zfock_os_build (params, g_zdens_ao, energies, g_zfock_ao)
      call ga_sync ()           !XXX needed?

C
C     Compute dipole interaction (updates field values inside), and add
C     to Fock matrix and add D.E(t).  We compute the contributions to
C     the fields acting on alpha and beta separately.
C

      call ga_zero (g_zscr(1))
      call ga_zero (g_zscr(2))
      call rt_tddft_calc_excite (params, tt, g_zscr)
      call ga_add (z1, g_zfock_ao(1), zn1, g_zscr(1), g_zfock_ao(1))
      call ga_add (z1, g_zfock_ao(2), zn1, g_zscr(2), g_zfock_ao(2))

      if (.not.ga_destroy(g_zscr(1)))
     $     call errquit ("failed to destroy zscr1", 0, GA_ERR)
      if (.not.ga_destroy(g_zscr(2)))
     $     call errquit ("failed to destroy zscr2", 0, GA_ERR)

      rt_tddft_os_tdfock = .true.
      end function
c $Id$
