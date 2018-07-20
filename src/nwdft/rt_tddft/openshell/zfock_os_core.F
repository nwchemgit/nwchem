C
C     zfock_os_core.F
C
C     Calculates 1 electron part (Hcore) of a open shell complex Fock
C     matrix.  Note that input is the real part of the density matrix,
C     and output is real part of the Fock matrix.
C
C     Of course, we can need only compute this once, but I do it every
C     time step since it is cheap (XXX do once and store).
C
      subroutine zfock_os_core (params, g_densre, g_fcore)
      implicit none

#include "rtdb.fh"
#include "bas.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "cdft.fh"
#include "global.fh"
#include "rt_tddft.fh"


CXXX  [KAL]: gives slightly different Ecore than from SCF???


C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_densre(2)

      
C     == Outputs ==
      integer, intent(in)           :: g_fcore(2) ! complex fock mat; (nbf x nbf)

      
C     == Parameters ==
      character(*), parameter :: pname = "zfock_os_core: "

      
C     == Variables ==
      double precision elapsed


C     xxx for SO
C      call rt_tddft_os_confirm (params,'junk')

      
      if (params%prof) call prof_start (elapsed)


C
C     Both alpha and beta spins have same core part: compute alpha and
C     copy to beta later.
C


C     (compute on the fly)
c$$$      call ga_zero (g_fcore(1))
c$$$      call ga_zero (g_fcore(2))
c$$$
c$$$      call int_1e_ga(ao_bas_han,ao_bas_han,g_fcore(1),'kinetic',oskel)
c$$$      call int_1e_ga(ao_bas_han,ao_bas_han,g_fcore(1),'potential',oskel)


C     (use pre-computed/static)
      call ga_zero (g_fcore(1))
      call ga_copy (params%g_fock_ao_core, g_fcore(1))

      call ga_copy (g_fcore(1), g_fcore(2))

      if (params%prof) call prof_end (elapsed, "Fock OS core")

      end subroutine
c $Id$
