C
C     zfock_cs_core.F
C
C     Calculates 1 electron part (Hcore) of a closed shell complex Fock
C     matrix.  Note that input is the real part of the density matrix,
C     and output is real part of the Fock matrix.
C
C     Of course, we can need only compute this once, but I do it every
C     time step since it is cheap (XXX do once and store).
C
      subroutine zfock_cs_core (params, g_densre, g_fcore)
      implicit none

#include "rtdb.fh"
#include "bas.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "cdft.fh"
#include "global.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_densre

      
C     == Outputs ==
      integer, intent(in)           :: g_fcore 

      
C     == Parameters ==
      character(*), parameter  :: pname = "zfock_cs_core: "

      
C     == Variables ==
      double precision elapsed
      integer g_tmp
      double precision junk

      call rt_tddft_cs_confirm (params,'zfock_cs_core.F')

      if (params%prof) call prof_start (elapsed)

C     (compute on the fly)
c$$$      call ga_zero (g_fcore)
c$$$      call int_1e_ga(ao_bas_han, ao_bas_han, g_fcore,'kinetic', oskel)
c$$$      call int_1e_ga(ao_bas_han, ao_bas_han, g_fcore,'potential', oskel)
c$$$C      call int_1e_ga(ao_bas_han, ao_bas_han, g_fcore,'potential0',oskel)  #0 means dont include bq part

C     (use pre-computed/static)
      call ga_zero (g_fcore)
      call ga_copy (params%g_fock_ao_core, g_fcore)
      
      
      if (params%prof) call prof_end (elapsed, "Fock CS core")

      end subroutine
c $Id$
