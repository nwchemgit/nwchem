C
C     zfock_cs_core.F
C
C     Calculates 1 electron part (Hcore) of a closed shell complex Fock
C     matrix.  Note that input is the real part of the density matrix,
C     and output is real part of the Fock matrix.
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

      call rt_tddft_cs_confirm (params)

      if (params%prof) call prof_start (elapsed)

      call ga_zero (g_fcore)

      call int_1e_ga(ao_bas_han, ao_bas_han, g_fcore,'kinetic', oskel)
      call int_1e_ga(ao_bas_han, ao_bas_han, g_fcore,'potential', oskel)
C      call int_1e_ga(ao_bas_han, ao_bas_han, g_fcore,'potential0',oskel)  #0 means dont include bq part


C
C     Testing Bq
C
c$$$      if (.not. ga_duplicate (g_fcore, g_tmp, "tmp"))
c$$$     $     call errquit (pname//"dup failed", 0, GA_ERR)
c$$$      call ga_zero (g_tmp)
c$$$      call int_1e_ga(ao_bas_han, ao_bas_han, g_tmp,'bq_pot', oskel)
c$$$      if (ga_nodeid().eq.0) write (luout, *) "junk = ",junk
c$$$      call ga_print (g_tmp)
c$$$      call ga_dadd (1d0, g_tmp, 1d0, g_fcore, g_fcore)
c$$$      if (.not. ga_destroy (g_tmp))
c$$$     $     call errquit (pname//"destroy failed", 0, GA_ERR)


      if (params%prof) call prof_end (elapsed, "Fock CS core")

      end subroutine
