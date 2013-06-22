C
C     Re-initialize fock_2e, specifically how to cache integrals.  This
C     was cleaned up already in dft_scf() so we need to do it again so
C     we don't mess with the SCF routine.
C
      subroutine rt_tddft_init_fock2e (params)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
C#include "cdft.fh"
#include "rt_tddft.fh"




C     == Inputs ==
      type(rt_params_t), intent(inout) :: params
      

C     == Parameters ==
      character(len=*), parameter :: pname = "rt_tddft_init_fock2e: "


C     == Variables ==
      integer itol2e


C
C     Note that we have tighter 2e integral tolerances for the imaginary
C     matricies (e.g,. when we have exact exchange), which are typically
C     orders of magnitude smaller than the real matricies.  This
C     magnitude issue seems to cause some problems (especially running
C     direct) so we tighten it up.
C
      if (.not.rtdb_get(params%rtdb,'dft:itol2e',mt_int,1,itol2e))
     $     call errquit(pname//"failed to read itol2e",0, RTDB_ERR)

      params%tol2e_re = 10.0d0**(-itol2e)
      params%tol2e_im = 0.01d0*params%tol2e_re

      if (params%tol2e_re .gt. 1.1d-8)
     $     call rt_tddft_print_warning
     $     ("Dangerously loose tol2e, be careful")

      call scf_get_fock_param (params%rtdb, params%tol2e_re)

      end subroutine rt_tddft_init_fock2e


c $Id$
