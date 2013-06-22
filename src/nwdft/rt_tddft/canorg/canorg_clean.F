C
C     canorg_clean.F
C
C     Clean up memory for canonical orthogonalization routines.
C
      subroutine canorg_clean (params)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "geom.fh"
C#include "cdft.fh"
#include "rtdb.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in) :: params


C     == Parameters ==
      character, parameter  :: pname = "canorg_clean: "

      if (.not. ga_destroy (params%g_zX))
     $     call errquit (pname//"failed to destroy g_zX", 0, 0)

      if (.not. ga_destroy (params%g_zY))
     $     call errquit (pname//"failed to destroy g_zY", 0, 0)

      if (.not. ga_destroy (params%g_s))    !allocd in rt_tddft_init_overlap_canorg.F
     $     call errquit (pname//"failed to destroy g_s", 0, 0)
      
c$$$      if (.not. ga_destroy (params%g_s12))  !allocd in rt_tddft_init_overlap_canorg.F
c$$$     $     call errquit (pname//"failed to destroy g_s12", 0, 0)
c$$$
c$$$      if (.not. ga_destroy (params%g_sn12)) !allocd in rt_tddft_init_overlap_canorg.F
c$$$     $     call errquit (pname//"failed to destroy g_sn12", 0, 0)

      end subroutine
c $Id$
