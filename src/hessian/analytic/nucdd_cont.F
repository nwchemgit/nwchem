C
C $Id$
C
      subroutine nucdd_cont(rtdb)
c
c  Add the nuclear contributions the hessian matrix.
c  Start with RHF and work our way up.
c
      implicit none
c
#include "bas.fh"
#include "nwc_const.fh"
#include "geomP.fh"
#include "global.fh"
#include "hess_info.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
c
      integer rtdb  ! [input] RTDB handle
c
      call nucdd_calc(coords(1,1,geom),charge(1,geom),ncenter(geom),
     1                dbl_mb(k_exy))
c
c Print info
c
      if (util_print('nucdd_cont',print_debug) .and.
     $    (ga_nodeid().eq.0)) then
        write(LuOut,*)
        write(LuOut,*) 'Hessian after nuclear contribution'
        call hess_hssout(dbl_mb(k_exy),n3xyz,n3xyz,n3xyz)
      endif
c
      return
      end
