C
C $Id: nucdd_cont.F,v 1.3 2000-05-30 14:47:30 windus Exp $
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
c
      integer rtdb  ! [input] RTDB handle
c
      call nucdd_calc(coords(1,1,geom),charge(1,geom),ncenter(geom),
     1                dbl_mb(k_exy))
c
      return
      end
