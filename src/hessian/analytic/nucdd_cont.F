C
C $Id: nucdd_cont.F,v 1.1 1999-05-06 17:29:35 windus Exp $
C
      subroutine nucdd_cont(rtdb)
c
c  Add the nuclear contributions the hessian matrix.
c  Start with RHF and work our way up.
c
      implicit none
c
#include "bas.fh"
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
      end
