      logical function argos_prop_rdmri(lfnmri,ndec,mropt)
c
      implicit none
c
#include "argos_prop.fh"
#include "mafdecls.fh"
c
      logical argos_prop_rmri
      external argos_prop_rmri
c
      integer lfnmri,ndec,mropt
c
      argos_prop_rdmri=argos_prop_rmri(lfnmri,ndec,mropt,
     + dbl_mb(i_dfr),dbl_mb(i_dfrm))
c
      return
      end
c $Id$
