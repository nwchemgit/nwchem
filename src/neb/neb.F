*
* $Id: neb.F,v 1.1 2004-02-11 02:37:09 marat Exp $
*
      logical function neb(rtdb)
      implicit none
      integer rtdb
      logical nebnew

#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "neb_data.fh"

      neb=.true.
      call neb_init(rtdb)
      call runall_bead_list(neb_list)
      
      
      return
      end
c

