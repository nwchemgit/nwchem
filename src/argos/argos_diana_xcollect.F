      subroutine argos_diana_xcollect(x,ifrom,nx,mx)
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      integer ifrom,nx,mx
      real*8 x(mx,3)
c
      integer i,j
c
      if(np.gt.1) then
c
      do 1 j=1,3
      do 2 i=1,ifrom-1
      x(i,j)=zero
    2 continue
    1 continue
c
      do 3 j=1,3
      do 4 i=ifrom+nx+1,mx
      x(i,j)=zero
    4 continue
    3 continue
c
      call ga_dgop(mag_d08,x,3*mx,'+')
c
      endif
c
      return
      end
c $Id$
