      subroutine argos_pme_getgrd(node,nsl,lsldo,maxsl,grid)
c
      implicit none
c
#include "argos_pme_common.fh"
#include "global.fh"
c
      integer node,nsl
      logical lsldo(np,mgz)
      integer maxsl(4,np,mgz)
      real*8 grid(mgrx,mgry,*)
c
      integer k,ilsl,ihsl,jlsl,jhsl
      integer minx,maxx,miny,maxy
c
      do 1 k=1,nsl
      if(lsldo(node+1,k)) then
      minx=maxsl(1,node+1,k)
      maxx=maxsl(2,node+1,k)
      miny=maxsl(3,node+1,k)
      maxy=maxsl(4,node+1,k)
c      minx=1
c      maxx=ngx
c      miny=1
c      maxy=ngy
      call ga_distribution(lslab(k),node,ilsl,ihsl,jlsl,jhsl)
      call ga_get(lslab(k),ilsl+minx-1,ilsl+maxx-1,jlsl+miny-1,
     + jlsl+maxy-1,grid(minx,miny,k),mgrx)
      endif
    1 continue
c
      return
      end
c $Id$
