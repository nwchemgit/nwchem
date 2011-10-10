      subroutine argos_pme_putgrd(gridc,grid)
c
      implicit none
c
#include "argos_pme_common.fh"
#include "global.fh"
c
      real*8 gridc(2,mgx,mgy,*),grid(mgrx,mgry)
c
      integer i,j,k,m,n,ilsl,ihsl,jlsl,jhsl
c
      do 1 k=1,ngzloc
      do 2 i=1,ngx
      do 3 j=1,ngy
      grid(i,j)=gridc(1,i,j,k)
    3 continue
    2 continue
      do 4 m=1,morder
      do 5 i=1,ngx
      grid(i,ngy+m)=grid(i,m)
    5 continue
      do 6 j=1,ngy
      grid(ngx+m,j)=grid(m,j)
    6 continue
      do 7 n=1,morder
      grid(ngx+m,ngy+n)=grid(m,n)
    7 continue
    4 continue
      call ga_distribution(lslab(k),me,ilsl,ihsl,jlsl,jhsl)
      call ga_put(lslab(k),ilsl,ilsl+ngrx-1,jlsl,jlsl+ngry-1,
     + grid(1,1),mgrx)
    1 continue
c
      return
      end
c $Id$
