      subroutine ana_xaver(x,xa,nx,mx,nave)
c
c $Id: ana_xaver.F,v 1.1 1999-12-10 23:54:47 d3j191 Exp $
c
      implicit none
c
#include "msgids.fh"
#include "global.fh"
c
      real*8 zero
      parameter(zero=0.0d0)
c
      integer nx,mx,nave
      real*8 x(mx,3),xa(mx,3)
c
      integer i,j
c
      if(nave.eq.0) then
      do 1 j=1,3
      do 2 i=1,nx
      xa(i,j)=x(i,j)
    2 continue
    1 continue
      else
      do 3 j=1,3
      do 4 i=1,nx
      xa(i,j)=xa(i,j)+x(i,j)
    4 continue
    3 continue
      endif
c
      nave=nave+1
c
      return
      end
      subroutine ana_xmean(xa,nx,mx,nave)
c
      implicit none
c
#include "msgids.fh"
#include "global.fh"
c
      integer nx,mx,nave
      real*8 xa(mx,3)
      integer i,j
      real*8 fact
c
      fact=1.0d0/dble(nave)
c
      do 1 j=1,3
      do 2 i=1,nx
      xa(i,j)=xa(i,j)*fact
    2 continue
    1 continue
c
      return
      end
