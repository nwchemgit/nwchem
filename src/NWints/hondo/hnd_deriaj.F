C
C $Id$
c
      subroutine hnd_deriaj(dxdi,dydi,dzdi,dxdj,dydj,dzdj,
     1                      x,y,z,lit,ljt,ai,aj)
c
c This is a routine taken from HONDO.  It was originally
c hnd_deri, but was extended to also do the derivative wrt j.
c
      implicit none
c
      double precision dxdi,dydi,dzdi,dxdj,dydj,dzdj
      double precision x,y,z,ai,aj
      integer lit, ljt
      integer j, i
c
      dimension x(lit+1,ljt+1),y(lit+1,ljt+1),z(lit+1,ljt+1)
      dimension dxdi(lit,ljt),dydi(lit,ljt),dzdi(lit,ljt)
      dimension dxdj(lit,ljt),dydj(lit,ljt),dzdj(lit,ljt)
c
c     ----- dxdi ... -----
c
      do 110 j=1,ljt
      dxdi(1,j)=x(2,j)*(ai+ai)
      dydi(1,j)=y(2,j)*(ai+ai)
      dzdi(1,j)=z(2,j)*(ai+ai)
  110 continue
c
      if(lit.eq.1) go to 130
c
      do 120 i=2,lit
      do 120 j=1,ljt
      dxdi(i,j)=x(i+1,j)*(ai+ai)-x(i-1,j)*dble(i-1)
      dydi(i,j)=y(i+1,j)*(ai+ai)-y(i-1,j)*dble(i-1)
      dzdi(i,j)=z(i+1,j)*(ai+ai)-z(i-1,j)*dble(i-1)
  120 continue
c
  130 continue
c
c     ----- dxdj ... -----
c
      do 140 i=1,lit
      dxdj(i,1)=x(i,2)*(aj+aj) 
      dydj(i,1)=y(i,2)*(aj+aj) 
      dzdj(i,1)=z(i,2)*(aj+aj) 
  140 continue
c
      if (ljt.eq.1) return
c
      do 150 j=2,ljt
      do 150 i=1,lit
      dxdj(i,j)=x(i,j+1)*(aj+aj)-x(i,j-1)*dble(j-1)
      dydj(i,j)=y(i,j+1)*(aj+aj)-y(i,j-1)*dble(j-1)
      dzdj(i,j)=z(i,j+1)*(aj+aj)-z(i,j-1)*dble(j-1)
  150 continue
c
      return
      end
