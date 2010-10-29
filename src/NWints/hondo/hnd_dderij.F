c $Id$
c
c taken from hondo.
c
      subroutine hnd_dderij2(ddxdii,ddydii,ddzdii,ddxdij,ddydij,ddzdij,
     1                  ddxdjj,ddydjj,ddzdjj,
     2                  dxdi,dydi,dzdi,dxdj,dydj,dzdj,
     3                  x,y,z,lit,ljt,ai,aj,nder)
      implicit none
      double precision x, y, z, dxdi, dydi, dzdi, dxdj, dydj, dzdj
      double precision ddxdii, ddydii, ddzdii 
      double precision ddxdjj, ddydjj, ddzdjj 
      double precision ddxdij, ddydij, ddzdij
      double precision ai, aj
      integer lit, ljt, j, i, nder
      dimension      x(lit+nder,ljt+nder)
      dimension      y(lit+nder,ljt+nder)
      dimension      z(lit+nder,ljt+nder)
      dimension   dxdi(lit,ljt),  dydi(lit,ljt),  dzdi(lit,ljt)
      dimension   dxdj(lit,ljt),  dydj(lit,ljt),  dzdj(lit,ljt)
      dimension ddxdii(lit,ljt),ddydii(lit,ljt),ddzdii(lit,ljt)
      dimension ddxdjj(lit,ljt),ddydjj(lit,ljt),ddzdjj(lit,ljt)
      dimension ddxdij(lit,ljt),ddydij(lit,ljt),ddzdij(lit,ljt)
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
      if(ljt.eq.1) go to 200
c
      do 150 j=2,ljt
      do 150 i=1,lit
      dxdj(i,j)=x(i,j+1)*(aj+aj)-x(i,j-1)*dble(j-1)
      dydj(i,j)=y(i,j+1)*(aj+aj)-y(i,j-1)*dble(j-1)
      dzdj(i,j)=z(i,j+1)*(aj+aj)-z(i,j-1)*dble(j-1)
  150 continue
c
  200 continue
      if (nder.lt.2) return
c
c     ----- ddxdii ... -----
c
      do 210 j=1,ljt
      ddxdii(1,j)=(x(3,j)*(ai+ai)-x(1,j))*(ai+ai)
      ddydii(1,j)=(y(3,j)*(ai+ai)-y(1,j))*(ai+ai)
      ddzdii(1,j)=(z(3,j)*(ai+ai)-z(1,j))*(ai+ai)
  210 continue
c
      if(lit.eq.1) go to 300
c
      do 220 j=1,ljt
      ddxdii(2,j)=(x(4,j)*(ai+ai)-x(2,j)*dble(2+2-1))*(ai+ai)
      ddydii(2,j)=(y(4,j)*(ai+ai)-y(2,j)*dble(2+2-1))*(ai+ai)
      ddzdii(2,j)=(z(4,j)*(ai+ai)-z(2,j)*dble(2+2-1))*(ai+ai)
  220 continue
c
      if(lit.eq.2) go to 300
c
      do 230 i=3,lit
      do 230 j=1,ljt
      ddxdii(i,j)=(x(i+2,j)*(ai+ai)-x(i,j)*dble(i+i-1))*(ai+ai)
     1            +x(i-2,j)*dble((i-1)*(i-2))
      ddydii(i,j)=(y(i+2,j)*(ai+ai)-y(i,j)*dble(i+i-1))*(ai+ai)
     1            +y(i-2,j)*dble((i-1)*(i-2))
      ddzdii(i,j)=(z(i+2,j)*(ai+ai)-z(i,j)*dble(i+i-1))*(ai+ai)
     1            +z(i-2,j)*dble((i-1)*(i-2))
  230 continue
c
  300 continue
c
c     ----- ddxdij ... -----
c
      ddxdij(1,1)= x(2,2)*((aj+aj)*(ai+ai))
      ddydij(1,1)= y(2,2)*((aj+aj)*(ai+ai))
      ddzdij(1,1)= z(2,2)*((aj+aj)*(ai+ai))
c
      if(ljt.eq.1) go to 330
      do 320 j=2,ljt
      ddxdij(1,j)=(x(2,j+1)*(aj+aj)-x(2,j-1)*dble(j-1))*(ai+ai)
      ddydij(1,j)=(y(2,j+1)*(aj+aj)-y(2,j-1)*dble(j-1))*(ai+ai)
      ddzdij(1,j)=(z(2,j+1)*(aj+aj)-z(2,j-1)*dble(j-1))*(ai+ai)
  320 continue
  330 if(lit.eq.1) go to 400
      do 340 i=2,lit
      ddxdij(i,1)=(x(i+1,2)*(ai+ai)-x(i-1,2)*dble(i-1))*(aj+aj)
      ddydij(i,1)=(y(i+1,2)*(ai+ai)-y(i-1,2)*dble(i-1))*(aj+aj)
      ddzdij(i,1)=(z(i+1,2)*(ai+ai)-z(i-1,2)*dble(i-1))*(aj+aj)
  340 continue
      if(ljt.eq.1) go to 400
      do 350 i=2,lit
      do 350 j=2,ljt
      ddxdij(i,j)=
     1      (x(i+1,j+1)*(aj+aj)-x(i+1,j-1)*dble(j-1))*(ai+ai)
     2     -(x(i-1,j+1)*(aj+aj)-x(i-1,j-1)*dble(j-1))*dble(i-1)
      ddydij(i,j)=
     1      (y(i+1,j+1)*(aj+aj)-y(i+1,j-1)*dble(j-1))*(ai+ai)
     2     -(y(i-1,j+1)*(aj+aj)-y(i-1,j-1)*dble(j-1))*dble(i-1)
      ddzdij(i,j)=
     1      (z(i+1,j+1)*(aj+aj)-z(i+1,j-1)*dble(j-1))*(ai+ai)
     2     -(z(i-1,j+1)*(aj+aj)-z(i-1,j-1)*dble(j-1))*dble(i-1)
  350 continue
c
  400 continue
c
c     ----- ddxdjj ... -----
c
      do 410 i=1,lit
      ddxdjj(i,1)=(x(i,3)*(aj+aj)-x(i,1))*(aj+aj)
      ddydjj(i,1)=(y(i,3)*(aj+aj)-y(i,1))*(aj+aj)
      ddzdjj(i,1)=(z(i,3)*(aj+aj)-z(i,1))*(aj+aj)
  410 continue
c
      if(ljt.eq.1) return
c
      do 420 i=1,lit
      ddxdjj(i,2)=(x(i,4)*(aj+aj)-x(i,2)*dble(2+2-1))*(aj+aj)
      ddydjj(i,2)=(y(i,4)*(aj+aj)-y(i,2)*dble(2+2-1))*(aj+aj)
      ddzdjj(i,2)=(z(i,4)*(aj+aj)-z(i,2)*dble(2+2-1))*(aj+aj)
  420 continue
c
      if(ljt.eq.2) return
c
      do 430 i=1,lit
      do 430 j=3,ljt
      ddxdjj(i,j)=(x(i,j+2)*(aj+aj)-x(i,j)*dble(j+j-1))*(aj+aj)
     1            +x(i,j-2)*dble((j-1)*(j-2))
      ddydjj(i,j)=(y(i,j+2)*(aj+aj)-y(i,j)*dble(j+j-1))*(aj+aj)
     1            +y(i,j-2)*dble((j-1)*(j-2))
      ddzdjj(i,j)=(z(i,j+2)*(aj+aj)-z(i,j)*dble(j+j-1))*(aj+aj)
     1            +z(i,j-2)*dble((j-1)*(j-2))
  430 continue
      return
      end
