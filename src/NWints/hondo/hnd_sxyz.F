c
c $Id: hnd_sxyz.F,v 1.2 1999-04-10 00:15:17 windus Exp $
c
c Taken from the property/prp.F file
c
      subroutine hnd_sxyz
      implicit none
c
c     ----- gauss-hermite quadrature using minimum point formula -----
c
#include "hnd_whermt.fh"
c
      common/hnd_xyzder/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,
     &                  ni,nj,cx,cy,cz
      double precision xint, yint, zint, t, x0, y0, z0, xi, yi, zi
      double precision xj, yj, zj, zero, dum, px, py, pz
      double precision cx, cy, cz
      double precision ptx, pty, ptz, ax, ay, az, bx, by, bz
      integer ni, nj, min, max, npts, imin, imax, i
      dimension min(7),max(7)
      data min  /1,2,4, 7,11,16,22/
      data max  /1,3,6,10,15,21,28/
      data zero /0.0d+00/
c
      xint=zero
      yint=zero
      zint=zero
      npts=(ni+nj-2)/2+1
      imin=min(npts)
      imax=max(npts)
      do 16 i=imin,imax
      dum=w(i)
      px=dum
      py=dum
      pz=dum
      dum=h(i)*t
      ptx=dum+x0
      pty=dum+y0
      ptz=dum+z0
      ax=ptx-xi
      ay=pty-yi
      az=ptz-zi
      bx=ptx-xj
      by=pty-yj
      bz=ptz-zj
      go to (7,6,5,4,3,2,1),ni
    1 px=px*ax
      py=py*ay
      pz=pz*az
    2 px=px*ax
      py=py*ay
      pz=pz*az
    3 px=px*ax
      py=py*ay
      pz=pz*az
    4 px=px*ax
      py=py*ay
      pz=pz*az
    5 px=px*ax
      py=py*ay
      pz=pz*az
    6 px=px*ax
      py=py*ay
      pz=pz*az
    7 go to (15,14,13,12,11,10,9,8),nj
    8 px=px*bx
      py=py*by
      pz=pz*bz
    9 px=px*bx
      py=py*by
      pz=pz*bz
   10 px=px*bx
      py=py*by
      pz=pz*bz
   11 px=px*bx
      py=py*by
      pz=pz*bz
   12 px=px*bx
      py=py*by
      pz=pz*bz
   13 px=px*bx
      py=py*by
      pz=pz*bz
   14 px=px*bx
      py=py*by
      pz=pz*bz
   15 continue
      xint=xint+px
      yint=yint+py
      zint=zint+pz
   16 continue
      return
      end
