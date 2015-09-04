c
c $Id: hnd_sxyz.F 19696 2010-10-29 16:53:42Z d3y133 $
c
c Taken from the property/prp.F file
c
      subroutine dim_sxyz2
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
      integer ni, nj, npts, imin, imax, i, ii, jj
      data zero /0.0d+00/
c
      xint=zero
      yint=zero
      zint=zero
      npts=(ni+nj-2)/2+1
      imin=hermin(npts)
      imax=hermax(npts)
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
         do ii = 1, ni-1
            px=px*ax
            py=py*ay
            pz=pz*az
         enddo
         do jj = 1, nj-1
            px=px*bx
            py=py*by
            pz=pz*bz
         enddo
         xint=xint+px
         yint=yint+py
         zint=zint+pz
   16 continue
      return
      end
