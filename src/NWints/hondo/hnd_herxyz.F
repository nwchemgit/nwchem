c
c $Id$
c
      subroutine hnd_hermxyz
      implicit none
c
c     ----- gauss-hermite quadrature using minimum point formula -----
c
#include "hnd_whermt.fh"
c
      common/hnd_xyzher/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj,
     &                  xd,yd,zd,ni,nj,xc,yc,zc,nc,nd
      double precision xint, yint, zint, t, x0, y0, z0, xi, yi, zi
      double precision xj, yj, zj, zero, dum, px, py, pz, xd, yd, zd
      double precision cx, cy, cz, dx, dy, dz, xc, yc, zc
      double precision ptx, pty, ptz, ax, ay, az, bx, by, bz
      integer ni, nj, npts, imin, imax, i, ii, jj, nc, nd, kk, ll
      data zero /0.0d+00/
c
      xint=zero
      yint=zero
      zint=zero
      npts=(ni+nj+nc+nd-4)/2+1
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
         cx=ptx-xc
         cy=pty-yc
         cz=ptz-zc
         dx=ptx-xd
         dy=pty-yd
         dz=ptz-zd
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
         do kk = 1, nc-1
            px=px*cx
            py=py*cy
            pz=pz*cz
         enddo
         do ll = 1, nd-1
            px=px*dx
            py=py*dy
            pz=pz*dz
         enddo
         xint=xint+px
         yint=yint+py
         zint=zint+pz
   16 continue
      return
      end
