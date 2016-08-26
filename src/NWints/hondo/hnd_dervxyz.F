c
c $Id$
c
c Taken and modified from HONDO
c
      subroutine hnd_dervxyz(mder)
      implicit none
c
c     Routine can do ddvxyz and dvxyz by setting mder to 2 or 1
c
c     ----- gauss-hermite quadrature using minimum point formula -----
c
#include "hnd_whermt.fh"
c
      double precision xint, yint, zint, t, x0, y0, z0, xi, yi, zi
      double precision xj, yj, zj, cx, cy, cz, zero, ptx, pty, ptz
      double precision dum, px, py, pz, ax, ay, az, bx, by, bz
      integer ni, nj, npts, imin, imax, i, mder, ii, jj
      common/hnd_xyzder/xint,yint,zint,t,x0,y0,z0,xi,yi,zi,xj,yj,zj
     1                  ,ni,nj,cx,cy,cz
      data zero /0.0d+00/
      xint = zero
      yint = zero
      zint = zero
      npts = (ni+nj+mder-2)/2+1
      imin = hermin(npts)
      imax = hermax(npts)
      do 11 i = imin,imax
         dum = h(i)*t
         ptx = dum+x0
         pty = dum+y0
         ptz = dum+z0
         px = ptx-cx
         py = pty-cy
         pz = ptz-cz
         px = px**mder
         py = py**mder
         pz = pz**mder
         ax = ptx-xi
         ay = pty-yi
         az = ptz-zi
         bx = ptx-xj
         by = pty-yj
         bz = ptz-zj
         if(ni.eq.nj) then
         do ii = 1, ni-1
            px=px*ax
            py=py*ay
            pz=pz*az
            px=px*bx
            py=py*by
            pz=pz*bz
         enddo
         else
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
         endif
         dum = w(i)
         xint = xint+dum*px
         yint = yint+dum*py
         zint = zint+dum*pz
   11 continue
      return
      end
