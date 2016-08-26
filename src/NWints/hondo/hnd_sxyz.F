c
c $Id$
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
      integer ni, nj, npts, imin, imax, i, ii, jj
      data zero /0.0d+00/
c
      xint=zero
      yint=zero
      zint=zero
      npts=(ni+nj-2)/2+1
      imin=hermin(npts)
      imax=hermax(npts)
!DEC$ LOOP COUNT MAX=9, MIN=1
      do 16 i=imin,imax
         dum=w(i)
         px=dum
         py=dum
         pz=dum
         dum=h(i)*t
         ptx=dum+x0
         pty=dum+y0
         ptz=dum+z0
         if(ni.gt.1) then
            ax=ptx-xi
            ay=pty-yi
            az=ptz-zi
         endif
         if(nj.gt.1) then
            bx=ptx-xj
            by=pty-yj
            bz=ptz-zj
         endif
         if(ni.eq.nj) then
!DEC$ LOOP COUNT MAX=9, MIN=1
            do ii = 1, ni-1
               px=px*ax*bx
               py=py*ay*by
               pz=pz*az*bz
            enddo
         else
!DEC$ LOOP COUNT MAX=9, MIN=1
         do ii = 1, ni-1
            px=px*ax
            py=py*ay
            pz=pz*az
         enddo
!DEC$ LOOP COUNT MAX=9, MIN=1
         do jj = 1, nj-1
            px=px*bx
            py=py*by
            pz=pz*bz
         enddo
      endif
         xint=xint+px
         yint=yint+py
         zint=zint+pz
   16 continue
      return
      end
