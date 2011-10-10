      subroutine hnd_denxyz
      implicit double precision (a-h,o-z)
      common/hnd_xyzder/xint,yint,zint,tx,x0,y0,z0,xi,yi,zi,xj,yj,zj,
     1                  ni,nj,cx,cy,cz
      data one /1.0d+00/
c
      xint=one
      yint=one
      zint=one
      ptxi=x0-xi
      ptyi=y0-yi
      ptzi=z0-zi
      ptxj=x0-xj
      ptyj=y0-yj
      ptzj=z0-zj
      do ii=1, ni-1
         xint=xint*ptxi
         yint=yint*ptyi
         zint=zint*ptzi
      enddo
      do jj=1, nj-1
         xint=xint*ptxj
         yint=yint*ptyj
         zint=zint*ptzj
      enddo
      return
      end
c $Id$
