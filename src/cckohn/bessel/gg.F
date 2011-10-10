      function gg(x)
c
      implicit none
      double precision x,gg
      integer nsmall,lbig
      parameter (nsmall=4015,lbig=6)
      double precision xstart,ystart,ypstart,al,mu,h,xmax
      common /integ/ xstart,ystart,ypstart,al,mu,h,xmax
      double precision e,eta,ak
      common /energy/ e,eta,ak
      double precision rd,r,w,cs,rd26
      integer l
      common/spl/rd,r(nsmall),w(nsmall,0:lbig),cs(nsmall),rd26,l
      double precision klo,a,b
      integer iklo
c
      klo=(x-xstart)/rd
      klo=klo+1
      iklo=klo
      a=(r(iklo+1)-x)/rd
      b=(x-r(iklo))/rd
      gg=a*w(iklo,l)+b*w(iklo+1,l)+(a*(a*a-1)*cs(iklo)+b*
     1 (b*b-1)*cs(iklo+1))*rd26
      return
      end
c $Id$
