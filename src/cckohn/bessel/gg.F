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
c
      klo=(x-xstart)/rd
      klo=klo+1
      a=(r(klo+1)-x)/rd
      b=(x-r(klo))/rd
      gg=a*w(klo,l)+b*w(klo+1,l)+(a*(a*a-1.)*cs(klo)+b*
     1 (b*b-1.)*cs(klo+1))*rd26
      return
      end
