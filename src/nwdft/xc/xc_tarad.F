      Subroutine xc_tarad(znumber,r,w,nr)

C$Id: xc_tarad.F,v 1.1 1999-01-13 03:07:13 edo Exp $
      implicit none

      double precision r(*) ! grid pts coord [output]
      double precision w(*) ! grid pts weights [output]
      integer znumber !   atomic number [input]
c
      integer nr,i
      double precision xi,wi
      double precision ri,wi,zetad(36),zeta
      double precision alpha,wsum
      parameter (alpha=0.6d0)
      data zetad/ 0.8d0, 0.9d0,
c Li-N
     *     1.8d0,1.4d0,1.3d0,1.1d0,0.9d0,0.9d0,0.9d0,0.9d0,
C Na-Ar
     *     1.4d0,1.3d0,1.3d0,1.2d0,1.1d0,1.0d0,1.0d0,1.0d0,
C K-Kr
     *     1.5d0,1.4d0,1.3d0,1.2d0,1.2d0,1.2d0,1.2d0,1.2d0,1.2d0,1.1d0,
C Cu-Kr
     *     1.1d0,1.1d0,1.1d0,1.0d0,0.9d0,0.9d0,0.9d0,0.9d0/
c
c
c Definition of Gauss-Chebyshev numerical quadrature points and weights
c for radial integrals.
c Transformation from 0<r<infty to -1<x<1 according to
c O Treutler and R Alrhichs, J Chem Phys 102, 346 (1995)
c
c***************************************************************************
      zeta=1.d0
      wsum=0.d0
      if(znumber.lt.37) zeta = zetad(znumber)

      do  i = 1,nr

        call gausscheb(wi,xi,i,nr)
        ri = zeta/log(2.D0)*(1.d0+xi)**alpha*log(2d0/(1d0-xi))
c        ri=(1d0+xi)/(1d0-xi)
c        write(6,*) znumber,xi,ri,wi

        r(i) = ri
c        w(i) = wi*ri*ri/(1d0-xi)**2*2d0
        w(i)=wi*ri*ri*ri*(alpha/(1d0+xi)+
     +       1d0/(1d0-xi)/log(2d0/(1d0-xi)))
c
c       check normalization
c
        wsum=wsum+wi
      enddo
      write(6,*) ' sum is ',wsum
      return
      end
      subroutine gausscheb(wi,xi,i_in,n)
      implicit none
      double precision wi,xi ! [output]
      integer i_in,n ! [input]
c
      double precision pi
      integer i
c
      pi=acos(-1d0)
      i=n-i_in+1
      wi=sin(pi*i/(1d0+n))**4/(3d0*(1d0+n))*16d0
      xi=(1.d0+n-i-i)/(1d0+n)+
     +     2d0/pi*(1d0+2d0/3d0*sin(pi*i/(1d0+n))**2)*
     *     cos(pi*i/(1d0+n))*sin(pi*i/(1d0+n))
      return
      end
