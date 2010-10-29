      subroutine xc_tarad(znumber,r,w,nr,origin,which,rlast)

C$Id$
      implicit none

      double precision r(*) ! grid pts coord [output]
      double precision w(*) ! grid pts weights [output]
      double precision origin ! origin of the grid [input]
      integer znumber !   atomic number [input]
c
      integer nr,i
      double precision xi,wi
      double precision ri,zetad(36),zeta,rlast
      character*8 which
      double precision alpha,rsc(400),wsc(400)
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
      if(znumber.lt.37) zeta = zetad(znumber)

      if(which.eq.'chebyshe') then
      do  i = 1,nr
        call grid_gausscheb(wi,xi,i,nr)
        ri = zeta/log(2.D0)*(1.d0+xi)**alpha*log(2d0/(1d0-xi))
        r(i) = ri+origin
        w(i)=wi*r(i)*r(i)*ri*(alpha/(1d0+xi)+
     +       1d0/(1d0-xi)/log(2d0/(1d0-xi)))
c
c       becke transf
c        ri=(1d0+xi)/(1d0-xi)
c        w(i) = wi*ri*ri/(1d0-xi)**2*2d0
      enddo
      elseif(which.eq.'chebyshr') then
c
c     get zeta from rmax
c
         call grid_gausscheb(wi,xi,nr,nr)
c M4
c         zeta=rlast*log(2d0)/
c     /        ((1.d0+xi)**alpha*log(2d0/(1d0-xi)))
c        ri = zeta/log(2.D0)*(1.d0+xi)**alpha*log(2d0/(1d0-xi))
c M3
         zeta=(rlast-origin)*log(2d0)/
     /        (log(2d0/(1d0-xi)))
        ri = zeta/log(2.D0)*log(2d0/(1d0-xi))
        r(nr) = ri+origin
        w(nr)=wi*r(nr)*r(nr)*ri*(alpha/(1d0+xi)+
     +       1d0/(1d0-xi)/log(2d0/(1d0-xi)))
         
      do  i = 1,nr-1
         call grid_gausscheb(wi,xi,i,nr)
        ri = zeta/log(2.D0)*(1.d0+xi)**alpha*log(2d0/(1d0-xi))
        r(i) = ri+origin
        w(i)=wi*r(i)*r(i)*ri*(alpha/(1d0+xi)+
     +       1d0/(1d0-xi)/log(2d0/(1d0-xi)))
c
c       becke transf
c        ri=(1d0+xi)/(1d0-xi)
c        w(i) = wi*ri*ri/(1d0-xi)**2*2d0
      enddo
      elseif(which.eq.'legendre') then
        call grid_gaussleg(nr, rsc, 
     &       wsc)
        do i=1,nr
          xi=rsc(nr-i+1)
          ri = zeta/log(2.D0)*(1.d0+xi)**alpha*log(2d0/(1d0-xi))
          r(i) = ri+origin
          r(i) = ri+origin
          wi = wsc(nr-i+1)
          w(i)=wi*r(i)*r(i)*ri*(alpha/(1d0+xi)+
     +         1d0/(1d0-xi)/log(2d0/(1d0-xi)))
        enddo
      endif
      return
      end
      subroutine xc_interv(rmin,rmax,r,w,nr,which)
      implicit none
#include "errquit.fh"
c
c     compute Gauss-Chebyshev quadrature for the 
c     interval rmin<r<rmax 
c     
      double precision rmin ! interval minimum [input]
      double precision rmax ! interval maximum [input]
      double precision r(*) ! grid pts coord [output]
      double precision w(*) ! grid pts weights [output]
      integer nr            ! [input]
      character*8 which
c
      double precision slope,intercept
      double precision wi,xi,ri,wsc(300),rsc(300)
      integer i
c
      if(nr.gt.300) call errquit(
     (    ' xc_interv: too many radial pts ',nr, UNKNOWN_ERR)
      slope=(rmax-rmin)*.5d0
      intercept=(rmin+rmax)*.5d0
      if(which.eq.'chebyshe') then
        do i=1,nr
          call grid_gausscheb(wi,xi,i,nr)
          ri = slope*xi + intercept
          r(i) = ri
          w(i) = wi*ri*ri*slope
        enddo
      elseif(which.eq.'legendre') then
        call grid_gaussleg(nr, rsc, 
     &       wsc)
        do i=1,nr
          ri = slope*rsc(nr-i+1) + intercept
          r(i) = ri
          w(i) = wsc(nr-i+1)*ri*ri*slope
        enddo
          else
        call errquit('unrecognized quadr. ',999, INPUT_ERR)
      endif
      return
      end

      subroutine grid_gausscheb(wi,xi,i_in,n)
      implicit none
      double precision wi,xi ! [output]
      integer i_in,n ! [input]
c
c     Gauss Chebyshev quadrature scheme
c     JM Perez-Yorda and E San-Fabian, CPC 77, 46 (1993)
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
