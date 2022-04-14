      Subroutine xc_mitani_dmr(r,w,nr,range,lreset, alpha_de2)
c double exponential DE2 grid      
c M. Mitani
cTheor Chem Acc (2011) 130:645–669
c DOI 10.1007/s00214-011-0985-x      
C$Id$
      implicit none

      double precision r(*) ! grid pts coord [output]
      double precision w(*) ! grid pts weights [output]
c
      integer nr                 ! no. grid points [input]
      double precision range     ! max extent of basis f
      double precision alpha_de2 ! alpha [input]
c
      logical lreset
#include "global.fh"
#include "util.fh"
c      
      integer i,maxiter
      double precision h ! mesh size
      double precision xi,ri,wi,alpha
c
      integer n_lo, n_hi
      integer attempts,max_attempts
      logical oprint
      double precision alpha_old, delta
      double precision lambertW
      external lambertW
      integer max_loops
      parameter(max_loops=10)
      integer n_param
      parameter (n_param=3)
      double precision paramn(n_param)
      double precision wrmin,wrmax,lrs,lrmin,lrmax
      double precision eps, alpha_0,rmin, invalpha
      parameter(eps=1d-6,alpha_0=1.d0,rmin=1.0d-7)
      parameter(maxiter=100)
c     Mitani 2011                                          alpha: 1
c     Mitani-Yoshioka 2012  DOI: 10.1007/s00214-012-1169-z alpha: 1.9
c     Dasgupta and Herbert  DOI: 10.1002/jcc.24761c       2.2 <alpha < 3.2
c     alpha can be set from input with the directive 'set dft:alpha_de2 1.9'
c***************************************************************************
c     print options
      oprint = util_print('quadrature', print_high).and.
     A     (ga_nodeid().eq.0)
c
      if(alpha_de2.ne.0d0) then
         alpha=alpha_de2
      else
         alpha=alpha_0
      endif

      invalpha = 1d0/alpha

      wrmin = lambertW(1d0/(alpha*rmin**invalpha))
      wrmax = lambertW(1d0/(alpha*range**invalpha))
      lrmax = log(range)
      lrmin = log(rmin)
      lrs = log(range/rmin)

      n_hi=nint((lrmax + alpha*wrmax)*(nr-1)/(lrs+alpha*(wrmax-wrmin)))
      n_lo=n_hi+1-nr

      ! we will adjust h to fix rmin = 1d-7
      h = (lrs + alpha*(wrmax - wrmin))/(alpha*(nr-1))
      do i=1,maxiter
        delta = (alpha*n_lo*h - exp(-n_lo*h) - lrmin)/
     $          (alpha*n_lo + n_lo*exp(-n_lo*h))
        if (abs(delta/h).lt.1d-16) exit
        h = h - delta
      enddo

      if(oprint) then
         write(6,*) ' Mitani: R_max ',range
         write(6,*) ' Mitani: alpha ',alpha
         write(6,*) ' Mitani: h ',h
         write(6,*) ' Mitani: n_lo ',n_lo
         write(6,*) ' Mitani: n_hi ',n_hi
         write(6,*) ' Mitani: R_max ',exp(alpha*n_hi*h-exp(-n_hi*h))
         write(6,*) ' Mitani: R_min ',exp(alpha*n_lo*h-exp(-n_lo*h))
      endif

      do i = n_lo,n_hi

         xi = dble(i)*h
         ri = exp(alpha*xi - exp(-xi))
         wi = h*exp(3*alpha*xi - 3*exp(-xi))*(alpha + exp(-xi))
         
         r(i-n_lo+1) = ri
         w(i-n_lo+1) = wi

      enddo

      end
      double precision function lambertW(z)
      implicit none
      double precision z,delta
      integer maxiter,i
      parameter (maxiter = 100)
      if (z.eq.1d0) then
        lambertW = 0.56714329040978387300d0
        return
      elseif (z.lt.1d0) then
        lambertW = min(0.5671d0,z)
      elseif (z.lt.2.72d0) then
        lambertW = min(max(0.5671d0,z),1d0)
      else
        lambertW = log(z)
      endif

      do i=1,maxiter
        delta=(lambertW*exp(lambertW)-z)/(exp(lambertW)*(1d0+lambertW))
        if (abs(delta/lambertW).lt.1d-16) exit
        lambertW = lambertW - delta
      enddo
      if (i.gt.maxiter) write(6,*) 'something wrong'


      return
      end
