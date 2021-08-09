      Subroutine xc_mitani(r,w,nr,range,lreset, alpha_de2)
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
      integer i
      double precision h ! mesh size
      double precision xi,ri,wi,alpha
c
      integer n_lo, n_hi
      integer attempts,max_attempts
      logical oprint
      double precision alpha_old
      double precision xc_de2_h
      external xc_de2_h
      integer max_loops
      parameter(max_loops=10)
      integer n_param
      parameter (n_param=3)
      double precision paramn(n_param)
      double precision eps, alpha_0
      parameter(eps=1d-14,alpha_0=1.d0)
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
c     odd and even number of grid points
      if(mod(nr,2).eq.1) then
         n_lo=-(nr-1)/2
      else
         n_lo=-nr/2
      endif
      n_hi = nr + n_lo -1
c     approximate value of h
      h = log(range)/(alpha*n_hi)
      paramn(1)=alpha
      paramn(2)=n_hi
      paramn(3)=log(range)
c     find h value
      call xc_bisect(h, xc_de2_h, 3, paramn)
      if(oprint) then
         write(6,*) 'first root search '
      write(6,*) ' h is ',h
      write(6,*) ' residual: ',xc_de2_h(h, 3, paramn)
      write(6,*) ' Mitani: n_lo n_hi ',n_lo,n_hi
      write(6,*) ' Mitani: alpha ',alpha
      endif
c     
      alpha_old=alpha
      if(lreset) then
         h = log(range)/(alpha*n_hi)
         xi=n_hi*h
         ri = exp(alpha*xi - exp(-xi))
         do i=1,max_loops
            alpha=(log(range) + exp(-n_hi*h))/(n_hi*h)
            h = log(range)/(alpha*n_hi)
            paramn(1)=alpha
            call xc_bisect(h, xc_de2_h, 3, paramn)
            if(dabs(alpha_old-alpha).lt.eps) goto 10
            alpha_old=alpha
         enddo
      endif
 10   continue
      if(oprint) then
         write(6,*) 'second root search '
         write(6,*) ' Mitani: R_max ',range
         write(6,*) ' Mitani: alpha ',alpha
         write(6,*) ' Mitani: h ',h
         write(6,*) ' residual: ',xc_de2_h(h, 3, paramn)
         if(alpha_old.ne.alpha) then
            write(6,*) 'reset loop i',i
            write(6,*) ' reset: rnhi i   ', ri
            write(6,*) ' delta alpha',dabs(alpha_old-alpha)
         endif
      endif

      do i = n_lo,n_hi

         xi = dble(i)*h
         ri = exp(alpha*xi - exp(-xi))
         wi = h*exp(3*alpha*xi - 3*exp(-xi))*(alpha + exp(-xi))
         
         r(i-n_lo+1) = ri
         w(i-n_lo+1) = wi

      enddo

      end
      double precision function xc_de2_h(h,n_param,paramn)
      implicit none
      double precision h
      integer n_param
      double precision paramn(n_param)
c
      double precision alpha,lnr
      integer n
c
      if(n_param.ne.3) call errquit(' incomp func ',0,0)
      alpha =  paramn(1)
      n     =  paramn(2)
      lnr   =  paramn(3)        
      xc_de2_h=alpha*n*h - exp(-n*h)  - lnR
      return
      end
      subroutine xc_bisect(x,func3,n_param,paramn)
      implicit none
      double precision x ! [in:guess,out:solution]
      integer n_param
      double precision  paramn(n_param)
c      
      double precision func3
      external func3
c
      double precision x0,x1,fmed,fmin,fmax
      double precision step
      integer attempts
      integer i,max_loops,mx_attempts
      parameter(max_loops=10,mx_attempts=20)
      double precision eps
      parameter(eps=1d-14)
c
      step=x/100d0
      do i=1,max_loops
         x0 = x - step
         x1 = x + step
         fmin = func3(x0, n_param, paramn)
         fmax = func3(x1, n_param, paramn)
         fmed = func3(x , n_param, paramn)
         if((fmin*fmax).gt.0d0)then
            step = step*2d0
         else
            goto 8
         endif
      enddo
      write(6,*) ' xc_bisect: sign fmin = fmax'
      call errquit('xc_bisect: err1 ',0,0)
 8    continue
      attempts=0
 1    continue
      attempts=attempts+1
      if(attempts.gt.mx_attempts) then
         write(6,*) ' xc_bisect did not cvg in ',attempts,' attempts'
         call errquit('xc_bisect: err2 ',0,0)
         goto 1999
       endif
      if(fmax*fmed.gt.0) then
c	zero between emin and efermi
         x1=x
         fmax=fmed
      else	
c     zero between fmed and fmax
         x0=x
         fmin=fmed
      endif
      if(abs(fmax-fmin).le.0.d0) then
         write(0,*) ' WARNING: fmax=fmin'
      endif
      if(attempts.lt.mx_attempts.or.abs(fmax-fmin).lt.0.d0) then
      x=0.5d0*(x0+x1)
      else
         x=-fmin*(x1-x0)/(fmax-fmin)+x0
      endif
      
      fmed = func3(x , n_param, paramn)
      if(abs(fmed).gt.eps.and.abs(x1-x0).gt.eps*1d-1) goto 1
 1999 continue
#ifdef DEBUG      
      write(6,*) ' x is ',x, ' in ',attempts,' attempts'
#endif      
      return
      end
      double precision function xc_de2alphasg2(znumber)
      implicit none
      integer znumber ! [input]
#include "errquit.fh"
      integer maxz
      parameter(maxz=17)
c     alpha SG-2 from Dasgupta and Herbert
      double precision alphad(maxz)
cH
      data alphad(1)/2.6d0/
cHe fit
      data alphad(2)/1.6d0/
cLi 	
      data alphad(3)/3.2d0/
cBe 	
      data alphad(4)/2.4d0/
cB 	
      data alphad(5)/2.4d0/
cC 	
      data alphad(6)/2.2d0/
cN 	
      data alphad(7)/2.2d0/
cO 
      data alphad(8)/2.2d0/
cF 	
      data alphad(9)/2.2d0/
cNa 	
      data alphad(11)/3.2d0/
cNe
      data alphad(10)/2.3d0/
cMg 	
      data alphad(12)/2.4d0/
cAl 	
      data alphad(13)/2.5d0/
cSi 	
      data alphad(14)/2.3d0/
cP 	
      data alphad(15)/2.5d0/
cS 	
      data alphad(16)/2.5d0/
cCl 	
      data alphad(17)/2.5d0/
      if(znumber.gt.maxz) call errquit(
     c     ' xc_de2alphasg2 z gt ',maxz,0)
      xc_de2alphasg2=alphad(znumber)
      if(xc_de2alphasg2.eq.0d0) call errquit(
     c     ' alpha not available for  z= ',znumber,0)
      return
      end
      double precision function xc_de2alphasg3(znumber)
      implicit none
      integer znumber ! [input]
#include "errquit.fh"
      integer maxz
      parameter(maxz=17)
c     alphaS SG-3 from Dasgupta and Herbert
      double precision alphad(maxz)
c     H
      data alphad(1)/2.7d0/
cHe fit
      data alphad(2)/1.6d0/
cLi 	
      data alphad(3)/3.0d0/
cBe 	
      data alphad(4)/2.4d0/
cB 	
      data alphad(5)/2.4d0/
cC 	
      data alphad(6)/2.4d0/
cN 	
      data alphad(7)/2.4d0/
cO 
      data alphad(8)/2.6d0/
cF 	
      data alphad(9)/2.1d0/
cNa 	
      data alphad(11)/3.2d0/
cNe fit
      data alphad(10)/2.2d0/
cMg 	
      data alphad(12)/2.6d0/
cAl 	
      data alphad(13)/2.6d0/
cSi 	
      data alphad(14)/2.8d0/
cP 	
      data alphad(15)/2.4d0/
cS 	
      data alphad(16)/2.4d0/
cCl 	
      data alphad(17)/2.6d0/
      if(znumber.gt.maxz) call errquit(
     c     ' xc_de2alphasg3 z gt ',maxz,0)
      xc_de2alphasg3=alphad(znumber)
      if(xc_de2alphasg3.eq.0d0) call errquit(
     c     ' alpha not available for  z= ',znumber,0)
      return
      end
