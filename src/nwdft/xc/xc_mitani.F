      Subroutine xc_mitani(r,w,nr,range,lreset)

C$Id$
      implicit none

      double precision r(*) ! grid pts coord [output]
      double precision w(*) ! grid pts weights [output]
c
      double precision rm    
      parameter (rm=3d0)
c
      integer nr
      double precision range ! max extent of basis f
c
      logical lreset
c      
      integer i
      double precision h ! mesh size
      double precision xi,ri,wi,alpha
c
      integer n_lo, n_hi
      integer attempts,max_attempts
      double precision alpha_old
      double precision h0,h1,h_step,fmin,fmax,fmed
      double precision xc_de2_h
      external xc_de2_h
      integer max_loops
      parameter(max_loops=10)
      double precision eps
      parameter(eps=1d-5)
c double exponential DE2 grid      
c M. Mitani
cTheor Chem Acc (2011) 130:645–669
c DOI 10.1007/s00214-011-0985-x      
c***************************************************************************
c
c     set mesh size
c     Mitani: r= 10R for GTO; N*h = 10R
c
c     h = 10d*R/nr
c     h = range/(nr-1)
c     h = log(range)/(nr-n0)
c     odd and even
      alpha=1d0
c     Dasgupta and Herbert        DOI: 10.1002/jcc.24761
c     2.2 <alpha < 3.2
      if(mod(nr,2).eq.1) then
         n_lo=-(nr-1)/2
      else
         n_lo=-nr/2
      endif
      n_hi = nr + n_lo -1
      h = log(range)/(alpha*n_hi)
      h_step=h*1d-2
      do i=1,max_loops
         write(6,*) 'h loop i',i
         h0 = h - h_step
         h1 = h + h_step
         fmin = xc_de2_h(h - h_step,alpha,n_hi,log(range))
         fmax = xc_de2_h(h + h_step,alpha,n_hi,log(range))
         fmed = xc_de2_h(h ,alpha,n_hi,log(range))
         write(6,*) ' fmin ',fmin
         write(6,*) ' fmed ',fmed
         write(6,*) ' fmax ',fmax
         if((fmin*fmax).gt.0d0)then
            h_step = h_step*2d0
         else
            goto 8
         endif
      enddo
 8    continue
      write(6,*) ' out of loop '
      attempts=0
 1    continue
      attempts=attempts+1
      if(attempts.gt.20) then
#ifdef DEBUG
         if(ga_nodeid().eq.0) then
         write(LuOut,*) ' efermi emax emin ',
     ,        efermi,emax,emin
         write(LuOut,*) ' fmed fmax fmin ',fmed,fmax,fmin
         write(LuOut,*) 'ZERO could no reach required accuracy'
         endif
#endif
         write(6,*) ' did not cvg '
         goto 1999
       endif
      if(fmax*fmed.gt.0) then
c	zero between emin and efermi
         h1=h
         fmax=fmed
      else	
c     zero between fmed and fmax
         h0=h
         fmin=fmed
      endif
      if(abs(fmax-fmin).le.0.d0) then
         write(0,*) ' WARNING: fmax=fmin'
      endif
      if(attempts.lt.15.or.abs(fmax-fmin).lt.0.d0) then
      h=0.5d0*(h0+h1)
      else
         h=-fmin*(h1-h0)/(fmax-fmin)+h0
      endif
      
      fmed = xc_de2_h(h ,alpha,n_hi,log(range))
      if(abs(fmed).gt.eps.and.abs(h1-h0).gt.eps*1d-1) goto 1
1999  continue
      write(6,*) ' h is ',h, ' in ',attempts,' attempts'
      write(6,*) ' Mitani: n_lo n_hi ',n_lo,n_hi
c     
      if(lreset) then
         alpha=1d0
         alpha_old=alpha
         h = log(range)/(alpha*n_hi)
         write(6,*) ' before reset: h ',h
         xi=n_hi*h
         ri = exp(alpha*xi - exp(-xi))
         write(6,*) ' before rest: rnhi  ', ri
         do i=1,max_loops
            write(6,*) 'reset loop i',i
            alpha=(log(range) + exp(-n_hi*h))/(n_hi*h)
            write(6,*) ' reset: alpha is ',alpha
            h = log(range)/(alpha*n_hi)
            write(6,*) ' reset: h is     ', h
            xi=n_hi*h
            ri = exp(alpha*xi - exp(-xi))
            write(6,*) ' reset: rnhi i   ', ri
            write(6,*) ' delta alpha',dabs(alpha_old-alpha)
            alpha_old=alpha
         enddo
      endif
 10   continue
      write(6,*) ' Mitani: R_max ',range
      write(6,*) ' Mitani: alpha ',alpha
      write(6,*) ' Mitani: h ',h

      do i = n_lo,n_hi

         xi = dble(i)*h
         ri = exp(alpha*xi - exp(-xi))
          write(6,*) ' Mitani: i xi Rrrr ',i,xi,ri,i-n_lo+1
         wi = h*exp(3*alpha*xi - 3*exp(-xi))*(alpha + exp(-xi))
         
         r(i-n_lo+1) = ri
         w(i-n_lo+1) = wi

      enddo
      write(6,*) ' Mitani: r0 ',r(1)
      write(6,*) ' Mitani: rn ',r(nr)

      end
      double precision function xc_de2_h(h,alpha,n,lnr)
      implicit none
      double precision h,alpha,lnr
      integer n
      xc_de2_h=alpha*n*h - exp(-n*h)  - lnR
      return
      end
