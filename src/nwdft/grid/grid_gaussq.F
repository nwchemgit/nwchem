        subroutine grid_gausslag(nn,x,realwt,idum)
c
c$Id$
c
c
c this subroutine calculates the Gauss-Laguerre points and weights for
c  associated Laguerre polynomials using subroutines provided in
c  Stroud and Secrest.  These routines are modified to allow continuous
c  scaling of the coefficient relationship, thereby extending the
c  range to very high order.  This is done in lgrecr.
c  The routine works easily up to order 400.
c
c  Input:  nn=order of the integration
c  Output: x(i),i=1,..,nn are the integration points
c          realwt(i),i=1,..,nn are the integration weights for the integral
c               int(0,inf) f(x) x**2 dx.
c          This is changed to
c               int(0,inf) exp(x)x**(2-alf)f(x) x**alf exp(-x) dx.
c          The real weights are thus a(i)*exp(x)*x**(2-alf) where a(i)
c          are the original Gauss-Laguerre weights.
c
C$Id$
      implicit none
#include "cdft.fh"      
#include "stdio.fh"      
c
c     The following arg was added to match the calls to this routine,
c     all of which have 4 arguments.  This is necessary to link under
c     WIN32.  Someone should fix this so that all 4 args are properly
c     used. BGJ (9/99)
c
      integer idum ! [input, !!!Dummy argument!!!]
      integer nn ! [input]
      double precision x(1),realwt(1)
C
      integer n
      double precision a(MXrad),b(MXrad),c(MXrad)
      double precision alf,eps
      double precision csx,tsx,tsa,csa
      parameter(alf=1.5d0,eps=1.d-10)
c
      do n=1,nn
        b(n)=alf+2*n-1
        c(n)=(n-1.d0)*(alf+n-1.d0)
      enddo
        call grid_laguer(nn,x,a,realwt,alf,b,c,eps,csx,csa,tsx,tsa)
c       write(LuOut,100)(x(i),realwt(i),i=1,nn)
100     format(2D18.10)
        return
        end
c
c
c
c
c
        subroutine grid_laguer(nn,x,a,realwt,alf,b,c,eps,csx,csa,
     ,     tsx,tsa)
        implicit none
        
        double precision x(*),realwt(*),eps,csx,csa,tsx,tsa,alf
        double precision a(*),b(*),c(*)
        integer nn
C
        double precision fn,cc,xt,fi,r1,r2,ratio,dpn,scale,pn1
        integer i,j
        double precision grid_gamaed
        external grid_gamaed
C
        fn=nn
        csx=0.d0
        csa=0.d0
        cc=log(grid_gamaed(alf+1.d0))
        tsx=fn*(fn+alf)
        tsa=cc
        do 1 j=2,nn
1       cc=cc+log(c(j))
        do 7 i=1,nn
c       write(LuOut,*)i
        if (i.eq.1) then
c smallest zero
2       xt=(1.d0+alf)*(3.d0+.92d0*alf)/(1.d0+2.4d0*fn+1.8d0*alf)
        elseif (i.eq.2) then
c second zero
4       xt=xt+(15.d0+6.25d0*alf)/(1.d0+2.5d0*fn+.9d0*alf)
        elseif(i.gt.2) then
c all other zeros
5       fi=i-2
        r1=(1.d0+2.55d0*fi)/(1.9d0*fi)
        r2=1.26d0*fi*alf/(1.d0+3.5d0*fi)
        ratio=(r1+r2)/(1.d0+.3d0*alf)
        xt=xt+ratio*(xt-x(i-2))
        endif
c
6       call grid_lgroot(xt,nn,alf,dpn,pn1,b,c,eps,scale)
c       write(LuOut,*)'back from lgroot'
        x(i)=xt
c       a(i)=sign(1.d0,dpn)*(exp(cc-log(abs(dpn))))/pn1
c       write(LuOut,20)alf,nn,i,xt,a(i)
        realwt(i)=sign(1.d0,dpn)*sign(1.d0,pn1)*
     &   exp(xt+(2.d0-alf)*log(xt)+cc-
     &   2*scale-log(abs(dpn))-log(abs(pn1)))
        csx=csx+xt
7       csa=csa+a(i)
c       write(LuOut,20)alf,nn,i,csx,csa,tsx,tsa
20      format(f6.2,2i3,2(1x,d14.8),2x,2(1x,d14.8))
        return
        end
c
c
c
c
        subroutine grid_lgroot(x,nn,alf,dpn,pn1,b,c,eps,scale)
        implicit none
        integer nn
        double precision b(*),c(*),x,alf,dpn,pn1,eps,scale
c
        integer iter
        double precision d,p,dp
c
        iter=0
1       iter=iter+1
        call grid_lgrecr(p,dp,pn1,x,nn,alf,b,c,scale)
        d=p/dp
        x=x-d
        if(abs(d/x)-eps)  3,3,2
2       if(iter.lt.10)  goto 1
3       dpn=dp
        return
        end
c
c
c
c
        subroutine grid_lgrecr(pn,dpn,pn1,x,nn,alf,b,c,scale)
        implicit none
        double precision b(*),c(*),pn,dpn,pn1,x,alf,scale
        integer nn
c
        double precision p1,p,dp1,dp,q,dq,scl1
        integer j,jmax,i,jmin
c
        p1=1.d0
        p=x-alf-1.d0
        dp1=0.d0
        dp=1.d0
        scale=0.d0
        jmax=1
c do all scaling first
        do 10 i=1,10000
        jmin=jmax+1
        jmax=min(nn,jmin+20)
          do 2 j=jmin,jmax
          q=(x-b(j))*p -c(j)*p1
          dq=(x-b(j))*dp + p -c(j)*dp1
          p1=p
          p=q
          dp1=dp
2         dp=dq
c scale everything down by the scale factor
          scl1=max(abs(p),abs(p1),abs(dp),abs(dp1))
          p=p/scl1
          p1=p1/scl1
          dp=dp/scl1
          dp1=dp1/scl1
          scale=scale+log(scl1)
        if(jmax.eq.nn)goto 12
10      continue
12      pn=p
        dpn=dp
        pn1=p1
        return
        end
c
c
c
c
        double precision function grid_gamaed(x)
        implicit none
#include "stdio.fh"      
        
        double precision x
c
        double precision z,za
        double precision y,gam
c
        gam(y)=(((((((.035868343d0*y - .193527818d0)*y + .482199394d0)*y
     1    - .756704078d0)*y + .918206857d0)*y - .897056937d0)*y
     2    + .988205891d0)*y - .577191652d0)*y + 1.
        z=x
        if(z)  1,1,4
1       grid_gamaed=0.d0
        write(LuOut,*)'error arg for gamma=',z
        go to 14
4       if(z-70.d0) 6,1,1
6       if(z-1.d0)  8,7,9
7       grid_gamaed=1.d0
        go to 14
8       grid_gamaed=gam(z)/z
        go to 14
9       za=1.d0
10      z=z-1.d0
        if(z-1.d0)  13,11,12
11      grid_gamaed=za
        go to 14
12      za=za*z
        go to 10
13      grid_gamaed=za*gam(z)
14      return
        end
c
      subroutine grid_gaussleg(nn,x,a)
c This subroutine calculates the Gauss-Legendre points and weights for
c Legendre polynomials using subroutines provided in Gaussian Quadrature
c Formulas by A. H. Stroud and D. Secrest.  These routines are for Jacobi
c polynomials specialized to the case of alf=bta=0.  The routine works easily
c up to order 200.
c
c  Input:  nn=order of the integration
c  Output: x(i),i=1,..,nn are the integration points
c          a(i),i=1,..,nn are the integration weights for the integral
c               int(-1,1) f(x) dx.
c
      implicit none
#include "cdft.fh"      
#include "stdio.fh"      
      integer nn
      double precision c(MXrad),x(*),a(*)
c
      integer n
      double precision eps,csx,csa,tsx,tsa
c
      eps=2.D-15
      do 10 n=1,nn
       c(n)=(n-1.D0)**2/((2.D0*n-1.D0)*(2.D0*n-3.D0))
   10 continue
      call grid_legendre(nn,x,a,c,eps,csx,csa,tsx,tsa)
c     write(LuOut,*) nn
c     write(LuOut,'(2E24.16)') (x(i),a(i), i=1,nn)
c     write(LuOut,*) '(csx,csa)=',csx,csa
c     write(LuOut,*) '(tsx,tsa)=',tsx,tsa
      end
c
      subroutine grid_legendre(nn,x,a,c,eps,csx,csa,tsx,tsa)
      implicit none
      integer nn
      double precision c(nn),x(nn),a(nn),
     ,     csx,csa,tsx,tsa,eps
C
      double precision fn,beta,cc,xt,dpn,pn1
      integer j,i
C
      fn=nn
      csx=0.D0
      csa=0.D0
      beta=1.D0
      cc=2.D0
      tsx=0.D0
      tsa=cc
      do 1 j=2,nn
    1  cc=cc*c(j)
      do i=1,nn
c       if (i-1)  12,2,3
       if (i.eq.1) then
c largest zero
    2  xt=1.D0-2.78D0/(4.D0+fn*fn)
       elseif(i.gt.1) then
c    3  if (i-2)  12,4,5
         if(i.eq.2) then
c second zero
    4  xt=xt-4.1D0*(1.D0+.06D0*(fn-8.D0)/fn)*(1.D0-xt)
         elseif(i.gt.2) then
c    5  if (i-3)  12,6,7
           if(i.eq.3) then
c third zero
    6  xt=xt-1.67D0*(1.D0+.22D0*(fn-8.D0)/fn)*(x(1)-xt)
           elseif(i.gt.3) then
cold    7  if (nn-i-1)  10,9,8
    7  if (nn.gt.i-1)  then
c middle zeros
    8  xt=3.D0*x(i-1)-3.D0*x(i-2)+x(i-3)
       elseif(nn.eq.i-1) then
c second last zero
    9  xt=xt+(xt-x(i-2))/
     &       (0.766D0*(1.D0+.639D0*(fn-4.D0)/(1.D0+.71D0*(fn-4.D0))))
       else
c last zero
   10  xt=xt+(xt-x(i-2))/
     &       (1.67D0*(1.D0+.22D0*(fn-8.D0)/fn))
       endif
c
   11  call grid_root(xt,nn,dpn,pn1,c,eps)
       x(i)=xt
       a(i)=cc/(dpn*pn1)
       csx=csx+x(i)
       csa=csa+a(i)
        endif
      endif
      endif
       enddo
      end
c
      Subroutine grid_root(x,nn,dpn,pn1,c,eps)
c Improves the approximate root "x".  In addition the following quantities
c are obtained:
c
c     dpn = derivative of P(n) at x
c     pn1 = value of P(n-1) at x
c
      implicit none
#include "stdio.fh"      
      integer nn
      double precision c(nn)
      double precision pn,dpn,pn1,eps,x
c
      integer i
      double precision d
c
      do 10 i=1,10
       call grid_recur(pn,dpn,pn1,x,nn,c)
       d=pn/dpn
       x=x-d
       if(abs(d).le.eps) return
   10 continue
      write(LuOut,*)'WARNING:  FAILED TO CONVERGE A ROOT OF THE '//
     &                    'LEGENDRE POLYNOMIAL TO THE DESIRED'
      write(LuOut,*)'          PRECISION (SEE GAUSSLEG.F).'
      write(LuOut,*)'            DESIRED PRECISION:  ',eps
      write(LuOut,*)'            OBTAINED PRECISION:  ',abs(d)
      end
c
      subroutine grid_recur(pn,dpn,pn1,x,nn,c)
      implicit none

      integer nn
      double precision c(nn),pn,dpn,pn1,x
c
      double precision p1,p,q,dq,dp1,dp
      integer j
c
      p1=1.D0
      p=x
      dp1=0.D0
      dp=1.D0
      do 1 j=2,nn
       q=x*p-c(j)*p1
       dq=x*dp+p-c(j)*dp1
       p1=p
       p=q
       dp1=dp
       dp=dq
    1 continue
      pn=p
      dpn=dp
      pn1=p1
      end
