#if defined(FUJITSU_VPP)
!ocl scalar
#endif
#ifndef SECOND_DERIV
      Subroutine xc_att_xc(rho,ipol,Ex,Amat,Cmat)
#else
      Subroutine xc_att_xc_d2(rho,ipol,Ex,Amat,Cmat,Amat2,Cmat2,Cmat3)
#endif
c
C$Id: xc_att_xc.F,v 1.5 2009-02-18 06:12:48 niri Exp $
c
      implicit none
c      
#include "stdio.fh"
c
c
ccase...start
#include "case.fh"
ccase...end
c
      double precision rho, Ex, Amat, Cmat
      integer ipol
#ifdef SECOND_DERIV
c
c     Second Derivatives of the Exchange Energy Functional
c
      double precision Amat2, Cmat2, Cmat3
#endif
c
c
c References:
c
c
c***************************************************************************
c
      double precision a, b, c, btmp,bfactor
c
      double precision a_first,a2_first,btmp_first, btmp1
c
      double precision sqrt_pi,t1,t2,t3,t4,t5,t6,t7
      double precision alpha,beta, DERF
      double precision f10, f01, b_first

      Parameter (sqrt_pi = 1.77245385090552d0)
      Parameter (t7 = 2.666666666666666667d0)
#ifdef SECOND_DERIV
c
c     Second Derivatives of the Exchange Energy Functional
c
      double precision a_second, a2_second, f20
      double precision b_second, btmp_second, t8
      double precision a3_second
 
#endif


c calculate the a_sigma parameter

c         write(luout,*) 'alpha',alpha
c         write(luout,*) 'beta',beta
c         write(luout,*) 'mu',mu
c           
          if (ipol.eq.1) then
             Ex = Ex/2d0
             rho = rho/2d0
          endif
          a = cam_omega*sqrt(-2d0*Ex)/(6d0*sqrt_pi*rho)
          alpha = cam_alpha
          beta = cam_beta
c 
          f10 = Amat/(2d0*Ex) -1d0/rho
          a_first = f10*a 
          f01 = Cmat/(2d0*Ex)
          a2_first = f01*a 
#ifdef SECOND_DERIV
          f20 = Amat2/(2d0*Ex) - Amat*Amat/(2d0*Ex*Ex) 
     &        + 1d0/(rho*rho)

          a_second = a*(f10*f10 + f20) 

          a2_second = a*(f10*f01 + Cmat2/(2d0*Ex)  
     &              - Amat*Cmat/(2d0*Ex*Ex))

          a3_second = a*(Cmat3/(2d0*Ex) - Cmat*Cmat/(4d0*Ex*Ex))
#endif


          if (a .lt. 0.14d0) then
c             write(luout,*) 'a is small'
           a = 2d0*a
           btmp = 1d0-(4d0/3d0)*sqrt_pi*a + 2d0*a*a - (2d0/3d0)*a*a*a*a
              btmp = 1d0-btmp

              btmp_first = (4d0/3d0)*(-sqrt_pi + 3d0*a +
     &                   (2d0*exp(-1/(a*a)) - 2d0)*a*a*a)  
              btmp_first = 2d0*btmp_first
             a = a /2d0
          else if (a .lt. 4.25d0) then
c            write(luout,*) 'a is medium'
c             stop
             b = exp(-1d0/(4d0*a*a)) - 1d0
             c = 2d0*a*a*b + 0.5d0
             btmp = (8d0/3d0)*a*(sqrt_pi*DERF(1/(2d0*a)) + 2d0*a*(b-c))
             t1 = 1/a
             t2 = a*a
             t3 = 1/t2
             t4 = exp(-0.25d0*t3)
             t5 = t4 -1d0
             t6 = t4 -2d0*t2*t5 - 1.5d0
             btmp_first = -t7*a *
     &       (2*a*(t4/(2*a**3) - 4d0*a*t5 - t1*t4) + 2d0*t6 -t3*t4) -
     &         t7*(2*a*t6 + sqrt_pi*DERF(0.5d0*t1))
          else
c            write(luout,*) 'a is large'
c            stop
             a = 2d0*a
             btmp = 1d0 - 1d0/(9d0*a*a) + 1d0/(60d0*a**4d0) -
     &           1d0/(420d0*a**6d0) + 1d0/(3240d0*a**8d0) -
     &           1d0/(27720d0*a**10d0)

             btmp_first = -1d0/(4.5d0*a**3) + 1d0/(15d0*a**5d0) -
     &                  1d0/(70d0*a**7d0) + 1d0/(405d0*a**9d0) 
             btmp_first = btmp_first*2d0
             a = a /2d0
          end if
#ifdef SECOND_DERIV
          if (abs(a) .lt. 1d-40) then 
             btmp_second = 16d0
          else if (a .ge. 5d0)  then
             btmp_second = 1d0/(6d0*a**4d0) - 1d0/(48d0*a**6d0) +
     &                  1d0/(640d0*a**8d0) - 1d0/(11520d0*a**10d0)
             
          else 
             t1 = a*a
             t2 = 1d0/t1
             t3 = exp(-0.25d0*t2)
             t4 = 1d0/(a*a*a)
             t5 = t3 - 1d0
             t6 = -t2*t3
             t8 = -t3/a + 0.5d0*t4*t3 - 4d0*a*t5
             btmp_second = -(8d0*a*(2d0*a*(t3/(4*a**6d0) -
     &       2d0*t3/(a**4d0) +t6 - 4d0*t5) -t3/(2*a**5d0) +
     &       4d0*t8 + 2d0*t4*t3)/3d0 + 16d0*(2d0*a*t8 +
     &       2d0*(t3 - 2d0*t1*t5-1.5d0) + t6)/3d0) 
          end if
#endif
          bfactor = 1d0 - alpha - beta*btmp
          b_first = beta*btmp_first
#ifdef SECOND_DERIV
          b_second = beta*btmp_second
          Amat2 = bfactor*Amat2 + 2d0*Amat*b_first*a_first 
     &          + Ex*b_second*a_first*a_first  
     &          + Ex*b_first*a_second

          Cmat2 = bfactor*Cmat2 + Amat*b_first*a2_first
     &          + Cmat*b_first*a_first
     &          + Ex*b_second*a_first*a2_first
     &          + Ex*b_first*a2_second
          
          Cmat3 = bfactor*Cmat3 + 2d0*Cmat*b_first*a2_first
     &          + Ex*b_second*a2_first*a2_first
     &          + Ex*b_first*a3_second
          
#endif
          Amat = bfactor*Amat + Ex*b_first*a_first
          Cmat = bfactor*Cmat + Ex*b_first*a2_first
          Ex = Ex*bfactor

          if (ipol.eq.1) then
             Ex = 2d0*Ex
             rho = 2d0*rho
          endif

c         write(luout,*) 'bfactor',bfactor
c
c

c
      return
      end
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "xc_att_xc.F"
#endif
