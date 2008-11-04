#if defined(FUJITSU_VPP)
!ocl scalar
#endif
      Subroutine xc_att_xc(rho,ipol,Ex,Amat,Cmat)
c
C$Id: xc_att_xc.F,v 1.1 2008-11-04 16:45:54 niri Exp $
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
c      Parameter (alpha = 0d0)
c      Parameter (beta = 1d0)
      Parameter (sqrt_pi = 1.77245385090552d0)
      Parameter (t7 = 2.666666666666666667d0)

c calculate the a_sigma parameter

c         write(luout,*) 'alpha',alpha
c         write(luout,*) 'beta',beta
c         write(luout,*) 'mu',mu
          a = cam_omega*sqrt(-4d0*Ex/ipol)/(6d0*sqrt_pi*rho)
c
          alpha = cam_alpha
          beta = cam_beta
c
c         write(luout,*) 'accessing xc_att_xc'
c         write(luout,*) 'alpha',alpha
c         write(luout,*) 'beta',beta
c         write(luout,*) 'mu',cam_omega
c 
          a_first = (Amat/(2d0*Ex) - 1d0/rho)*a 
          a2_first = (Cmat/(2d0*Ex))*a 

          if (a < 0.14d0) then
c             write(luout,*) 'a is small'
           a = 2d0*a
           btmp = 1d0-(4d0/3d0)*sqrt_pi*a + 2d0*a*a - (2d0/3d0)*a*a*a*a
              btmp = 1d0-btmp

              btmp_first = (4d0/3d0)*(-sqrt_pi + 3d0*a +
     &                   (2d0*exp(-1/(a*a)) - 2d0)*a*a*a)  
              btmp_first = 2d0*btmp_first
          else if (a <4.25d0) then
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
          end if

          bfactor = 1d0 - alpha - beta*btmp
          Amat = bfactor*Amat + Ex*beta*btmp_first*a_first
          Cmat = bfactor*Cmat + Ex*beta*btmp_first*a2_first
          Ex = Ex*bfactor

c         write(luout,*) 'bfactor',bfactor
c
c

c
      return
      end
