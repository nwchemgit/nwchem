      Subroutine xc_dirac(tol_rho,fac,lfac,nlfac,rho,Amat,nq,ipol,
     &                    Ex,qwght,ldew,func)

C$Id: xc_dirac.F,v 1.6 1997-02-15 01:34:15 gutowski Exp $
      Implicit none

      integer nq,ipol
      double precision fac,Ex
      logical ldew, lfac, nlfac
      double precision func(*)  ! value of the functional [output]

c--> Charge Density

      double precision rho(ipol*(ipol+1)/2,nq)

c--> Quadrature Weights

      double precision qwght(nq)

c--> Partial Derivatives of the Exchange Energy Functional

      double precision Amat(ipol,nq)
c
c Compute the partial derivatives of the exchange functional of Dirac.
c
      double precision P1,P2,P3,P4,tol_rho
c--> P1 =       -(3/PI)**(1/3)
c--> P2 = -(3/4)*(3/PI)**(1/3)
c--> P3 =       -(6/PI)**(1/3)
c--> P4 = -(3/4)*(6/PI)**(1/3)

      Parameter (P1 = -0.9847450218426959D+00)
      Parameter (P2 = -0.7385587663820219D+00)
      Parameter (P3 = -0.1240700981798799D+01)
      Parameter (P4 = -0.9305257363490993D+00)
      double precision rho13,rho32,rho33,Exp
      integer n
c***************************************************************************
        Exp = 0.d0

      if( ipol.eq.1 )then

c                      ======> SPIN-RESTRICTED <======
       do 10 n = 1,nq
         if(abs(rho(1,n)).gt.tol_rho) then
           rho13=rho(1,n)**(1.d0/3.d0)
           Amat(1,n) = Amat(1,n) + rho13
           Exp      = Exp + rho(1,n)*rho13*qwght(n)
           if(ldew) func(n) = func(n) + rho(1,n)*rho13*fac*P2
         endif
   10  continue
        call dscal(nq,P1*fac,Amat(1,1),1)
        Ex = Ex + P2*Exp*fac

      else

c                      ======> SPIN-UNRESTRICTED <======

       do 20 n = 1,nq
         if(abs(rho(1,n)).lt.tol_rho) goto 20
        rho32=rho(2,n)**(1.d0/3.d0)
        rho33=rho(3,n)**(1.d0/3.d0)
        Amat(1,n) = Amat(1,n) + P3*rho32*fac
        Amat(2,n) = Amat(2,n) + P3*rho33*fac
        if(ldew) func(n) = 
     &       func(n) + (rho32*rho(2,n)+rho33*rho(3,n))*P4*fac
        Exp      = Exp + (rho32*rho(2,n)+rho33*rho(3,n))*qwght(n)
   20  continue
        Ex = Ex + P4*Exp*fac

      end if
      return
      end
