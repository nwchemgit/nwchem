      Subroutine xc_dirac(rho,Amat,nq,ipol,Ex,qwght)

C$Id: xc_dirac.F,v 1.1 1995-04-26 03:40:11 og845 Exp $
      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

c--> P1 =       -(3/PI)**(1/3)
c--> P2 = -(3/4)*(3/PI)**(1/3)
c--> P3 =       -(6/PI)**(1/3)
c--> P4 = -(3/4)*(6/PI)**(1/3)

      Parameter (P1 = -0.9847450218426959D+00)
      Parameter (P2 = -0.7385587663820219D+00)
      Parameter (P3 = -0.1240700981798799D+01)
      Parameter (P4 = -0.9305257363490993D+00)
      parameter(TOLL=1.D-15)

c--> Charge Density

      Dimension rho(ipol*(ipol+1)/2,nq)

c--> Quadrature Weights

      Dimension qwght(nq)

c--> Partial Derivatives of the Exchange Energy Functional

      Dimension Amat(ipol,nq)
c
c Compute the partial derivatives of the exchange functional of Dirac.
c
c***************************************************************************
        Exp = 0.d0

      if( ipol.eq.1 )then

c                      ======> SPIN-RESTRICTED <======
       do 10 n = 1,nq
         if(rho(1,n).gt.TOLL) then
           rho13=rho(1,n)**(1.d0/3.d0)
           Amat(1,n) = Amat(1,n) + rho13
           Exp      = Exp + rho(1,n)*rho13*qwght(n)
         endif
   10  continue
        call dscal(nq,P1,Amat(1,1),1)
        Ex = Ex + P2*Exp

      else

c                      ======> SPIN-UNRESTRICTED <======

       do 20 n = 1,nq
         if(dabs(rho(1,n)).lt.TOLL) goto 20
        rho32=rho(2,n)**(1.d0/3.d0)
        rho33=rho(3,n)**(1.d0/3.d0)
        Amat(1,n) = Amat(1,n) + P3*rho32
        Amat(2,n) = Amat(2,n) + P3*rho33
        Exp      = Exp + (rho32*rho(2,n)+rho33*rho(3,n))*qwght(n)
   20  continue
        Ex = Ex + P4*Exp

      end if
      return
      end
