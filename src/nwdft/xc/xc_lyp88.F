      Subroutine xc_lyp88(fac,
     *     rho,delrho,Amat,Cmat,nq,ipol,Ec,qwght)

C$Id: xc_lyp88.F,v 1.4 1996-05-04 21:08:35 og845 Exp $
      implicit none 

      double precision fac ! [input]
      integer nq
      integer ipol
      double precision Ec
c--> Charge Density & Its Cube Root

      double precision rho(ipol*(ipol+1)/2,nq)

c--> Charge Density Gradient

      double precision delrho(ipol,3,nq)

c--> Quadrature Weights

      double precision qwght(nq)

c--> Sampling Matrices for the XC Potential & Energy

      double precision Amat(ipol,nq),Cmat(ipol,3,nq)
c
      double precision P13,P43,P83,P19,EPS,
     ,     P1,TOLL,A,B,C,D,
     ,     P2, P3

      Parameter (P13=1.D0/3.D0,P43=4.D0*P13,P83=8.D0*P13,P19=1.D0/9.D0)
      Parameter (EPS=1.D-32)

c--> P1 = 2**(11/3)*(3/10)*(3*PI**2)**(2/3)

      Parameter (P1=0.3646239897876487D+02)
      Parameter (TOLL=1.d-20)

c--> Colle-Salvetti Empirical Parameters

      Parameter (A = 0.04918D0)
      Parameter (B = 0.13200D0)
      Parameter (C = 0.25330D0)
      Parameter (D = 0.34900D0)

      Parameter (P2 = A*B/D, P3 = 4.D0*A/D)
c Compute the partial derivatives of the correlation functional of Lee, Yang
c and Parr.
c
c References:
c
c    Colle & Salvetti, Theor. Chim. Acta 37, 329 (1975)
c    Lee, Yang & Parr, Phys. Rev. B 37, 785 (1988)
c    Miehlich, Savin, Stoll & Preuss, Chem. Phys. Lett. 157, 200 (1989)
c    Johnson, Gill & Pople, J. Chem. Phys. 98, 5612 (1993)
c
      integer n
      double precision Cd,CDa,Cdb,CD13,CD13a,CD13b,CD43,CD83a,CD83b
      double precision gaa,Gbb,Gab,Ra,Rb,Rab,Raa,Rbb,Rab9
      double precision xc,xd,Rx,delta,omega
      double precision Faa,Fbb,Fab,dElypdGaa,dElypdGbb,dElypdGab
      double precision T1,T2a,T2b,T2,Elyp
      double precision deltap,Ro,dFaadCDa,dFbbdCDa,dFabdCDa,
     ,     d2ElypdCDadGaa,d2ElypdCDadGbb,d2ElypdCDadGab,dT1dCDa,
     ,  dT2dCDa, dElypdCDa
      double precision zeta
      double precision  dFaadCDb, dFbbdCDb, dFabdCDb,
     , d2ElypdCDbdGaa, d2ElypdCDbdGbb, d2ElypdCDbdGab,
     , dT1dCDb, dT2dCDb, dElypdCDb

c***************************************************************************

      if( ipol.eq.1 )then

c                      ======> SPIN-RESTRICTED <======

       do 10 n = 1,nq

        if(dabs(rho(1,n)).lt.TOLL) goto 10
        CD  = rho(1,n)
        CDa = 0.5D0*CD
        CDb = CDa

        CD13  = rho(1,n)**(1.D0/3.D0)
        CD13a = (0.5D0**P13)*CD13
        CD13b = CD13a

        Gaa = 0.25D0*(delrho(1,1,n)*delrho(1,1,n)+
     +       delrho(1,2,n)*delrho(1,2,n)+
     +       delrho(1,3,n)*delrho(1,3,n))
        Gbb = Gaa
        Gab = Gaa

        CD43  =  (CD13 **2)**2
        CD83a = ((CD13a**2)**2)**2
        CD83b =   CD83a

        Ra = 0.5D0
        Rb = 0.5D0

        Rab = 0.25D0
        Raa = 0.25D0
        Rbb = 0.25D0

        xc = C/CD13
        xd = D/CD13

        Rx = xd/(1.D0 + xd)

        delta = xc + Rx
        omega = P2*Rx*max(exp(-xc),EPS)/CD43

        Rab9 = P19*Rab

        Faa = Rab9*(1.D0 - 3.D0*delta - (delta - 11.D0)*Ra)
        Fbb = Faa
        Fab = Rab9*(47.D0 - 7.D0*delta)

        dElypdGaa = -omega*(Faa - Rbb)
        dElypdGbb = dElypdGaa
        dElypdGab = -omega*(Fab - P43)

        T1 = -P3*Rx*Rab*CD43

        T2a = -(P1*omega*Rab)*CD83a
        T2b = T2a

        T2 = T2a + T2b

        Elyp = T1 + T2 + dElypdGaa*Gaa + dElypdGbb*Gbb + dElypdGab*Gab

        deltap = P13*(Rx*Rx - delta)/CD

        Ro = P13*(delta - 5.D0)/CD

        dFaadCDa = -Rab9*( (3.D0 + Ra)*deltap + (delta - 11.D0)*Rb/CD )
        dFbbdCDa = -Rab9*( (3.D0 + Rb)*deltap - (delta - 11.D0)*Rb/CD )

        dFabdCDa = -7.D0*Rab9*deltap

        d2ElypdCDadGaa = Ro*dElypdGaa - omega*(dFaadCDa + Rbb*(2.D0/CD))

        d2ElypdCDadGbb = Ro*dElypdGbb - omega*(dFbbdCDa - Rab*(2.D0/CD))

        d2ElypdCDadGab = Ro*dElypdGab - omega* dFabdCDa

        dT1dCDa = (T1/CD)*(P13*Rx + 1.D0)

        dT2dCDa = Ro*T2b + (Ro + (P83/CDa))*T2a

        dElypdCDa = dT1dCDa + dT2dCDa + d2ElypdCDadGaa*Gaa +
     &                                  d2ElypdCDadGbb*Gbb +
     &                                  d2ElypdCDadGab*Gab

c ...

        Amat(1,n) = Amat(1,n) + dElypdCDa*fac
        gaa=(dElypdGaa + 0.5d0*dElypdGab)*fac
        if(dabs(gaa).gt.TOLL) then
          Cmat(1,1,n)=Cmat(1,1,n)+gaa*delrho(1,1,n)
          Cmat(1,2,n)=Cmat(1,2,n)+gaa*delrho(1,2,n)
          Cmat(1,3,n)=Cmat(1,3,n)+gaa*delrho(1,3,n)
        endif
        Ec      =  Ec + Elyp*qwght(n)*fac

   10  continue

      else

c                      ======> SPIN-UNRESTRICTED <======

       do 20 n = 1,nq
         if(dabs(rho(1,n)).lt.TOLL) goto 20
        CD  = rho(1,n)
        CDa = rho(2,n)
        CDb = rho(3,n)

        CD13  = rho(1,n)**(1.d0/3.D0)
        CD13a = rho(2,n)**(1.d0/3.D0)
        CD13b = rho(3,n)**(1.d0/3.D0)

        Gaa = delrho(1,1,n)*delrho(1,1,n)+
     +       delrho(1,2,n)*delrho(1,2,n)+
     +       delrho(1,3,n)*delrho(1,3,n)
        Gbb = delrho(2,1,n)*delrho(2,1,n)+
     +       delrho(2,2,n)*delrho(2,2,n)+
     +       delrho(2,3,n)*delrho(2,3,n)
        Gab = delrho(1,1,n)*delrho(2,1,n)+
     +       delrho(1,2,n)*delrho(2,2,n)+
     +       delrho(1,3,n)*delrho(2,3,n)

        CD43  =  (CD13 **2)**2
        CD83a = ((CD13a**2)**2)**2
        CD83b = ((CD13b**2)**2)**2

        Ra = CDa/CD
        Rb = CDb/CD

        Rab = Ra*Rb
        Raa = Ra*Ra
        Rbb = Rb*Rb

        zeta = Ra - Rb

        xc = C/CD13
        xd = D/CD13

        Rx = xd/(1.D0 + xd)

        delta = xc + Rx
        omega = P2*Rx*max(exp(-xc),EPS)/CD43

        Rab9 = P19*Rab

        Faa = Rab9*(1.D0 - 3.D0*delta - (delta - 11.D0)*Ra)
        Fbb = Rab9*(1.D0 - 3.D0*delta - (delta - 11.D0)*Rb)
        Fab = Rab9*(47.D0 - 7.D0*delta)

        dElypdGaa = -omega*(Faa - Rbb)
        dElypdGbb = -omega*(Fbb - Raa)
        dElypdGab = -omega*(Fab - P43)

        T1 = -P3*Rx*Rab*CD43

        T2a = -(P1*omega*Rab)*CD83a
        T2b = -(P1*omega*Rab)*CD83b

        T2 = T2a + T2b

        Elyp = T1 + T2 + dElypdGaa*Gaa + dElypdGbb*Gbb + dElypdGab*Gab

        deltap = P13*(Rx**2 - delta)/CD

        Ro = P13*(delta - 5.D0)/CD

        dFaadCDa = -(zeta/CDa)*Faa -
     &               Rab9*( (3.D0 + Ra)*deltap + (delta - 11.D0)*Rb/CD )
        dFaadCDb =  (zeta/CDb)*Faa - 
     &               Rab9*( (3.D0 + Ra)*deltap - (delta - 11.D0)*Ra/CD )

        dFbbdCDa = -(zeta/CDa)*Fbb - 
     &               Rab9*( (3.D0 + Rb)*deltap - (delta - 11.D0)*Rb/CD )
        dFbbdCDb =  (zeta/CDb)*Fbb - 
     &               Rab9*( (3.D0 + Rb)*deltap + (delta - 11.D0)*Ra/CD )

        dFabdCDa = -(zeta/CDa)*Fab - 7.D0*Rab9*deltap
        dFabdCDb =  (zeta/CDb)*Fab - 7.D0*Rab9*deltap

        d2ElypdCDadGaa = Ro*dElypdGaa - omega*(dFaadCDa + Rbb*(2.D0/CD))
        d2ElypdCDbdGaa = Ro*dElypdGaa - omega*(dFaadCDb - Rab*(2.D0/CD))

        d2ElypdCDadGbb = Ro*dElypdGbb - omega*(dFbbdCDa - Rab*(2.D0/CD))
        d2ElypdCDbdGbb = Ro*dElypdGbb - omega*(dFbbdCDb + Raa*(2.D0/CD))

        d2ElypdCDadGab = Ro*dElypdGab - omega* dFabdCDa
        d2ElypdCDbdGab = Ro*dElypdGab - omega* dFabdCDb

        dT1dCDa = (T1/CD)*(P13*Rx + Rb/Ra)
        dT1dCDb = (T1/CD)*(P13*Rx + Ra/Rb)

        dT2dCDa = (Ro - (zeta/CDa))*T2b + (Ro + (P83 - zeta)/CDa)*T2a
        dT2dCDb = (Ro + (zeta/CDb))*T2a + (Ro + (P83 + zeta)/CDb)*T2b

        dElypdCDa = dT1dCDa + dT2dCDa + d2ElypdCDadGaa*Gaa +
     &                                  d2ElypdCDadGbb*Gbb +
     &                                  d2ElypdCDadGab*Gab

        dElypdCDb = dT1dCDb + dT2dCDb + d2ElypdCDbdGaa*Gaa +
     &                                  d2ElypdCDbdGbb*Gbb +
     &                                  d2ElypdCDbdGab*Gab

c ...

        Amat(1,n) = Amat(1,n) + dElypdCDa*fac
        Amat(2,n) = Amat(2,n) + dElypdCDb*fac
        gaa = 2.d0*dElypdGaa*fac
        gbb = 2.d0*dElypdGbb*fac
        gab = dElypdGab*fac
        Cmat(1,1,n)=Cmat(1,1,n)+gaa*delrho(1,1,n)+gab*delrho(2,1,n)
        Cmat(2,1,n)=Cmat(2,1,n)+gbb*delrho(2,1,n)+gab*delrho(1,1,n)
        Cmat(1,2,n)=Cmat(1,2,n)+gaa*delrho(1,2,n)+gab*delrho(2,2,n)
        Cmat(2,2,n)=Cmat(2,2,n)+gbb*delrho(2,2,n)+gab*delrho(1,2,n)
        Cmat(1,3,n)=Cmat(1,3,n)+gaa*delrho(1,3,n)+gab*delrho(2,3,n)
        Cmat(2,3,n)=Cmat(2,3,n)+gbb*delrho(2,3,n)+gab*delrho(1,3,n)

        Ec      =  Ec + Elyp*qwght(n)*fac

   20  continue

      end if

      return
      end
