      Subroutine xc_lyp88(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                    Amat, Cmat, nq, ipol, Ec, qwght, ldew, func)
c
C$Id: xc_lyp88.F,v 1.13 2000-04-26 17:43:37 jolguin Exp $
c
      implicit none 
c
#include "dft2drv.fh"
c
      double precision fac ! [input]
      integer nq
      integer ipol
      double precision Ec
      logical lfac, nlfac, ldew
      double precision func(*)  ! value of the functional [output]
c
c     Charge Density & Its Cube Root
c
      double precision rho(nq,(ipol*(ipol+1))/2)
c
c     Charge Density Gradient
c
      double precision delrho(nq,3,ipol)
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(nq,ipol), Cmat(nq,*)
      double precision P13, P43, P83, P19, EPS, P1, tol_rho, A, 
     &                 B, C, D, P2,  P3
c
      Parameter (P13 = 1.D0/3.D0, P43 = 4.D0*P13, P83 = 8.D0*P13, 
     &           P19 = 1.D0/9.D0)
      Parameter (EPS = 1.D-32)
c
c     P1 = 2**(11/3)*(3/10)*(3*PI**2)**(2/3)
c
      Parameter (P1 = 0.3646239897876487D+02)
c
c     Colle-Salvetti Empirical Parameters
c
      Parameter (A = 0.04918D0)
      Parameter (B = 0.13200D0)
      Parameter (C = 0.25330D0)
      Parameter (D = 0.34900D0)
      Parameter (P2 = A*B/D,  P3 = 4.D0*A/D)
c
c     Compute the partial derivatives of the correlation functional of 
c     Lee, Yang and Parr.
c
c     References:
c
c     Colle & Salvetti, Theor. Chim. Acta 37, 329 (1975)
c     Lee, Yang & Parr, Phys. Rev. B 37, 785 (1988)
c     Miehlich, Savin, Stoll & Preuss, Chem. Phys. Lett. 157, 200 (1989)
c     Johnson, Gill & Pople, J. Chem. Phys. 98, 5612 (1993)
c
      integer n
      double precision Cd, CDa, Cdb, CD13, CD13a, CD13b, CD43, CD83a, 
     &                 CD83b
      double precision gaa, Gbb, Gab, Ra, Rb, Rab, Raa, Rbb, Rab9
      double precision xc, xd, Rx, delta, omega
      double precision Faa, Fbb, Fab, dElypdGaa, dElypdGbb, dElypdGab
      double precision T1, T2a, T2b, T2, Elyp
      double precision deltap, Ro, dFaadCDa, dFbbdCDa, dFabdCDa, 
     &                 d2ElypdCDadGaa, d2ElypdCDadGbb, d2ElypdCDadGab, 
     &                 dT1dCDa, dT2dCDa,  dElypdCDa
      double precision zeta
      double precision dFaadCDb,  dFbbdCDb,  dFabdCDb, d2ElypdCDbdGaa,  
     &                 d2ElypdCDbdGbb,  d2ElypdCDbdGab, dT1dCDb, 
     &                 dT2dCDb, dElypdCDb
c
      if (ipol.eq.1)then
c
c        ======> SPIN-RESTRICTED <======
c
         do 10 n = 1, nq
            if (abs(rho(n,1)).lt.tol_rho) goto 10
            CD  = rho(n,1)
            CDa = 0.5D0*CD
            CDb = CDa
            CD13  = rho(n,1)**(1.D0/3.D0)
            CD13a = (0.5D0**P13)*CD13
            CD13b = CD13a
            Gaa = 0.25D0*(delrho(n,1,1)*delrho(n,1,1) +
     &                    delrho(n,2,1)*delrho(n,2,1) +
     &                    delrho(n,3,1)*delrho(n,3,1))
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
            if(abs(xc).lt.eps) then
               omega=0d0
            else
               omega = P2*Rx*exp(-xc)/CD43
            endif
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
            Elyp = T1 + T2 + dElypdGaa*Gaa + 
     &                       dElypdGbb*Gbb + 
     &                       dElypdGab*Gab
            deltap = P13*(Rx*Rx - delta)/CD
            Ro = P13*(delta - 5.D0)/CD
            dFaadCDa = -Rab9*((3.D0 + Ra)*deltap + 
     &                        (delta - 11.D0)*Rb/CD)
            dFbbdCDa = -Rab9*((3.D0 + Rb)*deltap - 
     &                        (delta - 11.D0)*Rb/CD)
            dFabdCDa = -7.D0*Rab9*deltap
            d2ElypdCDadGaa = Ro*dElypdGaa - 
     &                       omega*(dFaadCDa + Rbb*(2.D0/CD))
            d2ElypdCDadGbb = Ro*dElypdGbb - 
     &                       omega*(dFbbdCDa - Rab*(2.D0/CD))
            d2ElypdCDadGab = Ro*dElypdGab - omega* dFabdCDa
            dT1dCDa = (T1/CD)*(P13*Rx + 1.D0)
            dT2dCDa = Ro*T2b + (Ro + (P83/CDa))*T2a
            dElypdCDa = dT1dCDa + dT2dCDa + d2ElypdCDadGaa*Gaa +
     &                                      d2ElypdCDadGbb*Gbb +
     &                                      d2ElypdCDadGab*Gab
c
            Amat(n,1) = Amat(n,1) + dElypdCDa*fac
            gaa = (dElypdGaa + 0.5d0*dElypdGab)*fac
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + dElypdGaa*fac
            Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + dElypdGab*fac
#if 0
            Cmat(n,1,1) = Cmat(n,1,1) + gaa*delrho(n,1,1)
            Cmat(n,2,1) = Cmat(n,2,1) + gaa*delrho(n,2,1)
            Cmat(n,3,1) = Cmat(n,3,1) + gaa*delrho(n,3,1)
#endif
            Ec = Ec + (Elyp*fac)*qwght(n)
            if (ldew) func(n) =func(n) + (Elyp*fac)
   10    continue
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1,nq
            if (abs(rho(n,1)).lt.tol_rho)goto 20
            CD  = rho(n,1)
            CDa = rho(n,2)
            CDb = rho(n,3)
            CD13  = rho(n,1)**(1.d0/3.D0)
            CD13a = rho(n,2)**(1.d0/3.D0)
            Gaa = delrho(n,1,1)*delrho(n,1,1) +
     &            delrho(n,2,1)*delrho(n,2,1) +
     &            delrho(n,3,1)*delrho(n,3,1)
            CD43  =  (CD13 **2)**2
            CD83a = ((CD13a**2)**2)**2
            Ra = CDa/CD
            Raa = Ra*Ra
            xc = C/CD13
!            xd = D/CD13
!            Rx = xd/(1.D0 + xd)
            rx=d/(d+cd13)
            delta = xc + Rx
            if(abs(xc).lt.eps) then
               omega=0d0
            else
               omega = P2*Rx*exp(-xc)/CD43
            endif
            if(cdb.gt.eps) then
               CD13b = rho(n,3)**(1.d0/3.D0)
               Gbb = delrho(n,1,2)*delrho(n,1,2) +
     &              delrho(n,2,2)*delrho(n,2,2) +
     &              delrho(n,3,2)*delrho(n,3,2)
               Gab = delrho(n,1,1)*delrho(n,1,2) +
     &              delrho(n,2,1)*delrho(n,2,2) +
     &              delrho(n,3,1)*delrho(n,3,2)
               CD83b = ((CD13b**2)**2)**2
               Rb = CDb/CD
               Rab = Ra*Rb
               Rbb = Rb*Rb
               Rab9 = P19*Rab
               zeta = Ra - Rb
               Faa = Rab9*(1.D0 - 3.D0*delta - (delta - 11.D0)*Ra)
               Fbb = Rab9*(1.D0 - 3.D0*delta - (delta - 11.D0)*Rb)
               Fab = Rab9*(47.D0 - 7.D0*delta)
               T1 = -P3*Rx*Rab*CD43
               T2a = -(P1*omega*Rab)*CD83a
               T2b = -(P1*omega*Rab)*CD83b
            else
               cd13b=0d0
               gbb=0d0
               gab=0d0
               cd83b=0d0
               rb=0d0
               rab=0d0
               rbb=0d0
               rab9=0d0
               zeta=ra
               faa=0d0
               fbb=0d0
               fab=0d0
               t1=0d0
               t2a=0d0
               t2b=0d0
            endif
            dElypdGaa = -omega*(Faa - Rbb)
            dElypdGbb = -omega*(Fbb - Raa)
            dElypdGab = -omega*(Fab - P43)
            T2 = T2a + T2b
            Elyp = T1 + T2 + dElypdGaa*Gaa + 
     &                       dElypdGbb*Gbb + 
     &                       dElypdGab*Gab
            deltap = P13*(Rx**2 - delta)/CD
            Ro = P13*(delta - 5.D0)/CD
            dFaadCDa = -(zeta/CDa)*Faa - Rab9*((3.D0 + Ra)*deltap + 
     &                  (delta - 11.D0)*Rb/CD)
            if(cdb.gt.eps) then
               dFaadCDb =  (zeta/CDb)*Faa - Rab9*((3.D0 + Ra)*deltap - 
     &                  (delta - 11.D0)*Ra/CD)
               dFbbdCDa=-(zeta/CDa)*Fbb - Rab9*((3.D0 + Rb)*deltap - 
     &              (delta - 11.D0)*Rb/CD)
               dFbbdCDb= (zeta/CDb)*Fbb - Rab9*((3.D0 + Rb)*deltap + 
     &              (delta - 11.D0)*Ra/CD)
               dFabdCDb =  (zeta/CDb)*Fab - 7.D0*Rab9*deltap
               dFabdCDa = -(zeta/CDa)*Fab - 7.D0*Rab9*deltap
            else
               dFaadCDb =0d0
               dFbbdCDa=0d0
               dFbbdCDb=  0d0
               dFabdCDb =0d0
               dFabdCDa = 0d0
            endif
            d2ElypdCDadGaa = Ro*dElypdGaa - 
     &                       omega*(dFaadCDa + Rbb*(2.D0/CD))
            d2ElypdCDbdGaa = Ro*dElypdGaa - 
     &                       omega*(dFaadCDb - Rab*(2.D0/CD))
            d2ElypdCDadGbb = Ro*dElypdGbb - 
     &                       omega*(dFbbdCDa - Rab*(2.D0/CD))
            d2ElypdCDbdGbb = Ro*dElypdGbb - 
     &                       omega*(dFbbdCDb + Raa*(2.D0/CD))
            d2ElypdCDadGab = Ro*dElypdGab - omega* dFabdCDa
            d2ElypdCDbdGab = Ro*dElypdGab - omega* dFabdCDb
            dT1dCDa = (T1/CD)*(P13*Rx + Rb/Ra)
            dT2dCDa = (Ro - (zeta/CDa))*T2b + 
     &                (Ro + (P83 - zeta)/CDa)*T2a
            if(cdb.gt.eps) then
               dT1dCDb = (T1/CD)*(P13*Rx + Ra/Rb)
               dT2dCDb = (Ro + (zeta/CDb))*T2a + 
     &              (Ro + (P83 + zeta)/CDb)*T2b
            else
               dT1dCDb = 0d0
               dT2dCDb = 0d0
            endif
            dElypdCDa = dT1dCDa + dT2dCDa + d2ElypdCDadGaa*Gaa +
     &                                      d2ElypdCDadGbb*Gbb +
     &                                      d2ElypdCDadGab*Gab
            dElypdCDb = dT1dCDb + dT2dCDb + d2ElypdCDbdGaa*Gaa +
     &                                      d2ElypdCDbdGbb*Gbb +
     &                                      d2ElypdCDbdGab*Gab
c
            Amat(n,1) = Amat(n,1) + dElypdCDa*fac
            Amat(n,2) = Amat(n,2) + dElypdCDb*fac
            gaa = 2.d0*dElypdGaa*fac
            gbb = 2.d0*dElypdGbb*fac
            gab = dElypdGab*fac
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + dElypdGaa*fac
            Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + dElypdGab*fac
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + dElypdGbb*fac
#if 0
            Cmat(n,1,1) = Cmat(n,1,1) + gaa*delrho(n,1,1) + 
     &                                  gab*delrho(n,1,2)
            Cmat(n,1,2) = Cmat(n,1,2) + gbb*delrho(n,1,2) + 
     &                                  gab*delrho(n,1,1)
            Cmat(n,2,1) = Cmat(n,2,1) + gaa*delrho(n,2,1) + 
     &                                  gab*delrho(n,2,2)
            Cmat(n,2,2) = Cmat(n,2,2) + gbb*delrho(n,2,2) + 
     &                                  gab*delrho(n,2,1)
            Cmat(n,3,1) = Cmat(n,3,1) + gaa*delrho(n,3,1) + 
     &                                  gab*delrho(n,3,2)
            Cmat(n,3,2) = Cmat(n,3,2) + gbb*delrho(n,3,2) + 
     &                                  gab*delrho(n,3,1)
#endif
            Ec = Ec + (Elyp*fac)*qwght(n)
            if (ldew) func(n) =func(n) + (Elyp*fac)
   20    continue
      endif
      return
      end
