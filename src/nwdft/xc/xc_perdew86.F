      Subroutine xc_perdew86(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                       Amat, Cmat, nq, ipol, Ec, qwght,
     ,     ldew, ffunc)
c
c$Id: xc_perdew86.F,v 1.16 1999-01-29 01:04:42 edo Exp $
c
      implicit none
c
      double precision tol_rho, fac ! [input]
      integer nq, ipol              ! [input]
      double precision Ec           ! [input/output]
      logical lfac, nlfac, ldew
      double precision ffunc(*)  ! value of the functional [output]
c
c     Charge Density 
c
      double precision rho(nq,ipol*(ipol+1)/2)
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
      double precision Amat(nq,ipol), Cmat(nq,3,ipol)
      double precision TOLL, EXPTOL, alpha, beta, pgamma, delta, 
     &                 beta10, ftilde, zzz, fff, pfff, CINF, ONE, 
     &                 ONETHIRD, THREE, FOURTHIRDS, SEV6, FIVE3, 
     &                 TWO3, FIVE6, pi
      Parameter (TOLL = 1.D-40, EXPTOL = 80.d0)
      Parameter (alpha = 0.023266D0, beta  =  7.389D-6, 
     &   pgamma = 8.723d0, delta = 0.472d0,  beta10 = 10000.d0*beta)
      parameter (ftilde = 0.11d0, zzz = 0.001667d0, fff = 0.002568d0)
      parameter(pfff = 1.745d0, CINF = zzz+fff)
      Parameter (ONE = 1.D0, ONETHIRD = 1.d0/3.d0, THREE = 3.d0)
      Parameter (FOURTHIRDS = 4.D0/3.D0, SEV6 = 7.d0/6.d0)
      parameter (FIVE3 = 5.d0/3.d0, TWO3 = 2.d0/3.d0, FIVE6 = 5.d0/6.d0)
      parameter (pi = 3.1415926535897932385d0)
c     
c     Mlynarski Salahub PRB 43, 1399 (1991)
c     
      integer n
      double precision rsfact, rs, rs2, rs3
      double precision rhoval, rho13, rho43, rho76, arho
      double precision drsdrho, gammaval, dsqgamma, dphidrho
      double precision anum, aden, CC, CCPRIME, expfac, phi, func
      double precision AAadd, g, zeta, add, adp, dzetada, dzetadb
      double precision gaa, gab, gbb
      double precision dzeta, dprime
      double precision x
      dzeta(x) = 2.d0**ONETHIRD*dsqrt(((ONE + x)*.5d0)**FIVE3 + 
     &                                ((ONE - x)*.5d0)**FIVE3   )
      dprime(x) = (4.d0**(-ONETHIRD)*FIVE6)*(((ONE + x)*.5d0)**TWO3-
     -              ((ONE - x)*.5d0)**TWO3 ) /
     &              dsqrt(((ONE + x)*.5d0)**FIVE3 +
     &                    ((ONE - x)*.5d0)**FIVE3)
c
      rsfact = (0.75d0/pi)**ONETHIRD
c
      if (ipol.eq.1 )then
c
c        ======> SPIN-RESTRICTED <======
c
         do 10 n = 1, nq
            rhoval = rho(n,1)
            if (abs(rhoval).lt.tol_rho) goto 10
            rho13  = rhoval**ONETHIRD
            rho43  = rhoval*rho13
            rho76  = rhoval**SEV6 
            rs = rsfact/rho13
            rs2 = rs*rs
            rs3 = rs2*rs
            drsdrho = -ONETHIRD*rs/rhoval
            gammaval = delrho(n,1,1)*delrho(n,1,1) +
     &                 delrho(n,2,1)*delrho(n,2,1) +
     &                 delrho(n,3,1)*delrho(n,3,1)
            dsqgamma = dsqrt(gammaval)
c         
c           C(n)
c         
            anum = fff+alpha*rs+beta*rs2
            aden = 1.d0/(1.d0+pgamma*rs+delta*rs2+beta10*rs3)
            CC = zzz+(anum*aden)
            CCPRIME = aden*((alpha+beta*(rs+rs)) -
     &                (anum*aden)*(pgamma+(rs+rs)*delta + 
     &                beta10*rs2*THREE))
            CCPRIME = CCPRIME*drsdrho
c         
c           phi(n,gradn)
c         
            expfac = 0.d0
            phi = (pfff*ftilde)*(CINF/CC)*dsqgamma/rho76
            if (phi.lt.EXPTOL)expfac = exp(-phi)
c         
c           functional
c         
            func = expfac*CC*gammaval/rho43*fac
            dphidrho =  -phi * ( SEV6/rhoval  +   (CCPRIME/CC)) 
            AAadd = ((CCPRIME/CC) - (dphidrho + 
     &              FOURTHIRDS/rhoval))*func
            Amat(n,1) = Amat(n,1) + AAadd*fac
            if (dsqgamma.gt.0.d0)then
               g = func/gammaval*(2.d0-phi)*fac
               Cmat(n,1,1) = Cmat(n,1,1) + g*delrho(n,1,1)
               Cmat(n,2,1) = Cmat(n,2,1) + g*delrho(n,2,1)
               Cmat(n,3,1) = Cmat(n,3,1) + g*delrho(n,3,1)
            endif
            Ec = Ec + func*qwght(n)*fac
            if (ldew) ffunc(n)=ffunc(n)+func*fac
   10    continue
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1, nq
            rhoval = rho(n,1)
            if (abs(rhoval).lt.tol_rho) goto 20
            arho=1.d0/rhoval
            rho13  = rhoval**ONETHIRD
            rho43  = rhoval*rho13
            rho76  = rhoval**SEV6 
            rs = rsfact/rho13
            rs2 = rs*rs
            rs3 = rs2*rs
            drsdrho = -ONETHIRD*rs*arho
            gammaval = delrho(n,1,1)*delrho(n,1,1) +
     &                 delrho(n,2,1)*delrho(n,2,1) +
     &                 delrho(n,3,1)*delrho(n,3,1) +
     &                 delrho(n,1,2)*delrho(n,1,2) +
     &                 delrho(n,2,2)*delrho(n,2,2) +
     &                 delrho(n,3,2)*delrho(n,3,2) +
     &           2.d0*(delrho(n,1,1)*delrho(n,1,2) +
     &                 delrho(n,2,1)*delrho(n,2,2) +
     &                 delrho(n,3,1)*delrho(n,3,2)   )
            dsqgamma = sqrt(gammaval)
            zeta = (rho(n,2) - rho(n,3))*arho
            add = 1.d0/dzeta(zeta)
            adp = dprime(zeta)
            dzetada =  (1.d0 - zeta)*arho*adp
            dzetadb = -(1.d0 + zeta)*arho*adp
c         
c           C(n)
c         
            anum = fff + alpha*rs + beta*rs2
            aden = 1.d0/(1.d0 + pgamma*rs + delta*rs2 + beta10*rs3)
            CC = zzz + (anum*aden)
            CCPRIME = aden*((alpha + beta*(rs + rs)) -
     &                (anum*aden)*(pgamma+(rs + rs)*delta + 
     &                beta10*rs2*THREE))
            CCPRIME = CCPRIME*drsdrho
c         
c           phi(n,gradn)
c         
            expfac = 0.d0
            phi = (pfff*ftilde)*(CINF/CC)*dsqgamma/rho76
            if (phi.lt.EXPTOL) expfac = exp(-phi)
c         
c           functional
c         
            func = expfac*CC*gammaval/rho43*add
            dphidrho =  -phi*(SEV6/rhoval + (CCPRIME/CC)) 
            AAadd = ((CCPRIME/CC) -
     &               (dphidrho + FOURTHIRDS/rhoval))*func
            Amat(n,1) = Amat(n,1) + (AAadd - (func*add)*dzetada)*fac
            Amat(n,2) = Amat(n,2) + (AAadd - (func*add)*dzetadb)*fac
            if (dsqgamma.gt.TOLL)then
               gaa = func/gammaval*(2.d0-phi)*fac
               gbb = gaa
               gab = gaa
               Cmat(n,1,1) = Cmat(n,1,1) + gaa*delrho(n,1,1) +
     &                                     gab*delrho(n,1,2)
               Cmat(n,2,1) = Cmat(n,2,1) + gaa*delrho(n,2,1) +
     &                                     gab*delrho(n,2,2)
               Cmat(n,3,1) = Cmat(n,3,1) + gaa*delrho(n,3,1) +
     &                                     gab*delrho(n,3,2)
               Cmat(n,1,2) = Cmat(n,1,2) + gbb*delrho(n,1,2) +
     &                                     gab*delrho(n,1,1)
               Cmat(n,2,2) = Cmat(n,2,2) + gbb*delrho(n,2,2) +
     &                                     gab*delrho(n,2,1)
               Cmat(n,3,2) = Cmat(n,3,2) + gbb*delrho(n,3,2) +
     &                                     gab*delrho(n,3,1)
            endif
            Ec = Ec + func*qwght(n)*fac
            if (ldew) ffunc(n)=ffunc(n)+func*fac
   20    continue
      endif
      return
      end

      Subroutine xc_p81(tol_rho, fac, lfac, nlfac, rho, Amat, nq, ipol, 
     &                  Ec, qwght, ldew, func)
c
c     Ceperley Alder LDA from Perdew Zunger PRB 23, 5048 (1981)
c
      implicit none
c
      integer nq, ipol
      logical lfac, nlfac, ldew
      double precision func(*)  ! value of the functional [output]
      double precision Ec, fac
c
c     Charge Density 
c
      double precision rho(nq,ipol*(ipol+1)/2)
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(nq,ipol)
      double precision AU, BU, CU, DU, GU, B1U, B2U, AP, BP, CP, DP, GP, 
     &                 B1P, B2P, pi, tol_rho, ONE3, SEV6, FOUR3, TWO3
      parameter(AU = 0.0311d0, BU = -0.048d0, CU = 0.0020d0, 
     &          DU = -0.0116d0, GU = -.1423d0, B1U = 1.0529d0, 
     &          B2U = 0.3334d0, AP = 0.01555d0, BP = -0.0269d0, 
     &          CP = 0.0007d0, DP = -0.0048d0, GP = -.0843d0, 
     &          B1P = 1.3981d0, B2P = 0.2611d0)
      parameter (pi = 3.1415926535897932385d0)
      Parameter (ONE3 = 1.d0/3.d0, SEV6 = 7.d0/6.d0)
      Parameter (FOUR3 = 4.D0/3.D0, TWO3 = 2.d0/3.d0)
      integer n
      double precision rhoval, rs, alnrs, eU, vU, drs, adU, ePol, 
     &                 fz, zeta, eP, vP, adP, AAadd
      double precision x, fzeta, fzzprime, rsfact
      fzeta(x) = ((1.d0+x)**FOUR3 + 
     &            (1.d0-x)**FOUR3 - 2.d0) / (2.d0**FOUR3-2.d0)
      fzzprime(x) = FOUR3*((1.d0+x)**ONE3 - 
     &                     (1.d0-x)**ONE3) / (2.d0**FOUR3-2.d0)
c
      rsfact = (0.75d0/pi)**ONE3
c
      if (ipol.eq.1 )then
c
c        ======> SPIN-RESTRICTED <======
c
         do n = 1, nq
            if (abs(rho(n,1)).gt.tol_rho)then
               rhoval = rho(n,1)
               rs = rsfact/rhoval**ONE3
               if (rs.lt.1.d0)then
                  alnrs = log(rs)
                  eU = AU*alnrs+BU+CU*rs*alnrs+DU*rs
                  vU = AU*alnrs + (BU-ONE3*AU) + TWO3*CU*rs*alnrs + 
     &                 ONE3*(DU+DU-CU)*rs
               else
                  drs = sqrt(rs)
                  adU = 1.d0/(1.d0+B1U*drs+B2U*rs)
                  eU = GU*adU
                  vU = eU*(1.d0+SEV6*B1U*drs+FOUR3*B2U*rs)*adU
               endif
               Ec = Ec+eU*qwght(n)*rhoval*fac
               if (ldew) func(n) = func(n) +eU*rhoval*fac
               Amat(n,1) = Amat(n,1) + vU*fac
            endif
         enddo
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do n = 1, nq
            if (abs(rho(n,1)).gt.tol_rho)then
               rhoval = rho(n,1)
               zeta = (rho(n,2)-rho(n,3))/rhoval
               fz = fzeta(zeta)
               rs = rsfact/rhoval**ONE3
               if (rs.lt.1.d0)then
                  alnrs = log(rs)
                  eU = AU*alnrs+BU+CU*rs*alnrs+DU*rs
                  eP = AP*alnrs+BP+CP*rs*alnrs+DP*rs
                  vU = AU*alnrs + (BU-ONE3*AU) + TWO3*CU*rs*alnrs + 
     &                 ONE3*(DU+DU-CU)*rs
                  vP = AP*alnrs + (BP-ONE3*AP) + TWO3*CP*rs*alnrs + 
     &                 ONE3*(DP+DP-CP)*rs
               else
                  drs = sqrt(rs)
                  adU = 1.d0/(1.d0+B1U*drs+B2U*rs)
                  adP = 1.d0/(1.d0+B1P*drs+B2P*rs)
                  eU = GU*adU
                  eP = GP*adP
                  vU = eU*(1.d0+SEV6*B1U*drs+FOUR3*B2U*rs)*adU
                  vP = eP*(1.d0+SEV6*B1P*drs+FOUR3*B2P*rs)*adP
               endif
               AAadd = vU + fz*(vP-vU) 
               ePOL = (eP-eU)*fzzprime(zeta)
               Ec = Ec+(eU+fz*(eP-eU))*qwght(n)*rhoval*fac
               if (ldew) func(n) = func(n) +
     *              (eU+fz*(eP-eU))*rhoval*fac
               Amat(n,1) = Amat(n,1) + (AAadd + (1.d0-zeta)*ePOL)*fac
               Amat(n,2) = Amat(n,2) + (AAadd - (1.d0+zeta)*ePOL)*fac
            endif
         enddo
      endif
      return
      end

