      Subroutine xc_csogga(tol_rho, cfac, lcfac, nlcfac, rho, delrho, 
     &                     nq, ipol, Ec, qwght, ldew, func,
     &                     Amat, Cmat,ijzy)     

c
c$Id$
c
c
c**********************************************************************c
c                                                                      c
c  xc_csogga evaluates the corelation part of the SOGGA11 and          c
c  SOGGA11-X functionals on the grid.                                  c
c                                                                      c
c     a) Peverati, Zhao and Truhlar, J.Phys.Chem.Lett, 2, 1991 (2011)  c
c     b) Peverati and Truhlar, J.Chem.Phys, 135, 191102 (2011)         c
c                                                                      c
c      ijzy = 1 - SOGGA11 functional (a)                               c
c      ijzy = 2 - SOGGA11-X functional (b)                             c
c                                                                      c
c Coded by Roberto Peverati (12/11)                                    c
c                                                                      c
c**********************************************************************c        
c
      implicit none
c
#include "errquit.fh"
#include "dft2drv.fh"
c
c     Input and other parameters
c
      integer ipol, nq
      double precision dummy

      double precision cfac
      logical lcfac, nlcfac
      double precision tol_rho
c
c     Correlation energy
c
      double precision Ec
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
      logical ldew
      double precision func(*)
c
c     Sampling Matrices for the XC Potential
c
      double precision Amat(nq,ipol), Cmat(nq,*)
c
c     Intermediate derivative results, etc.
c
      integer n, ijzy
c
      double precision CcA(0:5),CcB(0:5)
      double precision BETA, CC0
      double precision DELOCDA, DELOCDB, DFA1DG, DFA1DR, DFA1DZ
      double precision DFA2DG, DFA2DR, DFA2DZ, DFA3DG, DFA3DR, DFA3DZ
      double precision DFA4DG, DFA4DR, DFA4DZ, DFA5DG, DFA5DR, DFA5DZ
      double precision DFB1DG, DFB1DR, DFB1DZ, DFB2DG, DFB2DR, DFB2DZ
      double precision DFB3DG, DFB3DR, DFB3DZ, DFB4DG, DFB4DR, DFB4DZ
      double precision DFB5DG, DFB5DR, DFB5DZ
      double precision DFEXPDPON, DFFRACDPON, DFGGACDA, DFGGACDB
      double precision DFGGACDG, DFGGACDGX, DFGGACDR, DFGGACDZ
      double precision DGDGX, DGDZ, DLDA, DLDB, DLDS, DLDZ
      double precision DPONDG, DPONDR, DPONDZ, DSDR
      double precision DT2DG, DT2DR, DT2DZ, DZDA, DZDB, ELOC
      double precision F1O3, F1O6, F7O6 
      double precision FA0, FA1, FA2, FA3, FA4, FA5
      double precision FB0, FB1, FB2, FB3, FB4, FB5
      double precision FEXP, FFRAC, FGGAC, FKFAC
      double precision GAMMAVAL, G, G2, G3, GRHO
      double precision PI, PI34
      double precision PON, POTLC
      double precision RHOA, RHOB, RHOVAL, RHO16, RHO76, RS
      double precision SCFAC, SKFAC, T, T2, XNU, ZETA

 
      double precision F1, F2, F3, F4, F5, F6, F7
      Save F1, F2, F3, F4, F5, F6, F7
      DATA F1/1.0D+00/,  F2/2.0D+00/,  F3/3.0D+00/,
     $     F4/4.0D+00/,  F5/5.0D+00/,  F6/6.0D+00/,
     $     F7/7.0D+00/



      
c
c     ======> BOTH SPIN-RESTRICETED AND UNRESTRICTED <======
c

      if (ijzy.eq.1) then
       CcA(0) =  0.50000D+00
       CcA(1) = -4.62334D+00
       CcA(2) =  8.00410D+00
       CcA(3) = -130.226D+00
       CcA(4) =  38.2685D+00
       CcA(5) =  69.5599D+00
       CcB(0) =  0.50000D+00
       CcB(1) =  3.62334D+00
       CcB(2) =  9.36393D+00
       CcB(3) =  34.5114D+00
       CcB(4) = -18.5684D+00
       CcB(5) = -0.16519D+00
      elseif (ijzy.eq.2) then
       CcA(0) =  5.00000d-01
       CcA(1) =  7.82439d+01
       CcA(2) =  2.57211d+01
       CcA(3) = -1.38830d+01
       CcA(4) = -9.87375d+00
       CcA(5) = -1.41357d+01
       CcB(0) =  5.00000d-01
       CcB(1) = -7.92439d+01
       CcB(2) =  1.63725d+01
       CcB(3) =  2.08129d+00
       CcB(4) =  7.50769d+00
       CcB(5) = -1.01861d+01
      else
       call errquit("xc_csogga: illegal value of ijzy",ijzy,UERR)
      endif

      Pi = ACos(-F1)
      Pi34 = F3/(F4*Pi)
      F1o3 = F1/F3
      F1o6 = F1/F6
      F7o6 = F7/F6
      XNu = 15.75592D0
      CC0 = 0.004235D0
      beta= XNu*CC0
      SCfac=F1      


c
c     ======> BOTH SPIN-RESTRICETED AND UNRESTRICTED <======
c
      do 20 n = 1, nq

      if (ipol.eq.1) then
        RHOA = rho(n,1)/F2
        RHOB = RHOA
        gammaval =(delrho(n,1,1)*delrho(n,1,1) +
     &             delrho(n,2,1)*delrho(n,2,1) +
     &             delrho(n,3,1)*delrho(n,3,1))
      else
        RhoA = rho(n,2)
        RhoB = rho(n,3)
        gammaval = delrho(n,1,1)*delrho(n,1,1) +
     &             delrho(n,1,2)*delrho(n,1,2) +
     &             delrho(n,2,1)*delrho(n,2,1) +
     &             delrho(n,2,2)*delrho(n,2,2) +
     &             delrho(n,3,1)*delrho(n,3,1) +
     &             delrho(n,3,2)*delrho(n,3,2) +
     &       2.d0*(delrho(n,1,1)*delrho(n,1,2) +
     &             delrho(n,2,1)*delrho(n,2,2) +
     &             delrho(n,3,1)*delrho(n,3,2))
      endif
      Rhoval  = RhoA + RhoB
      GRho = max(dsqrt(gammaval),tol_rho)
         if (dabs(rhoval).le.tol_rho) goto 20
         rS = (Pi34/Rhoval)**F1o3
         Zeta = (RhoA-RhoB)/Rhoval
         Call lsdac(tol_rho,
     $        RS,Zeta,PotLC,dLdS,dLdZ,dummy,dummy,dummy)
         Call GZeta(Zeta,G,dGdZ,dummy,dummy)
c
         Eloc = SCfac*Rhoval*PotLC
         FKFac = (F3*Pi*Pi)**F1o3
         SKFac = Sqrt(FKFac/Pi)*F2
         Rho16 = Rhoval**F1o6
         Rho76 = Rhoval*Rho16
         T  = GRho/(F2*SKFac*Rho76*G)
         T2 = T*T
         G2 = G*G
         G3 = G*G2
         PON = (G3*beta*T2)/PotLC
         Ffrac = F1-F1/(F1-PON)
         Fexp  = F1-exp(PON)
         fa0 = CcA(0)
         fa1 = CcA(1) *Ffrac
         fa2 = CcA(2) *Ffrac**F2
         fa3 = CcA(3) *Ffrac**F3
         fa4 = CcA(4) *Ffrac**F4
         fa5 = CcA(5) *Ffrac**F5
         fb0 = CcB(0)
         fb1 = CcB(1) *Fexp
         fb2 = CcB(2) *Fexp**F2
         fb3 = CcB(3) *Fexp**F3
         fb4 = CcB(4) *Fexp**F4
         fb5 = CcB(5) *Fexp**F5
         FggaC = fa0+fa1+fa2+fa3+fa4+fa5+
     $           fb0+fb1+fb2+fb3+fb4+fb5
C
C     1st derivatives.
C
         dSdR  = -(F1o3*rS/Rhoval)
         dZdA = (F1-Zeta)/Rhoval
         dZdB = (-F1-Zeta)/Rhoval
         dLdA = dLdS*dSdR + dLdz*dZdA
         dLdB = dLdS*dSdR + dLdz*dZdB
         dElocdA = SCfac*PotLC + SCfac*Rhoval*dLdA
         dElocdB = SCfac*PotLC + SCfac*Rhoval*dLdB
         dT2dR = -F2*F7o6*T2/Rhoval
         dT2dG =  F2*T2/GRho
         dT2dZ = -F2*dGdZ*T2/G
         dPONdR = -(dLdS*dSdR*PON/PotLC)+dT2dR*PON/T2
         dPONdG = dT2dG*PON/T2
         dPONdZ = F3*dGdZ*PON/G-dLdZ*PON/PotLC+dT2dZ*PON/T2
         dFfracdPON = -F1/((F1-PON)**F2)
         dFexpdPON  = -exp(PON)
         dfa1dR = CcA(1)                 *dFfracdPON*dPONdR
         dfa2dR = CcA(2) *(F2 *Ffrac)    *dFfracdPON*dPONdR
         dfa3dR = CcA(3) *(F3 *Ffrac**F2)*dFfracdPON*dPONdR
         dfa4dR = CcA(4) *(F4 *Ffrac**F3)*dFfracdPON*dPONdR
         dfa5dR = CcA(5) *(F5 *Ffrac**F4)*dFfracdPON*dPONdR
         dfa1dG = CcA(1)                 *dFfracdPON*dPONdG
         dfa2dG = CcA(2) *(F2 *Ffrac)    *dFfracdPON*dPONdG
         dfa3dG = CcA(3) *(F3 *Ffrac**F2)*dFfracdPON*dPONdG
         dfa4dG = CcA(4) *(F4 *Ffrac**F3)*dFfracdPON*dPONdG
         dfa5dG = CcA(5) *(F5 *Ffrac**F4)*dFfracdPON*dPONdG
         dfa1dZ = CcA(1)                 *dFfracdPON*dPONdZ
         dfa2dZ = CcA(2) *(F2 *Ffrac)    *dFfracdPON*dPONdZ
         dfa3dZ = CcA(3) *(F3 *Ffrac**F2)*dFfracdPON*dPONdZ
         dfa4dZ = CcA(4) *(F4 *Ffrac**F3)*dFfracdPON*dPONdZ
         dfa5dZ = CcA(5) *(F5 *Ffrac**F4)*dFfracdPON*dPONdZ
         dfb1dR = CcB(1)                 *dFexpdPON*dPONdR
         dfb2dR = CcB(2) *(F2 *Fexp)     *dFexpdPON*dPONdR
         dfb3dR = CcB(3) *(F3 *Fexp**F2) *dFexpdPON*dPONdR
         dfb4dR = CcB(4) *(F4 *Fexp**F3) *dFexpdPON*dPONdR
         dfb5dR = CcB(5) *(F5 *Fexp**F4) *dFexpdPON*dPONdR
         dfb1dG = CcB(1)                 *dFexpdPON*dPONdG
         dfb2dG = CcB(2) *(F2 *Fexp)     *dFexpdPON*dPONdG
         dfb3dG = CcB(3) *(F3 *Fexp**F2) *dFexpdPON*dPONdG
         dfb4dG = CcB(4) *(F4 *Fexp**F3) *dFexpdPON*dPONdG
         dfb5dG = CcB(5) *(F5 *Fexp**F4) *dFexpdPON*dPONdG
         dfb1dZ = CcB(1)                 *dFexpdPON*dPONdZ
         dfb2dZ = CcB(2) *(F2 *Fexp)     *dFexpdPON*dPONdZ
         dfb3dZ = CcB(3) *(F3 *Fexp**F2) *dFexpdPON*dPONdZ
         dfb4dZ = CcB(4) *(F4 *Fexp**F3) *dFexpdPON*dPONdZ
         dfb5dZ = CcB(5) *(F5 *Fexp**F4) *dFexpdPON*dPONdZ
c
         dFggaCdR = dfa1dR+dfa2dR+dfa3dR+dfa4dR+dfa5dR +
     $              dfb1dR+dfb2dR+dfb3dR+dfb4dR+dfb5dR 
         dFggaCdG = dfa1dG+dfa2dG+dfa3dG+dfa4dG+dfa5dG +
     $              dfb1dG+dfb2dG+dfb3dG+dfb4dG+dfb5dG 
         dFggaCdZ = dfa1dZ+dfa2dZ+dfa3dZ+dfa4dZ+dfa5dZ +
     $              dfb1dZ+dfb2dZ+dfb3dZ+dfb4dZ+dfb5dZ 
c
         dFggaCdA = dFggaCdR + dFggaCdZ*dZdA
         dFggaCdB = dFggaCdR + dFggaCdZ*dZdB
         dGdGx = F1 / (F2*GRho)
         dFggaCdGx = dFggaCdG*dGdGx
c
         Ec = Ec+ (Eloc*FggaC)*qwght(n)
         if(ldew) func(n) = func(n)+Eloc*FggaC
c         
         Amat(n,1) = Amat(n,1) + dElocdA*FggaC + 
     $               Eloc*dFggaCdA
         if (ipol.eq.2) Amat(n,2) = Amat(n,2) + dElocdB*FggaC +
     $               Eloc*dFggaCdB
         Cmat(n,D1_GAA) = Cmat(n,D1_GAA) +   Eloc*dFggaCdGx
         Cmat(n,D1_GAB) = Cmat(n,D1_GAB) +F2*Eloc*dFggaCdGx
         if (ipol.eq.2) Cmat(n,D1_GBB) = Cmat(n,D1_GBB)+ Eloc*dFggaCdGx
   20 continue
      end
c
      Subroutine xc_csogga_d2()
      implicit none
      call errquit(' not coded ',0,0)
      return
      end


      Subroutine GZeta(Zeta,GZet,dGZdz,d2GZdz,d3GZdz)
      Implicit Real*8(A-H,O-Z)
C
C     Evaluate G(Zeta) and its derivatives for DFT.
C
      Save F0,F1,F2,F3,F4,F9,F27
      Data F0/0.0d0/,F1/1.0d0/,F2/2.0d0/,F3/3.0d0/,F4/4.0d0/,
     $     F9/9.0d0/,F27/27.0D0/

      F1o3 = F1/F3
      F1o9 = F1/F9
      F4o27= F4/F27
c
      OMZ3   = F0
      OPZ3   = F0
      GZet   = F0
      dGZdz  = F0
      d2GZdz = F0
      d3GZdz = F0
c
      OMZ = F1-Zeta
      OPZ = F1+Zeta
      OMZ2 = OMZ**2
      OPZ2 = OPZ**2
      if (OMZ.gt.0.0d0) then 
        OMZ3 = OMZ**(-F1o3)
        d2GZdz = d2GZdz-OMZ3/OMZ
        d3GZdz = d3GZdz-OMZ3/OMZ2
      endif
      GZet = GZet+OMZ*OMZ3
      dGZdz = dGZdz-OMZ3
      if (OPZ.gt.0.0d0) then 
        OPZ3 = OPZ**(-F1o3)
        d2GZdz = d2GZdz-OPZ3/OPZ
        d3GZdz = d3GZdz+OPZ3/OPZ2
      endif
      GZet = GZet+OPZ*OPZ3
      dGZdz = dGZdz+OPZ3
c
      GZet   = GZet/F2
      dGZdz  = dGZdz*F1o3
      d2GZdz = d2GZdz*F1o9
      d3GZdz = d3GZdz*F4o27
      Return
      End
