#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_c_sogga.F
C> The SOGGA11 and SOGGA-X correlation functionals
C>
C> @}
#endif
#endif
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief Evaluate the SOGGA11 and SOGGA11-X correlation functionals
C>
C> The SOGGA11 and SOGGA11-X functionals are GGA functionals that are
C> formulated and optimized such that the coefficient of the gradient
C> correction term matches that of the appropriate expansion of the
C> electronic energy [1,2].
C>
C> ### References ###
C>
C> [1] R. Peverati, Y. Zhao, D.G. Truhlar, 
C>     "Generalized gradient approximation that recovers the
C>     second-order density-gradient expansion with optimized
C>     across-the-board performance", J. Phys. Chem. Lett. <b>2</b>
C>     (2011) 1991-1997, DOI:
C>     <a href="https://doi.org/10.1021/jz200616w">
C>     10.1021/jz200616w</a>.
C>
C> [2] R. Peverati, D.G. Truhlar, 
C>     "Communication: A global hybrid generalized gradient
C>     approximation to the exchange-correlation functional that
C>     satisfies the second-order density-gradient constraint and has
C>     broad applicability in chemistry", J. Chem. Phys. <b>135</b>
C>     (2011) 191102, DOI:
C>     <a href="https://doi.org/10.1063/1.3663871">
C>     10.1063/1.3663871</a>.
C>
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      subroutine nwxc_c_sogga_p(param, tol_rho, ipol, nq, wght, rho, 
     &                          rgamma, ffunc)
#else
      subroutine nwxc_c_sogga(param, tol_rho, ipol, nq, wght, rho, 
     &                        rgamma, ffunc)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      subroutine nwxc_c_sogga_d2(param, tol_rho, ipol, nq, wght, rho, 
     &                        rgamma, ffunc)
#else
      subroutine nwxc_c_sogga_d3(param, tol_rho, ipol, nq, wght, rho, 
     &                        rgamma, ffunc)
#endif
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
#include "nwad.fh"
      implicit none
c
#include "intf_nwxc_c_lsda.fh"
#include "intf_nwxc_GZeta.fh"
c
#include "nwxc_param.fh"
c
c     Input and other parameters
c
#if defined(NWAD_PRINT)
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::param(*) !< [Input] Parameters of the functional,
      type(nwad_dble)::CcA(0:5),CcB(0:5)
      type(nwad_dble)::FA0, FB0
#else
      double precision param(*) !< [Input] Parameters of the functional,
      double precision CcA(0:5),CcB(0:5)
      double precision FA0, FB0
#endif
#else
      double precision param(*)!< [Input] Parameters of the functional,
                               !< see Table 1 of [1] and Table 1 of [2].
                               !< - param(1): \f$ a_0 \f$
                               !< - param(2): \f$ a_1 \f$
                               !< - param(3): \f$ a_2 \f$
                               !< - param(4): \f$ a_3 \f$
                               !< - param(5): \f$ a_4 \f$
                               !< - param(6): \f$ a_5 \f$
                               !< - param(7): \f$ b_0 \f$
                               !< - param(8): \f$ b_1 \f$
                               !< - param(9): \f$ b_2 \f$
                               !< - param(10): \f$ b_3 \f$
                               !< - param(11): \f$ b_4 \f$
                               !< - param(12): \f$ b_5 \f$
      double precision CcA(0:5),CcB(0:5)
      double precision FA0, FB0
#endif
      double precision tol_rho !< [Input] The lower limit on the density
      integer ipol             !< [Input] The number of spin channels
      integer nq               !< [Input] The number of points
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density 
c
      type(nwad_dble)::rho(nq,*)    !< [Input] The density
c
c     Charge Density Gradient
c
      type(nwad_dble)::rgamma(nq,*) !< [Input] The norm of the density gradients
c
c     Sampling Matrices for the XC Potential
c
      type(nwad_dble)::ffunc(nq)    !< [Output] The value of the functional
c     double precision Amat(nq,*)   !< [Output] The derivative wrt rho
c     double precision Cmat(nq,*)   !< [Output] The derivative wrt rgamma
c
c     Intermediate derivative results, etc.
c
      integer n
c
      type(nwad_dble)::RHOA,RHOB,rhoval,gammaval,RS,Zeta
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
      double precision DT2DG, DT2DR, DT2DZ, DZDA, DZDB
      double precision F1O3, F1O6, F7O6 
      type(nwad_dble)::FA1, FA2, FA3, FA4, FA5
      type(nwad_dble)::FB1, FB2, FB3, FB4, FB5
      type(nwad_dble)::FEXP, FFRAC
      type(nwad_dble)::FGGAC
      double precision FKFAC
      type(nwad_dble)::G, G2, G3
      double precision PI, PI34
c     type(nwad_dble)::PON, T, T2
      type(nwad_dble)::PON, T2
      type(nwad_dble)::RHO16, RHO76, GRHO, POTLC, ELOC
      double precision SCFAC, SKFAC, XNU
      double precision dummy

 
      double precision F1, F2, F3, F4, F5, F6, F7
      Save F1, F2, F3, F4, F5, F6, F7
      DATA F1/1.0D+00/,  F2/2.0D+00/,  F3/3.0D+00/,
     $     F4/4.0D+00/,  F5/5.0D+00/,  F6/6.0D+00/,
     $     F7/7.0D+00/
c
c     ======> BOTH SPIN-RESTRICETED AND UNRESTRICTED <======
c
c     if (ijzy.eq.1) then
c      CcA(0) =  0.50000D+00
c      CcA(1) = -4.62334D+00
c      CcA(2) =  8.00410D+00
c      CcA(3) = -130.226D+00
c      CcA(4) =  38.2685D+00
c      CcA(5) =  69.5599D+00
c      CcB(0) =  0.50000D+00
c      CcB(1) =  3.62334D+00
c      CcB(2) =  9.36393D+00
c      CcB(3) =  34.5114D+00
c      CcB(4) = -18.5684D+00
c      CcB(5) = -0.16519D+00
c     elseif (ijzy.eq.2) then
c      CcA(0) =  5.00000d-01
c      CcA(1) =  7.82439d+01
c      CcA(2) =  2.57211d+01
c      CcA(3) = -1.38830d+01
c      CcA(4) = -9.87375d+00
c      CcA(5) = -1.41357d+01
c      CcB(0) =  5.00000d-01
c      CcB(1) = -7.92439d+01
c      CcB(2) =  1.63725d+01
c      CcB(3) =  2.08129d+00
c      CcB(4) =  7.50769d+00
c      CcB(5) = -1.01861d+01
c     endif
      CcA(0) = param(1)
      CcA(1) = param(2)
      CcA(2) = param(3)
      CcA(3) = param(4)
      CcA(4) = param(5)
      CcA(5) = param(6)
      CcB(0) = param(7)
      CcB(1) = param(8)
      CcB(2) = param(9)
      CcB(3) = param(10)
      CcB(4) = param(11)
      CcB(5) = param(12)

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
           RHOA = rho(n,R_T)/F2
           RHOB = RHOA
           gammaval = rgamma(n,G_TT)
c          gammaval =(delrho(n,1,1)*delrho(n,1,1) +
c    &                delrho(n,2,1)*delrho(n,2,1) +
c    &                delrho(n,3,1)*delrho(n,3,1))
         else
           RhoA = 0.0d0
           RhoB = 0.0d0
           gammaval = 0.0d0
           if (rho(n,R_A).gt.0.5d0*tol_rho) then
             RhoA = rho(n,R_A)
             gammaval = gammaval + rgamma(n,G_AA)
           endif
           if (rho(n,R_B).gt.0.5d0*tol_rho) then
             RhoB = rho(n,R_B)
             gammaval = gammaval + rgamma(n,G_BB)
             if (rho(n,R_A).gt.0.5d0*tol_rho) then
               gammaval = gammaval + 2.0d0*rgamma(n,G_AB)
             endif
           endif
c          gammaval = delrho(n,1,1)*delrho(n,1,1) +
c    &                delrho(n,1,2)*delrho(n,1,2) +
c    &                delrho(n,2,1)*delrho(n,2,1) +
c    &                delrho(n,2,2)*delrho(n,2,2) +
c    &                delrho(n,3,1)*delrho(n,3,1) +
c    &                delrho(n,3,2)*delrho(n,3,2) +
c    &          2.d0*(delrho(n,1,1)*delrho(n,1,2) +
c    &                delrho(n,2,1)*delrho(n,2,2) +
c    &                delrho(n,3,1)*delrho(n,3,2))
         endif
         Rhoval  = RhoA + RhoB
c        GRho = max(dsqrt(gammaval),tol_rho)
c        GRho = sqrt(gammaval)
         GRho = gammaval
         if (rhoval.le.tol_rho) goto 20
         rS = (Pi34/Rhoval)**F1o3
         Zeta = (RhoA-RhoB)/Rhoval
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
         Call nwxc_c_lsda_p(tol_rho,
     $        RS,Zeta,PotLC,dLdS,dLdZ,dummy,dummy,dummy)
         Call nwxc_GZeta_p(Zeta,G,dGdZ,dummy,dummy)
#else
         Call nwxc_c_lsda(tol_rho,
     $        RS,Zeta,PotLC,dLdS,dLdZ,dummy,dummy,dummy)
         Call nwxc_GZeta(Zeta,G,dGdZ,dummy,dummy)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
         Call nwxc_c_lsda_d2(tol_rho,
     $        RS,Zeta,PotLC,dLdS,dLdZ,dummy,dummy,dummy)
         Call nwxc_GZeta_d2(Zeta,G,dGdZ,dummy,dummy)
#else
         Call nwxc_c_lsda_d3(tol_rho,
     $        RS,Zeta,PotLC,dLdS,dLdZ,dummy,dummy,dummy)
         Call nwxc_GZeta_d3(Zeta,G,dGdZ,dummy,dummy)
#endif
c
         Eloc = SCfac*Rhoval*PotLC
         FKFac = (F3*Pi*Pi)**F1o3
         SKFac = Sqrt(FKFac/Pi)*F2
         Rho16 = Rhoval**F1o6
         Rho76 = Rhoval*Rho16
c        T  = GRho/(F2*SKFac*Rho76*G)
c        T2 = T*T
         T2 = GRho/((F2*SKFac*Rho76*G)**2.0d0)
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
c        dSdR  = -(F1o3*rS/Rhoval)
c        dZdA = (F1-Zeta)/Rhoval
c        dZdB = (-F1-Zeta)/Rhoval
c        dLdA = dLdS*dSdR + dLdz*dZdA
c        dLdB = dLdS*dSdR + dLdz*dZdB
c        dElocdA = SCfac*PotLC + SCfac*Rhoval*dLdA
c        dElocdB = SCfac*PotLC + SCfac*Rhoval*dLdB
c        dT2dR = -F2*F7o6*T2/Rhoval
c        dT2dG =  F2*T2/GRho
c        dT2dZ = -F2*dGdZ*T2/G
c        dPONdR = -(dLdS*dSdR*PON/PotLC)+dT2dR*PON/T2
c        dPONdG = dT2dG*PON/T2
c        dPONdZ = F3*dGdZ*PON/G-dLdZ*PON/PotLC+dT2dZ*PON/T2
c        dFfracdPON = -F1/((F1-PON)**F2)
c        dFexpdPON  = -exp(PON)
c        dfa1dR = CcA(1)                 *dFfracdPON*dPONdR
c        dfa2dR = CcA(2) *(F2 *Ffrac)    *dFfracdPON*dPONdR
c        dfa3dR = CcA(3) *(F3 *Ffrac**F2)*dFfracdPON*dPONdR
c        dfa4dR = CcA(4) *(F4 *Ffrac**F3)*dFfracdPON*dPONdR
c        dfa5dR = CcA(5) *(F5 *Ffrac**F4)*dFfracdPON*dPONdR
c        dfa1dG = CcA(1)                 *dFfracdPON*dPONdG
c        dfa2dG = CcA(2) *(F2 *Ffrac)    *dFfracdPON*dPONdG
c        dfa3dG = CcA(3) *(F3 *Ffrac**F2)*dFfracdPON*dPONdG
c        dfa4dG = CcA(4) *(F4 *Ffrac**F3)*dFfracdPON*dPONdG
c        dfa5dG = CcA(5) *(F5 *Ffrac**F4)*dFfracdPON*dPONdG
c        dfa1dZ = CcA(1)                 *dFfracdPON*dPONdZ
c        dfa2dZ = CcA(2) *(F2 *Ffrac)    *dFfracdPON*dPONdZ
c        dfa3dZ = CcA(3) *(F3 *Ffrac**F2)*dFfracdPON*dPONdZ
c        dfa4dZ = CcA(4) *(F4 *Ffrac**F3)*dFfracdPON*dPONdZ
c        dfa5dZ = CcA(5) *(F5 *Ffrac**F4)*dFfracdPON*dPONdZ
c        dfb1dR = CcB(1)                 *dFexpdPON*dPONdR
c        dfb2dR = CcB(2) *(F2 *Fexp)     *dFexpdPON*dPONdR
c        dfb3dR = CcB(3) *(F3 *Fexp**F2) *dFexpdPON*dPONdR
c        dfb4dR = CcB(4) *(F4 *Fexp**F3) *dFexpdPON*dPONdR
c        dfb5dR = CcB(5) *(F5 *Fexp**F4) *dFexpdPON*dPONdR
c        dfb1dG = CcB(1)                 *dFexpdPON*dPONdG
c        dfb2dG = CcB(2) *(F2 *Fexp)     *dFexpdPON*dPONdG
c        dfb3dG = CcB(3) *(F3 *Fexp**F2) *dFexpdPON*dPONdG
c        dfb4dG = CcB(4) *(F4 *Fexp**F3) *dFexpdPON*dPONdG
c        dfb5dG = CcB(5) *(F5 *Fexp**F4) *dFexpdPON*dPONdG
c        dfb1dZ = CcB(1)                 *dFexpdPON*dPONdZ
c        dfb2dZ = CcB(2) *(F2 *Fexp)     *dFexpdPON*dPONdZ
c        dfb3dZ = CcB(3) *(F3 *Fexp**F2) *dFexpdPON*dPONdZ
c        dfb4dZ = CcB(4) *(F4 *Fexp**F3) *dFexpdPON*dPONdZ
c        dfb5dZ = CcB(5) *(F5 *Fexp**F4) *dFexpdPON*dPONdZ
c
c        dFggaCdR = dfa1dR+dfa2dR+dfa3dR+dfa4dR+dfa5dR +
c    $              dfb1dR+dfb2dR+dfb3dR+dfb4dR+dfb5dR 
c        dFggaCdG = dfa1dG+dfa2dG+dfa3dG+dfa4dG+dfa5dG +
c    $              dfb1dG+dfb2dG+dfb3dG+dfb4dG+dfb5dG 
c        dFggaCdZ = dfa1dZ+dfa2dZ+dfa3dZ+dfa4dZ+dfa5dZ +
c    $              dfb1dZ+dfb2dZ+dfb3dZ+dfb4dZ+dfb5dZ 
c
c        dFggaCdA = dFggaCdR + dFggaCdZ*dZdA
c        dFggaCdB = dFggaCdR + dFggaCdZ*dZdB
c        dGdGx = F1 / (F2*GRho)
c        dFggaCdGx = dFggaCdG*dGdGx
c
c        Ec = Ec+ (Eloc*FggaC)*qwght(n)
         ffunc(n) = ffunc(n)+ (Eloc*FggaC)*wght
c         
c        Amat(n,D1_RA) = Amat(n,D1_RA) + dElocdA*FggaC*wght +
c    $                                   Eloc*dFggaCdA*wght
c        if (ipol.eq.2) then
c          Amat(n,D1_RB) = Amat(n,D1_RB) + dElocdB*FggaC*wght +
c    $                                     Eloc*dFggaCdB*wght
c        endif
c        if (ipol.eq.1) then
c          Cmat(n,D1_GAA) = Cmat(n,D1_GAA) +   Eloc*dFggaCdGx*wght
c          Cmat(n,D1_GAB) = Cmat(n,D1_GAB) +F2*Eloc*dFggaCdGx*wght
c        else
c          Cmat(n,D1_GAA) = Cmat(n,D1_GAA) +   Eloc*dFggaCdGx*wght
c          Cmat(n,D1_GAB) = Cmat(n,D1_GAB) +F2*Eloc*dFggaCdGx*wght
c          Cmat(n,D1_GBB) = Cmat(n,D1_GBB) +   Eloc*dFggaCdGx*wght
c        endif
   20 continue
      end
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_GZeta_p(Zeta,GZet,dGZdz,d2GZdz,d3GZdz)
#else
      Subroutine nwxc_GZeta(Zeta,GZet,dGZdz,d2GZdz,d3GZdz)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_GZeta_d2(Zeta,GZet,dGZdz,d2GZdz,d3GZdz)
#else
      Subroutine nwxc_GZeta_d3(Zeta,GZet,dGZdz,d2GZdz,d3GZdz)
#endif
#include "nwad.fh"
      Implicit Real*8(A-H,O-Z)
C
C     Evaluate G(Zeta) and its derivatives for DFT.
C
      type(nwad_dble)::Zeta,GZet
      type(nwad_dble)::OMZ,OMZ2,OMZ3,OPZ,OPZ2,OPZ3
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
c     dGZdz  = F0
c     d2GZdz = F0
c     d3GZdz = F0
c
      OMZ = F1-Zeta
      OPZ = F1+Zeta
      OMZ2 = OMZ**2.0d0
      OPZ2 = OPZ**2.0d0
      if (OMZ.gt.F0) then
        OMZ3 = OMZ**(-F1o3)
c       d2GZdz = d2GZdz-OMZ3/OMZ
c       d3GZdz = d3GZdz-OMZ3/OMZ2
      endif
      GZet = GZet+OMZ*OMZ3
c     dGZdz = dGZdz-OMZ3
      if (OPZ.gt.F0) then
        OPZ3 = OPZ**(-F1o3)
c       d2GZdz = d2GZdz-OPZ3/OPZ
c       d3GZdz = d3GZdz+OPZ3/OPZ2
      endif
      GZet = GZet+OPZ*OPZ3
c     dGZdz = dGZdz+OPZ3
c
      GZet   = GZet/F2
c     dGZdz  = dGZdz*F1o3
c     d2GZdz = d2GZdz*F1o9
c     d3GZdz = d3GZdz*F4o27
      Return
      End
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_c_sogga.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_c_sogga.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_c_sogga.F"
#endif
#undef NWAD_PRINT
C>
C> @}
