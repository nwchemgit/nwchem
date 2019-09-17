c
c$Id$
c
C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_x_ft97.F
C> The Filatov, Thiel exchange functional
C>
C> @}
C>
C> \ingroup nwpwxc_priv
C> @{
C>
C> \brief Evaluate the Filatov-Thiel exchange functional
C>
C> This subroutine calculates the Filatov-Thiel 97 exchange
C> functional [1,2] recommended variant B.
C> Also the derivatives with respect to the density components and 
C> the dot product of the gradients are computed. 
C>
C> This implementation includes the LDA exchange part [3] of the
C> exchange functional as well as the gradient correction part.
C>
C> The original code was provided by Prof. Walter Thiel.
C>
C> ### References ###
C>
C> [1] M. Filatov, W. Thiel,
C>     "A nonlocal correlation energy density functional from a
C>     Coulomb hole model", Int. J. Quant. Chem. <b>62</b> (1997)
C>     603-616, DOI:
C>     <a href="https://doi.org/10.1002/(SICI)1097-461X(1997)62:6<603::AID-QUA4>3.0.CO;2-%23">
C>     10.1002/(SICI)1097-461X(1997)62:6<603::AID-QUA4>3.0.CO;2-#</a>.
C>
C> [2] M. Filatov, W. Thiel,
C>     "A new gradient-corrected exchange-correlation 
C>     density functional", Mol. Phys. <b>91</b> (1997) 847-859, DOI:
C>     <a href="https://doi.org/10.1080/002689797170950">
C>     10.1080/002689797170950</a>.
C>
C> [3] P.A.M. Dirac, "Note on exchange phenomena in the Thomas atom",
C>     Math. Proc. Cambridge Philos. Soc. <b>26</b> (1930) 376-385, DOI:
C>     <a href="https://doi.org/10.1017/S0305004100016108">
C>     10.1017/S0305004100016108</a>.
C>
      Subroutine nwpwxc_x_ft97(tol_rho, ipol, nq, wght,
     &                       rho, rgamma, func, Amat, Cmat)
      implicit none
#include "nwpwxc_param.fh"
c
c     Input and other parameters
c
      double precision tol_rho !< [Input] The lower limit on the density
      integer ipol             !< [Input] The number of spin channels
      integer nq               !< [Input] The number of points
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density 
c
      double precision rho(nq,*)    !< [Input] The density
c
c     Charge Density Gradient
c
      double precision rgamma(nq,*) !< [Input] The norm of the density gradients
c
c     Sampling Matrices for the XC Potential
c
      double precision func(nq)     !< [Output] The value of the functional
      double precision Amat(nq,*)   !< [Output] The derivative wrt rho
      double precision Cmat(nq,*)   !< [Output] The derivative wrt rgamma
c
      integer n
      double precision gammaval
c to hcth
      double precision rhoa 
      double precision rhob 
      double precision za   
      double precision zb   
      double precision dfdza,dfdzb
c
      double precision fx_ft97,dfdrax,dfdgaax,dfdrbx,dfdgbbx
c
      if(ipol.eq.1) then
        do n=1,nq
          if(rho(n,R_T).gt.tol_rho) then 
            gammaval=rgamma(n,G_TT)
c           gammaval=delrho(n,1,1)*delrho(n,1,1) +
c    &           delrho(n,2,1)*delrho(n,2,1) +
c    &           delrho(n,3,1)*delrho(n,3,1)
            call nwpwxc_rks_x_ft97(rho(n,R_T),gammaval,
     *           fx_ft97,dfdrax,dfdgaax,tol_rho)
            func(n)=func(n)+fx_ft97*wght
            Amat(n,D1_RA) = Amat(n,D1_RA)+dfdrax*wght
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + dfdgaax*wght
          endif
        enddo
      else
        do n=1,nq
          if(rho(n,R_A)+rho(n,R_B).gt.tol_rho) then
            rhoa=rho(n,R_A)
            rhob=rho(n,R_B)
            za=rgamma(n,G_AA)
            zb=rgamma(n,G_BB)
c           za=delrho(n,1,1)*delrho(n,1,1) +
c    &           delrho(n,2,1)*delrho(n,2,1) +
c    &           delrho(n,3,1)*delrho(n,3,1)
c           zb=delrho(n,1,2)*delrho(n,1,2) +
c    &           delrho(n,2,2)*delrho(n,2,2) +
c    &           delrho(n,3,2)*delrho(n,3,2) 
            call nwpwxc_uks_x_ft97(tol_rho,
     (           rhoa,rhob,za,zb,
     *           fx_ft97,dfdrax,dfdrbx,dfdgaax,dfdgbbx)
            func(n)=func(n)+fx_ft97*wght 
            Amat(n,D1_RA) = Amat(n,D1_RA)+dfdrax*wght
            Amat(n,D1_RB) = Amat(n,D1_RB)+dfdrbx*wght
            dfdza=dfdgaax*wght
            dfdzb=dfdgbbx*wght
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + dfdza
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + dfdzb
          endif
        enddo
      endif
      return
      end
cDFT functional repository: xc_ft97 fortran77 source
CDFT repository Quantum Chemistry Group
C
C CCLRC Density functional repository Copyright notice
C This database was prepared as a result of work undertaken by CCLRC. 
C Users may view, print and download the content for personal use only
C and the content must not be used for any commercial purpose without CCLRC
C prior written permission
C

c-----------------------------------------------------------------------
      subroutine nwpwxc_rks_x_ft97(r,g,fx,dfdrax,dfdgaax,
     ,     tol_rho)
c
c     This subroutine calculates the Filatov-Thiel 97 
c     exchange-correlation functional [1,2] for the closed shell case. 
c     This functional is taken to be the correlation functional [1] plus
c     the recommended variant B of the exchange functional [2].
c     Also the derivatives with respect to the density components and 
c     the dot product of the gradients are computed. 
c
c     This implementation includes the LDA exchange part [3] of the
c     exchange functional as well as the gradient correction part.
c
c     The original code was provide by Prof. Walter Thiel.
c
c     [1] M. Filatov, and Walter Thiel,
c         "A nonlocal correlation energy density functional from a
c          Coulomb hole model",
c         Int.J.Quant.Chem. 62 (1997) 603-616.
c
c     [2] M. Filatov, and Walter Thiel,
c         "A new gradient-corrected exchange-correlation 
c          density functional",
c         Mol.Phys. 91 (1997) 847-859.
c
c     [3] P.A.M. Dirac, Proceedings of the Cambridge Philosophical
c         Society, Vol. 26 (1930) 376.
c
c     Parameters:
c
c     r      the total electron density
c     g      the dot product of the total density gradient with itself
c     f      On return the functional value
c     dfdra  On return the derivative of f with respect to the
c            alpha-electron density
c     dfdgaa On return the derivative of f with respect to the dot 
c            product of the alpha-electron density gradient with itself
c
      implicit none
c
      double precision r, g
      double precision fc ,dfdrac ,dfdrb ,dfdgaac ,dfdgbb
      double precision fx,dfdrax,       dfdgaax
c
      double precision rhoa, rhob, rhoa13, rhob13
      double precision gama, gamb
c
c...  Parameters
c
      double precision r13,tol_rho
      parameter (r13=1.0d0/3.0d0)
c
      rhoa   = 0.5d0*r
      rhoa13 = rhoa**r13
      gama   = 0.25d0*g
      rhob   = rhoa 
      rhob13 = rhoa13
      gamb   = gama
c     call FT97_ECFUN(rhoa,rhob,rhoa13,rhob13,gama,gamb,
c    +                fc,dfdrac,dfdrb,dfdgaac,dfdgbb,tol_rho)
      call nwpwxc_FT97_EXFUN(rhoa,rhob,rhoa13,rhob13,
     +                gama,gamb,fx,.false.,.false.,tol_rho)
      call nwpwxc_FT97_EXGRAD(rhoa,rhoa13,gama,dfdrax,dfdgaax,.false.)
      end
c-----------------------------------------------------------------------
      subroutine nwpwxc_uks_x_ft97(tol_rho,ra,rb,gaa,gbb,
     ,     fx,dfdrax,dfdrbx,dfdgaax,dfdgbbx)
c
c     This subroutine calculates the Filatov-Thiel 97 
c     exchange-correlation functional [1,2] for the spin polarised case.
c     This functional is taken to be the correlation functional [1] plus
c     the recommended variant B of the exchange functional [2].
c     Also the derivatives with respect to the density components and 
c     the dot product of the gradients are computed. 
c
c     This implementation includes the LDA exchange part [3] of the
c     exchange functional as well as the gradient correction part.
c
c     The original code was provide by Prof. Walter Thiel.
c
c     [1] M. Filatov, and Walter Thiel,
c         "A nonlocal correlation energy density functional from a
c          Coulomb hole model",
c         Int.J.Quant.Chem. 62 (1997) 603-616.
c
c     [2] M. Filatov, and Walter Thiel,
c         "A new gradient-corrected exchange-correlation 
c          density functional",
c         Mol.Phys. 91 (1997) 847-859.
c
c     [3] P.A.M. Dirac, Proceedings of the Cambridge Philosophical
c         Society, Vol. 26 (1930) 376.
c
c     Parameters:
c
c     ra     the alpha-electron density
c     rb     the beta-electron density
c     gaa    the dot product of the alpha density gradient with itself
c     gbb    the dot product of the beta density gradient with itself
c            the beta density
c     f      On return the functional value
c     dfdra  On return the derivative of f with respect to ra
c     dfdrb  On return the derivative of f with respect to rb
c     dfdgaa On return the derivative of f with respect to gaa
c     dfdgbb On return the derivative of f with respect to gbb
c
      implicit none
c
      double precision ra, rb, gaa, gbb
      double precision fc ,dfdrac ,dfdrbc ,dfdgaac ,dfdgbbc
      double precision fx,dfdrax,dfdrbx,dfdgaax,dfdgbbx
c
      double precision rhoa13, rhob13
c
      double precision r13,tol_rho
      parameter (r13=1.0d0/3.0d0)
c
      rhoa13=0d0
      if(ra.gt.tol_rho**2) rhoa13 = ra**r13
      rhob13=0d0
      if(rb.gt.tol_rho**2) rhob13 = rb**r13
c     call FT97_ECFUN(ra,rb,rhoa13,rhob13,gaa,gbb,
c    +                fc,dfdrac,dfdrbc,dfdgaac,dfdgbbc,tol_rho)
      call nwpwxc_FT97_EXFUN(ra,rb,rhoa13,rhob13,
     +           gaa,gbb,fx,.true.,.false.,tol_rho)
      if(ra.gt.tol_rho**2) 
     +     call nwpwxc_FT97_EXGRAD(ra,rhoa13,gaa,dfdrax,dfdgaax,.false.)
      if(rb.gt.tol_rho**2) 
     +     call nwpwxc_FT97_EXGRAD(rb,rhob13,gbb,dfdrbx,dfdgbbx,.false.)
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      SUBROUTINE nwpwxc_FT97_EXFUN (RHOA,RHOB,RHOA13,RHOB13,AMA,AMB,
     >                            FT97,
     +                            UHF,VARIA,tol_rho)
C     *
C     VALUE OF THE FT97 EXCHANGE xc_ft97 AT A GIVEN POINT IN SPACE.
C     *
C     REFERENCE: M.FILATOV AND W.THIEL, MOL.PHYS. 91, 847 (1997).
C     *
C     ARGUMENT LIST. I=INPUT, O=OUTPUT.
C     RHOA     DENSITY rho.alpha                                 (I)
C     RHOB     DENSITY rho.beta                                  (I)
C     RHOA13   RHOA**(1.0/3.0), CUBIC ROOT OF rho.alpha          (I)
C     RHOB13   RHOB**(1.0/3.0), CUBIC ROOT OF rho.beta           (I)
C     AMA      NORM OF THE GRADIENT OF RHOA WRT X,Y,Z SQUARED    (I)
C     AMB      NORM OF THE GRADIENT OF RHOB WRT X,Y,Z SQUARED    (I)
C     FT97     VALUE OF EXCHANGE xc_ft97                      (O)
C     UHF      LOGICAL FLAG (.TRUE. FOR UHF)                     (I)
C     VARIA    LOGICAL FLAG FOR CHOICE OF xc_ft97             (I)
C              =.TRUE.  VARIANT A                                (I)
C              =.FALSE. VARIANT B (RECOMMENDED)                  (I)
C     *
      IMPLICIT double precision (A-H,O-Z)
      LOGICAL UHF,VARIA
      double precision tol_rho
C     Variants A and B, see eqs.(15) and (16).
      PARAMETER (BETAA=2.930000D-3)
      PARAMETER (AX=9.474169D-4, BX=2.501149D03, BETA0=2.913644D-3)
C     BETA   : Beta(sigma), eq.(16a) for variant B.
      IF(VARIA) THEN
        BETA = BETAA
      ELSE
        BETA = BETA0+AX*AMA/(BX*BX+AMA)
      ENDIF
C     XALPHA : Reduced density gradient xi(sigma), below eq.(11).
      XALPHA = (AMA/(RHOA*RHOA13)**2)
C     DENOMA : Denominator in f(xi), eqs. (15) and (16).
      SINHIA = LOG(XALPHA+SQRT(XALPHA*XALPHA+1.0D0))
      DENOMA = SQRT(1.0D0+9.D0*(BETA**2)*XALPHA*(SINHIA**2))
C     CORRAI : Nonlocal correction to exchange energy density ex.
C     CORRAI : See eqs.(9),(11),(15) and (16).
C     CORRAI : Part of integrand in eq.(11) beyond ex(LDA).
      CORRAI =-RHOA*RHOA13*BETA*XALPHA/DENOMA
      IF(UHF) THEN
         IF(RHOB.LT.tol_rho) THEN
            FT97   = CORRAI
         ELSE
            IF(.NOT.VARIA) BETA = BETA0+AX*AMB/(BX*BX+AMB)
            XBETA  = (AMB/(RHOB*RHOB13)**2)
            SINHIB = LOG(XBETA+SQRT(XBETA*XBETA+1.0D0))
            DENOMB = SQRT(1.0D0+9.D0*BETA**2*XBETA*SINHIB**2)
            CORRBI =-BETA*RHOB*RHOB13*XBETA/DENOMB
            FT97   = CORRAI + CORRBI
         ENDIF
      ELSE
         FT97 = 2.0D0*CORRAI
      ENDIF
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE nwpwxc_FT97_EXGRAD (RHO,RHO13,AM,TERM1,U,VARIA)
C     *
C     DERIVATIVES OF THE FT97 EXCHANGE xc_ft97 WITH RESPECT TO THE
C     DENSITY AT A GIVEN POINT IN SPACE (EXCHANGE POTENTIAL).
C     *
C     REFERENCE: M.FILATOV AND W.THIEL, MOL.PHYS. 91, 847 (1997).
C     *
C     ARGUMENT LIST. I=INPUT, O=OUTPUT.
C     RHO      DENSITY rho                                       (I)
C     RHO13    RHO**(1.0/3.0), CUBIC ROOT OF rho                 (I)
C     AM       NORM OF THE GRADIENT OF RHO WRT X,Y,Z SQUARED     (I)
C     TERM1    DERIVATIVE d(Ex)/d(rho)                           (O)
C     U        DERIVATIVE d(Ex)/d(grad(rho)^2)                   (O)
C     VARIA    LOGICAL FLAG FOR CHOICE OF xc_ft97             (I)
C              =.TRUE.  VARIANT A                                (I)
C              =.FALSE. VARIANT B (RECOMMENDED)                  (I)
C
      IMPLICIT double precision (A-H,O-Z)
      LOGICAL VARIA
C     Variants A and B, see eqs.(15) and (16).
      PARAMETER (BETAA=2.930000D-3)
      PARAMETER (AX=9.474169D-4, BX=2.501149D03, BETA0=2.913644D-3)
      PARAMETER (FT=4.0D0/3.0D0)
C     BETA   : Beta(sigma), eq.(16a) for variant B.
C     DBDG   : Derivative of BETA w.r.t. square of density gradient.
      IF(VARIA) THEN
        BETA = BETAA
        dbdg=0d0
      ELSE
        BETA = BETA0+AX*AM/(BX*BX+AM)
        DBDG = AX*BX*BX/(BX*BX+AM)**2
      ENDIF
      XL     = 1.0D0/(RHO*RHO13)
C     ZT     : Reduced density gradient xi(sigma), below eq.(11).
      ZT     = AM*(XL)**2
C     DENOM  : Denominator in f(xi), eqs. (15) and (16).
      SQZT   = SQRT(ZT*ZT+1.D0)
      SH     = LOG(ZT+SQZT)
      DEN    = 1.D0+9.D0*BETA**2*ZT*SH**2
      DENOM  = SQRT(DEN)
C     F      : Last factor in eqs.(15) and (16).
      F      = BETA*ZT/DENOM
C     FZ     : Derivative of F with respect to ZT (see above).
      FZ     = -4.5D0*BETA**3*ZT*(2.D0*ZT*SH/SQRT(1.D0+ZT**2)+SH**2)/
     1         (DEN*DENOM) + BETA/DENOM
C     TERM1  : Derivative d(Ex)/d(rho).
      TERM1  = FT*RHO13*(2.D0*ZT*FZ - F)
C     U      : Derivative d(Ex)/d(grad(rho)^2), variants A and B.
      IF(VARIA) THEN
        U    =-0.5d0*(FZ+FZ)*XL
      ELSE
C       EB   : Derivative of exchange energy with respect to BETA.
        EB   =-RHO*RHO13*ZT/(DEN*DENOM)
        U    =-(FZ*XL-EB*DBDG)
      ENDIF
      RETURN
      END
c-----------------------------------------------------------------------
C> @}
