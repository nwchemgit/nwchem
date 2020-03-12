c
c$Id$
c
C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_c_ft97.F
C> The Filatov, Thiel correlation functional
C>
C> @}
C>
C> \ingroup nwpwxc_priv
C> @{
C>
C> \brief Evaluate the Filatov-Thiel correlation functional
C>
C> This subroutine calculates the Filatov-Thiel 97 correlation
C> functional [1,2].
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
      Subroutine nwpwxc_c_ft97(tol_rho, ipol, nq, wght, rho, rgamma,
     &                       func, Amat, Cmat) 
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
      integer n
      double precision gammaval
c to hcth
      double precision rhoa 
      double precision rhob 
      double precision za   
      double precision zb   
      double precision dfdza,dfdzb
c
      double precision fc_ft97,dfdrac,dfdgaac,dfdrbc,dfdgbbc
c
      if(ipol.eq.1) then
        do n=1,nq
          if(rho(n,R_T).gt.tol_rho) then 
            gammaval=rgamma(n,G_TT)
c           gammaval=delrho(n,1,1)*delrho(n,1,1) +
c    &           delrho(n,2,1)*delrho(n,2,1) +
c    &           delrho(n,3,1)*delrho(n,3,1)
            call nwpwxc_rks_c_ft97(rho(n,R_T),gammaval,
     *           fc_ft97,dfdrac,dfdgaac,tol_rho)
            func(n)=func(n)+fc_ft97*wght
            Amat(n,D1_RA) = Amat(n,D1_RA)+dfdrac*wght
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + dfdgaac*wght
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
            call nwpwxc_uks_c_ft97(tol_rho,
     (           rhoa,rhob,za,zb,
     *           fc_ft97,dfdrac,dfdrbc,dfdgaac,dfdgbbc)
            func(n)=func(n)+fc_ft97*wght
            Amat(n,D1_RA) = Amat(n,D1_RA)+dfdrac*wght
            Amat(n,D1_RB) = Amat(n,D1_RB)+dfdrbc*wght
            dfdza=dfdgaac*wght
            dfdzb=dfdgbbc*wght
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
      subroutine nwpwxc_rks_c_ft97(r,g,fc,dfdrac,dfdgaac,tol_rho)
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
      call nwpwxc_FT97_ECFUN(rhoa,rhob,rhoa13,rhob13,gama,gamb,
     +                fc,dfdrac,dfdrb,dfdgaac,dfdgbb,tol_rho)
c     call FT97_EXFUN(rhoa,rhob,rhoa13,rhob13,
c    +                gama,gamb,fx,.false.,.false.,tol_rho)
c     call FT97_EXGRAD(rhoa,rhoa13,gama,dfdrax,dfdgaax,.false.)
      end
c-----------------------------------------------------------------------
      subroutine nwpwxc_uks_c_ft97(tol_rho,ra,rb,gaa,gbb,
     ,     fc,dfdrac,dfdrbc,dfdgaac,dfdgbbc)
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
      call nwpwxc_FT97_ECFUN(ra,rb,rhoa13,rhob13,gaa,gbb,
     +                fc,dfdrac,dfdrbc,dfdgaac,dfdgbbc,tol_rho)
c     call FT97_EXFUN(ra,rb,rhoa13,rhob13,
c    +           gaa,gbb,fx,.true.,.false.,tol_rho)
c     if(ra.gt.tol_rho**2) 
c    +     call FT97_EXGRAD(ra,rhoa13,gaa,dfdrax,dfdgaax,.false.)
c     if(rb.gt.tol_rho**2) 
c    +     call FT97_EXGRAD(rb,rhob13,gbb,dfdrbx,dfdgbbx,.false.)
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      SUBROUTINE nwpwxc_FT97_ECFUN (RHOA,RHOB,RHOA13,RHOB13,AMA,AMB,EC,
     1                       V1,V2,V3,V4,tol_rho)
C     *
C     VALUE OF THE CORRELATION xc_ft97 (EC) AND ITS DERIVATIVES
C     (V1,V2,V3,V4) WITH REGARD TO THE DENSITY AT A GIVEN POINT
C     IN SPACE WITH CARTESIAN COORDINATES X,Y,Z.
C     *
C     REFERENCE:
C     M.FILATOV AND W.THIEL, 
C     "A nonlocal correlation energy density functional from a Coulomb
C      hole model", INT.J.QUANTUM CHEM. Vol. 62, (1997) 603-616.
C     *
C     ARGUMENT LIST. I=INPUT, O=OUTPUT.
C     RHOA     DENSITY rho.alpha                                 (I)
C     RHOB     DENSITY rho.beta                                  (I)
C     RHOA13   RHOA**(1.0/3.0), CUBIC ROOT OF rho.alpha          (I)
C     RHOB13   RHOB**(1.0/3.0), CUBIC ROOT OF rho.beta           (I)
C     AMA      NORM OF THE GRADIENT OF RHOA WRT X,Y,Z SQUARED    (I)
C     AMB      NORM OF THE GRADIENT OF RHOB WRT X,Y,Z SQUARED    (I)
C     EC       CORRELATION ENERGY                                (O)
C     V1       DERIVATIVE d(Ec)/d(rho.alpha)                     (I,O)
C     V2       DERIVATIVE d(Ec)/d(rho.beta)                      (I,O)
C     V3       DERIVATIVE d(Ec)/d(grad(rho.alpha)^2)             (I,O)
C     V4       DERIVATIVE d(Ec)/d(grad(rho.beta)^2)              (I,O)
C     UHF      LOGICAL FLAG (.TRUE. FOR UHF)                     (I)
C     *
      IMPLICIT double precision (A-H,O-Z)
      EC     = 0.D0
      V1     = 0.D0
      V2     = 0.D0
      V3     = 0.D0
      V4     = 0.D0
      CALL nwpwxc_FT97_VCOUR (RHOA,RHOB,RHOA13,RHOB13,AMA,AMB,EAA,EBB,
     1            EAB,EBA,
     2            DEAADR,DEBBDR,DEABDR,DEBADR,
     3            DEAADG,DEBBDG,DEABDG,DEBADG,tol_rho)
C     EC IS THE CONTRIBUTION TO THE CORRELATION ENERGY AT (X,Y,Z).
C     EC IS DEFINED IN EQUATION (4) OF THE REFERENCE.
C     Add parallel-spin contribution to the energy.
      EC    = EC+RHOA*EAA+RHOB*EBB
C     Add opposite-spin contribution to the energy.
      EC    = EC+RHOA*EAB+RHOB*EBA
C     Add parallel-spin contribution to the derivatives.
      V1    = V1 +EAA+RHOA*DEAADR
      V2    = V2 +EBB+RHOB*DEBBDR
      V3    = V3 +RHOA*DEAADG
      V4    = V4 +RHOB*DEBBDG
C     Add opposite-spin contribution to the derivatives.
      V1    = V1 +EAB+RHOB*DEBADR
      V2    = V2 +EBA+RHOA*DEABDR
      V3    = V3 +RHOB*DEBADG
      V4    = V4 +RHOA*DEABDG
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE nwpwxc_FT97_VCOUR (RHOA,RHOB,RHOA13,RHOB13,AMA,AMB,
     1                       EAA,EBB,EAB,EBA,
     2                       DEAADR,DEBBDR,DEABDR,DEBADR,
     3                       DEAADG,DEBBDG,DEABDG,DEBADG,tol_rho)
C     *
C     COMPUTATION OF TERMS IN THE CORRELATION xc_ft97.
C     *
C     REFERENCE:
C     M.FILATOV AND W.THIEL, INT.J.QUANTUM CHEM. 62, 603 (1997).
C     *
C     NOTATION FOR ARGUMENTS.
C     RHOA   - rho(alpha), density, eq.(2)
C     RHOB   - rho(beta )
C     RHOA13 - rho(alpha)**1/3, cubic root of density
C     RHOB13 - rho(beta )**1/3 
C     AMA    - grad(rho(alpha))**2, gradient norm
C     AMB    - grad(rho(beta ))**2
C     EAA    - epsilon(alpha,alpha), correlation energy density, eq.(5)
C     EBB    - epsilon(beta ,beta )
C     EAB    - epsilon(alpha,beta )
C     EBA    - epsilon(beta ,alpha)
C     DEDRAA - d(Ec)/d(rho.alpha)
C     DEDGAA - d(Ec)/d(grad(rho.alpha)^2)
C     *
      IMPLICIT double precision (A-H,O-Z)
      double precision KSSP0,KSS0,MUAB,MUBA,MUBB,MUAA,
     *     KSSPBIG,big,kss0big
      double precision nwpwxc_FT97_EXPEI
C     NOTATION FOR CONSTANTS.
C     C0 = (1-ln2)/(2*Pi**2), see eq.(9)
C     C1 = 2*C0/3           , see eq.(13,(28),(33)
C     C2 = (3/(4*Pi))**1/3  , see eq.(8)
C     C3 = C2/3
      PARAMETER (C0=0.1554534543483D-1)
      PARAMETER (C1=0.2072712724644D-1)
      PARAMETER (C2=0.6203504908994D00)
      PARAMETER (C3=0.2067834969665D00)
C     OPTIMIZED PARAMETERS IN GRADIENT-CORRECTED CORRELATION xc_ft97.
C     See formulas for FSSP, eq.(45), and FSS, eq.(46).
      PARAMETER (A1=1.622118767D0)
      PARAMETER (A2=0.489958076D0)
      PARAMETER (A3=1.379021941D0)
      PARAMETER (A4=4.946281353D0)
      PARAMETER (A5=3.600612059D0)
C     NUMERICAL CUTOFF FOR MUAA, MUAB, MUBA, MUBB, see eqs.(13) and (33).
      PARAMETER (CUTOFF=1.0D7)
      parameter(ksspbig=1.291551074D0 - 0.349064173D0,big=1d4)
      parameter(KSS0BIG  = 1.200801774D0 + 0.859614445D0-
     -           0.812904345D0)
C
C *** DEFINITION OF FORMULA FUNCTIONS.
C
C     R        : Density parameter Rs (Wigner radius), see eq.(8).
C     KSSP0    : See eq.(39) for k(sigma,sigma').
      KSSP0(R) = 1.291551074D0 - 0.349064173D0*(1.D0-
     .           EXP(-0.083275880D0*R**0.8D0))
C     KSS0     : See eq.(40) for k(sigma,sigma).
      KSS0(R)  = 1.200801774D0 + 0.859614445D0*(1.D0-
     .           EXP(-1.089338848*SQRT(R)))-
     .           0.812904345D0*(1.D0-EXP(-0.655638823D0*R**0.4D0))
C
C     FSSP     : See eq.(45) for F(sigma,sigma').
      FSSP(R,GR) = (1.D0+A1*GR+(A2*GR)**2)*
     .             EXP(-(A2*GR)**2)/SQRT(1.D0+A3*GR/R)
C     USSP1    : Derivative of ln(k(sigma,sigma')) with respect to Rs.
      USSP1(R) = -4.65098019D-2*EXP(-0.083275880D0*R**0.8D0)*R**0.8D0 
C     USSP2    : Derivative of ln(F(sigma,sigma')) with respect to Rs.
      USSP2(R,GR) = A3*GR/(R+A3*GR)
C     WSSP     : Derivative of ln(F(sigma,sigma')) w.r.t. Nabla(Rs)^2
      WSSP(R,GR) = ((A1+A1)-A1*(A2+A2)**2*GR**2-(A2*(A2+A2))**2*GR**3)
     .             /(1.D0+A1*GR+(A2*GR)**2) - A3/(R+A3*GR)
C
C     FSS      : See eq.(46) for F(sigma,sigma).
      FSS(R,GR)= (1.D0+(A4*GR)**2)*EXP(-(A4*GR)**2)/
     .           SQRT(1.D0+A5*GR/R)
C     FACTOR   : See eq.(34).
      FACTOR(R)= EXP(-(R/(0.939016D0*SQRT(R)+1.733170D0*R))**2)
C     USS1     : Derivative of ln(k(sigma,sigma)) with respect to Rs.
      USS1(R)  =-4.26377318D-1*R**0.4D0*EXP(-0.655638823D0*R**0.4D0)+
     .           0.93641140924D0*SQRT(R)*EXP(-1.089338848*SQRT(R))
C     USS2      : Derivative of ln(F(sigma,sigma)) with respect to Rs.
      USS2(R,GR)= A5*GR/(R+A5*GR)
C     WSS       : Derivative of ln(F(sigma,sigma)) w.r.t. Nabla(Rs)^2
      WSS(R,GR) =-(A4*(A4+A4))**2*GR**3/(1.D0+(A4*GR)**2)-A5/(R+A5*GR)
C     YSS       : Derivative of FACTOR (eq.(34)) with respect to Rs.
      YSS(R)    = 0.939016D0*SQRT(R)*R*R/(0.939016D0*SQRT(R)+
     .            1.733170D0*R)**3
C
C     FF        : Used in approximate normalization, eq.(15).
      FF(T)     = (3.D0+2.D0*(SQRT(T)+T))/
     .            (3.D0+6.D0*(SQRT(T)+T))
C     DFFDMU    : Derivative of FF with respect to the argument T.
      DFFDMU(T) =-(2.D0*(SQRT(T)+T+T))/
     .            (3.D0*(1.D0+2.D0*(SQRT(T)+T))**2)
C
C *** RSA - Rs(alpha), GA - grad(Rs(alpha))
C
      IF(RHOA.GT.tol_rho) THEN
C         Wigner radius Rs(alpha), eq.(8)
          RSA = C2/RHOA13
          CFT = C3/(RHOA*RHOA13)
C         GA  : Nabla(Rs(alpha))
C         GA2 : Square of Nabla(Rs(alpha))
          GA2 = CFT*CFT*AMA
C         mu(beta,alpha), eq.(13)
C         MUBA = C1*RSA/(KSSP0(RSA)*FSSP(RSA,GA2))**2
          EBA = 0.D0
          ECMUC0 = 0.D0
          DEBADR = 0.D0
          DEBADG = 0.D0
          if((ga2*a2)**2.gt.big) then
             DENOM  = 0d0
          else
             if(rsa.gt.big) then
                DENOM  = (KSSPBIG*FSSP(RSA,GA2))**2
             else
                DENOM  = (KSSP0(RSA)*FSSP(RSA,GA2))**2
             endif
          endif
          IF(C1*RSA .LE. DENOM*CUTOFF) THEN
             MUBA = C1*RSA/DENOM
             EEI = nwpwxc_FT97_EXPEI(-MUBA)
             EEI1 = MUBA*EEI+1.D0
C            EBA : Correlation energy density, eq.(19)
             EBA = C0*(EEI+2.D0*FF(MUBA)*EEI1)
C            ECMUC0 : Derivative of EBA w.r.t. MUBA
             ECMUC0 = C0*(EEI1*(1.D0+2.D0*(DFFDMU(MUBA)+MUBA*
     .                FF(MUBA)))+2.D0*FF(MUBA)*MUBA*EEI)
C            d(Ec)/d(Nabla(rho)^2)  =  -CFT^2*(mu*dEcdmu)*Wss'
             DEBADG = -WSSP(RSA,GA2)*ECMUC0*CFT*CFT 
C            d(Ec)/d(rho)
             DEBADR = -ECMUC0*(1.D0-USSP1(RSA)/KSSP0(RSA)-
     .                 USSP2(RSA,GA2) -
     .                 8.D0*WSSP(RSA,GA2)*GA2)
     .                 / (3.D0*RHOA)
         ENDIF
C        mu(alpha,alpha), eq.(33)
C        MUAA = C1*RSA/(KSS0(RSA)*FSS(RSA,GA2))**2
         FA   = FACTOR(RSA)
         EAA  = 0.D0
         ECMUC0 = 0.D0
         DEAADR = 0.D0
         DEAADG = 0.D0
         if((ga2*a2)**2.gt.big) then
            DENOM  = 0d0
         else
            if(rsa.gt.big) then
              DENOM  = (KSS0big*FSS(RSA,GA2))**2
           else
              DENOM  = (KSS0(RSA)*FSS(RSA,GA2))**2
           endif
        endif
         IF(C1*RSA .LE. DENOM*CUTOFF) THEN
            MUAA = C1*RSA/DENOM
            EEI  = nwpwxc_FT97_EXPEI(-MUAA)
            EEI1 = MUAA*EEI+1.D0
C           EAA  : Correlation energy density, eqs.(19) and (31).
            EAA  = FA*C0*(EEI+2.D0*FF(MUAA)*EEI1)
C           ECMUC0 : Derivative of EAA w.r.t. MUAA
            ECMUC0 = FA*C0*(EEI1*(1.D0+2.D0*(DFFDMU(MUAA)+MUAA*
     .               FF(MUAA)))+2.D0*FF(MUAA)*MUAA*EEI)
C           d(Ec)/d(rho)
            DEAADR = (-ECMUC0*(1.D0-USS1(RSA)/KSS0(RSA)-USS2(RSA,GA2) -
     .               8.D0*WSS(RSA,GA2)*GA2) + YSS(RSA)*EAA)
     .               / (3.D0*RHOA)
C           d(Ec)/d(Nabla(rho)^2) = -Cft^2*(mu*dEcdmu*Fa)*Wss
            DEAADG = -WSS(RSA,GA2)*ECMUC0*CFT*CFT
         ENDIF
      ELSE
         EAA = 0.D0
         EBA = 0.D0
         DEAADR = 0.D0
         DEBADR = 0.D0
         DEAADG = 0.D0
         DEBADG = 0.D0
      ENDIF
C *** RSB - Rs(beta), GB - grad(Rs(beta))
      IF(RHOB.GT.tol_rho) THEN
C        Wigner radius Rs(beta), eq.(8)
         RSB = C2/RHOB13
         CFT = C3/(RHOB*RHOB13)
C        GB  : Nabla(Rs(beta))
C        GB2 : Square of Nabla(Rs(beta))
         GB2 = CFT*CFT*AMB
C        mu(alpha,beta), eq.(13)
C        MUAB = C1*RSB/(KSSP0(RSB)*FSSP(RSB,GB2))**2
         EAB  = 0.D0
         ECMUC0 = 0.D0
         DEABDR = 0.D0
         DEABDG = 0.D0
         DENOM  = (KSSP0(RSB)*FSSP(RSB,GB2))**2
         if((gb2*a2)**2.gt.big) then
            DENOM  = 0d0
         else
            if(rsb.gt.big) then
               DENOM  = (KSSPBIG*FSSP(RSB,GB2))**2
            else
               DENOM  = (KSSP0(RSB)*FSSP(RSB,GB2))**2
            endif
         endif
         IF(C1*RSB .LE. DENOM*CUTOFF) THEN
            MUAB = C1*RSB/DENOM
            EEI  = nwpwxc_FT97_EXPEI(-MUAB)
            EEI1 = MUAB*EEI+1.D0
C           EAB  : Correlation energy density, eqs.(19).
            EAB  = C0*(EEI+2.D0*FF(MUAB)*EEI1)
C           ECMUC0 : Derivative of EAB w.r.t. MUAB
            ECMUC0 = C0*(EEI1*(1.D0+2.D0*(DFFDMU(MUAB)+MUAB*
     .               FF(MUAB)))+2.D0*FF(MUAB)*MUAB*EEI)
C           d(Ec)/d(Nabla(rho)^2)  =  -CFT^2*(mu*dEcdmu)*Wss'
            DEABDG = -WSSP(RSB,GB2)*ECMUC0*CFT*CFT 
C           d(Ec)/d(rho)
            DEABDR = -ECMUC0*(1.D0-USSP1(RSB)/KSSP0(RSB)-
     .               USSP2(RSB,GB2) -
     .               8.D0*WSSP(RSB,GB2)*GB2)
     .               / (3.D0*RHOB)
         ENDIF
C        mu(beta,beta), eq.(33)
C        MUBB = C1*RSB/(KSS0(RSB)*FSS(RSB,GB2))**2
         FA   = FACTOR(RSB)
         EBB  = 0.D0
         ECMUC0 = 0.D0
         DEBBDR = 0.D0
         DEBBDG = 0.D0
         if((gb2*a2)**2.gt.big) then
            DENOM  = 0d0
         else
            if(rsb.gt.big) then
               DENOM  = (KSS0big*FSS(RSB,GB2))**2
            else
               DENOM  = (KSS0(RSB)*FSS(RSB,GB2))**2
            endif
         endif
         IF(C1*RSB .LE. DENOM*CUTOFF) THEN
            MUBB = C1*RSB/DENOM
            EEI  = nwpwxc_FT97_EXPEI(-MUBB)
            EEI1 = MUBB*EEI+1.D0
C           EBB  : Correlation energy density, eqs.(19) and (31).
            EBB  = FA*C0*(EEI+2.D0*FF(MUBB)*EEI1)
C           ECMUC0 : Derivative of EBB w.r.t. MUBB
            ECMUC0 = FA*C0*(EEI1*(1.D0+2.D0*(DFFDMU(MUBB)+MUBB*
     .                FF(MUBB)))+2.D0*FF(MUBB)*MUBB*EEI)
C           d(Ec)/d(rho)
            DEBBDR = (-ECMUC0*(1.D0-USS1(RSB)/KSS0(RSB)
     .                -USS2(RSB,GB2) -
     .               8.D0*WSS(RSB,GB2)*GB2) + YSS(RSB)*EBB)
     .               / (3.D0*RHOB)
C           d(Ec)/d(Nabla(rho)^2) =  -Cft^2*(mu*dEcdmu*Fa)*Wss
            DEBBDG = -WSS(RSB,GB2)*ECMUC0*CFT*CFT
         ENDIF
      ELSE
         EAB = 0.D0
         EBB = 0.D0
         DEABDR = 0.D0
         DEBBDR = 0.D0
         DEABDG = 0.D0
         DEBBDG = 0.D0
      ENDIF
      RETURN
      END
c-----------------------------------------------------------------------
      FUNCTION nwpwxc_FT97_EXPEI(X)
C
C This function program computes approximate values for the
C   function  exp(-x) * Ei(x), where  Ei(x)  is the exponential
C   integral, and  x  is real.
C
C  Author: W. J. Cody
C
C  Latest modification: January 12, 1988
C
      INTEGER INT
      double precision  nwpwxc_FT97_EXPEI, X, RESULT
      INT = 3
      CALL nwpwxc_FT97_CALCEI(X,RESULT,INT)
      nwpwxc_FT97_EXPEI = RESULT
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE nwpwxc_FT97_CALCEI(ARG,RESULT,INT)
C
C  This Fortran 77 packet computes the exponential integrals Ei(x),
C  E1(x), and  exp(-x)*Ei(x)  for real arguments  x  where
C
C           integral (from t=-infinity to t=x) (exp(t)/t),  x > 0,
C  Ei(x) =
C          -integral (from t=-x to t=infinity) (exp(t)/t),  x < 0,
C
C  and where the first integral is a principal value integral.
C  The packet contains three function type subprograms: EI, EONE,
C  and EXPEI;  and one subroutine type subprogram: CALCEI.  The
C  calling statements for the primary entries are
C
C                 Y = EI(X),            where  X .NE. 0,
C
C                 Y = EONE(X),          where  X .GT. 0,
C  and
C                 Y = EXPEI(X),         where  X .NE. 0,
C
C  and where the entry points correspond to the functions Ei(x),
C  E1(x), and exp(-x)*Ei(x), respectively.  The routine CALCEI
C  is intended for internal packet use only, all computations within
C  the packet being concentrated in this routine.  The function
C  subprograms invoke CALCEI with the Fortran statement
C         CALL CALCEI(ARG,RESULT,INT)
C  where the parameter usage is as follows
C
C     Function                  Parameters for CALCEI
C       Call                 ARG             RESULT         INT
C
C      EI(X)              X .NE. 0          Ei(X)            1
C      EONE(X)            X .GT. 0         -Ei(-X)           2
C      EXPEI(X)           X .NE. 0          exp(-X)*Ei(X)    3
C
C  The main computation involves evaluation of rational Chebyshev
C  approximations published in Math. Comp. 22, 641-649 (1968), and
C  Math. Comp. 23, 289-303 (1969) by Cody and Thacher.  This
C  transportable program is patterned after the machine-dependent
C  FUNPACK packet  NATSEI,  but cannot match that version for
C  efficiency or accuracy.  This version uses rational functions
C  that theoretically approximate the exponential integrals to
C  at least 18 significant decimal digits.  The accuracy achieved
C  depends on the arithmetic system, the compiler, the intrinsic
C  functions, and proper selection of the machine-dependent
C  constants.
C
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants
C
C   beta = radix for the floating-point system.
C   minexp = smallest representable power of beta.
C   maxexp = smallest power of beta that overflows.
C   XBIG = largest argument acceptable to EONE; solution to
C          equation:
C                     exp(-x)/x * (1 + 1/x) = beta ** minexp.
C   XINF = largest positive machine number; approximately
C                     beta ** maxexp
C   XMAX = largest argument acceptable to EI; solution to
C          equation:  exp(x)/x * (1 + 1/x) = beta ** maxexp.
C
C     Approximate values for some important machines are:
C
C                           beta      minexp      maxexp
C
C  CRAY-1        (S.P.)       2       -8193        8191
C  Cyber 180/185 
C    under NOS   (S.P.)       2        -975        1070
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)       2        -126         128
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)       2       -1022        1024
C  IBM 3033      (D.P.)      16         -65          63
C  VAX D-Format  (D.P.)       2        -128         127
C  VAX G-Format  (D.P.)       2       -1024        1023
C
C                           XBIG       XINF       XMAX
C
C  CRAY-1        (S.P.)    5670.31  5.45E+2465   5686.21
C  Cyber 180/185 
C    under NOS   (S.P.)     669.31  1.26E+322     748.28
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)      82.93  3.40E+38       93.24
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)     701.84  1.79D+308     716.35
C  IBM 3033      (D.P.)     175.05  7.23D+75      179.85
C  VAX D-Format  (D.P.)      84.30  1.70D+38       92.54
C  VAX G-Format  (D.P.)     703.22  8.98D+307     715.66
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  The following table shows the types of error that may be
C  encountered in this routine and the function value supplied
C  in each case.
C
C       Error       Argument         Function values for
C                    Range         EI      EXPEI     EONE
C
C     UNDERFLOW  (-)X .GT. XBIG     0        -         0
C     OVERFLOW      X .GE. XMAX    XINF      -         -
C     ILLEGAL X       X = 0       -XINF    -XINF     XINF
C     ILLEGAL X      X .LT. 0       -        -     USE ABS(X)
C
C Intrinsic functions required are:
C
C     ABS, SQRT, EXP
C
C
C  Author: W. J. Cody
C          Mathematics abd Computer Science Division
C          Argonne National Laboratory
C          Argonne, IL 60439
C
C  Latest modification: September 9, 1988
C
C----------------------------------------------------------------------
      INTEGER I,INT
      double precision 
     1       A,ARG,B,C,D,EXP40,E,EI,F,FOUR,FOURTY,FRAC,HALF,ONE,P,
     2       PLG,PX,P037,P1,P2,q,QLG,QX,Q1,Q2,R,RESULT,S,SIX,SUMP,
     3       SUMQ,T,THREE,TWELVE,TWO,TWO4,W,X,XBIG,XINF,XMAX,XMX0,
     4       X0,X01,X02,X11,Y,YSQ,ZERO
      DIMENSION  A(7),B(6),C(9),D(9),E(10),F(10),P(10),q(10),R(10),
     1   S(9),P1(10),Q1(9),P2(10),Q2(9),PLG(4),QLG(4),PX(10),QX(10)
C----------------------------------------------------------------------
C  Mathematical constants
C   EXP40 = exp(40)
C   X0 = zero of Ei
C   X01/X11 + X02 = zero of Ei to extra precision
C----------------------------------------------------------------------
      DATA ZERO,P037,HALF,ONE,TWO/0.0D0,0.037D0,0.5D0,1.0D0,2.0D0/,
     1     THREE,FOUR,SIX,TWELVE,TWO4/3.0D0,4.0D0,6.0D0,12.D0,24.0D0/,
     2     FOURTY,EXP40/40.0D0,2.3538526683701998541D17/,
     3     X01,X11,X02/381.5D0,1024.0D0,-5.1182968633365538008D-5/,
     4     X0/3.7250741078136663466D-1/
C----------------------------------------------------------------------
C Machine-dependent constants
C----------------------------------------------------------------------
CS    DATA XINF/3.40E+38/,XMAX/93.246E0/,XBIG/82.93E0/
      DATA XINF/1.79D+308/,XMAX/716.351D0/,XBIG/701.84D0/
C----------------------------------------------------------------------
C Coefficients  for -1.0 <= X < 0.0
C----------------------------------------------------------------------
      DATA A/1.1669552669734461083368D2, 2.1500672908092918123209D3,
     1       1.5924175980637303639884D4, 8.9904972007457256553251D4,
     2       1.5026059476436982420737D5,-1.4815102102575750838086D5,
     3       5.0196785185439843791020D0/
      DATA B/4.0205465640027706061433D1, 7.5043163907103936624165D2,
     1       8.1258035174768735759855D3, 5.2440529172056355429883D4,
     2       1.8434070063353677359298D5, 2.5666493484897117319268D5/
C----------------------------------------------------------------------
C Coefficients for -4.0 <= X < -1.0
C----------------------------------------------------------------------
      DATA C/3.828573121022477169108D-1, 1.107326627786831743809D+1,
     1       7.246689782858597021199D+1, 1.700632978311516129328D+2,
     2       1.698106763764238382705D+2, 7.633628843705946890896D+1,
     3       1.487967702840464066613D+1, 9.999989642347613068437D-1,
     4       1.737331760720576030932D-8/
      DATA D/8.258160008564488034698D-2, 4.344836335509282083360D+0,
     1       4.662179610356861756812D+1, 1.775728186717289799677D+2,
     2       2.953136335677908517423D+2, 2.342573504717625153053D+2,
     3       9.021658450529372642314D+1, 1.587964570758947927903D+1,
     4       1.000000000000000000000D+0/
C----------------------------------------------------------------------
C Coefficients for X < -4.0
C----------------------------------------------------------------------
      DATA E/1.3276881505637444622987D+2,3.5846198743996904308695D+4,
     1       1.7283375773777593926828D+5,2.6181454937205639647381D+5,
     2       1.7503273087497081314708D+5,5.9346841538837119172356D+4,
     3       1.0816852399095915622498D+4,1.0611777263550331766871D03,
     4       5.2199632588522572481039D+1,9.9999999999999999087819D-1/
      DATA F/3.9147856245556345627078D+4,2.5989762083608489777411D+5,
     1       5.5903756210022864003380D+5,5.4616842050691155735758D+5,
     2       2.7858134710520842139357D+5,7.9231787945279043698718D+4,
     3       1.2842808586627297365998D+4,1.1635769915320848035459D+3,
     4       5.4199632588522559414924D+1,1.0D0/
C----------------------------------------------------------------------
C  Coefficients for rational approximation to ln(x/a), |1-x/a| < .1
C----------------------------------------------------------------------
      DATA PLG/-2.4562334077563243311D+01,2.3642701335621505212D+02,
     1         -5.4989956895857911039D+02,3.5687548468071500413D+02/
      DATA QLG/-3.5553900764052419184D+01,1.9400230218539473193D+02,
     1         -3.3442903192607538956D+02,1.7843774234035750207D+02/
C----------------------------------------------------------------------
C Coefficients for  0.0 < X < 6.0,
C  ratio of Chebyshev polynomials
C----------------------------------------------------------------------
      DATA P/-1.2963702602474830028590D01,-1.2831220659262000678155D03,
     1       -1.4287072500197005777376D04,-1.4299841572091610380064D06,
     2       -3.1398660864247265862050D05,-3.5377809694431133484800D08,
     3        3.1984354235237738511048D08,-2.5301823984599019348858D10,
     4        1.2177698136199594677580D10,-2.0829040666802497120940D11/
      DATA q/ 7.6886718750000000000000D01,-5.5648470543369082846819D03,
     1        1.9418469440759880361415D05,-4.2648434812177161405483D06,
     2        6.4698830956576428587653D07,-7.0108568774215954065376D08,
     3        5.4229617984472955011862D09,-2.8986272696554495342658D10,
     4        9.8900934262481749439886D10,-8.9673749185755048616855D10/
C----------------------------------------------------------------------
C J-fraction coefficients for 6.0 <= X < 12.0
C----------------------------------------------------------------------
      DATA R/-2.645677793077147237806D00,-2.378372882815725244124D00,
     1       -2.421106956980653511550D01, 1.052976392459015155422D01,
     2        1.945603779539281810439D01,-3.015761863840593359165D01,
     3        1.120011024227297451523D01,-3.988850730390541057912D00,
     4        9.565134591978630774217D00, 9.981193787537396413219D-1/
      DATA S/ 1.598517957704779356479D-4, 4.644185932583286942650D00,
     1        3.697412299772985940785D02,-8.791401054875438925029D00,
     2        7.608194509086645763123D02, 2.852397548119248700147D01,
     3        4.731097187816050252967D02,-2.369210235636181001661D02,
     4        1.249884822712447891440D00/
C----------------------------------------------------------------------
C J-fraction coefficients for 12.0 <= X < 24.0
C----------------------------------------------------------------------
      DATA P1/-1.647721172463463140042D00,-1.860092121726437582253D01,
     1        -1.000641913989284829961D01,-2.105740799548040450394D01,
     2        -9.134835699998742552432D-1,-3.323612579343962284333D01,
     3         2.495487730402059440626D01, 2.652575818452799819855D01,
     4        -1.845086232391278674524D00, 9.999933106160568739091D-1/
      DATA Q1/ 9.792403599217290296840D01, 6.403800405352415551324D01,
     1         5.994932325667407355255D01, 2.538819315630708031713D02,
     2         4.429413178337928401161D01, 1.192832423968601006985D03,
     3         1.991004470817742470726D02,-1.093556195391091143924D01,
     4         1.001533852045342697818D00/
C----------------------------------------------------------------------
C J-fraction coefficients for  X .GE. 24.0
C----------------------------------------------------------------------
      DATA P2/ 1.75338801265465972390D02,-2.23127670777632409550D02,
     1        -1.81949664929868906455D01,-2.79798528624305389340D01,
     2        -7.63147701620253630855D00,-1.52856623636929636839D01,
     3        -7.06810977895029358836D00,-5.00006640413131002475D00,
     4        -3.00000000320981265753D00, 1.00000000000000485503D00/
      DATA Q2/ 3.97845977167414720840D04, 3.97277109100414518365D00,
     1         1.37790390235747998793D02, 1.17179220502086455287D02,
     2         7.04831847180424675988D01,-1.20187763547154743238D01,
     3        -7.99243595776339741065D00,-2.99999894040324959612D00,
     4         1.99999999999048104167D00/
C----------------------------------------------------------------------
      X = ARG
      IF (X .EQ. ZERO) THEN
            EI = -XINF
            IF (INT .EQ. 2) EI = -EI
         ELSE IF ((X .LT. ZERO) .OR. (INT .EQ. 2)) THEN 
C----------------------------------------------------------------------
C Calculate EI for negative argument or for E1.
C----------------------------------------------------------------------
            Y = ABS(X)
            IF (Y .LE. ONE) THEN
                  SUMP = A(7) * Y + A(1)
                  SUMQ = Y + B(1)
                  DO 110 I = 2, 6
                     SUMP = SUMP * Y + A(I)
                     SUMQ = SUMQ * Y + B(I)
  110             CONTINUE
                  EI = LOG(Y) - SUMP / SUMQ
                  IF (INT .EQ. 3) EI = EI * EXP(Y)
               ELSE IF (Y .LE. FOUR) THEN
                  W = ONE / Y
                  SUMP = C(1)
                  SUMQ = D(1)
                  DO 130 I = 2, 9
                     SUMP = SUMP * W + C(I)
                     SUMQ = SUMQ * W + D(I)
  130             CONTINUE
                  EI = - SUMP / SUMQ
                  IF (INT .NE. 3) EI = EI * EXP(-Y)
               ELSE
                  IF ((Y .GT. XBIG) .AND. (INT .LT. 3)) THEN
                        EI = ZERO
                     ELSE
                        W = ONE / Y
                        SUMP = E(1) 
                        SUMQ = F(1)
                        DO 150 I = 2, 10
                           SUMP = SUMP * W + E(I)
                           SUMQ = SUMQ * W + F(I)
  150                   CONTINUE
                        EI = -W * (ONE - W * SUMP / SUMQ )
                        IF (INT .NE. 3) EI = EI * EXP(-Y)
                  END IF
            END IF
            IF (INT .EQ. 2) EI = -EI
         ELSE IF (X .LT. SIX) THEN
C----------------------------------------------------------------------
C  To improve conditioning, rational approximations are expressed
C    in terms of Chebyshev polynomials for 0 <= X < 6, and in
C    continued fraction form for larger X.
C----------------------------------------------------------------------
            T = X + X
            T = T / THREE - TWO
            PX(1) = ZERO
            QX(1) = ZERO
            PX(2) = P(1)
            QX(2) = q(1)
            DO 210 I = 2, 9
               PX(I+1) = T * PX(I) - PX(I-1) + P(I)
               QX(I+1) = T * QX(I) - QX(I-1) + q(I)
  210       CONTINUE
            SUMP = HALF * T * PX(10) - PX(9) + P(10)
            SUMQ = HALF * T * QX(10) - QX(9) + q(10)
            FRAC = SUMP / SUMQ
            XMX0 = (X - X01/X11) - X02
            IF (ABS(XMX0) .GE. P037) THEN
                  EI = LOG(X/X0) + XMX0 * FRAC
                  IF (INT .EQ. 3) EI = EXP(-X) * EI
               ELSE
C----------------------------------------------------------------------
C Special approximation to  ln(X/X0)  for X close to X0
C----------------------------------------------------------------------
                  Y = XMX0 / (X + X0)
                  YSQ = Y*Y
                  SUMP = PLG(1)
                  SUMQ = YSQ + QLG(1)
                  DO 220 I = 2, 4
                     SUMP = SUMP*YSQ + PLG(I)
                     SUMQ = SUMQ*YSQ + QLG(I)
  220             CONTINUE
                  EI = (SUMP / (SUMQ*(X+X0)) + FRAC) * XMX0
                  IF (INT .EQ. 3) EI = EXP(-X) * EI
            END IF
         ELSE IF (X .LT. TWELVE) THEN
            FRAC = ZERO
            DO 230 I = 1, 9
               FRAC = S(I) / (R(I) + X + FRAC)
  230       CONTINUE
            EI = (R(10) + FRAC) / X
            IF (INT .NE. 3) EI = EI * EXP(X)
         ELSE IF (X .LE. TWO4) THEN
            FRAC = ZERO
            DO 240 I = 1, 9
               FRAC = Q1(I) / (P1(I) + X + FRAC)
  240       CONTINUE
            EI = (P1(10) + FRAC) / X
            IF (INT .NE. 3) EI = EI * EXP(X)
         ELSE
            IF ((X .GE. XMAX) .AND. (INT .LT. 3)) THEN
                  EI = XINF
               ELSE
                  Y = ONE / X
                  FRAC = ZERO
                  DO 250 I = 1, 9
                     FRAC = Q2(I) / (P2(I) + X + FRAC)
  250             CONTINUE
                  FRAC = P2(10) + FRAC
                  EI = Y + Y * Y * FRAC
                  IF (INT .NE. 3) THEN
                        IF (X .LE. XMAX-TWO4) THEN
                              EI = EI * EXP(X)
                           ELSE
C----------------------------------------------------------------------
C Calculation reformulated to avoid premature overflow
C----------------------------------------------------------------------
                              EI = (EI * EXP(X-FOURTY)) * EXP40
                        END IF
                  END IF
            END IF
      END IF
      RESULT = EI
      RETURN
      END
c-----------------------------------------------------------------------
C> @}
