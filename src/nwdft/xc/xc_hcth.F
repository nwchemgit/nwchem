      Subroutine xc_hcth(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                    Amat, Cmat, nq, ipol, Ec, qwght,
     &                    ldew,func)     
c
c$Id: xc_hcth.F,v 1.2 1998-10-13 00:28:21 og845 Exp $
c
      Implicit none
c
      logical ldew
      logical lfac, nlfac
      double precision func(*)
      double precision fac
c
      integer ipol, nq        ! [input]
      double precision tol_rho! [input]
      double precision Ec ! Correlation energy [input/output] 
      double precision rho(nq,ipol*(ipol+1)/2)! Charge Density [input] 
      double precision delrho(nq,3,ipol) ! Charge Density Gradient[input] 
      double precision qwght(nq) ! Quadrature Weights [input]
      double precision Amat(nq,ipol)  !Sampling Matrices for the XC 
      double precision Cmat(nq,3,ipol)!Potential & Energy [output]
c     
c References:
c    F.A.Hamprecht, A.J.Cohen, D.J.Tozer and N.C.Handy, 
c    J. Chem. Phys. 109, 6264-6271 (1998)
c    subroutine supplied by Fred Hamprecht fah@igc.phys.chem.ethz.ch
c
      integer n
      double precision gammaval
c to hcth
      double precision rhoa 
      double precision rhob 
      double precision za   
      double precision zb   
      double precision totalF_xc 
      double precision dfdra     
      double precision dfdrb     
      double precision dfdza     
      double precision dfdzb     
c
      if(ipol.eq.1) then
        do n=1,nq
          rhoa=0.5d0*rho(n,1)
          rhob=rhoa
          za=0.d0
          gammaval=0.25d0*(delrho(n,1,1)*delrho(n,1,1) +
     &         delrho(n,2,1)*delrho(n,2,1) +
     &         delrho(n,3,1)*delrho(n,3,1))
          if(gammaval.gt.tol_rho)
     &         za=sqrt(gammaval)
          zb=za
          call hcth(dfdra, dfdza, dfdrb, dfdzb, 
     1         rhoa, rhob, za, zb, .false., totalF_xc, tol_rho)
          if(ldew) func(n)=func(n)+totalF_xc*fac
          Ec=Ec+totalF_xc*qwght(n)*fac
          Amat(n,1) = Amat(n,1)+dfdra*fac
          dfdza=dfdza*fac*0.5d0
          Cmat(n,1,1) = Cmat(n,1,1) + dfdza*delrho(n,1,1)
          Cmat(n,2,1) = Cmat(n,2,1) + dfdza*delrho(n,2,1)
          Cmat(n,3,1) = Cmat(n,3,1) + dfdza*delrho(n,3,1)
        enddo
      else
        do n=1,nq
          rhoa=rho(n,2)
          rhob=rho(n,3)
          za=0.d0
          gammaval=delrho(n,1,1)*delrho(n,1,1) +
     &         delrho(n,2,1)*delrho(n,2,1) +
     &         delrho(n,3,1)*delrho(n,3,1)
          if(gammaval.gt.tol_rho)
     &         za=sqrt(gammaval)
          zb=0.d0
          gammaval=
     &         delrho(n,1,2)*delrho(n,1,2) +
     &         delrho(n,2,2)*delrho(n,2,2) +
     &         delrho(n,3,2)*delrho(n,3,2) 
          if(gammaval.gt.tol_rho)
     &   zb=sqrt(gammaval)
          call hcth(dfdra, dfdza, dfdrb, dfdzb, 
     1         rhoa, rhob, za, zb, .false., totalF_xc, tol_rho)
          if(ldew) func(n)=func(n)+totalF_xc*fac
          Ec=Ec+totalF_xc*qwght(n)*fac
          Amat(n,1) = Amat(n,1)+dfdra*fac
          Amat(n,2) = Amat(n,2)+dfdrb*fac
          dfdza=dfdza*fac
          dfdzb=dfdzb*fac
          Cmat(n,1,1) = Cmat(n,1,1) + dfdza*delrho(n,1,1)
          Cmat(n,2,1) = Cmat(n,2,1) + dfdza*delrho(n,2,1)
          Cmat(n,3,1) = Cmat(n,3,1) + dfdza*delrho(n,3,1)
          Cmat(n,1,2) = Cmat(n,1,2) + dfdzb*delrho(n,1,2)
          Cmat(n,2,2) = Cmat(n,2,2) + dfdzb*delrho(n,2,2)
          Cmat(n,3,2) = Cmat(n,3,2) + dfdzb*delrho(n,3,2)
        enddo
      endif
      return
      end

      SUBROUTINE hcth(dfdra, dfdza, dfdrb, dfdzb, 
     1                  rhoa, rhob, za, zb, energy, totalF_xc, 
     2     tol_rho)

C    SUPPLIED TO THE ROUTINE:
C    
C    rhoa   -- value of rhoalpha at a given grid point 
C    rhob   -- value of rhobeta at a given grid point
C    za     -- zeta_alpha, as defined in the TH1 paper (JCP 108 2545), 
C              that is mod(grad(rhoalpha)), a scalar quantity.
C    zb     -- mod(grad(rhobeta)) 
C    zab    -- zeta_{alpha beta} as defined in the TH1 paper, that is
C              grad(rhoalpha).grad(rhobeta) 
C    energy -- a boolean variable deciding whether to compute the energy 
C              contribution at the point in space (true) or the
C              appropriate derivatives (false) needed for the KS matrix
C              _and_ the energy contribution.

C    RETURNED FROM THE ROUTINE:

C    totalF_xc --   the contribution to the energy at this point in space.
C    dfdra  -- partial functional derivative of F_xc with respect to 
C              rhoalpha
C    dfdrb  -- partial functional derivative of F_xc with respect to   
C              rhobeta
C    dfdza  -- partial functional derivative of F_xc with respect to   
C              mod(grad(rhoalpha)), divided by za !!!!!!!!!
C              i.e.  1    d f 
C                   --- * ---- 
C                   za    d za
C              This is a consequence of the Cadpac implementation
C    dfdzb  -- partial functional derivative of F_xc with respect to   
C              mod(grad(rhobeta)), divided by zb !!!!!!!!!

C please note that a cutoff of 1.D-10 for the density is hard-wired in the 
C code in several different places; if the density at a point in space is 
C lower than that value, its contribution is neglected. 

      implicit none
      
      double precision rhoa ![input]
      double precision rhob ![input]
      double precision za   ![input]
      double precision zb   ![input]
      double precision totalF_xc ![output]
      double precision tol_rho ![input]
      double precision dfdra     ![output]
      double precision dfdrb     ![output]
      double precision dfdza     ![output]
      double precision dfdzb     ![output]
      integer max_pow_u,n
      PARAMETER(max_pow_u = 4)
      double precision pi
      PARAMETER (PI=3.1415926535898D0)

Cfah max_pow_u is equivalent to "m" in the Becke V paper, that is, the greatest 
Cfah power of u appearing in the power expansion. 
Cfah The "Becke V" paper is: 
Cfah Becke A. D.  Density-functional thermochemistry. V.
Cfah Systematic optimization of exchange-correlation functionals, 
Cfah J. Chem. Phys., 1997, 107, 8554-8560

      double precision sol((max_pow_u+1)*3), F((max_pow_u+1)*3,5), 
     &          FF((max_pow_u+1)*3,5,5),
     &          F_xc((max_pow_u+1)*3)

Cfah sol -- contains the coefficients of the terms in F_xc
Cfah        convention: sol(1) = c_{x alpha, 0}, c_{x beta, 0}
Cfah                    sol(2) = c_{c alpha alpha, 0}, c_{c beta beta, 0} 
Cfah                    sol(3) = c_{c alpha beta, 0} 
Cfah                    sol(4) = c_{x alpha, 1}, c_{x beta, 1}
Cfah                    sol(5) = c_{c alpha alpha, 1}, c_{c beta beta, 1} 
Cfah                    sol(6) = c_{c alpha beta, 1} 
Cfah                           
Cfah                           etc.
Cfah 
Cfah f(5) -- contains the partial first functional derivatives of F_xc with 
Cfah respect to 
Cfah the four quantities (IN THIS ORDER): ra, rb, za, zb
Cfah 
Cfah ff(5,5) contains the second derivatives with
Cfah respect to the same five quantities

Cfah F_xa -- contains the alpha exchange bit containing the various powers 
Cfah         of u_{x alpha} (eq. (18) of Becke V paper) 
Cfah F_xb --              beta       
Cfah            u_{x beta} 
Cfah F_caa -- contains the alpha parallel spin correlation bit with the powers
Cfah          of u_{c alpha alpha} 
Cfah F_cbb --              beta 
Cfah             u_{c beta beta} 
Cfah F_cab -- contains the anti-parallel spin correlation bit with the powers 
Cfah          of u_{c alpha beta} 

Cfah these transformed variables u will be defined and given short-cut names 
Cfah below. 

      LOGICAL energy

C     Initialise

      dfdra  = 0.D0
      dfdrb  = 0.D0
      dfdza  = 0.D0
      dfdzb  = 0.D0
      totalF_xc = 0.D0

      IF ((rhoa + rhob) .LT. 1.D-10) RETURN
Cfah numerical cutoff: if the density is too low, its contribution is 
Cfah neglectable. 

C please refer to these coeffs as THCH1/iterate-e750-g500-v1-m4-n4

      sol( 1) =     0.109320D+01
      sol( 2) =     0.222601D+00
      sol( 3) =     0.729974D+00
      sol( 4) =    -0.744056D+00
      sol( 5) =    -0.338622D-01
      sol( 6) =     0.335287D+01
      sol( 7) =     0.559920D+01
      sol( 8) =    -0.125170D-01
      sol( 9) =    -0.115430D+02
      sol(10) =    -0.678549D+01
      sol(11) =    -0.802496D+00
      sol(12) =     0.808564D+01
      sol(13) =     0.449357D+01
      sol(14) =     0.155396D+01
      sol(15) =    -0.447857D+01

      IF (energy) THEN
        CALL deriv(F,FF,
     &             F_xc,
     &             .FALSE.,
Cfah this logical indicates whether or not the
Cfah derivatives are required
     &             rhoa,rhob,za,zb,
     &                 max_pow_u,
     &             sol,tol_rho)
 

      ELSE

        CALL deriv(F,FF,
     &             F_xc,
     &             .TRUE.,
Cfah this logical indicates whether or not the
Cfah derivatives are required 
     &             rhoa,rhob,za,zb,
     &                 max_pow_u,
     &             sol,tol_rho)


        DO n = 1, (max_pow_u+1)*3 
          dfdra = dfdra + F(n,1) 
          dfdrb = dfdrb + F(n,2) 
          if(za.gt.tol_rho) dfdza = dfdza + F(n,3) / za
          if(zb.gt.tol_rho) dfdzb = dfdzb + F(n,4) / zb  
        ENDDO
Cfah big thanks to NCH: cadpac requires df/(za * dza), NOT 
Cfah                                    df/dza 
        
      ENDIF
        DO n = 0, max_pow_u 
          totalF_xc = totalF_xc + 
     &                F_xc (n*3 + 1) +
     &                F_xc (n*3 + 2) +
     &                F_xc (n*3 + 3)
        ENDDO

      RETURN
      END
Cfah------------------------------------------------------------ 
Cfah the functions below were generated by Mathematica and were not 
Cfah intended to be human-readable. 

      double precision FUNCTION F_xs(
     (     n, c_xs, rhos, zs,tol_rho)
      IMPLICIT none
      double precision rhos, zs,tol_rho
Cfah cheating about the array size here to avoid introduction of 
Cfah another variable: 
      double precision c_xs(1),pi
      integer n
      PARAMETER (PI=3.1415926535898D0)
Cfah this is one term of the exchange bit of F_xc; n is the power of u in that
Cfah term of the power expansion.
      F_xs = (-3.D0*c_xs(1)*
     -    (3.D0/Pi)**
     -     (1.D0/3.D0)*
     -    rhos**(4.D0/3.D0)*
     -    ((0.004D0*zs**2.D0)/
     -       (0.004D0*zs**2.D0 +
     -         rhos**
     -          (8.D0/3.D0)))
     -      **n)/
     -  (2.D0*2.D0**(2.D0/3.D0))

      IF (rhos.lt.tol_rho)  F_xs = 0.D0
Cfah to make safe for the case of beta part of Hydrogen.

      END 

      double precision FUNCTION dF_xs_by_drhos (n, c_xs, rhos, zs)
      double precision c_xs(1)
      IMPLICIT none
      double precision pi, rhos, zs
      PARAMETER (PI=3.1415926535898D0)
      integer n
Cfah  computes the derivative of the term with u^n of the exchange part of
Cfah  F_xc with respect to rho of the same spin.
Cfah  n     -- the power of u involved in this term
Cfah  c_xs  -- the coefficient c_xs(n) of the term of spin s with the
Cfah           power n of u; is NOT passed over as an array.
Cfah  rhos -- rhosigma, that is, either rhoalpha or rhobeta
Cfah  zs    -- mod(grad(rhosigma)), again for alpha or beta
       dF_xs_by_drhos =  -(c_xs(1)*(6.D0/Pi)**
     -      (1.D0/3.D0)*
     -     rhos**(1.D0/3.D0)*
     -     ((0.004D0*zs**2.D0)/
     -        (rhos**
     -           (8.D0/3.D0)
     -           + 0.004D0*zs**2.D0))**n)
     -    + (2.D0*c_xs(1)*0.004D0*n*
     -     (6.D0/Pi)**
     -      (1.D0/3.D0)*
     -     rhos**3.D0*zs**2.D0*
     -     ((0.004D0*zs**2.D0)/
     -        (rhos**
     -           (8.D0/3.D0)
     -           + 0.004D0*zs**2.D0))**
     -      (-1.D0 + n))/
     -   (rhos**(8.D0/3.D0)+
     -      0.004D0*zs**2.D0)**2.D0
      END

      double precision FUNCTION dF_xs_by_dzs (n, c_xs, rhos, zs)
      IMPLICIT none
      double precision c_xs(1),pi,rhos,zs    
      PARAMETER (PI=3.1415926535898D0)
      integer n
Cfah  idem, but derivation with respect to zs
Cfah  see above (function dF_xs_by_drhos) for definition of the
Cfah  other variables
      dF_xs_by_dzs = (-3*c_xs(1)*n*
     -    (3.D0/Pi)**
     -     (1.D0/3.D0)*
     -    rhos**(4.D0/3.D0)*
     -    ((0.004D0*zs**2.D0)/
     -       (rhos**
     -          (8.D0/3.D0)+
     -         0.004D0*zs**2.D0))**
     -     (-1.D0 + n)*
     -    ((-2.D0*0.004D0**2.D0*zs**3.D0)/
     -       (rhos**
     -           (8.D0/3.D0)
     -           + 0.004D0*zs**2.D0)**2.D0
     -       + (2.D0*0.004D0*zs)/
     -       (rhos**
     -          (8.D0/3.D0)+
     -         0.004D0*zs**2.D0)))/
     -  (2.D0*2.D0**(2.D0/3.D0))
      END

      double precision FUNCTION d2F_xs_by_drhos_drhos(
     (     n, c_xs, rhos, zs)
      IMPLICIT none
      double precision c_xs(1),rhos,zs,pi
      PARAMETER (PI=3.1415926535898D0)
      integer n
Cfah  this does the second derivative as specified by the name of the
Cfah  function.
      d2F_xs_by_drhos_drhos = -((c_xs(1)*(2.D0/Pi)**
     -       (1.D0/3.D0)*
     -      ((0.004D0*zs**2.D0)/
     -         (rhos**
     -            (8.D0/3.D0)
     -            + 0.004D0*zs**2.D0))**n
     -       *((1.D0 - 10.D0*n + 16.D0*n**2.D0)*
     -         rhos**
     -          (16.D0/3.D0) +
     -        2*0.004D0*(1.D0 - 13.D0*n)*
     -         rhos**
     -          (8.D0/3.D0)*
     -         zs**2 +
     -        0.004D0**2*zs**4))/
     -    (3**(2.D0/3.D0)*
     -      rhos**
     -       (2.D0/3.D0)*
     -      (rhos**
     -          (8.D0/3.D0) +
     -         0.004D0*zs**2)**2))
      END

      double precision FUNCTION d2F_xs_by_drhos_dzs(
     (     n, c_xs, rhos, zs)
      IMPLICIT none
      double precision c_xs(1),rhos,zs,pi
      PARAMETER (PI=3.1415926535898D0)
      integer n
Cfah  this does the second derivative as specified by the name of the
Cfah  function.
      d2F_xs_by_drhos_dzs = (2.D0*c_xs(1)*n*
     -    (6.D0/Pi)**
     -     (1.D0/3.D0)*
     -    rhos**3.D0*
     -    ((-1.D0 + 2.D0*n)*
     -       rhos**(8.D0/3.D0)
     -        - 3.D0*0.004D0*zs**2.D0)*
     -    ((0.004D0*zs**2.D0)/
     -       (rhos**
     -          (8.D0/3.D0) +
     -         0.004D0*zs**2.D0))**n)/
     -  (zs*(rhos**
     -        (8.D0/3.D0) +
     -       0.004D0*zs**2.D0)**2.D0)
      END

      double precision FUNCTION d2F_xs_by_dzs_dzs(
     (     n, c_xs, rhos, zs)
      implicit none
      double precision c_xs(1),rhos,zs,pi
      PARAMETER (PI=3.1415926535898D0)
      integer n
Cfah  this does the second derivative as specified by the name of the
Cfah  function.
      d2F_xs_by_dzs_dzs = (-3.D0*c_xs(1)*n*
     -    (3.D0/Pi)**
     -     (1.D0/3.D0)*
     -    rhos**4.D0*
     -    ((-1.D0 + 2.D0*n)*
     -       rhos**(8.D0/3.D0)
     -        - 3.D0*0.004D0*zs**2.D0)*
     -    ((0.004D0*zs**2.D0)/
     -       (rhos**
     -          (8.D0/3.D0) +
     -         0.004D0*zs**2.D0))**n)/
     -  (2.D0**(2.D0/3.D0)*
     -    zs**2.D0*
     -    (rhos**
     -        (8.D0/3.D0) +
     -       0.004D0*zs**2.D0)**2.D0)
      END

Cfah------------------------------------------------------------ 

      SUBROUTINE deriv(F,FF,
     &             F_xc,
     &             derivestuff,
Cfah this logical indicates whether or not the 
Cfah derivatives are required (important for cadpac use only)
     &             rhoa,rhob,za,zb,
     &                 max_pow_u,
     &             sol,tol_rho)

      IMPLICIT NONE

      INTEGER max_pow_u
      double precision f_xc((max_pow_u + 1)*3)
      double precision f((max_pow_u+1)*3,5), ff((max_pow_u+1)*3,5,5)
      double precision rhoa, rhob, za, zb,tol_rho

Cfah  COMMON/special/h_atom
Cfah  LOGICAL h_atom
      LOGICAL derivestuff 
      DOUBLE PRECISION sol((max_pow_u+1)*3)
    
      
      DOUBLE PRECISION dF_xa(4)
      DOUBLE PRECISION dF_xb(4)
      DOUBLE PRECISION dF_caa(4)
      DOUBLE PRECISION dF_cbb(4)
      DOUBLE PRECISION dF_cab(4)
Cfah these are the first derivatives of the terms of F_xc with respect to 
Cfah the 4 quantities. the index
Cfah runs over the particular partial derivatives of each term.  
Cfah More explicitly: these are the partial functional derivatives of 
Cfah F_??? with respect to rhoa, rhob, za and zb. 

      DOUBLE PRECISION d2F_xa(4,4)
      DOUBLE PRECISION d2F_xb(4,4)
      DOUBLE PRECISION d2F_caa(4,4)
      DOUBLE PRECISION d2F_cbb(4,4)
      DOUBLE PRECISION d2F_cab(4,4)

Cfah these are the first derivatives of the different transformed variables 
Cfah u with respect to rhoa, rhob, za and zb. These different derivatives 
Cfah with respect to these 4 quantities named above are stored in these 
Cfah arrays.

      DOUBLE PRECISION Pi 
      PARAMETER (Pi = 3.1415926535898D0)
      double precision rho
      DOUBLE PRECISION s_a2, s_b2, s_avg2, u_caa, u_cbb, u_cab
      DOUBLE PRECISION du_caa_by_drhoa, du_caa_by_dza, du_cbb_by_drhob 
      DOUBLE PRECISION du_cbb_by_dzb, du_cab_by_drhoa, du_cab_by_drhob 
      DOUBLE PRECISION du_cab_by_dza, du_cab_by_dzb, du_caa_by_drhoa_dza 
      DOUBLE PRECISION du_caa_by_dza_dza, du_cbb_by_dzb_dzb
      DOUBLE PRECISION du_cbb_by_drhob_dzb, du_cab_by_drhoa_dza 
      DOUBLE PRECISION du_cab_by_drhoa_dzb, du_cab_by_drhob_dza 
      DOUBLE PRECISION du_cab_by_drhob_dzb, du_cab_by_dza_dza, 
     ,du_cab_by_dza_dzb 
      DOUBLE PRECISION du_cab_by_dzb_dzb 
      DOUBLE PRECISION rsa, rsa12, rsa32, rsa21, rsb, 
     ,rsb12, rsb32, rsb21 
      DOUBLE PRECISION rsab, rsab12, rsab32, rsab21 
      DOUBLE PRECISION drsa_by_drhoa, drsb_by_drhob, drsab_by_drhoa
      DOUBLE PRECISION drsab_by_drhob 
      DOUBLE PRECISION zeta, dzeta_by_drhoa, dzeta_by_drhob 
      DOUBLE PRECISION fzeta, dfzeta_by_dzeta, 
     ,     e_crsa1, e_crsb1
      DOUBLE PRECISION e_crsab1, e_crsab0, a_crsab 
      DOUBLE PRECISION e_crsabzeta, de_crsa1_by_drsa, de_crsb1_by_drsb 
      DOUBLE PRECISION da_crsab_by_drsab, de_crsab0_by_drsab 
      DOUBLE PRECISION de_crsab1_by_drsab, de_crsabzeta_by_drsab 
      DOUBLE PRECISION de_crsabzeta_by_dzeta, e_caa, e_cbb, e_cab, 
     & de_caa_by_drhoa, de_cbb_by_drhob, de_cab_by_drhoa, 
     & de_cab_by_drhob,
     & c_naa, c_nbb, c_nab
      DOUBLE PRECISION F_xs ! this is a function which is called. 
      DOUBLE PRECISION dF_xs_by_drhos, dF_xs_by_dzs, 
     & d2F_xs_by_drhos_dzs, d2F_xs_by_dzs_dzs


      INTEGER i, j, k, n


      DO j = 1, 4
        DO n = 1, (max_pow_u+1)*3
          F(n,j) = 0.D0
Cfah  later on, n has a different meaning: n as power of u, not 
Cfah  as number of the coefficient. 
        ENDDO
        dF_xa(j) = 0.D0
        dF_xb(j) = 0.D0
        dF_caa(j) = 0.D0
        dF_cbb(j) = 0.D0
        dF_cab(j) = 0.D0
        DO k = 1, 4
          DO n = 1, (max_pow_u+1)*3
            FF(n,j,k) = 0.D0
          ENDDO
          d2F_xa(j,k) = 0.D0
          d2F_xb(j,k) = 0.D0
          d2F_caa(j,k) = 0.D0
          d2F_cbb(j,k) = 0.D0
          d2F_cab(j,k) = 0.D0
        ENDDO
      ENDDO
      DO j = 1, (max_pow_u+1)*3
        F_xc(j) = 0.D0
      ENDDO 

Cfah --------------------------------------------------------------

Cfah call the expensive correlation parts here just once, and store their
Cfah values in a temporary variable. Then compute the actual F_c derivatives
Cfah with the various powers of u.  

      rho = rhoa + rhob

      s_a2=0.d0
      if(za.gt.tol_rho) s_a2 = za**2.D0 / rhoa**(8.D0/3.D0)
      s_b2=0.d0
      if(zb.gt.tol_rho) s_b2 = zb**2.D0 / rhob**(8.D0/3.D0)
      s_avg2 = 0.5D0*(s_a2 + s_b2)

      u_caa = 0.2D0*s_a2/(1.D0+0.2D0*s_a2) 
      u_cbb = 0.2D0*s_b2/(1.D0+0.2D0*s_b2) 
      u_cab = 0.006D0*s_avg2/(1.d0+0.006D0*s_avg2)

      rsa = ((3.d0/Pi)**(1.D0/3.D0)*
     -    (1.d0/rhoa)**(1.D0/3.D0))/
     -  2**(2.D0/3.D0)
      rsa12 = rsa**(1.D0/2.D0)
      rsa32 = rsa**(3.D0/2.D0)
      rsa21 = rsa**2.D0

      rsb = ((3.d0/Pi)**(1.D0/3.D0)*
     -    (1.d0/rhob)**(1.D0/3.D0))/
     -  2**(2.D0/3.D0)
      rsb12 = rsb**(1.D0/2.D0)
      rsb32 = rsb**(3.D0/2.D0)
      rsb21 = rsb**2.D0

      rsab = ((3.d0/Pi)**(1.D0/3.D0)*
     -    (1.d0/rho)**(1.D0/3.D0))/
     -  2**(2.D0/3.D0)
      rsab12 = rsab**(1.D0/2.D0)
      rsab32 = rsab**(3.D0/2.D0)
      rsab21 = rsab**2.D0

      zeta = (rhoa-rhob)/rho

      fzeta = (-2.d0 + (1.d0 - zeta)**
     -     (4.D0/3.D0) +
     -    (1.d0 + zeta)**(4.D0/3.D0))/
     -  (-2.d0 + 2*2**(1.D0/3.D0))

      e_crsa1 = -0.03108999999999999d0*
     -  dlog(1 + 32.16468317787069d0/
     -     (14.1189d0*rsa12 +
     -       6.1977d0*rsa + 3.3662d0*rsa32 +
     -       0.6251699999999999d0*rsa21))*
     -  (1.d0 + 0.20548d0*rsa)

      e_crsb1 = -0.03108999999999999d0*
     -  dlog(1.d0 + 32.16468317787069d0/
     -     (14.1189d0*rsb12 +
     -       6.1977d0*rsb + 3.3662d0*rsb32 +
     -       0.6251699999999999d0*rsb21))*
     -  (1.d0 + 0.20548d0*rsb)

      e_crsab1 = -0.03108999999999999d0*
     -  dlog(1.d0 + 32.16468317787069d0/
     -     (14.1189d0*rsab12 +
     -       6.1977d0*rsab + 3.3662d0*rsab32 +
     -       0.6251699999999999d0*rsab21))*
     -  (1.d0 + 0.20548d0*rsab)

      e_crsab0 = -0.062182d0*dlog(1.d0 +
     -    16.0818243221511d0/
     -     (7.595699999999999d0*rsab12 +
     -       3.5876d0*rsab +
     -       1.6382d0*rsab32 +
     -       0.49294d0*rsab21))*
     -  (1.d0 + 0.2137d0*rsab)

      a_crsab = 0.03377399999999999d0*
     -  dlog(1.d0 + 29.60857464321667/
     -     (10.35699999999999d0*rsab12 +
     -       3.623099999999999d0*rsab +
     -       0.88026d0*rsab32 +
     -       0.49671d0*rsab21))*
     -  (1.d0 + 0.11125d0*rsab)

      e_crsabzeta = e_crsab0+a_crsab*fzeta*(1.d0-zeta**4)/1.709921D0+
     -  (e_crsab1-e_crsab0)*fzeta*zeta**4
  
Cfah       print*, 'deriv2', e_crsab0, a_crsab, fzeta, zeta, e_crsab1

      e_caa = rhoa*e_crsa1
      e_cbb = rhob*e_crsb1
      e_cab = rho*e_crsabzeta - rhoa*e_crsa1 - rhob*e_crsb1

Cfah       print*, 'deriv1:',rho, e_crsabzeta, e_crsa1, e_crsb1, rhoa, rhob


Cfah derive all this stuff only if we need it, not if only F_c itself 
Cfah is required. 
      IF (derivestuff) THEN  

      du_caa_by_drhoa = (-8*0.2D0*za**2*rhoa**(5.D0/3.D0))/
     -  (3.*(0.2D0*za**2 + 
     -       rhoa**(8.D0/3.D0))**2)
      du_caa_by_dza = (2*0.2D0*za*rhoa**(8.D0/3.D0))/
     -  (0.2D0*za**2 + rhoa**(8.D0/3.D0))**2

      du_cbb_by_drhob = (-8*0.2D0*zb**2*rhob**(5.D0/3.D0))/
     -  (3.d0*(0.2D0*zb**2 + 
     -       rhob**(8.D0/3.D0))**2)
      du_cbb_by_dzb = (2*0.2D0*zb*rhob**(8.D0/3.D0))/
     -  (0.2D0*zb**2 + rhob**(8.D0/3.D0))**2

      du_cab_by_drhoa = (-16*0.006D0*za**2*rhoa**(5.D0/3.D0)*
     -    rhob**(16.D0/3.D0))/
     -  (3.d0*(0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
     -       0.006D0*za**2*rhob**(8.D0/3.D0) + 
     -       2.d0*rhoa**(8.D0/3.D0)*
     -        rhob**(8.D0/3.D0))**2) 
      du_cab_by_drhob = (-16*0.006D0*zb**2*rhob**(5.D0/3.D0)*
     -    rhoa**(16.D0/3.D0))/
     -  (3.d0*(0.006D0*za**2*rhob**(8.D0/3.D0) + 
     -       0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
     -       2.d0*rhob**(8.D0/3.D0)*
     -        rhoa**(8.D0/3.D0))**2) 
      du_cab_by_dza = (4*0.006D0*za*rhoa**(8.D0/3.D0)*
     -    rhob**(16.D0/3.D0))/
     -  (0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
     -     0.006D0*za**2*rhob**(8.D0/3.D0) + 
     -     2.d0*rhoa**(8.D0/3.D0)*
     -      rhob**(8.D0/3.D0))**2
      du_cab_by_dzb = (4*0.006D0*zb*rhoa**(16.D0/3.D0)*
     -    rhob**(8.D0/3.D0))/
     -  (0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
     -     0.006D0*za**2*rhob**(8.D0/3.D0) + 
     -     2.d0*rhoa**(8.D0/3.D0)*
     -      rhob**(8.D0/3.D0))**2

Cfah Second derivatives are not required by cadpac. 
Cfah   du_caa_by_drhoa_dza = (16*0.2D0*za*rhoa**(5.D0/3.D0)*
Cfah -    (0.2D0*za**2 - rhoa**(8.D0/3.D0)))/
Cfah -  (3.*(0.2D0*za**2 + 
Cfah -       rhoa**(8.D0/3.D0))**3)
Cfah   du_cbb_by_drhob_dzb = (16*0.2D0*zb*rhob**(5.D0/3.D0)*
Cfah -    (0.2D0*zb**2 - rhob**(8.D0/3.D0)))/
Cfah -  (3.*(0.2D0*zb**2 + 
Cfah -       rhob**(8.D0/3.D0))**3)
Cfah
Cfah   du_caa_by_dza_dza = (2*0.2D0*rhoa**(8.D0/3.D0)*
Cfah -    (-3*0.2D0*za**2 + rhoa**(8.D0/3.D0))
Cfah -    )/
Cfah -  (0.2D0*za**2 + rhoa**(8.D0/3.D0))**3
Cfah   du_cbb_by_dzb_dzb = (2*0.2D0*rhob**(8.D0/3.D0)*
Cfah -    (-3*0.2D0*zb**2 + rhob**(8.D0/3.D0))
Cfah -    )/
Cfah -  (0.2D0*zb**2 + rhob**(8.D0/3.D0))**3
Cfah
Cfah   du_cab_by_drhoa_dza = (-32*0.006D0*rhoa**(5.D0/3.D0)*
Cfah -    (0.006D0*za*zb**2*
Cfah -       rhoa**(8.D0/3.D0)*
Cfah -       rhob**(16.D0/3.D0) - 
Cfah -      0.006D0*za**3*rhob**8 + 
Cfah -      2*za*rhoa**(8.D0/3.D0)*rhob**8))/
Cfah -  (3.*(0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
Cfah -       0.006D0*za**2*rhob**(8.D0/3.D0) + 
Cfah -       2*rhoa**(8.D0/3.D0)*
Cfah -        rhob**(8.D0/3.D0))**3) 
Cfah   du_cab_by_drhoa_dzb = (64*0.006D0**2*za**2*zb*
Cfah -    rhoa**(13.D0/3.D0)*
Cfah -    rhob**(16.D0/3.D0))/
Cfah -  (3.*(0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
Cfah -       0.006D0*za**2*rhob**(8.D0/3.D0) + 
Cfah -       2*rhoa**(8.D0/3.D0)*
Cfah -        rhob**(8.D0/3.D0))**3) 
Cfah   du_cab_by_drhob_dza = (64*0.006D0**2*za*zb**2*
Cfah -    rhoa**(16.D0/3.D0)*
Cfah -    rhob**(13.D0/3.D0))/
Cfah -  (3.*(0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
Cfah -       0.006D0*za**2*rhob**(8.D0/3.D0) + 
Cfah -       2*rhoa**(8.D0/3.D0)*
Cfah -        rhob**(8.D0/3.D0))**3) 
Cfah   du_cab_by_drhob_dzb = (-32*0.006D0*rhob**(5.D0/3.D0)*
Cfah -    (-(0.006D0*zb**3*rhoa**8) + 
Cfah -      0.006D0*za**2*zb*
Cfah -       rhoa**(16.D0/3.D0)*
Cfah -       rhob**(8.D0/3.D0) + 
Cfah -      2*zb*rhoa**8*rhob**(8.D0/3.D0)))/
Cfah -  (3.*(0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
Cfah -       0.006D0*za**2*rhob**(8.D0/3.D0) + 
Cfah -       2*rhoa**(8.D0/3.D0)*
Cfah -        rhob**(8.D0/3.D0))**3) 
Cfah   du_cab_by_dza_dza = (4*0.006D0*(0.006D0*zb**2*
Cfah -       rhoa**(16.D0/3.D0)*
Cfah -       rhob**(16.D0/3.D0) - 
Cfah -      3*0.006D0*za**2*rhoa**(8.D0/3.D0)*
Cfah -       rhob**8 + 2*rhoa**(16.D0/3.D0)*rhob**8
Cfah -      ))/
Cfah -  (0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
Cfah -     0.006D0*za**2*rhob**(8.D0/3.D0) + 
Cfah -     2*rhoa**(8.D0/3.D0)*
Cfah -      rhob**(8.D0/3.D0))**3 
Cfah   du_cab_by_dza_dzb = (-16*0.006D0**2*za*zb*
Cfah -    rhoa**(16.D0/3.D0)*
Cfah -    rhob**(16.D0/3.D0))/
Cfah -  (0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
Cfah -     0.006D0*za**2*rhob**(8.D0/3.D0) + 
Cfah -     2*rhoa**(8.D0/3.D0)*
Cfah -      rhob**(8.D0/3.D0))**3 
Cfah   du_cab_by_dzb_dzb = (4*0.006D0*rhoa**(16.D0/3.D0)*
Cfah -    rhob**(8.D0/3.D0)*
Cfah -    (-3*0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
Cfah -      0.006D0*za**2*rhob**(8.D0/3.D0) + 
Cfah -      2*rhoa**(8.D0/3.D0)*
Cfah -       rhob**(8.D0/3.D0)))/
Cfah -  (0.006D0*zb**2*rhoa**(8.D0/3.D0) + 
Cfah -     0.006D0*za**2*rhob**(8.D0/3.D0) + 
Cfah -     2*rhoa**(8.D0/3.D0)*
Cfah -      rhob**(8.D0/3.D0))**3 

      drsa_by_drhoa = -((1.d0/rhoa)**(4.D0/3.D0)/
     -    (6.d0**(2.D0/3.D0)*Pi**(1.D0/3.D0)))
      drsb_by_drhob = -((1.d0/rhob)**(4.D0/3.D0)/
     -    (6.d0**(2.D0/3.D0)*Pi**(1.D0/3.D0)))
      drsab_by_drhoa = -((1.d0/rho)**(4.D0/3.D0)/
     -    (6.d0**(2.D0/3.D0)*Pi**(1.D0/3.D0))) 
      drsab_by_drhob = drsab_by_drhoa 

      dzeta_by_drhoa = 2.d0*rhob/rho**2
      dzeta_by_drhob = -2.d0*rhoa/rho**2

      dfzeta_by_dzeta = ((-4d0*(1.d0 - zeta)**(1.D0/3.D0))/
     -     3.d0 + (4d0*(1.d0 + zeta)**
     -        (1.D0/3.D0))/3.)/
     -  (-2.d0 + 2d0*2**(1.D0/3.D0))

      de_crsa1_by_drsa = (1.d0*(1.d0 + 0.20548d0*rsa)*
     -     (6.1977 + 7.05945/rsa12 + 1.25034d0*rsa + 
     -       5.0493d0*rsa12))/
     -   ((6.1977d0*rsa + 14.1189d0*rsa12 + 
     -        0.6251699999999999d0*rsa21 + 3.3662d0*rsa32)**
     -      2d0*(1.d0 + 32.16468317787069/
     -        (6.1977d0*rsa + 14.1189d0*rsa12 + 
     -          0.6251699999999999d0*rsa21 + 3.3662d0*rsa32)
     -       )) - 0.006388373199999999d0*
     -   dlog(1.d0 + 32.16468317787069/
     -      (6.1977d0*rsa + 14.1189d0*rsa12 + 
     -        0.6251699999999999d0*rsa21 + 3.3662d0*rsa32)) 

      de_crsb1_by_drsb = (1.d0*(1.d0 + 0.20548d0*rsb)*
     -     (6.1977 + 7.05945/rsb12 + 1.25034d0*rsb + 
     -       5.0493d0*rsb12))/
     -   ((6.1977d0*rsb + 14.1189d0*rsb12 + 
     -        0.6251699999999999d0*rsb21 + 3.3662d0*rsb32)**
     -      2d0*(1.d0 + 32.16468317787069d0/
     -        (6.1977d0*rsb + 14.1189d0*rsb12 + 
     -          0.6251699999999999d0*rsb21 + 3.3662d0*rsb32)
     -       )) - 0.006388373199999999d0*
     -   dlog(1.d0 + 32.16468317787069d0/
     -      (6.1977d0*rsb + 14.1189d0*rsb12 + 
     -        0.6251699999999999d0*rsb21 + 3.3662d0*rsb32)) 

      da_crsab_by_drsab = (-1.d0*(1.d0 + 0.11125d0*rsab)*
     -     (3.623099999999999d0 + 
     -       5.178499999999999d0/rsab12 + 
     -       0.99342d0*rsab + 1.32039d0*rsab12))/
     -   ((1.d0 + 29.60857464321667d0/
     -        (3.623099999999999d0*rsab + 
     -          10.35699999999999d0*rsab12 + 
     -          0.49671d0*rsab21 + 0.88026d0*rsab32))*
     -     (3.623099999999999d0*rsab + 
     -        10.35699999999999d0*rsab12 + 
     -        0.49671d0*rsab21 + 0.880260*rsab32)**2) + 
     -  0.003757357499999999d0*
     -   dlog(1.d0 + 29.60857464321667/
     -      (3.623099999999999d0*rsab + 
     -        10.35699999999999d0*rsab12 + 
     -        0.49671d0*rsab21 + 0.88026d0*rsab32)) 

      de_crsab0_by_drsab = (1.d0*(1.d0 + 0.2137d0*rsab)*
     -     (3.5876d0 + 3.797849999999999d0/rsab12 + 
     -       0.98588d0*rsab + 2.4573*rsab12))/
     -   ((3.5876d0*rsab + 7.595699999999999d0*rsab12 + 
     -        0.49294d0*rsab21 + 1.6382d0*rsab32)**2*
     -     (1.d0 + 16.0818243221511d0/
     -        (3.5876d0*rsab + 7.595699999999999d0*rsab12 + 
     -          0.49294d0*rsab21 + 1.6382d0*rsab32))) - 
     -  0.01328829339999999d0*
     -   dlog(1.d0 + 16.0818243221511d0/
     -      (3.5876*rsab + 7.595699999999999d0*rsab12 + 
     -        0.49294d0*rsab21 + 1.6382d0*rsab32))

      de_crsab1_by_drsab = (1.*(1.d0 + 0.20548d0*rsab)*
     -     (6.1977 + 7.05945d0/Sqrt(rsab) + 
     -       1.25034d0*rsab + 5.0493d0*rsab12))/
     -   ((6.1977*rsab + 14.1189d0*rsab12 + 
     -        0.6251699999999999d0*rsab21 + 3.3662d0*rsab32)
     -       **2*(1.d0 + 32.16468317787069d0/
     -        (6.1977*rsab + 14.1189*rsab12 + 
     -          0.6251699999999999d0*rsab21 + 
     -          3.3662d0*rsab32))) - 
     -  0.006388373199999999d0*
     -   dlog(1.d0 + 32.16468317787069d0/
     -      (6.1977d0*rsab + 14.1189d0*rsab12 + 
     -        0.6251699999999999d0*rsab21 + 3.3662*rsab32)
     -     )

      de_crsabzeta_by_drsab = 1.124999956683108D0*(1.d0 - zeta**4)*
     -   (-2.d0 + (1.d0 - zeta)**(4.D0/3.D0) + 
     -     (1.d0 + zeta)**(4.D0/3.D0))*
     -   da_crsab_by_drsab + 
     -  de_crsab0_by_drsab + 
     -  fzeta*zeta**4*
     -   (- de_crsab0_by_drsab + 
     -      de_crsab1_by_drsab  )

      de_crsabzeta_by_dzeta = 1.499999942244144D0*(-1.d0 + zeta**4)*
     -   ((1.d0 - zeta)**(1.D0/3.D0) - 
     -     (1.d0 + zeta)**(1.D0/3.D0))*a_crsab - 
     -  4.499999826732434D0*zeta**3*
     -   (-2.d0 + (1.d0 - zeta)**(4.D0/3.D0) + 
     -     (1.d0 + zeta)**(4.D0/3.D0))*a_crsab + 
     -  (2*zeta**4*((1.d0 - zeta)**(1.D0/3.D0) - 
     -       (1.d0 + zeta)**(1.D0/3.D0))*
     -     (e_crsab0 - e_crsab1))/
     -   (3.*(-1.d0 + 2**(1.D0/3.D0))) + 
     -  4*fzeta*(-e_crsab0 + 
     -     e_crsab1)*zeta**3

Cfah this is with application of the chain rule; I keep it that general
Cfah because this way, I only have to define one "G". 
      de_caa_by_drhoa = e_crsa1 + rhoa*de_crsa1_by_drsa*drsa_by_drhoa 
      de_cbb_by_drhob = e_crsb1 + rhob*de_crsb1_by_drsb*drsb_by_drhob 

      de_cab_by_drhoa = -e_crsa1 + 
     -  e_crsabzeta - 
     -  rhoa*de_crsa1_by_drsa* 
     -  drsa_by_drhoa + 
     -  rho*(de_crsabzeta_by_drsab*
     -  drsab_by_drhoa + 
     -  de_crsabzeta_by_dzeta*
     -  dzeta_by_drhoa)

      de_cab_by_drhob = -e_crsb1 + 
     -  e_crsabzeta - 
     -  rhob*de_crsb1_by_drsb* 
     -  drsb_by_drhob + 
     -  rho*(de_crsabzeta_by_drsab*
     -  drsab_by_drhob + 
     -  de_crsabzeta_by_dzeta*
     -  dzeta_by_drhob)

      ENDIF ! if (derivestuff) 


Cfah Here starts the big outer loop over the powers u 

      DO n = 0, max_pow_u 
        c_naa = sol((n*3) + 2)
        c_nbb = c_naa
        c_nab = sol((n*3) + 3) 

Cfah construction of the F_xc itself
Cfah -------------------------------
        IF (rhoa.GT.tol_rho) THEN
          F_xc(n*3+1) = F_xs (n, sol((n*3) + 1), rhoa, za, 
     ,         tol_rho)
          F_xc(n*3+2) = e_caa*u_caa**n*c_naa
        ENDIF

        IF (rhob.GT.tol_rho) THEN
          F_xc(n*3+1) = F_xc(n*3+1)+F_xs(n, sol((n*3) + 1), rhob, zb,
     ,         tol_rho)
          F_xc(n*3+2) = F_xc(n*3+2)+e_cbb*u_cbb**n*c_nbb
        ENDIF

        IF (rhoa.GT.tol_rho .AND. rhob.GT.tol_rho) THEN
          F_xc(n*3+3) = e_cab*u_cab**n*c_nab
        ENDIF

Cfah       print*, 'in deriv:', e_cab, u_cab, c_nab

Cfah    First Derivatives
Cfah ---------------------

        if(za.gt.tol_rho) 
     *       dF_xa(1) = dF_xs_by_drhos (n, sol((n*3) + 1), rhoa, za) 
        dF_xa(2) = 0 
        if(za.gt.tol_rho) 
     *        dF_xa(3) = dF_xs_by_dzs (n, sol((n*3) + 1), rhoa, za)
        dF_xa(4) = 0 

        dF_xb(1) = 0 
        if(zb.gt.tol_rho) 
     *        dF_xb(2) = dF_xs_by_drhos (n, sol((n*3) + 1), rhob, zb)
        dF_xb(3) = 0 
        if(zb.gt.tol_rho) 
     *        dF_xb(4) = dF_xs_by_dzs (n, sol((n*3) + 1), rhob, zb)

        dF_caa(1) = c_naa*u_caa**n*de_caa_by_drhoa +  
     -   c_naa*n*e_caa*u_caa**(-1 + n)*
     -   du_caa_by_drhoa
        dF_caa(2) = 0.D0
        dF_caa(3) = c_naa*n*e_caa*u_caa**(-1 + n)*
     -  du_caa_by_dza
        dF_caa(4) = 0.D0

        dF_cbb(1) = 0.D0
        dF_cbb(2) = c_nbb*u_cbb**n*de_cbb_by_drhob + 
     -   c_nbb*n*e_cbb*u_cbb**(-1 + n)*
     -   du_cbb_by_drhob
        dF_cbb(3) = 0.D0
        dF_cbb(4) = c_nbb*n*e_cbb*u_cbb**(-1 + n)*
     -  du_cbb_by_dzb


        dF_cab(1) = c_nab*u_cab**n*
     -   de_cab_by_drhoa + 
     -   c_nab*n*e_cab*
     -   u_cab**(-1 + n)*
     -   du_cab_by_drhoa 
        dF_cab(2) = c_nab*u_cab**n*
     -   de_cab_by_drhob +
     -   c_nab*n*e_cab*
     -   u_cab**(-1 + n)*
     -   du_cab_by_drhob
        dF_cab(3) = c_nab*n*e_cab*
     -  u_cab**(-1 + n)*
     -  du_cab_by_dza
        dF_cab(4) = c_nab*n*e_cab*
     -  u_cab**(-1 + n)*
     -  du_cab_by_dzb

Cfah Second Derivatives
Cfah ------------------

Cfah         d2F_xa(1,1) = d2F_xs_by_drhos_drhos (n, sol((n*3) + 1), 
Cfah      &                                         rhoa, za)
Cfah see comment below, for the (2,2) term. 
Cfah    d2F_xa(1,2) = 0 
Cfah    d2F_xa(1,3) = d2F_xs_by_drhos_dzs (n, sol((n*3) + 1), rhoa, za)
Cfah    d2F_xa(1,4) = 0 
Cfah    d2F_xa(2,2) = 0  
Cfah    d2F_xa(2,3) = 0  
Cfah    d2F_xa(2,4) = 0 
Cfah    d2F_xa(3,3) = d2F_xs_by_dzs_dzs (n, sol((n*3) + 1), rhoa, za)
Cfah    d2F_xa(3,4) = 0 
Cfah    d2F_xa(4,4) = 0 

Cfah for alpha spin, elements are non-zero when both indices are odd; 
Cfah for beta spin, elements are non-zero when both indices are even. 
Cfah the matrix is symmetric, and the upper triangle contains the 
Cfah 10 elements given above and below. 

Cfah    d2F_xb(1,1) = 0
Cfah    d2F_xb(1,2) = 0  
Cfah    d2F_xb(1,3) = 0           
Cfah    d2F_xb(1,4) = 0 
Cfah        d2F_xb(2,2) = d2F_xs_by_drhos_drhos (n, sol((n*3) + 1), 
Cfah     &                                         rhob, zb)
Cfah this term is NOT zero, but needs not be evaluated since we don't 
Cfah need it for the construction of v (cf. routine "va" in the fit 
Cfah program) 
Cfah    d2F_xb(2,3) = 0  
Cfah    d2F_xb(2,4) = d2F_xs_by_drhos_dzs (n, sol((n*3) + 1), rhob, zb)
Cfah    d2F_xb(3,3) = 0
Cfah    d2F_xb(3,4) = 0 
Cfah    d2F_xb(4,4) = d2F_xs_by_dzs_dzs (n, sol((n*3) + 1), rhob, zb)


Cfah    d2F_caa(1,1) = !=0, but not needed 
Cfah    d2F_caa(1,2) = 0.D0 (not needed)  
Cfah    d2F_caa(1,3) = c_naa*n*u_caa**(-1 + n)*
Cfah -   de_caa_by_drhoa*
Cfah -   du_caa_by_dza + 
Cfah -  c_naa*(-1 + n)*n*e_caa*
Cfah -   u_caa**(-2 + n)*
Cfah -   du_caa_by_dza*
Cfah -   du_caa_by_drhoa + 
Cfah -  c_naa*n*e_caa*u_caa**(-1 + n)*
Cfah -   du_caa_by_drhoa_dza
Cfah    d2F_caa(1,4) = 0.D0
Cfah    d2F_caa(2,2) = 0.D0 (not needed)
Cfah    d2F_caa(2,3) = 0.D0 
Cfah    d2F_caa(2,4) = 0.D0 
Cfah    d2F_caa(3,3) = c_naa*n*e_caa*u_caa**(-2 + n)*
Cfah -  ((-1 + n)*du_caa_by_dza**
Cfah -      2 + u_caa*
Cfah -     du_caa_by_dza_dza)
Cfah    d2F_caa(3,4) = 0.D0 
Cfah    d2F_caa(4,4) = 0.D0 


Cfah    d2F_cbb(1,1) = 0.D0 (not needed)
Cfah    d2F_cbb(1,2) = 0.D0 (not needed)
Cfah    d2F_cbb(1,3) = 0.D0
Cfah    d2F_cbb(1,4) = 0.D0
Cfah    d2F_cbb(2,2) = !=0, but not needed 
Cfah    d2F_cbb(2,3) = 0.D0
Cfah    d2F_cbb(2,4) = c_nbb*n*u_cbb**(-1 + n)*
Cfah -   de_cbb_by_drhob*
Cfah -   du_cbb_by_dzb +
Cfah -  c_nbb*(-1 + n)*n*e_cbb*
Cfah -   u_cbb**(-2 + n)*
Cfah -   du_cbb_by_dzb*
Cfah -   du_cbb_by_drhob +
Cfah -  c_nbb*n*e_cbb*u_cbb**(-1 + n)*
Cfah -   du_cbb_by_drhob_dzb
Cfah    d2F_cbb(3,3) = 0.D0
Cfah    d2F_cbb(3,4) = 0.D0
Cfah    d2F_cbb(4,4) =  c_nbb*n*e_cbb*u_cbb**(-2 + n)*
Cfah -  ((-1 + n)*du_cbb_by_dzb**
Cfah -      2 + u_cbb*
Cfah -     du_cbb_by_dzb_dzb)

Cfah    d2F_cab(1,1) = not needed
Cfah    d2F_cab(1,2) = not needed
Cfah    d2F_cab(1,3) = c_nab*n*u_cab**(-2 + n)*
Cfah -  ((-1 + n)*e_cab*
Cfah -     du_cab_by_dza
Cfah -      *du_cab_by_drhoa + 
Cfah -    u_cab*(de_cab_by_drhoa*
Cfah -        du_cab_by_dza + e_cab*
Cfah -        du_cab_by_drhoa_dza
Cfah -  ))
Cfah    d2F_cab(1,4) = c_nab*n*u_cab**(-2 + n)*
Cfah -  (  (-1 + n)*e_cab*
Cfah -     du_cab_by_dzb
Cfah -      *du_cab_by_drhoa + 
Cfah -    u_cab*
Cfah -     (de_cab_by_drhoa*
Cfah -        du_cab_by_dzb + 
Cfah -       e_cab*
Cfah -        du_cab_by_drhoa_dzb))
Cfah    d2F_cab(2,2) = not needed
Cfah    d2F_cab(2,3) = c_nab*n*u_cab**(-2 + n)*
Cfah -  ((-1 + n)*e_cab*
Cfah -     du_cab_by_dza
Cfah -      *du_cab_by_drhob +
Cfah -    u_cab*
Cfah -     (de_cab_by_drhob*
Cfah -        du_cab_by_dza + 
Cfah -       e_cab*
Cfah -        du_cab_by_drhob_dza))
Cfah    d2F_cab(2,4) = c_nab*n*u_cab**(-2 + n)*
Cfah -  ((-1 + n)*e_cab*
Cfah -     du_cab_by_dzb
Cfah -      *du_cab_by_drhob + 
Cfah -    u_cab*(de_cab_by_drhob*
Cfah -        du_cab_by_dzb + e_cab*
Cfah -        du_cab_by_drhob_dzb ))
Cfah    d2F_cab(3,3) = c_nab*n*e_cab*
Cfah -  u_cab**(-2 + n)*
Cfah -  ((-1 + n)*du_cab_by_dza**2 + 
Cfah -    u_cab*
Cfah -     du_cab_by_dza_dza)
Cfah    d2F_cab(3,4) = c_nab*n*e_cab*
Cfah -  u_cab**(-2 + n)*
Cfah -  ((-1 + n)*du_cab_by_dzb*
Cfah -     du_cab_by_dza
Cfah -       + u_cab*
Cfah -     du_cab_by_dza_dzb)
Cfah    d2F_cab(4,4) = c_nab*n*e_cab*
Cfah -  u_cab**(-2 + n)*
Cfah -  ((-1 + n)*du_cab_by_dzb**2 + 
Cfah -    u_cab*
Cfah -     du_cab_by_dzb_dzb)
Cfah
Cfah here, the second derivatives are completed (Schwartz's rule: 
Cfah df/(dadb) = df/(dbda) 
Cfah    DO i = 1, 4
Cfah      DO j = i, 4
Cfah        d2F_xa(j,i) = d2F_xa(i,j)  
Cfah        d2F_xb(j,i) = d2F_xb(i,j)
Cfah        d2F_caa(j,i) = d2F_caa(i,j) 
Cfah        d2F_cbb(j,i) = d2F_cbb(i,j) 
Cfah        d2F_cab(j,i) = d2F_cab(i,j) 
Cfah      ENDDO
Cfah    ENDDO

Cfah test for zero densities (as in beta part of H atom):
        IF (rhob.LT.tol_rho) THEN 
          DO i = 1, 4
            dF_xb(i) = 0.D0
            dF_cbb(i) = 0.D0
            dF_cab(i) = 0.D0
Cfah        DO j = 1, 4
Cfah          d2F_xb(i,j) = 0.D0
Cfah          d2F_cbb(i,j) = 0.D0
Cfah          d2F_cab(i,j) = 0.D0
Cfah        ENDDO 
          ENDDO 
        ENDIF

        IF (rhoa.LT.tol_rho) THEN
          DO i = 1, 4
            dF_xa(i) = 0.D0
            dF_caa(i) = 0.D0
            dF_cab(i) = 0.D0
Cfah        DO j = 1, 4
Cfah          d2F_xa(i,j) = 0.D0
Cfah          d2F_caa(i,j) = 0.D0
Cfah          d2F_cab(i,j) = 0.D0
Cfah        ENDDO
          ENDDO
        ENDIF


Cfah Sum up all the partial derivatives with respect to the same function
Cfah of terms containing different powers of u with the help of the big outer 
Cfah loop: 

Cfah have the partial derivative 

        DO i = 1, 4
          F(n*3+1,i) = dF_xa(i) + dF_xb(i) 
          F(n*3+2,i) = dF_caa(i) +dF_cbb(i) 
          F(n*3+3,i) = dF_cab(i) 
Cfah      DO j = 1, 4
Cfah        FF(n*3+1,i,j) = d2F_xa(i,j) + d2F_xb(i,j) 
Cfah        FF(n*3+2,i,j) = d2F_caa(i,j) + d2F_cbb(i,j) 
Cfah        FF(n*3+3,i,j) = d2F_cab(i,j)  
Cfah      ENDDO 
        ENDDO 

Cfah these partial derivatives have not been computed because they are
Cfah zero since we don't have a gradrhoagradrhob term in the Becke V functional
        F(n*3+1,5) = 0
        F(n*3+2,5) = 0
        F(n*3+3,5) = 0
Cfah    DO i = 1, 5
Cfah      FF(n*3+1,i,5) = 0
Cfah      FF(n*3+2,i,5) = 0
Cfah      FF(n*3+3,i,5) = 0
Cfah      FF(n*3+1,5,i) = 0
Cfah      FF(n*3+2,5,i) = 0
Cfah      FF(n*3+3,5,i) = 0
Cfah    ENDDO

      ENDDO

      RETURN

      END

Cfah-----------------------------------------------------------

