*     ************************************************
*     *                                              *
*     *              nwpw_tpss03_x                   *
*     *                                              *
*     ************************************************
      subroutine nwpw_tpss03_x(pi,thrd,twthrd,frthrd,fvthrd,etthrd,
     >                         C920,C1,C2,C3,kappa,mu,b,c,e,
     >                         Cx,P23,es,
     >                         n,agr,tau,
     >                         xe,dfdnx,dfdagrx,dfdtaux)
      implicit none
*     ***** input *****
      real*8 pi,thrd,twthrd,frthrd,fvthrd,etthrd,C920,C1,C2,C3
      real*8 kappa,mu,b,c,e
      real*8 Cx,P23,es
      real*8 n,agr,tau
*     ***** output *****
      real*8 xe,dfdnx,dfdagrx,dfdtaux
*     ***** local declarations *****
      real*8 n_13,n_53,n_83,inv_n,agr2,tauW,tauU
      real*8 p,z,fz,alpha,qb,p2,p3,z2,z3,qb2,fa,thresA,dalpha_dfa
      real*8 x1a,x1,x2,x3a,x3b,x3,x4,x5,x6,x7,xnum,x
      real*8 dp_dn,dz_dn,dalpha_dn,qba,qbb,dqba_dn,dqbb_dn,dqb_dn
      real*8 dx1_dn,dx2_dn,dx3_dn,dx4_dn,dx5_dn,dx6_dn,dx7_dn 
      real*8 dx3a_dn,dx3b_dn,dxnum_dn,dx_dn
      real*8 dp_dagr,dz_dagr,dalpha_dagr,dqba_dagr,dqbb_dagr,dqb_dagr
      real*8 dx1_dagr,dx2_dagr,dx3_dagr,dx4_dagr,dx5_dagr,dx6_dagr
      real*8 dx7_dagr,dx3a_dagr,dx3b_dagr,dxnum_dagr,dx_dagr
      real*8 dz_dtau,dalpha_dtau,dqb_dtau
      real*8 dx1_dtau,dx2_dtau,dx3_dtau,dx5_dtau,dxnum_dtau,dx_dtau
      real*8 Fx,ex0,nex0,dFx_dx,dFx_dn,dFx_dagr,dFx_dtau


      n_13  = n**thrd
      n_53  = n_13*n_13*n
      n_83  = n_53*n
      inv_n = 1.0d0/n
      agr2  = agr*agr
        
      p       =  agr2/(4.0d0*P23*n_83)
      dp_dn   = -etthrd*p*inv_n
      dp_dagr =  2.0d0*p/agr
c     dp_dtau =  0.0d0

      p2 = p*p
      p3 = p2*p

      tauW  = 0.125d0*agr2*inv_n
      tauU  = 0.3d0*P23*n_53
      fz    = tauW/tau

      if (fz .gt. 1.0d0) then
         z       = 1.0d0
         dz_dn   = 0.0d0
         dz_dagr = 0.0d0
         dz_dtau = 0.0d0
      else 
         z       =  fz
         dz_dn   = -z*inv_n
         dz_dagr =  2.0d0*z/agr
         dz_dtau = -z/tau
      end if

      z2 = z*z
      z3 = z2*z

      alpha = fvthrd*p*(1.0d0/z - 1.0d0)
c     alpha = (tau - tauW)/tauU

      if (alpha .le. 0.0d0) then
        alpha       = 0.0d0
        dalpha_dn   = 0.0d0
        dalpha_dagr = 0.0d0
        dalpha_dtau = 0.0d0
      else
        dalpha_dn   = fvthrd*(-p*dz_dn/z2 + dp_dn*(1.0d0/z - 1.0d0))
        dalpha_dagr = (alpha/p)*dp_dagr - fvthrd*(p/z2)*dz_dagr
        dalpha_dtau = 1.0d0/tauU
c       dalpha_dtau = fvthrd*p*(-1.0d0/z2)*dz_dtau
      end if

      qb = C920*(alpha - 1.0d0)/dsqrt(1.0d0+b*alpha*(alpha - 1.0d0))
     &   + twthrd*p

      qb2 = qb*qb

      x1a = (c*z2)/((1.0d0 + z2)**2.0d0)
      x1  = (C1 + x1a)*p
      x2  = C2*qb2
      x3a = C3*qb
      x3b = dsqrt(0.5d0*(0.36d0*z2 + p2))
      x3  = x3a*x3b
      x4  = C1*C1*p2/kappa
      x5  = 2.0d0*es*C1*0.36d0*z2
      x6  = e*mu*p3
      x7  = 1.0d0/((1.0d0 + es*p)**2.0d0)

      xnum = x1 + x2 + x3 + x4 + x5 + x6
      x    = xnum*x7

*     ***** dxdn *****
      qba     = alpha - 1.0d0
      qbb     = 1.0d0/dsqrt(1.0d0 + b*alpha*(alpha - 1.0d0))
      dqba_dn = dalpha_dn
      dqbb_dn = b*dalpha_dn*(2.0d0*alpha - 1.0d0)*
     &          (-0.5d0)*(qbb**3.0d0)
      dqb_dn  = C920*(qbb*dqba_dn + qba*dqbb_dn) + twthrd*dp_dn

      dx1_dn  = (x1/p)*dp_dn 
     &        + 2.0d0*c*p*z*dz_dn/((1.0d0 + z2)**3.0d0)*(1.0d0 - z2)
      dx2_dn  = 2.0d0*C2*qb*dqb_dn

      dx3a_dn = C3*dqb_dn
      dx3b_dn = 0.5d0/x3b*(0.36d0*z*dz_dn + p*dp_dn)
      dx3_dn  = x3b*dx3a_dn+x3a*dx3b_dn

      dx4_dn  = (2.0d0*x4/p)*dp_dn
      dx5_dn  = (2.0d0*x5/z)*dz_dn
      dx6_dn  = (3.0d0*x6/p)*dp_dn
      dx7_dn  = -2.0d0*es*dp_dn/(1.0d0 + es*p)**3.0d0

      dxnum_dn = dx1_dn + dx2_dn + dx3_dn + dx4_dn + dx5_dn + dx6_dn
      dx_dn    = x7*dxnum_dn + xnum*dx7_dn

*     ***** dxdagr *****
      dqba_dagr = dalpha_dagr
      dqbb_dagr = -0.5d0*dalpha_dagr*b*(2.0d0*alpha - 1.0d0)
     &            *qbb**3.0d0
      dqb_dagr  = C920*(qba*dqbb_dagr + qbb*dqba_dagr) 
     &          + twthrd*dp_dagr
      dx1_dagr  = (x1/p)*dp_dagr 
     &          + 2.0d0*c*p*z*dz_dagr/((1.0d0 + z2)**3.0d0)*(1.0d0 - z2)

      dx2_dagr  = C2*2.0d0*qb*dqb_dagr

      dx3a_dagr = C3*dqb_dagr
      dx3b_dagr = 0.5d0/x3b*( 0.36d0*z*dz_dagr + p*dp_dagr)
      dx3_dagr  = x3b*dx3a_dagr + x3a*dx3b_dagr

      dx4_dagr  = (2.0d0*x4/p)*dp_dagr
      dx5_dagr  = (2.0d0*x5/z)*dz_dagr
      dx6_dagr  = (3.0d0*x6/p)*dp_dagr

      dx7_dagr  = -2.0d0*es*dp_dagr/(1.0d0 + es*p)**3.0d0

      dxnum_dagr = dx1_dagr + dx2_dagr + dx3_dagr 
     &           + dx4_dagr + dx5_dagr + dx6_dagr
      dx_dagr    = x7*dxnum_dagr + xnum*dx7_dagr

*     ***** dxdtau *****

      dqb_dtau = C920*dalpha_dtau*qbb*(1.0d0
     &         - 0.5d0*b*qba*qbb*qbb*(2.0d0*alpha - 1.0d0))

      dx1_dtau = c*p*dz_dtau*2.0d0*z*(1.0d0 - z2)/
     &           ((1.0d0 + z2)**3.0d0)
      dx2_dtau = 2.0d0*c2*qb*dqb_dtau
      dx3_dtau = x3*(dqb_dtau/qb 
     &         + 0.5d0*0.36d0*z*dz_dtau/(x3b*x3b))
      dx5_dtau = 2.0d0*(x5/z)*dz_dtau

      dxnum_dtau = dx1_dtau + dx2_dtau + dx3_dtau + dx5_dtau
      dx_dtau    = x7*dxnum_dtau

      Fx = 1.0d0 + kappa - kappa/(1.0d0 + x/kappa)

      ex0  = Cx*n_13
      xe   = ex0*Fx
      nex0 = n*ex0

      dFx_dx = 1.0d0/(1.0d0 + x/kappa)**2.0d0
      dFx_dn = dFx_dx*dx_dn
      dfdnx  = frthrd*xe + nex0*dFx_dn

      dFx_dagr = dFx_dx*dx_dagr
      dfdagrx  = nex0*dFx_dagr

      dFx_dtau = dFx_dx*dx_dtau
      dfdtaux  = nex0*dFx_dtau
     
      return
      end

*     ************************************************
*     *                                              *
*     *              gen_TPSS03_restricted           *
*     *                                              *
*     ************************************************

*    This function returns the TPSS03 exchange-correlation
*  energy density, xce, and its derivatives with respect
*  to n, |grad n|, tau.

*
*   Entry - n2ft3d   : number of grid points
*           rho_in(*) :  density (nup+ndn)
*           agr_in(*): |grad rho_in|
*           tau_in(*): tau
*           x_parameter: scale parameter for exchange
*           c_parameter: scale parameter for correlation
*
*     Exit  - xce(n2ft3d) : TPSS03 exchange correlation energy density
*             fn(n2ft3d)  : d(n*xce)/dn
*             fdn(n2ft3d) : d(n*xce)/d|grad n|
*             fdtau(n2ft3d) : d(n*xce)/dtau
*

      subroutine gen_TPSS03_restricted(n2ft3d,rho_in,agr_in,tau_in,
     >                               x_parameter,c_parameter,
     >                               xce,fn,fdn,fdtau)
      implicit none
*     ***** input *****
      integer n2ft3d
      real*8 rho_in(*),agr_in(*),tau_in(*)     
      real*8 x_parameter,c_parameter
*     ***** output *****
      real*8 xce(*),fn(*),fdn(*),fdtau(*)      
*     ***** local declarations *****
      integer i
      real*8 n,agr,tau
      real*8 eup0c,fnup0c,fdnup0c
      real*8 e0c,fn0c,fdn0c
      real*8 Cx,P23,es
      real*8 ex,fnx,fdnx,fdtaux
      real*8 z,z2,z3,tmp1,tmp2,agr2,c00
      real*8 tauW,dz2_dn,dz2_dagr,dz2_dtau,dz3_dn,dz3_dagr,dz3_dtau
      real*8 PKZB1,dPKZB1_dn,dPKZB1_dagr,dPKZB1_dtau
      real*8 pbec,dpbec_dn,dpbec_dagr
      real*8 pbeupc,dpbeupc_dn,dpbeupc_dagr
      real*8 etil,detil_dn,detil_dagr
      real*8 PKZB2,dPKZB2_dn,dPKZB2_dagr,dPKZB2_dtau
      real*8 revPKZB,drev_dn,drev_dagr,drev_dtau
      real*8 revz3,drevz3_dn,drevz3_dagr,drevz3_dtau
      real*8 ec,fnc,fdnc,fdtauc
*     ***** constants *****
      real*8 pi,thrd,twthrd,frthrd,fvthrd,etthrd
      parameter (pi     =  3.14159265358979311599d0)
      parameter (thrd   =  1.0d0/3.0d0)
      parameter (twthrd =  2.0d0/3.0d0)
      parameter (frthrd =  4.0d0/3.0d0)
      parameter (fvthrd =  5.0d0/3.0d0)
      parameter (etthrd =  8.0d0/3.0d0)
*     ***** density cutoff parameters *****
      real*8 tol,ETA
      parameter (tol = 1.0d-18)
      parameter (ETA            =      1.0d-20)
*     ***** TPSS03 constants *****
      real*8 kappa,mu,b,c,d,e,C920,C1,C2,C3
      parameter (kappa =  0.8040d0)
      parameter (mu    =  0.21951d0)
      parameter (b     =  0.40d0)
      parameter (c     =  1.59096d0)
      parameter (d     =  2.8d0)
      parameter (e     =  1.537d0)
      parameter (C920  =  9.0d0/20.0d0)
      parameter (C1    =  10.0d0/81.0d0)
      parameter (C2    =  146.0d0/2025.0d0)
      parameter (C3    = -73.0d0/405.0d0)

      Cx  = (-0.75d0)*(3.0d0/pi)**thrd
      P23 = (3.0d0*pi**2.0d0)**twthrd
      es  = dsqrt(e)

      do i=1,n2ft3d
        n       =       rho_in(i) + ETA
        agr     =       agr_in(i) + ETA
        tau     = 2.0d0*tau_in(i) + ETA

*       ***** TPSS03 Exchange *****
        call nwpw_tpss03_x(pi,thrd,twthrd,frthrd,fvthrd,etthrd,
     >                     C920,C1,C2,C3,kappa,mu,b,c,e,
     >                     Cx,P23,es,
     >                     n,agr,tau,
     >                     ex,fnx,fdnx,fdtaux)

*       ***** TPSS03 Correlation *****
        call gen_PBE96_c_restricted(rho_in(i),agr_in(i),
     >                            e0c,fn0c,fdn0c)

        pbec       = e0c
        dpbec_dn   = (fn0c - pbec)/n
        dpbec_dagr = fdn0c/n

        call gen_PBE96_c_unrestricted(0.50d0*rho_in(i),0.50d0*agr_in(i),
     >                              eup0c,fnup0c,fdnup0c)

        pbeupc       = eup0c
        dpbeupc_dn   = (fnup0c - pbeupc)/n
        dpbeupc_dagr = fdnup0c/n

        c00  = 0.53d0

        agr2 = agr*agr

        tauW = 0.125d0*agr2/n
        z    = tauW/tau

        if (z .gt. 1.0d0) then
          z2       = 1.0d0
          dz2_dn   = 0.0d0
          dz2_dagr = 0.0d0
          dz2_dtau = 0.0d0

          z3       = 1.0d0
          dz3_dn   = 0.0d0
          dz3_dagr = 0.0d0
          dz3_dtau = 0.0d0
        else
          z2       =  z*z
          dz2_dn   = -2.0d0*z2/n
          dz2_dagr =  4.0d0*z2/agr
          dz2_dtau = -2.0d0*z2/tau

          z3       =  z2*z
          dz3_dn   = -3.0d0*z3/n
          dz3_dagr =  6.0d0*z3/agr
          dz3_dtau = -3.0d0*z3/tau 
        end if

        tmp1        = 1.0d0 + c00*z2
        PKZB1       = pbec*tmp1
        dPKZB1_dn   = dpbec_dn*tmp1     + pbec*c00*dz2_dn
        dPKZB1_dagr = dpbec_dagr*tmp1   + pbec*c00*dz2_dagr
        dPKZB1_dtau = pbec*c00*dz2_dtau

        if (pbeupc .lt. pbec) then
          etil       = pbec
          detil_dn   = dpbec_dn
          detil_dagr = dpbec_dagr
        else
          etil       = pbeupc
          detil_dn   = dpbeupc_dn
          detil_dagr = dpbeupc_dagr
        endif
  
        tmp1        = -(1.0d0 + c00)
        tmp2        =  etil 
        PKZB2       =  tmp1*z2*tmp2
        dPKZB2_dn   =  tmp1*(dz2_dn*tmp2   + z2*detil_dn)
        dPKZB2_dagr =  tmp1*(dz2_dagr*tmp2 + z2*detil_dagr)
        dPKZB2_dtau =  tmp1*dz2_dtau*tmp2

        revPKZB   = PKZB1       + PKZB2
        drev_dn   = dPKZB1_dn   + dPKZB2_dn
        drev_dagr = dPKZB1_dagr + dPKZB2_dagr
        drev_dtau = dPKZB1_dtau + dPKZB2_dtau

        revz3       = 1.0d0           + d*revPKZB*z3
        drevz3_dn   = d*(drev_dn*z3   + revPKZB*dz3_dn)
        drevz3_dagr = d*(drev_dagr*z3 + revPKZB*dz3_dagr)
        drevz3_dtau = d*(drev_dtau*z3 + revPKZB*dz3_dtau)
 
        ec     = revPKZB*revz3
        fnc    = n*(drev_dn*revz3   + revPKZB*drevz3_dn)    + ec
        fdnc   = n*(drev_dagr*revz3 + revPKZB*drevz3_dagr) 
        fdtauc = n*(drev_dtau*revz3 + revPKZB*drevz3_dtau)

        xce(i)   = x_parameter*ex     + c_parameter*ec
        fn(i)    = x_parameter*fnx    + c_parameter*fnc
        fdn(i)   = x_parameter*fdnx   + c_parameter*fdnc
        fdtau(i) = x_parameter*fdtaux + c_parameter*fdtauc

      end do
      return
      end

*     ************************************************
*     *                                              *
*     *            gen_TPSS03_unrestricted           *
*     *                                              *
*     ************************************************

*    This function returns the TPSS03 exchange-correlation
*  energy density, xce, and its derivatives with respect
*  to nup, ndn, |grad nup|, |grad ndn|, |grad n|, tauup, taudn.

*
*   Entry - n2ft3d   : number of grid points
*           rho_in(*,2) :  density (nup and ndn)
*           agr_in(*,3): |grad rho_in| (nup, ndn and n)
*           tau_in(*,2): tau (nup and ndn)
*           x_parameter: scale parameter for exchange
*           c_parameter: scale parameter for correlation
*
*     Exit  - xce(n2ft3d) : TPSS03 exchange correlation energy density
*             fn(n2ft3d,2)  : d(n*xce)/dnup, d(n*xce)/dndn
*             fdn(n2ft3d,3) : d(n*xce)/d|grad nup|, d(n*xce)/d|grad ndn|, d(n*xce)/d|grad n|
*             fdtau(n2ft3d,2) : d(n*xce)/dtauup, d(n*xce)/dtaudn
*

      subroutine gen_TPSS03_unrestricted(n2ft3d,rho_in,agr_in,tau_in,
     >                               x_parameter,c_parameter,
     >                               xce,fn,fdn,fdtau)
      implicit none
*     ***** input *****
      integer n2ft3d
      real*8 rho_in(n2ft3d,2),agr_in(n2ft3d,3),tau_in(n2ft3d,2)      
      real*8 x_parameter,c_parameter
*     ***** output *****
      real*8 xce(n2ft3d),fn(n2ft3d,2),fdn(n2ft3d,3),fdtau(n2ft3d,2)  
*     ***** local declarations *****
      integer i
      real*8 n,agr,tau
      real*8 nup,agrup,tauup
      real*8 ndn,agrdn,taudn
      real*8 Cx,P23,es
      real*8 z,z2,z3
      real*8 ex,eupx,fnupx,fdnupx,fdtauupx,ednx,fndnx,fdndnx,fdtaudnx
      real*8 e0c,fn0c1,fn0c2,fdn0c1,fdn0c2,fdn0c3
      real*8 eup0c,fnup0c,fdnup0c
      real*8 edn0c,fndn0c,fdndn0c
      real*8 zeta,dzeta_dnup,dzeta_dndn
      real*8 onepzeta,onemzeta,onepzeta2,onemzeta2
      real*8 zeta2,one2mzeta2
      real*8 n2,agr2,agraa,agrbb,agrab
      real*8 tmp1,tmp2,denxi,denxi2,gzeta2
      real*8 dgzeta2_dnup,dgzeta2_dndn
      real*8 dgzeta2_dagrup,dgzeta2_dagrdn,dgzeta2_dagr
      real*8 xi2,dxi2_dnup,dxi2_dndn
      real*8 dxi2_dagrup,dxi2_dagrdn,dxi2_dagr
      real*8 cz0,dcz0_dzeta,onezetap43,onezetam73
      real*8 denczx,ddenczx_dnup,ddenczx_dndn,denczx4,denczx5,denczx8
      real*8 czx,dczx_dnup,dczx_dndn
      real*8 dczx_dagrup,dczx_dagrdn,dczx_dagr 
      real*8 tauW,dz2_dn,dz2_dagr,dz2_dtau,dz3_dn,dz3_dagr,dz3_dtau
      real*8 PKZB1,dPKZB1_dnup,dPKZB1_dndn
      real*8 dPKZB1_dagrup,dPKZB1_dagrdn,dPKZB1_dagr
      real*8 dPKZB1_dtauup,dPKZB1_dtaudn
      real*8 pbec,dpbec_dnup,dpbec_dndn
      real*8 dpbec_dagrup,dpbec_dagrdn,dpbec_dagr
      real*8 pbeupc,dpbeupc_dnup,dpbeupc_dndn
      real*8 dpbeupc_dagrup,dpbeupc_dagrdn,dpbeupc_dagr
      real*8 pbednc,dpbednc_dnup,dpbednc_dndn
      real*8 dpbednc_dagrup,dpbednc_dagrdn,dpbednc_dagr
      real*8 etilup,detilup_dnup,detilup_dndn
      real*8 detilup_dagrup,detilup_dagrdn,detilup_dagr
      real*8 etildn,detildn_dnup,detildn_dndn
      real*8 detildn_dagrup,detildn_dagrdn,detildn_dagr
      real*8 fa,fb
      real*8 PKZB2,dPKZB2_dnup,dPKZB2_dndn
      real*8 dPKZB2_dagrup,dPKZB2_dagrdn,dPKZB2_dagr
      real*8 dPKZB2_dtauup,dPKZB2_dtaudn
      real*8 revPKZB,drev_dnup,drev_dndn
      real*8 drev_dagrup,drev_dagrdn,drev_dagr
      real*8 drev_dtauup,drev_dtaudn
      real*8 revz3,drevz3_dnup,drevz3_dndn
      real*8 drevz3_dagrup,drevz3_dagrdn,drevz3_dagr
      real*8 drevz3_dtauup,drevz3_dtaudn
      real*8 ec,fnupc,fndnc,fdnupc,fdndnc,fdnc,fdtauupc,fdtaudnc
*     ***** constants *****
      real*8 pi,thrd,twthrd,frthrd,fvthrd,svthrd,etthrd
      parameter (pi     = 3.14159265358979311599d0)
      parameter (thrd   = 1.0d0/3.0d0)
      parameter (twthrd = 2.0d0/3.0d0)
      parameter (frthrd = 4.0d0/3.0d0)
      parameter (fvthrd = 5.0d0/3.0d0)
      parameter (svthrd = 7.0d0/3.0d0)
      parameter (etthrd = 8.0d0/3.0d0)
*     ***** density cutoff parameters *****
      real*8 tol,ETA
      parameter (tol = 1.0d-18)
      parameter (ETA            =      1.0d-20)
*     ***** TPSS03 constants *****
      real*8 kappa,mu,b,c,d,e,C920,C1,C2,C3
      parameter (kappa =  0.8040d0)
      parameter (mu    =  0.21951d0)
      parameter (b     =  0.40d0)
      parameter (c     =  1.59096d0)
      parameter (d     =  2.8d0)
      parameter (e     =  1.537d0)
      parameter (C920  =  9.0d0/20.0d0)
      parameter (C1    =  10.0d0/81.0d0)
      parameter (C2    =  146.0d0/2025.0d0)
      parameter (C3    = -73.0d0/405.0d0)

      Cx  = (-0.75d0)*(3.0d0/pi)**thrd
      P23 = (3.0d0*pi**2.0d0)**twthrd
      es  = dsqrt(e)

      do i=1,n2ft3d
        nup       = rho_in(i,1) + ETA
        agrup     = agr_in(i,1) + ETA
        tauup     = tau_in(i,1) + ETA
        ndn       = rho_in(i,2) + ETA
        agrdn     = agr_in(i,2) + ETA
        taudn     = tau_in(i,2) + ETA

        
*       ***** TPSS03 Exchange *****
*       ***** UP *****
        n   = 2.0d0*nup
        agr = 2.0d0*agrup
        tau = 2.0d0*tauup

        call nwpw_tpss03_x(pi,thrd,twthrd,frthrd,fvthrd,etthrd,
     >                     C920,C1,C2,C3,kappa,mu,b,c,e,
     >                     Cx,P23,es,
     >                     n,agr,tau,eupx,fnupx,fdnupx,fdtauupx)

*       ***** DOWN *****
        n   = 2.0d0*ndn
        agr = 2.0d0*agrdn
        tau = 2.0d0*taudn

        call nwpw_tpss03_x(pi,thrd,twthrd,frthrd,fvthrd,etthrd,
     >                     C920,C1,C2,C3,kappa,mu,b,c,e,
     >                     Cx,P23,es,
     >                     n,agr,tau,ednx,fndnx,fdndnx,fdtaudnx)

        n  = nup + ndn

        ex = (eupx*nup + ednx*ndn)/n

*       ***** TPSS03 Correlation *****
        agr       = agr_in(i,3) + ETA

        tau  = tauup + taudn
        n2   = n*n
        agr2 = agr*agr

        call gen_PBE96_c_full_unrestricted(rho_in(i,1),rho_in(i,2),
     >                            agr_in(i,1),agr_in(i,2),agr_in(i,3),
     >                            e0c,fn0c1,fn0c2,fdn0c1,fdn0c2,fdn0c3) 

        pbec         = e0c
        dpbec_dnup   = (fn0c1 - pbec)/n
        dpbec_dndn   = (fn0c2 - pbec)/n
        dpbec_dagrup = fdn0c1/n
        dpbec_dagrdn = fdn0c2/n
        dpbec_dagr   = fdn0c3/n

        call gen_PBE96_c_unrestricted(rho_in(i,1),agr_in(i,1),
     >                            eup0c,fnup0c,fdnup0c)

        pbeupc         = eup0c
        dpbeupc_dnup   = (fnup0c - pbeupc)/nup  
        dpbeupc_dndn   = 0.0d0   
        dpbeupc_dagrup = fdnup0c/nup  
        dpbeupc_dagrdn = 0.0d0 
        dpbeupc_dagr   = 0.0d0  

        call gen_PBE96_c_unrestricted(rho_in(i,2),agr_in(i,2),
     >                            edn0c,fndn0c,fdndn0c) 

        pbednc         = edn0c
        dpbednc_dnup   = 0.0d0
        dpbednc_dndn   = (fndn0c - pbednc)/ndn  
        dpbednc_dagrup = 0.0d0
        dpbednc_dagrdn = fdndn0c/ndn  
        dpbednc_dagr   = 0.0d0  

c       if (dabs(nup - ndn) .lt. 1.0d-7) then
c         czx         = 0.53d0
c         dczx_dnup   = 0.0d0
c         dczx_dndn   = 0.0d0
c         dczx_dagrup = 0.0d0
c         dczx_dagrdn = 0.0d0
c         dczx_dagr   = 0.0d0
c       else
        zeta = (nup - ndn)/n

        onepzeta  = 1.0d0 + zeta
        onemzeta  = 1.0d0 - zeta
        onepzeta2 = onepzeta*onepzeta
        onemzeta2 = onemzeta*onemzeta

        zeta2      = zeta*zeta
        one2mzeta2 = 1.0d0 - zeta2

        dzeta_dnup =  onemzeta/n
        dzeta_dndn = -onepzeta/n

        agraa = agrup*agrup
        agrbb = agrdn*agrdn
        agrab = agraa + agrbb - agr2

        denxi  = 2.0d0*((3.0d0*pi*pi*n)**thrd)
        denxi2 = denxi*denxi

        gzeta2       = (agraa*onemzeta2 + agrbb*onepzeta2 
     &               + agrab*one2mzeta2)/n2

        tmp1         = (-agraa*onemzeta + agrbb*onepzeta 
     &               - agrab*zeta)*2.0d0*dzeta_dnup
        tmp1         = tmp1/n2
        tmp2         = -2.0d0*gzeta2/n
        dgzeta2_dnup = tmp1 + tmp2

        tmp1         = (-agraa*onemzeta + agrbb*onepzeta 
     &               - agrab*zeta)*2.0d0*dzeta_dndn
        tmp1         = tmp1/n2
        tmp2         = -2.0d0*gzeta2/n
        dgzeta2_dndn = tmp1 + tmp2

        dgzeta2_dagrup =  2.0d0*agrup*(onemzeta2 + one2mzeta2)/n2
        dgzeta2_dagrdn =  2.0d0*agrdn*(onepzeta2 + one2mzeta2)/n2
        dgzeta2_dagr   = -2.0d0*agr*one2mzeta2/n2

        xi2 = gzeta2/denxi2
  
        dxi2_dnup = (dgzeta2_dnup - twthrd*gzeta2/n)/denxi2
        dxi2_dndn = (dgzeta2_dndn - twthrd*gzeta2/n)/denxi2

        dxi2_dagrup = dgzeta2_dagrup/denxi2
        dxi2_dagrdn = dgzeta2_dagrdn/denxi2
        dxi2_dagr   = dgzeta2_dagr/denxi2
       
        cz0        = 0.53d0 + 0.87d0*zeta**2.0d0 
     &             + 0.50d0*zeta**4.0d0 + 2.26*zeta**6.0d0
        dcz0_dzeta = 1.74d0*zeta + 2.0d0*zeta**3.0d0 
     &             + 13.56d0*zeta**5.0d0
        onezetap43 = onepzeta**(-frthrd) + onemzeta**(-frthrd)
        onezetam73 = onemzeta**(-svthrd) - onepzeta**(-svthrd)

        denczx       = 1.0d0 + 0.50d0*xi2*onezetap43
        ddenczx_dnup = 2.0d0*(denczx**3.0d0)*(dxi2_dnup*onezetap43 
     &               + xi2*frthrd*onezetam73*dzeta_dnup)
        ddenczx_dndn = 2.0d0*(denczx**3.0d0)*(dxi2_dndn*onezetap43
     &               + xi2*frthrd*onezetam73*dzeta_dndn)

        denczx4 = denczx**4.0d0
        denczx5 = denczx4*denczx
        denczx8 = denczx4*denczx4

        czx         = cz0/denczx4
        dczx_dnup   = dcz0_dzeta*dzeta_dnup/denczx4
     &              - cz0*ddenczx_dnup/denczx8
        dczx_dndn   = dcz0_dzeta*dzeta_dndn/denczx4
     &              - cz0*ddenczx_dndn/denczx8
        dczx_dagrup = -2.0d0*cz0/denczx5*onezetap43*dxi2_dagrup 
        dczx_dagrdn = -2.0d0*cz0/denczx5*onezetap43*dxi2_dagrdn 
        dczx_dagr   = -2.0d0*cz0/denczx5*onezetap43*dxi2_dagr
c       end if

        tauW = 0.125d0*agr2/n
        z    = tauW/tau

        if (z .gt. 1.0d0) then
          z2       = 1.0d0
          dz2_dn   = 0.0d0
          dz2_dagr = 0.0d0
          dz2_dtau = 0.0d0

          z3       = 1.0d0
          dz3_dn   = 0.0d0
          dz3_dagr = 0.0d0
          dz3_dtau = 0.0d0
        else
          z2       =  z*z
          dz2_dn   = -2.0d0*z2/n
c         dz2_dnup == dz2_dndn == dz2_dn
          dz2_dagr =  4.0d0*z2/agr
          dz2_dtau = -2.0d0*z2/tau
c         dz2_dtauup == dz2_dtaudn == dz2_dtau

          z3       =  z2*z
          dz3_dn   = -3.0d0*z3/n
c         dz3_dnup == dz3_dndn == dz3_dn
          dz3_dagr =  6.0d0*z3/agr
          dz3_dtau = -3.0d0*z3/tau 
c         dz3_dtauup == dz3_dtaudn == dz3_dtau
        end if

        tmp1          = 1.0d0 + czx*z2
        PKZB1         = pbec*tmp1
        dPKZB1_dnup   = dpbec_dnup*tmp1
     &                + pbec*(dczx_dnup*z2 + czx*dz2_dn)
        dPKZB1_dndn   = dpbec_dndn*tmp1
     &                + pbec*(dczx_dndn*z2 + czx*dz2_dn)
        dPKZB1_dagrup = dpbec_dagrup*tmp1
     &                + pbec*dczx_dagrup*z2 
        dPKZB1_dagrdn = dpbec_dagrdn*tmp1
     &                + pbec*dczx_dagrdn*z2 
        dPKZB1_dagr   = dpbec_dagr*tmp1
     &                + pbec*(dczx_dagr*z2 + czx*dz2_dagr)
        dPKZB1_dtauup = pbec*czx*dz2_dtau
        dPKZB1_dtaudn = dPKZB1_dtauup

        if (pbeupc .lt. pbec) then
          etilup         = pbec
          detilup_dnup   = dpbec_dnup
          detilup_dndn   = dpbec_dndn
          detilup_dagrup = dpbec_dagrup
          detilup_dagrdn = dpbec_dagrdn
          detilup_dagr   = dpbec_dagr
        else
          etilup         = pbeupc
          detilup_dnup   = dpbeupc_dnup
          detilup_dndn   = dpbeupc_dndn 
          detilup_dagrup = dpbeupc_dagrup
          detilup_dagrdn = dpbeupc_dagrdn
          detilup_dagr   = dpbeupc_dagr
        endif
  
        if (pbednc .lt. pbec) then
          etildn         = pbec
          detildn_dnup   = dpbec_dnup 
          detildn_dndn   = dpbec_dndn
          detildn_dagrup = dpbec_dagrup
          detildn_dagrdn = dpbec_dagrdn
          detildn_dagr   = dpbec_dagr
        else
          etildn         = pbednc
          detildn_dnup   = dpbednc_dnup
          detildn_dndn   = dpbednc_dndn
          detildn_dagrup = dpbednc_dagrup
          detildn_dagrdn = dpbednc_dagrdn
          detildn_dagr   = dpbednc_dagr
        endif

        fa = nup/n
        fb = ndn/n

        tmp1          = -(1.0d0 + czx)
        tmp2          =  fa*etilup + fb*etildn
        PKZB2         =  tmp1*z2*tmp2
        dPKZB2_dnup   = -dczx_dnup*z2*tmp2 + tmp1*(dz2_dn*tmp2 
     &                +  z2*(ndn/n2*(etilup - etildn) 
     &                +  fa*detilup_dnup + fb*detildn_dnup))
        dPKZB2_dndn   = -dczx_dndn*z2*tmp2 + tmp1*(dz2_dn*tmp2
     &                +  z2*(nup/n2*(etildn - etilup)
     &                +  fb*detildn_dndn + fa*detilup_dndn))
        dPKZB2_dagrup = -dczx_dagrup*z2*tmp2
     &                +  tmp1*z2*(fa*detilup_dagrup + fb*detildn_dagrup)
        dPKZB2_dagrdn = -dczx_dagrdn*z2*tmp2
     &                +  tmp1*z2*(fb*detildn_dagrdn + fa*detilup_dagrdn)
        dPKZB2_dagr   = -dczx_dagr*z2*tmp2 + tmp1*(dz2_dagr*tmp2
     &                +  z2*(fa*detilup_dagr + fb*detildn_dagr))
        dPKZB2_dtauup =  tmp1*dz2_dtau*tmp2
        dPKZB2_dtaudn =  dPKZB2_dtauup

        revPKZB     = PKZB1         + PKZB2
        drev_dnup   = dPKZB1_dnup   + dPKZB2_dnup
        drev_dndn   = dPKZB1_dndn   + dPKZB2_dndn
        drev_dagrup = dPKZB1_dagrup + dPKZB2_dagrup
        drev_dagrdn = dPKZB1_dagrdn + dPKZB2_dagrdn
        drev_dagr   = dPKZB1_dagr   + dPKZB2_dagr
        drev_dtauup = dPKZB1_dtauup + dPKZB2_dtauup
        drev_dtaudn = dPKZB1_dtaudn + dPKZB2_dtaudn

        revz3         = 1.0d0             + d*revPKZB*z3
        drevz3_dnup   = d*(drev_dnup*z3   + revPKZB*dz3_dn)
        drevz3_dndn   = d*(drev_dndn*z3   + revPKZB*dz3_dn)
        drevz3_dagrup = d*drev_dagrup*z3
        drevz3_dagrdn = d*drev_dagrdn*z3
        drevz3_dagr   = d*(drev_dagr*z3   + revPKZB*dz3_dagr)
        drevz3_dtauup = d*(drev_dtauup*z3 + revPKZB*dz3_dtau)
        drevz3_dtaudn = d*(drev_dtaudn*z3 + revPKZB*dz3_dtau)
 
        ec       = revPKZB*revz3
        fnupc    = n*(drev_dnup*revz3   + revPKZB*drevz3_dnup)   + ec
        fndnc    = n*(drev_dndn*revz3   + revPKZB*drevz3_dndn)   + ec
        fdnupc   = n*(drev_dagrup*revz3 + revPKZB*drevz3_dagrup) 
        fdndnc   = n*(drev_dagrdn*revz3 + revPKZB*drevz3_dagrdn) 
        fdnc     = n*(drev_dagr*revz3   + revPKZB*drevz3_dagr) 
        fdtauupc = n*(drev_dtauup*revz3 + revPKZB*drevz3_dtauup)
        fdtaudnc = n*(drev_dtaudn*revz3 + revPKZB*drevz3_dtaudn)
 
      
        xce(i)     = x_parameter*ex       + c_parameter*ec
        fn(i,1)    = x_parameter*fnupx    + c_parameter*fnupc
        fn(i,2)    = x_parameter*fndnx    + c_parameter*fndnc
        fdn(i,1)   = x_parameter*fdnupx   + c_parameter*fdnupc
        fdn(i,2)   = x_parameter*fdndnx   + c_parameter*fdndnc
        fdn(i,3)   = c_parameter*fdnc
        fdtau(i,1) = x_parameter*fdtauupx + c_parameter*fdtauupc
        fdtau(i,2) = x_parameter*fdtaudnx + c_parameter*fdtaudnc

      end do 

      return
      end 
