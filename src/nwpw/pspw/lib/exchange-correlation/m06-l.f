*     ************************************************
*     *                                              *
*     *              gen_M06L_restricted             *
*     *                                              *
*     ************************************************

*    This function returns the M06-L exchange-correlation
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
*     Exit  - xce(n2ft3d) : M06-L exchange correlation energy density
*             fn(n2ft3d)  : d(n*xce)/dn
*             fdn(n2ft3d) : d(n*xce)/d|grad n|
*             fdtau(n2ft3d) : d(n*xce)/dtau
*
      subroutine gen_M06L_restricted(n2ft3d,rho_in,agr_in,tau_in,
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
      real*8 Cx,P23 
      real*8 nup,agrup,tauup
      real*8 n,agr,tau
      real*8 xr,zr,gamma
      real*8 inv_n,n_13,n_43,n_53,n_83,agr2
      real*8 x,dx_dn,dx_dagr,z,dz_dn,dz_dtau
      real*8 GG,dGdx,dGdz,dG_dn,dG_dagr,dG_dtau
      real*8 n_onethird,kf,ks,ss,P0,FF,Fs,ex_lda,fdnx_const
      real*8 pbex,pbefnx,pbefdnx
      real*8 tauU,t,dt_dn,dt_dtau
      real*8 w,dw_dt,w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11
      real*8 fw,dfw_dw,dfw_dn,dfw_dtau
      real*8 ex,fnx,fdnx,fdtaux
      real*8 ess0c,dess0c_drs,dess0c_dn
      real*8 U,dU_dx,U2,U3,U4
      real*8 gss,dgss_dU,dgss_dx,gopp,dgopp_dU,dgopp_dx
      real*8 hss,dhss_dx,dhss_dz,hopp,dhopp_dx,dhopp_dz
      real*8 ghss,ghopp,dgh_dx,dgh_dn,dgh_dagr,dgh_dtau
      real*8 eud0c,deud0c_drs
      real*8 eudc,eud1c,fud1c,dfud1c_dn
      real*8 dfudc_dn,dfudc_dagr,dfudc_dtau
      real*8 rs,drs_dn,dummy
      real*8 wt1,dwt1,aw,bw,cw,dw
      real*8 D,dD_dx,dD_dz,dD_dn,dD_dagr,dD_dtau
      real*8 Dt,dD_dDt,dDt_dx,dDt_dz
      real*8 ec,fnc,fdnc,fdtauc
*     ***** constants *****
      real*8 pi,cf,thrd,twthrd,frthrd,fvthrd,etthrd
      parameter (pi     = 3.14159265358979311599d0)
      parameter (cf     = 9.115599720d0)
c     cf = 3.0d0/5.0d0*(6.0d0*pi*pi)**(2.0d0/3.0d0))
      parameter (thrd   = 1.0d0/3.0d0)
      parameter (twthrd = 2.0d0/3.0d0)
      parameter (frthrd = 4.0d0/3.0d0)
      parameter (fvthrd = 5.0d0/3.0d0)
      parameter (etthrd = 8.0d0/3.0d0)
*     ***** density cutoff parametersi *****
      real*8 tol,ETA,t1,t2,thresd,thresx
      parameter (tol = 1.0d-10)
      parameter (ETA            =      1.0d-20)
      parameter (t1      =   1.0d7)
      parameter (t2      =   5.0d7)
      parameter (thresd  =   5.0d-3)
      parameter (thresx  =   1.0d8)
c     ***** PBE96 GGA exchange constants ******
      real*8 MU,KAPPA
      parameter (MU    = 0.2195149727645171d0)
      parameter (KAPPA = 0.8040000000000000d0)
*     ***** VS98 constants *****
      real*8 aax,bbx,ccx,ddx,eex,ffx,alphax
      real*8 aass,bbss,ccss,ddss,eess,ffss,alphass
      real*8 aaopp,bbopp,ccopp,ddopp,eeopp,ffopp,alphaopp
      real*8 clda
      parameter (clda     =  0.9305257363491d0)
      parameter (aax      =  6.012244d-1)
      parameter (bbx      =  4.748822d-3)
      parameter (ccx      = -8.635108d-3)
      parameter (ddx      = -9.308062d-6)
      parameter (eex      =  4.482811d-5)
      parameter (ffx      =  0.000000d0)
      parameter (alphax   =  0.00186726d0)
      parameter (aass     =  4.650534d-1)
      parameter (bbss     =  1.617589d-1)
      parameter (ccss     =  1.833657d-1)
      parameter (ddss     =  4.692100d-4)
      parameter (eess     = -4.990573d-3)
      parameter (ffss     =  0.000000d0)
      parameter (alphass  =  0.00515088d0)
      parameter (aaopp    =  3.957626d-1)
      parameter (bbopp    = -5.614546d-1)
      parameter (ccopp    =  1.403963d-2)
      parameter (ddopp    =  9.831442d-4)
      parameter (eeopp    = -3.577176d-3)
      parameter (ffopp    =  0.000000d0)
      parameter (alphaopp =  0.00304966d0)
*     ***** M06 constants *****
      real*8 ma0,ma1,ma2,ma3,ma4,ma5,ma6,ma7,ma8,ma9,ma10,ma11
      real*8 C1,C2
      real*8 css,ccss0,ccss1,ccss2,ccss3,ccss4
      real*8 copp,ccopp0,ccopp1,ccopp2,ccopp3,ccopp4
      parameter (ma0    =  3.987756d-1)
      parameter (ma1    =  2.548219d-1)
      parameter (ma2    =  3.923994d-1)
      parameter (ma3    = -2.103655d0)
      parameter (ma4    = -6.302147d0)
      parameter (ma5    =  1.097615d1)
      parameter (ma6    =  3.097273d1)
      parameter (ma7    = -2.318489d1)
      parameter (ma8    = -5.673480d1)
      parameter (ma9    =  2.160364d1)
      parameter (ma10   =  3.421814d1)
      parameter (ma11   = -9.049762d0)
      parameter (C1     =  3.36116d-3)
      parameter (C2     =  4.49267d-3)
      parameter (css    =  0.06d0)
      parameter (ccss0  =  5.349466d-1)
      parameter (ccss1  =  5.396620d-1)
      parameter (ccss2  = -3.161217d1)
      parameter (ccss3  =  5.149592d1)
      parameter (ccss4  = -2.919613d1)
      parameter (copp   =  0.0031d0)
      parameter (ccopp0 =  6.042374d-1)
      parameter (ccopp1 =  1.776783d2)
      parameter (ccopp2 = -2.513252d2)
      parameter (ccopp3 =  7.635173d1)
      parameter (ccopp4 = -1.255699d1)

      Cx         = -1.50d0*(0.75d0/pi)**thrd
      P23        =  0.60d0*(6.0d0*pi*pi)**twthrd
      fdnx_const = -3.0d0/(8.0d0*pi)

      do i=1,n2ft3d
        n        = rho_in(i) + ETA
        agr      = agr_in(i) + ETA
        tau      = 2.0d0*tau_in(i) + ETA

        n   = 0.50d0*n
        agr = 0.50d0*agr
        tau = 0.50d0*tau

        agr2  = agr*agr
        inv_n = 1.0d0/n
        n_13  = n**thrd
        n_43  = n_13*n
        n_53  = n_43*n_13
        n_83  = n_53*n

        x       =  agr2/n_83
        dx_dn   = -etthrd*x*inv_n
        dx_dagr =  2.0d0*agr/n_83
        z       =  tau/n_53 - cf
        dz_dn   = -fvthrd*tau/n_83
        dz_dtau =  1.0d0/n_53

*       ***** VS98 Exchange *****
        gamma = 1.0d0 + alphax*(x + z)
        xr    = x/gamma
        zr    = z/gamma

        call nwpw_GVT4(aax,bbx,ccx,ddx,eex,ffx,alphax,alphax,
     >              xr,zr,gamma,GG,dGdx,dGdz)

        dG_dn   = dGdx*dx_dn + dGdz*dz_dn
        dG_dagr = dGdx*dx_dagr
        dG_dtau = dGdz*dz_dtau

        ex     = -clda*n_13*GG
        fnx    =  frthrd*ex - clda*n_43*dG_dn
        fdnx   = -clda*n_43*dG_dagr
        fdtaux = -clda*n_43*dG_dtau

*       ***** M06-L Enhancement Factor *****
        tauU    =  P23*n_53
        t       =  tauU/tau
        dt_dn   =  fvthrd*t/n
        dt_dtau = -t/tau

        w     = 1.0d0
        dw_dt = 0.0d0

        if(t .le. t1) then
          w     = (t - 1.0d0)/(t + 1.0d0)
          dw_dt = 2.0d0/((1.0d0 + t)**2.0d0)
        else if(t .gt. t1 .and. t .lt. t2) then
          wt1 = (t1 - 1.0d0)/(t1 + 1.0d0)
          dwt1 = 2.0d0/((1.0d0 + t1)**2.0d0)
          aw = (3.0d0*t1*t2*t2 - t2**3.0d0)*wt1/(t1 - t2)
     &       + (t1**3.0d0 - 3.0d0*t1*t1*t2)/(t1 - t2)
     &       - t1*t2*t2*dwt1
          aw = aw/((t1 - t2)*(t1 - t2))
          bw = -6.0d0*t1*t2*wt1/(t1 - t2) + 6.0d0*t1*t2/(t1 - t2)
     &       + (t2*t2 + 2.0d0*t1*t2)*dwt1
          bw = bw/((t1 - t2)*(t1 - t2))
          cw = 3.0d0*(t1 + t2)*wt1/(t1 - t2)
     &       - 3.0d0*(t1 + t2)/(t1 - t2) - (t1 + 2.0d0*t2)*dwt1
          cw = cw/((t1 - t2)*(t1 - t2))
          dw = -2.0d0*wt1/(t1 - t2) + 2.0d0/(t1 - t2) + dwt1
          dw = dw/((t1 - t2)*(t1 - t2))

          w     = aw + bw*t + cw*t*t + dw*t*t*t
          dw_dt = bw + 2.0d0*cw*t + 3.0d0*dw*t*t
        else if(t .ge. t2) then
          w     = 1.0d0
          dw_dt = 0.0d0
        end if

        w1  = w
        w2  = w1*w
        w3  = w2*w
        w4  = w3*w
        w5  = w4*w
        w6  = w5*w
        w7  = w6*w
        w8  = w7*w
        w9  = w8*w
        w10 = w9*w
        w11 = w10*w

        fw       = ma0    + ma1*w1 + ma2*w2 + ma3*w3 + ma4*w4   
     &           + ma5*w5 + ma6*w6 + ma7*w7 + ma8*w8 + ma9*w9 
     &           + ma10*w10 + ma11*w11
        dfw_dw   = ma1            + 2.0d0*ma2*w1     + 3.0d0*ma3*w2
     &           + 4.0d0*ma4*w3   + 5.0d0*ma5*w4     + 6.0d0*ma6*w5
     &           + 7.0d0*ma7*w6   + 8.0d0*ma8*w7     + 9.0d0*ma9*w8
     &           + 10.0d0*ma10*w9 + 11.0d0*ma11*w10
        dfw_dn   = dfw_dw*dw_dt*dt_dn
        dfw_dtau = dfw_dw*dw_dt*dt_dtau

*       ***** PBE96 Exchange *****
        n_onethird = (3.0d0*n/pi)**thrd
        ex_lda     = -0.75d0*n_onethird

        kf = (3.0d0*pi*pi*n)**thrd
        ss = agr/(2.0d0*kf*n)
        P0 = 1.0d0 + (MU/KAPPA)*ss*ss

        FF  = (1.0d0 + KAPPA - KAPPA/P0)
        Fs  = 2.0d0*MU/(P0*P0)*ss

        pbex    = ex_lda*FF
        pbefnx  = frthrd*(pbex - ex_lda*Fs*ss)
        pbefdnx = fdnx_const*Fs

        ex     = ex      + pbex*fw
        fnx    = fnx     + pbefnx*fw        + n*pbex*dfw_dn
        fdnx   = fdnx    + pbefdnx*fw
        fdtaux = fdtaux  + n*pbex*dfw_dtau


*       ***** VS98 Correlation *****
*       ***** Same-Spin (alpha-alpha/beta-beta) *****
        rs     = (0.75d0/(pi*n))**thrd
        drs_dn = -thrd*rs/n

        call gen_PW91_c_rz(tol,rs,1.0d0,ess0c,dess0c_drs,dummy)
        dess0c_dn = dess0c_drs*drs_dn

        Dt      = 1.0d0 - 0.25d0*x/(z + cf)
        D       = 0.0d0
        dD_dn   = 0.0d0
        dD_dagr = 0.0d0
        dD_dtau = 0.0d0

        if (Dt .le. 0.0d0) then
           D       = 0.0d0
           dD_dn   = 0.0d0
           dD_dagr = 0.0d0
           dD_dtau = 0.0d0
        else if (Dt .gt. 0.0d0 .and. Dt .lt. thresd) then
           D       =  2.0d0*Dt*Dt/thresd - Dt*Dt*Dt/(thresd*thresd)
           dD_dDt  =  4.0d0*Dt/thresd - 3.0d0*Dt*Dt/(thresd*thresd)
           dDt_dx  = -0.25d0/(z + cf)
           dDt_dz  =  0.25d0*x/((z + cf)*(z + cf))
           dD_dn   =  dD_dDt*(dDt_dx*dx_dn + dDt_dz*dz_dn)
           dD_dagr =  dD_dDt*dDt_dx*dx_dagr
           dD_dtau =  dD_dDt*dDt_dz*dz_dtau
        else if(Dt .ge. thresd) then
           D       =  Dt
           dD_dx   = -0.25d0/(z + cf)
           dD_dz   =  0.25d0*x/((z + cf)*(z + cf))
           dD_dn   =  dD_dx*dx_dn + dD_dz*dz_dn
           dD_dagr =  dD_dx*dx_dagr
           dD_dtau =  dD_dz*dz_dtau
        end if
 
        gamma = 1.0d0 + alphass*(x + z)
        xr    = x/gamma
        zr    = z/gamma

        if (x .ge. thresx) then
           U     = 1.0d0
           dU_dx = 0.0d0
        else
           U     = css*x/(1.0d0 + css*x)
           dU_dx = css/((1.0d0 + css*x)*(1.0d0 + css*x))
        end if

        U2 = U*U
        U3 = U2*U
        U4 = U3*U

        gss     = ccss0 + ccss1*U + ccss2*U2 + ccss3*U3 + ccss4*U4
        dgss_dU = ccss1 + 2.0d0*ccss2*U + 3.0d0*ccss3*U2
     &          + 4.0d0*ccss4*U3
        dgss_dx = dgss_dU*dU_dx

        call nwpw_GVT4(aass,bbss,ccss,ddss,eess,ffss,alphass,alphass,
     >                 xr,zr,gamma,hss,dhss_dx,dhss_dz)

        ghss     = gss + hss
        dgh_dx   = dgss_dx + dhss_dx
        dgh_dn   = dgh_dx*dx_dn + dhss_dz*dz_dn
        dgh_dagr = dgh_dx*dx_dagr
        dgh_dtau = dhss_dz*dz_dtau

        ec     = ess0c*ghss*D
        fnc    = ec + n*dess0c_drs*drs_dn*ghss*D
     &         + n*ess0c*dgh_dn*D 
     &         + n*ess0c*ghss*dD_dn
        fdnc   = n*ess0c*(dgh_dagr*D + ghss*dD_dagr)
        fdtauc = n*ess0c*(dgh_dtau*D + ghss*dD_dtau)

*       ***** Opposite-Spin (alpha-beta) *****
        n = 2.0d0*n

        rs         =  (0.75d0/(pi*n))**thrd
        drs_dn     = -thrd*rs/n

        call gen_PW91_c_rz(tol,rs,0.0d0,eud0c,deud0c_drs,dummy)
        eud1c      = eud0c - ess0c
        fud1c      = n*eud1c
        dfud1c_dn  = eud0c + n*deud0c_drs*drs_dn
     &             - ess0c - 0.50d0*n*dess0c_dn

        x   = 2.0d0*x
        z   = 2.0d0*z

        gamma = 1.0d0 + alphaopp*(x + z)
        xr    = x/gamma
        zr    = z/gamma

        if (x .ge. thresx) then
           U     = 1.0d0
           dU_dx = 0.0d0
        else
           U     = copp*x/(1.0d0 + copp*x)
           dU_dx = copp/((1.0d0 + copp*x)*(1.0d0 + copp*x))
        end if

        U2 = U*U
        U3 = U2*U
        U4 = U3*U

        gopp     = ccopp0 + ccopp1*U + ccopp2*U2 + ccopp3*U3 + ccopp4*U4
        dgopp_dU = ccopp1          + 2.0d0*ccopp2*U
     &           + 3.0d0*ccopp3*U2 + 4.0d0*ccopp4*U3
        dgopp_dx = dgopp_dU*dU_dx

        call nwpw_GVT4(aaopp,bbopp,ccopp,ddopp,eeopp,ffopp,
     >                 alphaopp,alphaopp,
     >                 xr,zr,gamma,hopp,dhopp_dx,dhopp_dz)

        ghopp    = gopp + hopp
        dgh_dx   = dgopp_dx + dhopp_dx
        dgh_dn   = dgh_dx*dx_dn + dhopp_dz*dz_dn
        dgh_dagr = dgh_dx*dx_dagr
        dgh_dtau = dhopp_dz*dz_dtau

        eudc       = eud1c*ghopp
        dfudc_dn   = dfud1c_dn*ghopp + fud1c*dgh_dn
        dfudc_dagr = fud1c*dgh_dagr
        dfudc_dtau = fud1c*dgh_dtau

        ec     = ec     + eudc
        fnc    = fnc    + dfudc_dn
        fdnc   = fdnc   + dfudc_dagr
        fdtauc = fdtauc + dfudc_dtau

        xce(i)   = x_parameter*ex     + c_parameter*ec
        fn(i)    = x_parameter*fnx    + c_parameter*fnc
        fdn(i)   = x_parameter*fdnx   + c_parameter*fdnc
        fdtau(i) = x_parameter*fdtaux + c_parameter*fdtauc
      end do

      return
      end
*     ************************************************
*     *                                              *
*     *              gen_M06L_unrestricted           *
*     *                                              *
*     ************************************************

*    This function returns the M06-L exchange-correlation
*  energy density, xce, and its derivatives with respect
*  to nup, ndn, |grad nup|, |grad ndn|, tauup, taudn.

*
*   Entry - n2ft3d   : number of grid points
*           rho_in(*,2) :  density (nup and ndn)
*           agr_in(*,3): |grad rho_in| (nup, ndn and n)
*           tau_in(*,2): tau (nup and ndn)
*           x_parameter: scale parameter for exchange
*           c_parameter: scale parameter for correlation
*
*     Exit  - xce(n2ft3d) : M06-L exchange correlation energy density
*             fn(n2ft3d,2)  : d(n*xce)/dnup, d(n*xce)/dndn
*             fdn(n2ft3d,3) : d(n*xce)/d|grad nup|, d(n*xce)/d|grad ndn|, d(n*xce)/d|grad n|
*             fdtau(n2ft3d,2) : d(n*xce)/dtauup, d(n*xce)/dtaudn
*
      subroutine gen_M06L_unrestricted(n2ft3d,rho_in,agr_in,tau_in,
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
      real*8 Cx,P23 
      real*8 nup,agrup,tauup
      real*8 ndn,agrdn,taudn
      real*8 n,agr,tau
      real*8 x,dxdn,dxdagr,z,dzdn,dzdtau
      real*8 xr,zr,gamma
      real*8 inv_nup,nup_13,nup_43,nup_53,nup_83,agrup2
      real*8 inv_ndn,ndn_13,ndn_43,ndn_53,ndn_83,agrdn2
      real*8 xup,dxup_dnup,dxup_dagrup,zup,dzup_dnup,dzup_dtauup
      real*8 xdn,dxdn_dndn,dxdn_dagrdn,zdn,dzdn_dndn,dzdn_dtaudn
      real*8 GG,dGdx,dGdz
      real*8 n_onethird,kf,ks,ss,P0,FF,Fs,ex_lda,fdnx_const
      real*8 pbex,pbefnx,pbefdnx
      real*8 tauU,t,dt_dn,dt_dtau
      real*8 w,dw_dt,w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11
      real*8 fw,dfw_dw,dfw_dn,dfw_dtau
      real*8 eupx,fnupx,fdnupx,fdtauupx
      real*8 ednx,fndnx,fdndnx,fdtaudnx
      real*8 euu0c,deuu0c_dnup
      real*8 U,dU_dx,U2,U3,U4,gss,dgss_dU
      real*8 hss,dhss_dx,dhss_dz,ghss,dgh_dx,dgh_dn,dgh_dagr,dgh_dtau
      real*8 euuc,dfuuc_dnup,dfuuc_dagrup,dfuuc_dtauup
      real*8 edd0c,dedd0c_dndn
      real*8 eddc,dfddc_dndn,dfddc_dagrdn,dfddc_dtaudn
      real*8 eud0c,deud0c_drs,deud0c_dzeta
      real*8 eudc,eud1c,fud1c,dfud1c_dnup,dfud1c_dndn
      real*8 gopp,dgopp_dU,hopp,dhopp_dx,dhopp_dz
      real*8 ghopp,dgh_dnup,dgh_dndn
      real*8 dgh_dagrup,dgh_dagrdn,dgh_dtauup,dgh_dtaudn
      real*8 dfudc_dnup,dfudc_dagrup,dfudc_dtauup
      real*8 dfudc_dndn,dfudc_dagrdn,dfudc_dtaudn
      real*8 rs,drs_dn,zeta,dzeta_dnup,dzeta_dndn,deuu0c_drs,dedd0c_drs
      real*8 wt1,dwt1,aw,bw,cw,dw
      real*8 D,dD_dx,dD_dz,dD_dn,dD_dagr,dD_dtau,dummy
      real*8 Dt,dD_dDt,dDt_dx,dDt_dz
      real*8 eupc,fnupc,fdnupc,fdtauupc
      real*8 ednc,fndnc,fdndnc,fdtaudnc
      real*8 ex,ec
      real*8 xbar,xb2,xb3,xb4,xb5,dgss_dx,dgopp_dx
*     ***** constants *****
      real*8 pi,cf,thrd,twthrd,frthrd,fvthrd,etthrd
      parameter (pi     = 3.14159265358979311599d0)
      parameter (cf     = 9.115599720d0)
c     cf = 3.0d0/5.0d0*(6.0d0*pi*pi)**(2.0d0/3.0d0))
      parameter (thrd   = 1.0d0/3.0d0)
      parameter (twthrd = 2.0d0/3.0d0)
      parameter (frthrd = 4.0d0/3.0d0)
      parameter (fvthrd = 5.0d0/3.0d0)
      parameter (etthrd = 8.0d0/3.0d0)
*     ***** density cutoff parametersi *****
      real*8 tol,ETA,t1,t2,thresd,thresx
      parameter (tol = 1.0d-10)
      parameter (ETA            =      1.0d-20)
      parameter (t1      =   1.0d7)
      parameter (t2      =   5.0d7)
      parameter (thresd  =   5.0d-3)
      parameter (thresx  =   1.0d8)
c     ***** PBE96 GGA exchange constants ******
      real*8 MU,KAPPA
      parameter (MU    = 0.2195149727645171d0)
      parameter (KAPPA = 0.8040000000000000d0)
*     ***** VS98 constants *****
      real*8 aax,bbx,ccx,ddx,eex,ffx,alphax
      real*8 aass,bbss,ccss,ddss,eess,ffss,alphass
      real*8 aaopp,bbopp,ccopp,ddopp,eeopp,ffopp,alphaopp
      real*8 clda
      parameter (clda     = 0.9305257363491d0)
      parameter (aax      = 6.012244d-1)
      parameter (bbx      = 4.748822d-3)
      parameter (ccx      = -8.635108d-3)
      parameter (ddx      = -9.308062d-6)
      parameter (eex      = 4.482811d-5)
      parameter (ffx      = 0.000000d0)
      parameter (alphax   = 0.00186726d0)
      parameter (aass     = 4.650534d-1)
      parameter (bbss     = 1.617589d-1)
      parameter (ccss     = 1.833657d-1)
      parameter (ddss     = 4.692100d-4)
      parameter (eess     = -4.990573d-3)
      parameter (ffss     = 0.000000d0)
      parameter (alphass  = 0.00515088d0)
      parameter (aaopp    = 3.957626d-1)
      parameter (bbopp    = -5.614546d-1)
      parameter (ccopp    = 1.403963d-2)
      parameter (ddopp    = 9.831442d-4)
      parameter (eeopp    = -3.577176d-3)
      parameter (ffopp    =  0.000000d0)
      parameter (alphaopp =  0.00304966d0)
*     ***** M06 constants *****
      real*8 ma0,ma1,ma2,ma3,ma4,ma5,ma6,ma7,ma8,ma9,ma10,ma11
      real*8 C1,C2
      real*8 css,ccss0,ccss1,ccss2,ccss3,ccss4
      real*8 copp,ccopp0,ccopp1,ccopp2,ccopp3,ccopp4
      parameter (ma0    =  3.987756d-1)
      parameter (ma1    =  2.548219d-1)
      parameter (ma2    =  3.923994d-1)
      parameter (ma3    = -2.103655d0)
      parameter (ma4    = -6.302147d0)
      parameter (ma5    =  1.097615d1)
      parameter (ma6    =  3.097273d1)
      parameter (ma7    = -2.318489d1)
      parameter (ma8    = -5.673480d1)
      parameter (ma9    =  2.160364d1)
      parameter (ma10   =  3.421814d1)
      parameter (ma11   = -9.049762d0)
      parameter (C1     =  3.36116d-3)
      parameter (C2     =  4.49267d-3)
      parameter (css    =  0.06d0)
      parameter (ccss0  =  5.349466d-1)
      parameter (ccss1  =  5.396620d-1)
      parameter (ccss2  = -3.161217d1)
      parameter (ccss3  =  5.149592d1)
      parameter (ccss4  = -2.919613d1)
      parameter (copp   =  0.0031d0)
      parameter (ccopp0 =  6.042374d-1)
      parameter (ccopp1 =  1.776783d2)
      parameter (ccopp2 = -2.513252d2)
      parameter (ccopp3 =  7.635173d1)
      parameter (ccopp4 = -1.255699d1)

      Cx         = -1.50d0*(0.75d0/pi)**thrd
      P23        =  0.60d0*(6.0d0*pi*pi)**twthrd
      fdnx_const = -3.0d0/(8.0d0*pi)

      do i=1,n2ft3d
        nup        = rho_in(i,1) + ETA
        agrup      = agr_in(i,1) + ETA
        tauup      = tau_in(i,1) + ETA
        ndn        = rho_in(i,2) + ETA
        agrdn      = agr_in(i,2) + ETA
        taudn      = tau_in(i,2) + ETA

        n = nup + ndn

*       ***** M06-L Exchange *****
*       ***** Up *****
        agrup2  = agrup*agrup
        inv_nup = 1.0d0/nup
        nup_13  = nup**thrd
        nup_43  = nup_13*nup
        nup_53  = nup_43*nup_13
        nup_83  = nup_53*nup

        xup         =  agrup2/nup_83
        dxup_dnup   = -etthrd*xup*inv_nup
        dxup_dagrup =  2.0d0*agrup/nup_83
        zup         =  tauup/nup_53 - cf
        dzup_dnup   = -fvthrd*tauup/nup_83
        dzup_dtauup =  1.0d0/nup_53

*       ***** VS98 Exchange *****
        gamma = 1.0d0 + alphax*(xup + zup)
        xr    = xup/gamma
        zr    = zup/gamma

        call nwpw_GVT4(aax,bbx,ccx,ddx,eex,ffx,alphax,alphax,
     >              xr,zr,gamma,GG,dGdx,dGdz)

        eupx     = -clda*nup_13*GG
        fnupx    =  frthrd*eupx - clda*nup_43*(dGdx*dxup_dnup 
     &           +  dGdz*dzup_dnup)
        fdnupx   = -clda*nup_43*(dGdx*dxup_dagrup)
        fdtauupx = -clda*nup_43*(dGdz*dzup_dtauup)

*       ***** M06-L Enhancement Factor *****
        tauU    =  P23*nup_53
        t       =  tauU/tauup
        dt_dn   =  fvthrd*t/nup
        dt_dtau = -t/tauup

        w     = 1.0d0
        dw_dt = 0.0d0

        if(t .le. t1) then
          w     = (t - 1.0d0)/(t + 1.0d0)
          dw_dt = 2.0d0/((1.0d0 + t)**2.0d0)
        else if(t .gt. t1 .and. t .lt. t2) then
          wt1 = (t1 - 1.0d0)/(t1 + 1.0d0)
          dwt1 = 2.0d0/((1.0d0 + t1)**2.0d0)
          aw = (3.0d0*t1*t2*t2 - t2**3.0d0)*wt1/(t1 - t2)
     &       + (t1**3.0d0 - 3.0d0*t1*t1*t2)/(t1 - t2)
     &       - t1*t2*t2*dwt1
          aw = aw/((t1 - t2)*(t1 - t2))
          bw = -6.0d0*t1*t2*wt1/(t1 - t2) + 6.0d0*t1*t2/(t1 - t2)
     &       + (t2*t2 + 2.0d0*t1*t2)*dwt1
          bw = bw/((t1 - t2)*(t1 - t2))
          cw = 3.0d0*(t1 + t2)*wt1/(t1 - t2)
     &       - 3.0d0*(t1 + t2)/(t1 - t2) - (t1 + 2.0d0*t2)*dwt1
          cw = cw/((t1 - t2)*(t1 - t2))
          dw = -2.0d0*wt1/(t1 - t2) + 2.0d0/(t1 - t2) + dwt1
          dw = dw/((t1 - t2)*(t1 - t2))

          w     = aw + bw*t + cw*t*t + dw*t*t*t
          dw_dt = bw + 2.0d0*cw*t + 3.0d0*dw*t*t
        else if(t .ge. t2) then
          w     = 1.0d0
          dw_dt = 0.0d0
        end if

        w1  = w
        w2  = w1*w
        w3  = w2*w
        w4  = w3*w
        w5  = w4*w
        w6  = w5*w
        w7  = w6*w
        w8  = w7*w
        w9  = w8*w
        w10 = w9*w
        w11 = w10*w

        fw       = ma0    + ma1*w1 + ma2*w2 + ma3*w3 + ma4*w4   
     &           + ma5*w5 + ma6*w6 + ma7*w7 + ma8*w8 + ma9*w9 
     &           + ma10*w10 + ma11*w11
        dfw_dw   = ma1            + 2.0d0*ma2*w1     + 3.0d0*ma3*w2
     &           + 4.0d0*ma4*w3   + 5.0d0*ma5*w4     + 6.0d0*ma6*w5
     &           + 7.0d0*ma7*w6   + 8.0d0*ma8*w7     + 9.0d0*ma9*w8
     &           + 10.0d0*ma10*w9 + 11.0d0*ma11*w10
        dfw_dn   = dfw_dw*dw_dt*dt_dn
        dfw_dtau = dfw_dw*dw_dt*dt_dtau

*       ***** PBE96 Exchange *****
        n_onethird = (3.0d0*nup/pi)**thrd
        ex_lda     = -0.75d0*n_onethird

        kf = (3.0d0*pi*pi*nup)**thrd
        ss = agrup/(2.0d0*kf*nup)
        P0 = 1.0d0 + (MU/KAPPA)*ss*ss

        FF  = (1.0d0 + KAPPA - KAPPA/P0)
        Fs  = 2.0d0*MU/(P0*P0)*ss

        pbex    = ex_lda*FF
        pbefnx  = frthrd*(pbex - ex_lda*Fs*ss)
        pbefdnx = fdnx_const*Fs

        eupx     = eupx      + pbex*fw
        fnupx    = fnupx     + pbefnx*fw          + nup*pbex*dfw_dn
        fdnupx   = fdnupx    + pbefdnx*fw
        fdtauupx = fdtauupx  + nup*pbex*dfw_dtau


*       ***** Down *****
        agrdn2  = agrdn*agrdn
        inv_ndn = 1.0d0/ndn
        ndn_13  = ndn**thrd
        ndn_43  = ndn_13*ndn
        ndn_53  = ndn_43*ndn_13
        ndn_83  = ndn_53*ndn

        xdn         =  agrdn2/ndn_83
        dxdn_dndn   = -etthrd*xdn*inv_ndn
        dxdn_dagrdn =  2.0d0*agrdn/ndn_83
        zdn         =  taudn/ndn_53 - cf
        dzdn_dndn   = -fvthrd*taudn/ndn_83
        dzdn_dtaudn =  1.0d0/ndn_53

        gamma = 1.0d0 + alphax*(xdn + zdn)
        xr    = xdn/gamma
        zr    = zdn/gamma

        call nwpw_GVT4(aax,bbx,ccx,ddx,eex,ffx,alphax,alphax,
     >              xr,zr,gamma,GG,dGdx,dGdz)

        ednx     = -clda*ndn_13*GG
        fndnx    =  frthrd*ednx - clda*ndn_43*(dGdx*dxdn_dndn 
     &           +  dGdz*dzdn_dndn)
        fdndnx   = -clda*ndn_43*(dGdx*dxdn_dagrdn)
        fdtaudnx = -clda*ndn_43*(dGdz*dzdn_dtaudn)

*       ***** M06-L Enhancement Factor *****
        tauU    =  P23*ndn_53
        t       =  tauU/taudn
        dt_dn   =  fvthrd*t/ndn
        dt_dtau = -t/taudn

        w     = 1.0d0
        dw_dt = 0.0d0

        if(t .le. t1) then
          w     = (t - 1.0d0)/(t + 1.0d0)
          dw_dt = 2.0d0/((1.0d0 + t)**2.0d0)
        else if(t .gt. t1 .and. t .lt. t2) then
          wt1 = (t1 - 1.0d0)/(t1 + 1.0d0)
          dwt1 = 2.0d0/((1.0d0 + t1)**2.0d0)
          aw = (3.0d0*t1*t2*t2 - t2**3.0d0)*wt1/(t1 - t2)
     &       + (t1**3.0d0 - 3.0d0*t1*t1*t2)/(t1 - t2)
     &       - t1*t2*t2*dwt1
          aw = aw/((t1 - t2)*(t1 - t2))
          bw = -6.0d0*t1*t2*wt1/(t1 - t2) + 6.0d0*t1*t2/(t1 - t2)
     &       + (t2*t2 + 2.0d0*t1*t2)*dwt1
          bw = bw/((t1 - t2)*(t1 - t2))
          cw = 3.0d0*(t1 + t2)*wt1/(t1 - t2)
     &       - 3.0d0*(t1 + t2)/(t1 - t2) - (t1 + 2.0d0*t2)*dwt1
          cw = cw/((t1 - t2)*(t1 - t2))
          dw = -2.0d0*wt1/(t1 - t2) + 2.0d0/(t1 - t2) + dwt1
          dw = dw/((t1 - t2)*(t1 - t2))

          w     = aw + bw*t + cw*t*t + dw*t*t*t
          dw_dt = bw + 2.0d0*cw*t + 3.0d0*dw*t*t
        else if(t .ge. t2) then
          w     = 1.0d0
          dw_dt = 0.0d0
        end if

        w1  = w
        w2  = w1*w
        w3  = w2*w
        w4  = w3*w
        w5  = w4*w
        w6  = w5*w
        w7  = w6*w
        w8  = w7*w
        w9  = w8*w
        w10 = w9*w
        w11 = w10*w

        fw       = ma0    + ma1*w1 + ma2*w2 + ma3*w3 + ma4*w4   
     &           + ma5*w5 + ma6*w6 + ma7*w7 + ma8*w8 + ma9*w9 
     &           + ma10*w10 + ma11*w11
        dfw_dw   = ma1            + 2.0d0*ma2*w1     + 3.0d0*ma3*w2
     &           + 4.0d0*ma4*w3   + 5.0d0*ma5*w4     + 6.0d0*ma6*w5
     &           + 7.0d0*ma7*w6   + 8.0d0*ma8*w7     + 9.0d0*ma9*w8
     &           + 10.0d0*ma10*w9 + 11.0d0*ma11*w10
        dfw_dn   = dfw_dw*dw_dt*dt_dn
        dfw_dtau = dfw_dw*dw_dt*dt_dtau

*       ***** PBE96 Exchange *****
        n_onethird = (3.0d0*ndn/pi)**thrd
        ex_lda     = -0.75d0*n_onethird

        kf = (3.0d0*pi*pi*ndn)**thrd
        ss = agrdn/(2.0d0*kf*ndn)
        P0 = 1.0d0 + (MU/KAPPA)*ss*ss

        FF  = (1.0d0 + KAPPA - KAPPA/P0)
        Fs  = 2.0d0*MU/(P0*P0)*ss

        pbex    = ex_lda*FF
        pbefnx  = frthrd*(pbex - ex_lda*Fs*ss)
        pbefdnx = fdnx_const*Fs

        ednx     = ednx      + pbex*fw
        fndnx    = fndnx     + pbefnx*fw          + ndn*pbex*dfw_dn
        fdndnx   = fdndnx    + pbefdnx*fw
        fdtaudnx = fdtaudnx  + ndn*pbex*dfw_dtau

        ex = (eupx*nup + ednx*ndn)/n

*       ***** M06-L Correlation *****
*       ***** Same-Spin *****
*       ***** alpha-alpha *****
        rs     = (0.75d0/(pi*nup))**thrd
        drs_dn = -thrd*rs/nup

        call gen_PW91_c_rz(tol,rs,1.0d0,euu0c,deuu0c_drs,dummy)
        deuu0c_dnup = deuu0c_drs*drs_dn

        Dt      = 1.0d0 - 0.25d0*xup/(zup + cf)
        D       = 0.0d0
        dD_dn   = 0.0d0
        dD_dagr = 0.0d0
        dD_dtau = 0.0d0

        if (Dt .le. 0.0d0) then
           D       = 0.0d0
           dD_dn   = 0.0d0
           dD_dagr = 0.0d0
           dD_dtau = 0.0d0
        else if (Dt .gt. 0.0d0 .and. Dt .lt. thresd) then
           D       =  2.0d0*Dt*Dt/thresd - Dt*Dt*Dt/(thresd*thresd)
           dD_dDt  =  4.0d0*Dt/thresd - 3.0d0*Dt*Dt/(thresd*thresd)
           dDt_dx  = -0.25d0/(zup + cf)
           dDt_dz  =  0.25d0*xup/((zup + cf)*(zup + cf))
           dD_dn   =  dD_dDt*(dDt_dx*dxup_dnup + dDt_dz*dzup_dnup)
           dD_dagr =  dD_dDt*dDt_dx*dxup_dagrup
           dD_dtau =  dD_dDt*dDt_dz*dzup_dtauup
        else if(Dt .ge. thresd) then
           D       =  Dt
           dD_dx   = -0.25d0/(zup + cf)
           dD_dz   =  0.25d0*xup/((zup + cf)*(zup + cf))
           dD_dn   =  dD_dx*dxup_dnup + dD_dz*dzup_dnup
           dD_dagr =  dD_dx*dxup_dagrup
           dD_dtau =  dD_dz*dzup_dtauup
        end if

        gamma = 1.0d0 + alphass*(xup + zup)
        xr    = xup/gamma
        zr    = zup/gamma

        if (xup .ge. thresx) then
           U     = 1.0d0
           dU_dx = 0.0d0
        else
           U     = css*xup/(1.0d0 + css*xup)
           dU_dx = css/((1.0d0 + css*xup)*(1.0d0 + css*xup))
        end if

        U2 = U*U
        U3 = U2*U
        U4 = U3*U

        gss     = ccss0 + ccss1*U + ccss2*U2 + ccss3*U3 + ccss4*U4
        dgss_dU = ccss1 + 2.0d0*ccss2*U + 3.0d0*ccss3*U2 
     &          + 4.0d0*ccss4*U3
        dgss_dx = dgss_dU*dU_dx

        call nwpw_GVT4(aass,bbss,ccss,ddss,eess,ffss,alphass,alphass,
     >                 xr,zr,gamma,hss,dhss_dx,dhss_dz)

        ghss     = gss + hss
        dgh_dx   = dgss_dx + dhss_dx
        dgh_dn   = dgh_dx*dxup_dnup + dhss_dz*dzup_dnup
        dgh_dagr = dgh_dx*dxup_dagrup 
        dgh_dtau = dhss_dz*dzup_dtauup

        eupc     = euu0c*ghss*D
        fnupc    = eupc + nup*deuu0c_drs*drs_dn*ghss*D
     &           + nup*euu0c*dgh_dn*D 
     &           + nup*euu0c*ghss*dD_dn
        fdnupc   = nup*euu0c*(dgh_dagr*D + ghss*dD_dagr)
        fdtauupc = nup*euu0c*(dgh_dtau*D + ghss*dD_dtau)

*       ***** beta-beta *****
        rs     = (0.75d0/(pi*ndn))**thrd
        drs_dn = -thrd*rs/ndn

        call gen_PW91_c_rz(tol,rs,1.0d0,edd0c,dedd0c_drs,dummy)
        dedd0c_dndn = dedd0c_drs*drs_dn

        Dt      = 1.0d0 - 0.25d0*xdn/(zdn + cf)
        D       = 0.0d0
        dD_dn   = 0.0d0
        dD_dagr = 0.0d0
        dD_dtau = 0.0d0

        if (Dt .le. 0.0d0) then
           D       = 0.0d0
           dD_dn   = 0.0d0
           dD_dagr = 0.0d0
           dD_dtau = 0.0d0
        else if (Dt .gt. 0.0d0 .and. Dt .lt. thresd) then
           D       =  2.0d0*Dt*Dt/thresd - Dt*Dt*Dt/(thresd*thresd)
           dD_dDt  =  4.0d0*Dt/thresd - 3.0d0*Dt*Dt/(thresd*thresd)
           dDt_dx  = -0.25d0/(zdn + cf)
           dDt_dz  =  0.25d0*xdn/((zdn + cf)*(zdn + cf))
           dD_dn   =  dD_dDt*(dDt_dx*dxdn_dndn + dDt_dz*dzdn_dndn)
           dD_dagr =  dD_dDt*dDt_dx*dxdn_dagrdn
           dD_dtau =  dD_dDt*dDt_dz*dzdn_dtaudn
        else if(Dt .ge. thresd) then
           D       =  Dt
           dD_dx   = -0.25d0/(zdn + cf)
           dD_dz   =  0.25d0*xdn/((zdn + cf)*(zdn + cf))
           dD_dn   =  dD_dx*dxdn_dndn + dD_dz*dzdn_dndn
           dD_dagr =  dD_dx*dxdn_dagrdn
           dD_dtau =  dD_dz*dzdn_dtaudn
        end if

        gamma = 1.0d0 + alphass*(xdn + zdn)
        xr    = xdn/gamma
        zr    = zdn/gamma

        if (xdn .ge. thresx) then
           U     = 1.0d0
           dU_dx = 0.0d0
        else
           U     = css*xdn/(1.0d0 + css*xdn)
           dU_dx = css/((1.0d0 + css*xdn)*(1.0d0 + css*xdn))
        end if

        U2 = U*U
        U3 = U2*U
        U4 = U3*U
 
        gss     = ccss0 + ccss1*U + ccss2*U2 + ccss3*U3 + ccss4*U4
        dgss_dU = ccss1 + 2.0d0*ccss2*U + 3.0d0*ccss3*U2 
     &          + 4.0d0*ccss4*U3
        dgss_dx = dgss_dU*dU_dx

        call nwpw_GVT4(aass,bbss,ccss,ddss,eess,ffss,alphass,alphass,
     >                 xr,zr,gamma,hss,dhss_dx,dhss_dz)

        ghss     = gss + hss
        dgh_dx   = dgss_dx + dhss_dx
        dgh_dn   = dgh_dx*dxdn_dndn + dhss_dz*dzdn_dndn
        dgh_dagr = dgh_dx*dxdn_dagrdn 
        dgh_dtau = dhss_dz*dzdn_dtaudn

        ednc     = edd0c*ghss*D
        fndnc    = ednc + ndn*dedd0c_drs*drs_dn*ghss*D
     &           + ndn*edd0c*dgh_dn*D 
     &           + ndn*edd0c*ghss*dD_dn
        fdndnc   = ndn*edd0c*(dgh_dagr*D + ghss*dD_dagr)
        fdtaudnc = ndn*edd0c*(dgh_dtau*D + ghss*dD_dtau)

        ec = (eupc*nup + ednc*ndn)/n

*       ***** Opposite-Spin (alpha-beta) *****
        rs         =  (0.75d0/(pi*n))**thrd
        drs_dn     = -thrd*rs/n
        zeta       =  (nup - ndn)/n
        dzeta_dnup =  ( 1.0d0 - zeta)/n
        dzeta_dndn =  (-1.0d0 - zeta)/n

        call gen_PW91_c_rz(tol,rs,zeta,eud0c,deud0c_drs,deud0c_dzeta)
        eud1c       = eud0c - (euu0c*nup + edd0c*ndn)/n
        fud1c       = n*eud1c
        dfud1c_dnup = eud0c + n*(deud0c_drs*drs_dn
     &              + deud0c_dzeta*dzeta_dnup) - euu0c
     &              - nup*deuu0c_dnup
        dfud1c_dndn = eud0c + n*(deud0c_drs*drs_dn
     &              + deud0c_dzeta*dzeta_dndn) - edd0c
     &              - ndn*dedd0c_dndn

        x   = xup + xdn
        z   = zup + zdn

        gamma = 1.0d0 + alphaopp*(x + z)
        xr    = x/gamma
        zr    = z/gamma

        if (x .ge. thresx) then
           U     = 1.0d0
           dU_dx = 0.0d0
        else
           U     = copp*x/(1.0d0 + copp*x)
           dU_dx = copp/((1.0d0 + copp*x)*(1.0d0 + copp*x))
        end if

        U2 = U*U
        U3 = U2*U
        U4 = U3*U
 
        gopp     = ccopp0 + ccopp1*U + ccopp2*U2 + ccopp3*U3 + ccopp4*U4
        dgopp_dU = ccopp1          + 2.0d0*ccopp2*U
     &           + 3.0d0*ccopp3*U2 + 4.0d0*ccopp4*U3
        dgopp_dx = dgopp_dU*dU_dx

        call nwpw_GVT4(aaopp,bbopp,ccopp,ddopp,eeopp,ffopp,
     >                 alphaopp,alphaopp,
     >                 xr,zr,gamma,hopp,dhopp_dx,dhopp_dz)

        ghopp      = gopp + hopp
        dgh_dx     = dgopp_dx + dhopp_dx
        dgh_dnup   = dgh_dx*dxup_dnup + dhopp_dz*dzup_dnup
        dgh_dndn   = dgh_dx*dxdn_dndn + dhopp_dz*dzdn_dndn
        dgh_dagrup = dgh_dx*dxup_dagrup
        dgh_dagrdn = dgh_dx*dxdn_dagrdn
        dgh_dtauup = dhopp_dz*dzup_dtauup
        dgh_dtaudn = dhopp_dz*dzdn_dtaudn
        
        eudc         = eud1c*ghopp
        dfudc_dnup   = dfud1c_dnup*ghopp + fud1c*dgh_dnup
        dfudc_dndn   = dfud1c_dndn*ghopp + fud1c*dgh_dndn
        dfudc_dagrup = fud1c*dgh_dagrup
        dfudc_dagrdn = fud1c*dgh_dagrdn
        dfudc_dtauup = fud1c*dgh_dtauup
        dfudc_dtaudn = fud1c*dgh_dtaudn

        ec       = ec       + eudc
        fnupc    = fnupc    + dfudc_dnup
        fndnc    = fndnc    + dfudc_dndn
        fdnupc   = fdnupc   + dfudc_dagrup
        fdndnc   = fdndnc   + dfudc_dagrdn
        fdtauupc = fdtauupc + dfudc_dtauup
        fdtaudnc = fdtaudnc + dfudc_dtaudn

        xce(i)     = x_parameter*ex       + c_parameter*ec
        fn(i,1)    = x_parameter*fnupx    + c_parameter*fnupc
        fn(i,2)    = x_parameter*fndnx    + c_parameter*fndnc
        fdn(i,1)   = x_parameter*fdnupx   + c_parameter*fdnupc
        fdn(i,2)   = x_parameter*fdndnx   + c_parameter*fdndnc
        fdn(i,3)   = 0.0d0
        fdtau(i,1) = x_parameter*fdtauupx + c_parameter*fdtauupc
        fdtau(i,2) = x_parameter*fdtaudnx + c_parameter*fdtaudnc
      end do

      return
      end

