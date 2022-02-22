*
* $Id$
*

*     ************************************************
*     *                                              *
*     *              gen_SCAN_restricted             *
*     *                                              *
*     ************************************************

*    This function returns the SCAN exchange-correlation
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
*     Exit  - xce(n2ft3d) : SCAN exchange correlation energy density
*             fn(n2ft3d)  : d(n*xce)/dn
*             fdn(n2ft3d) : d(n*xce)/d|grad n|
*             fdtau(n2ft3d) : d(n*xce)/dtau
*

      subroutine gen_SCAN_restricted(n2ft3d,rho_in,agr_in,tau_in,
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
      real*8 Cx,P23,P13,P23t
      real*8 ex,fnx,fdnx,fdtaux
      real*8 n_13,n_53,n_83,agr2,tauW,tauU
      real*8 p,p_14,dp_dn,dp_dagr
      real*8 z,z2,fz,dz_dn,dz_dagr,dz_dtau
      real*8 alpha,dalpha_dn,dalpha_dagr,dalpha_dtau
      real*8 thr0,dthr0_dagr,dthr0_dtau,delta2,diff,diff2,diff3,dod,dod2
      real*8 tauUn,falpha,dfalpha_dn,dfalpha_dagr,dfalpha_dtau
      real*8 oma,oma2
      real*8 exp1,exp2,exp3,exp4,exp5
      real*8 x1,x2,x,dx1_dp,dx2_dp,dx_dp,dx_dalpha
      real*8 denh1x,numh1x,h1x,dh1x_dx,dh1x_dp,dh1x_dalpha
      real*8 gx,dgx_dp
      real*8 fxa,dfxa_dalpha
      real*8 Fx,dFx_dp,dFx_dalpha,dFx_dn,dFx_dagr,dFx_dtau
      real*8 ex0,nex0
      real*8 rs,drs_dn,rs_12
      real*8 zeta
      real*8 ecLDA1,decLDA1_drs,decLDA1_dzeta
      real*8 beta,dbeta_drs,w1fac,expw1,w1,dw1_drs
      real*8 A,dA_drs,t2,dt2_drs,dt2_dp
      real*8 At2,dengAt2,gAt2,tmp1,dtmp1_drs,dtmp1_dp
      real*8 H1,dH1_drs,dH1_dp
      real*8 ec1,dec1_drs,dec1_dp
      real*8 ecLDA0,decLDA0_drs
      real*8 w0fac,expw0,w0,dw0_drs,ginf,dginf_dp
      real*8 tmp0,dtmp0_drs,dtmp0_dp
      real*8 H0,dH0_drs,dH0_dp
      real*8 ec0,dec0_drs,dec0_dp
      real*8 exp6,exp7,fca,dfca_dalpha
      real*8 dec1_dn,dec1_dagr
      real*8 dec0_dn,dec0_dagr
      real*8 dfca_dn,dfca_dagr,dfca_dtau
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

*     ***** SCAN constants *****
      real*8 a1,b1,b2,b3,b4,c1x,c2x,dx,muAK,K1,h0x 
      real*8 gamma,beta0,beta1,beta2,b1c,b2c,b3c,c1c,c2c,dc,dxc,xi
      real*8 thr1,thr2,delta
      parameter (a1    = 4.9479d0)
      parameter (b3    = 0.5d0)
      parameter (c1x   = 0.667d0)
      parameter (c2x   = 0.8d0)
      parameter (dx    = 1.24d0)
      parameter (muAK  = 10.0d0/81.0d0)
      parameter (K1    = 0.065d0)
      parameter (h0x   = 1.174d0)
      parameter (gamma = 0.03109069086965489503494086371273d0)
      parameter (beta0 = 0.06672455060314922d0)
      parameter (beta1 = 0.1d0)
      parameter (beta2 = 0.1778d0)
      parameter (b1c   = 0.0285764d0)
      parameter (b2c   = 0.0889d0)
      parameter (b3c   = 0.125541d0)
      parameter (c1c   = 0.64d0)
      parameter (c2c   = 1.5d0)
      parameter (dc    = 0.7d0)
      parameter (dxc   = 2.3631d0)
      parameter (xi    = 0.12802585262625815d0)
      parameter (thr1  = 0.996d0)
      parameter (thr2  = 1.004d0)
      parameter (delta = 5.0d-5)

      b2 = dsqrt(5913.0d0/405000.0d0)
      b1 = (511.0d0/13500.0d0)/(2.0d0*b2)
      b4 = muAK*muAK/K1 - 1606.0d0/18225.0d0 - b1*b1

      Cx   = (-0.75d0)*(3.0d0/pi)**thrd
      P23  = (3.0d0*pi**2.0d0)**twthrd
      P13  = (3.0d0/(4.0d0*pi))**thrd
      P23t = (3.0d0*pi*pi/16.0d0)**twthrd

      do i=1,n2ft3d
        n       =       rho_in(i) + ETA
        agr     =       agr_in(i) + ETA
        tau     = 2.0d0*tau_in(i) + ETA

        n_13  = n**thrd
        n_53  = n_13*n_13*n
        n_83  = n_53*n
        agr2  = agr*agr

        p       =  agr2/(4.0d0*P23*n_83)
        p_14    =  dsqrt(dsqrt(p))
        dp_dn   = -etthrd*p/n
        dp_dagr =  2.0d0*p/agr
c       dp_dtau =  0.0d0

        tauW  = 0.125d0*agr2/n
        tauU  = 0.3d0*P23*n_53
        z    = tauW/tau
        dz_dn   = -z/n
        dz_dagr =  2.0d0*z/agr
        dz_dtau = -z/tau

        z2 = z*z

        alpha       = 0.0d0
        dalpha_dn   = 0.0d0
        dalpha_dagr = 0.0d0
        dalpha_dtau = 0.0d0

        thr0 = agr2/(8.0d0*tau)
        dthr0_dagr = 2.0d0*thr0/agr
        dthr0_dtau = -thr0/tau
        delta2 = delta*delta
        diff = n - thr0
        diff2 = diff*diff
        diff3 = diff2*diff
        dod = diff/delta
        dod2 = dod*dod
        tauUn = tauU*n

        if (n .le. thr0) then
          alpha       = 0.0d0
          dalpha_dn   = 0.0d0
          dalpha_dagr = 0.0d0
          dalpha_dtau = 0.0d0
        else if ((n .gt. thr0) .and. (n .lt. (thr0 + delta))) then
          falpha = tau*(-diff3/delta2 + 2.0d0*diff2/delta)
          dfalpha_dn =  tau*(-3.0d0*dod2 + 4.0d0*dod)
          dfalpha_dagr = tau*(3.0d0*dod2 - 4.0d0*dod)*dthr0_dagr
          dfalpha_dtau = falpha/tau
     &                 + tau*(3.0d0*dod2 - 4.0d0*dod)*dthr0_dtau
          alpha = falpha/tauUn
          dalpha_dn = -etthrd*falpha/(tauUn*n) + dfalpha_dn/tauUn
          dalpha_dagr = dfalpha_dagr/tauUn
          dalpha_dtau = dfalpha_dtau/tauUn
        else if (n .gt. (thr0 + delta)) then
          falpha = tau*n - agr2/8.0d0
          dfalpha_dn = tau
          dfalpha_dagr = -agr/4.0d0
          dfalpha_dtau = n
          alpha = falpha/tauUn
          dalpha_dn = -etthrd*falpha/(tauUn*n) + dfalpha_dn/tauUn
          dalpha_dagr = dfalpha_dagr/tauUn
          dalpha_dtau = dfalpha_dtau/tauUn
        end if

        oma  = 1.0d0 - alpha
        oma2 = oma*oma

*       ***** SCAN Exchange *****
        exp1 = dexp(-b4*p/muAK)
        exp2 = dexp(-b3*oma2)

        x1 = muAK*p*(1.0d0 + (b4*p/muAK)*exp1)
        x2 = b1*p + b2*oma*exp2
        x = x1 + x2*x2

        denh1x = K1 + x
        numh1x = denh1x + K1*x
        h1x    = numh1x/denh1x

        dx1_dp    = muAK + b4*p*exp1*(2.0d0 - p*b4/muAK)
        dx2_dp    = b1
        dx_dp     = dx1_dp + 2.0d0*x2*dx2_dp
        dx_dalpha = 2.0d0*b2*exp2*x2*(2.0d0*b3*oma2 - 1.0d0)

        dh1x_dx     = (K1/denh1x)**2.0d0
        dh1x_dp     = dh1x_dx*dx_dp
        dh1x_dalpha = dh1x_dx*dx_dalpha

        gx     = 1.0d0
        dgx_dp =  0.0d0
        if (p_14 .lt. 0.002d0) then
           exp3   = 0.0d0
           gx     = 1.0d0
           dgx_dp = 0.0d0
        else
           exp3   = dexp(-a1/p_14)
           gx     = 1.0d0 - exp3
           dgx_dp = -0.25d0*a1*exp3/(p*p_14)
        endif

        fxa = 0.d0
        dfxa_dalpha = 0.0d0

        if (alpha .le. thr1) then
           exp4 = dexp(-c1x*alpha/oma)
           exp5 = 0.0d0
           fxa  = exp4
           dfxa_dalpha = -c1x*exp4/oma2
        else if ((alpha .gt. thr1) .and. (alpha .lt. thr2)) then
           exp4 = 0.0d0
           exp5 = 0.0d0
           fxa = 0.0d0
           dfxa_dalpha = 0.0d0
        else if (alpha .ge. thr2) then
           exp4 = 0.0d0
           exp5 = dexp(c2x/oma)
           fxa = -dx*exp5
           dfxa_dalpha = -dx*c2x*exp5/oma2
        end if

        Fx = (h1x + fxa*(h0x - h1x))*gx

        dFx_dp     = dgx_dp*(h1x + fxa*(h0x - h1x)) 
     &             + gx*dh1x_dp*(1.0d0 - fxa)
        dFx_dalpha = gx*(dh1x_dalpha + dfxa_dalpha*(h0x - h1x) 
     &             - fxa*dh1x_dalpha)

        dFx_dn   = dFx_dalpha*dalpha_dn   + dFx_dp*dp_dn
        dFx_dagr = dFx_dalpha*dalpha_dagr + dFx_dp*dp_dagr
        dFx_dtau = dFx_dalpha*dalpha_dtau

        ex0  = Cx*n_13
        nex0 = n*ex0

        ex      = ex0*Fx
        fnx     = nex0*dFx_dn     + frthrd*ex
        fdnx    = nex0*dFx_dagr
        fdtaux  = nex0*dFx_dtau

*       ***** SCAN Correlation *****
        zeta = 0.0d0

        rs     =  P13/n_13
        drs_dn = -thrd*rs/n
        rs_12  =  dsqrt(rs)

        beta      = beta0*(1.0d0 + beta1*rs)/(1.0d0 + beta2*rs)
        dbeta_drs = beta0*(beta1 - beta2)/((1.0d0 + beta2*rs)**2.0d0)

        call gen_PW91_c_rz(tol,rs,zeta,ecLDA1,decLDA1_drs,decLDA1_dzeta)

        w1fac   =  ecLDA1/gamma
        expw1   =  dexp(-w1fac)
        w1      =  expw1 - 1.0d0
        dw1_drs = -expw1*decLDA1_drs/gamma

        A      = beta/(gamma*w1)
        dA_drs = dbeta_drs/(gamma*w1) - A*dw1_drs/w1

        t2      =  P23t*p/rs
        dt2_drs = -t2/rs
        dt2_dp  =  P23t/rs

        At2     = A*t2
        dengAt2 = 1.0d0 + 4.0d0*At2
        gAt2    = 1.0d0/(dsqrt(dsqrt(dengAt2)))

        tmp1      = 1.0d0 + w1*(1.0d0 - gAt2)
        dtmp1_drs = dw1_drs*(1.0d0 - gAt2) +
     &              gAt2*w1*(t2*dA_drs + A*dt2_drs)/dengAt2
        dtmp1_dp  = gAt2*w1*A*dt2_dp/dengAt2


        H1      = gamma*dlog(tmp1)
        dH1_drs = gamma*dtmp1_drs/tmp1
        dH1_dp  = gamma*dtmp1_dp/tmp1

        ec1      = ecLDA1 + H1
        dec1_drs = decLDA1_drs + dH1_drs
        dec1_dp  = dH1_dp

        ecLDA0      = -b1c/(1.0d0 + b2c*rs_12 + b3c*rs)
        decLDA0_drs =  b1c*(b3c + 0.5d0*b2c/rs_12)/
     &                 ((1.0d0 + b2c*rs_12 + b3c*rs)**2.0d0)

        w0fac   =  ecLDA0/b1c
        expw0   =  dexp(-w0fac)
        w0      =  expw0 - 1.0d0
        dw0_drs = -decLDA0_drs*expw0/b1c

        ginf     =  1.0d0/(dsqrt(dsqrt(1.0d0 + 4.0d0*xi*p)))
        dginf_dp = -xi*ginf/(1.0d0 + 4.0d0*xi*p)

        tmp0      =  1.0d0 + w0*(1.0d0 - ginf)
        dtmp0_drs =  dw0_drs*(1.0d0 - ginf)
        dtmp0_dp  = -w0*dginf_dp

        H0      = b1c*dlog(tmp0)
        dH0_drs = b1c*dtmp0_drs/tmp0
        dH0_dp  = b1c*dtmp0_dp/tmp0

        ec0      = ecLDA0 + H0
        dec0_drs = decLDA0_drs + dH0_drs
        dec0_dp  = dH0_dp

c       alpha = fvthrd*p*(1.0d0/z - 1.0d0)
c
c       if (alpha .le. 0.0d0) then
c         alpha       = 0.0d0
c         dalpha_dn   = 0.0d0
c         dalpha_dagr = 0.0d0
c         dalpha_dtau = 0.0d0
c       else
c         dalpha_dn   = fvthrd*(-p*dz_dn/z2 + dp_dn*(1.0d0/z - 1.0d0))
c         dalpha_dagr = (alpha/p)*dp_dagr - fvthrd*(p/z2)*dz_dagr
c         dalpha_dtau = 1.0d0/tauU
c         dalpha_dtau = fvthrd*p*(-1.0d0/z2)*dz_dtau
c       end if
c
c       oma  = 1.0d0 - alpha
c       oma2 = oma*oma

        fca = 0.d0
        dfca_dalpha = 0.0d0

        if (alpha .le. thr1) then
           exp6 = dexp(-c1c*alpha/oma)
           exp7 = 0.0d0
           fca  = exp6
           dfca_dalpha = -c1c*exp6/oma2
        else if ((alpha .ge. thr1) .and. (alpha .le. thr2)) then
           exp6 = 0.0d0
           exp7 = 0.0d0
           fca  = 0.0d0
           dfca_dalpha = 0.0d0
        else if (alpha .gt. thr2) then
           exp6 = 0.0d0
           exp7 = dexp(c2c/oma)
           fca = -dc*exp7
           dfca_dalpha = -dc*c2c*exp7/oma2
        end if
 
        dec1_dn   = dec1_drs*drs_dn + dec1_dp*dp_dn
        dec1_dagr = dec1_dp*dp_dagr
        
        dec0_dn   = dec0_drs*drs_dn + dec0_dp*dp_dn
        dec0_dagr = dec0_dp*dp_dagr

        dfca_dn   = dfca_dalpha*dalpha_dn
        dfca_dagr = dfca_dalpha*dalpha_dagr
        dfca_dtau = dfca_dalpha*dalpha_dtau
          
        ec     = ec1 + fca*(ec0 - ec1)
        fnc    = n*(dec1_dn + dfca_dn*(ec0 - ec1) 
     &         + fca*(dec0_dn - dec1_dn)) + ec
        fdnc   = n*(dec1_dagr + dfca_dagr*(ec0 - ec1) 
     &         + fca*(dec0_dagr - dec1_dagr)) 
        fdtauc = n*dfca_dtau*(ec0 - ec1)

        
        xce(i)   = x_parameter*ex     + c_parameter*ec
        fn(i)    = x_parameter*fnx    + c_parameter*fnc
        fdn(i)   = x_parameter*fdnx   + c_parameter*fdnc
        fdtau(i) = x_parameter*fdtaux + c_parameter*fdtauc

        end do

      return
      end

*     ************************************************
*     *                                              *
*     *            gen_SCAN_unrestricted             *
*     *                                              *
*     ************************************************

*    This function returns the SCAN exchange-correlation
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
*     Exit  - xce(n2ft3d) : SCAN exchange correlation energy density
*             fn(n2ft3d,2)  : d(n*xce)/dnup, d(n*xce)/dndn
*             fdn(n2ft3d,3) : d(n*xce)/d|grad nup|, d(n*xce)/d|grad ndn|, d(n*xce)/d|grad n|
*             fdtau(n2ft3d,2) : d(n*xce)/dtauup, d(n*xce)/dtaudn
*

      subroutine gen_SCAN_unrestricted(n2ft3d,rho_in,agr_in,tau_in,
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
      real*8 Cx,P23,P13,P23t
      real*8 n_13,n_53,n_83,n2,agr2,tauW,tauU
      real*8 p,p_14,dp_dn,dp_dagr
      real*8 z,z2,fz,dz_dn,dz_dagr,dz_dtau
      real*8 alpha,dalpha_dn,dalpha_dagr,dalpha_dtau
      real*8 thr0,dthr0_dagr,dthr0_dtau,delta2,diff,diff2,diff3,dod,dod2
      real*8 tauUn,falpha,dfalpha_dn,dfalpha_dagr,dfalpha_dtau
      real*8 oma,oma2
      real*8 exp1,exp2,exp3,exp4,exp5
      real*8 x1,x2,x,dx1_dp,dx2_dp,dx_dp,dx_dalpha
      real*8 denh1x,numh1x,h1x,dh1x_dx,dh1x_dp,dh1x_dalpha
      real*8 gx,dgx_dp
      real*8 fxa,dfxa_dalpha
      real*8 Fx,dFx_dp,dFx_dalpha,dFx_dn,dFx_dagr,dFx_dtau
      real*8 ex0,nex0
      real*8 ex,eupx,fnupx,fdnupx,fdtauupx,ednx,fndnx,fdndnx,fdtaudnx
      real*8 rs,drs_dn,rs_12
      real*8 zeta,dzeta_dnup,dzeta_dndn
      real*8 opz,omz,opz_23,omz_23,phi,phi2,phi3,dphi_dzeta
      real*8 ecLDA1,decLDA1_drs,decLDA1_dzeta
      real*8 beta,dbeta_drs,w1fac,expw1,w1,dw1_drs,dw1_dzeta
      real*8 A,dA_drs,dA_dzeta,t2,dt2_drs,dt2_dzeta,dt2_dp
      real*8 At2,dengAt2,gAt2,tmp1,dtmp1_drs,dtmp1_dzeta,dtmp1_dp
      real*8 H1,dH1_drs,dH1_dzeta,dH1_dp
      real*8 ec1,dec1_dzeta,dec1_drs,dec1_dp
      real*8 zeta12,omz12,zeta11
      real*8 ecLDA0,decLDA0_drs,dxz,ddxz_dzeta,gc,dgc_dzeta
      real*8 w0fac,expw0,w0,dw0_drs,ginf,dginf_dp
      real*8 tmp0,dtmp0_drs,dtmp0_dp
      real*8 H0,dH0_drs,dH0_dp
      real*8 ec0,dec0_drs,dec0_dzeta,dec0_dp
      real*8 ds,dds_dzeta
      real*8 dalpha_dzeta,tmpa1,tmpa2,dalpha_dnup,dalpha_dndn
      real*8 dalpha_dtauup,dalpha_dtaudn
      real*8 exp6,exp7,fca,dfca_dalpha
      real*8 dec1_dnup,dec1_dndn,dec1_dagr
      real*8 dec0_dnup,dec0_dndn,dec0_dagr
      real*8 dfca_dnup,dfca_dndn,dfca_dagr,dfca_dtauup,dfca_dtaudn
      real*8 ec,fnupc,fndnc,fdnupc,fdndnc,fdnc,fdtauupc,fdtaudnc
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
*     ***** SCAN constants *****
      real*8 a1,b1,b2,b3,b4,c1x,c2x,dx,muAK,K1,h0x 
      real*8 gamma,beta0,beta1,beta2,b1c,b2c,b3c,c1c,c2c,dc,dxc,xi
      real*8 thr1,thr2,delta
      parameter (a1    = 4.9479d0)
      parameter (b3    = 0.5d0)
      parameter (c1x   = 0.667d0)
      parameter (c2x   = 0.8d0)
      parameter (dx    = 1.24d0)
      parameter (muAK  = 10.0d0/81.0d0)
      parameter (K1    = 0.065d0)
      parameter (h0x   = 1.174d0)
      parameter (gamma = 0.03109069086965489503494086371273d0)
      parameter (beta0 = 0.06672455060314922d0)
      parameter (beta1 = 0.1d0)
      parameter (beta2 = 0.1778d0)
      parameter (b1c   = 0.0285764d0)
      parameter (b2c   = 0.0889d0)
      parameter (b3c   = 0.125541d0)
      parameter (c1c   = 0.64d0)
      parameter (c2c   = 1.5d0)
      parameter (dc    = 0.7d0)
      parameter (dxc   = 2.3631d0)
      parameter (xi    = 0.12802585262625815d0)
      parameter (thr1  = 0.996d0)
      parameter (thr2  = 1.004d0)
      parameter (delta = 1.0d-3)

      b2 = dsqrt(5913.0d0/405000.0d0)
      b1 = (511.0d0/13500.0d0)/(2.0d0*b2)
      b4 = muAK*muAK/K1 - 1606.0d0/18225.0d0 - b1*b1

      Cx   = (-0.75d0)*(3.0d0/pi)**thrd
      P23  = (3.0d0*pi**2.0d0)**twthrd
      P13  = (3.0d0/(4.0d0*pi))**thrd
      P23t = (3.0d0*pi*pi/16.0d0)**twthrd

      do i=1,n2ft3d
        nup       = rho_in(i,1) + ETA
        agrup     = agr_in(i,1) + ETA
        tauup     = tau_in(i,1) + ETA
        ndn       = rho_in(i,2) + ETA
        agrdn     = agr_in(i,2) + ETA
        taudn     = tau_in(i,2) + ETA

*       ***** SCAN Exchange *****
*       ***** UP *****
        n   = 2.0d0*nup
        agr = 2.0d0*agrup
        tau = 2.0d0*tauup

        n_13  = n**thrd
        n_53  = n_13*n_13*n
        n_83  = n_53*n
        agr2  = agr*agr

        p       =  agr2/(4.0d0*P23*n_83)
        p_14    =  dsqrt(dsqrt(p))
        dp_dn   = -etthrd*p/n
        dp_dagr =  2.0d0*p/agr
c       dp_dtau =  0.0d0

        tauW  = 0.125d0*agr2/n
        tauU  = 0.3d0*P23*n_53
        z    = tauW/tau
        dz_dn   = -z/n
        dz_dagr =  2.0d0*z/agr
        dz_dtau = -z/tau

        z2 = z*z

        alpha       = 0.0d0
        dalpha_dn   = 0.0d0
        dalpha_dagr = 0.0d0
        dalpha_dtau = 0.0d0

        thr0 = agr2/(8.0d0*tau)
        dthr0_dagr = 2.0d0*thr0/agr
        dthr0_dtau = -thr0/tau
        delta2 = delta*delta
        diff = n - thr0
        diff2 = diff*diff
        diff3 = diff2*diff
        dod = diff/delta
        dod2 = dod*dod
        tauUn = tauU*n

c       if (n .le. thr0) then
c         alpha       = 0.0d0
c         dalpha_dn   = 0.0d0
c         dalpha_dagr = 0.0d0
c         dalpha_dtau = 0.0d0
c       else if ((n .gt. thr0) .and. (n .lt. (thr0 + delta))) then
c         falpha = tau*(-diff3/delta2 + 2.0d0*diff2/delta)
c         dfalpha_dn =  tau*(-3.0d0*dod2 + 4.0d0*dod)
c         dfalpha_dagr = tau*(3.0d0*dod2 - 4.0d0*dod)*dthr0_dagr
c         dfalpha_dtau = falpha/tau
c    &                 + tau*(3.0d0*dod2 - 4.0d0*dod)*dthr0_dtau
c         alpha = falpha/tauUn
c         dalpha_dn = -etthrd*falpha/(tauUn*n) + dfalpha_dn/tauUn
c         dalpha_dagr = dfalpha_dagr/tauUn
c         dalpha_dtau = dfalpha_dtau/tauUn
c       else if (n .ge. (thr0 + delta)) then
c         falpha = tau*n - agr2/8.0d0
c         dfalpha_dn = tau
c         dfalpha_dagr = -agr/4.0d0
c         dfalpha_dtau = n
c         alpha = falpha/tauUn
c         dalpha_dn = -etthrd*falpha/(tauUn*n) + dfalpha_dn/tauUn
c         dalpha_dagr = dfalpha_dagr/tauUn
c         dalpha_dtau = dfalpha_dtau/tauUn
c       end if

        alpha = (tau - tauW)/tauU

        if (alpha .lt. 0.0d0) then
          alpha       = 0.0d0
          dalpha_dn   = 0.0d0
          dalpha_dagr = 0.0d0
          dalpha_dtau = 0.0d0
        else
          falpha = tau*n - agr2/8.0d0
          dfalpha_dn = tau
          dfalpha_dagr = -agr/4.0d0
          dfalpha_dtau = n
          alpha = falpha/tauUn
          dalpha_dn = -etthrd*falpha/(tauUn*n) + dfalpha_dn/tauUn
          dalpha_dagr = dfalpha_dagr/tauUn
          dalpha_dtau = dfalpha_dtau/tauUn
        end if

        oma  = 1.0d0 - alpha
        oma2 = oma*oma

        exp1 = dexp(-b4*p/muAK)
        exp2 = dexp(-b3*oma2)

        x1 = muAK*p*(1.0d0 + (b4*p/muAK)*exp1)
        x2 = b1*p + b2*oma*exp2
        x = x1 + x2*x2

        denh1x = K1 + x
        numh1x = denh1x + K1*x
        h1x    = numh1x/denh1x

        dx1_dp    = muAK + b4*p*exp1*(2.0d0 - p*b4/muAK)
        dx2_dp    = b1
        dx_dp     = dx1_dp + 2.0d0*x2*dx2_dp
        dx_dalpha = 2.0d0*b2*exp2*x2*(2.0d0*b3*oma2 - 1.0d0)

        dh1x_dx     = (K1/denh1x)**2.0d0
        dh1x_dp     = dh1x_dx*dx_dp
        dh1x_dalpha = dh1x_dx*dx_dalpha

        gx     = 1.0d0
        dgx_dp =  0.0d0
        if (p_14 .lt. 0.002d0) then
           exp3   = 0.0d0
           gx     = 1.0d0
           dgx_dp = 0.0d0
        else
           exp3   = dexp(-a1/p_14)
           gx     = 1.0d0 - exp3
           dgx_dp = -0.25d0*a1*exp3/(p*p_14)
        endif

        fxa = 0.d0
        dfxa_dalpha = 0.0d0

        if (alpha .le. thr1) then
           exp4 = dexp(-c1x*alpha/oma)
           exp5 = 0.0d0
           fxa  = exp4
           dfxa_dalpha = -c1x*exp4/oma2
        else if ((alpha .gt. thr1) .and. (alpha .lt. thr2)) then
           exp4 = 0.0d0
           exp5 = 0.0d0
           fxa = 0.0d0
           dfxa_dalpha = 0.0d0
        else if (alpha .ge. thr2) then
           exp4 = 0.0d0
           exp5 = dexp(c2x/oma)
           fxa = -dx*exp5
           dfxa_dalpha = -dx*c2x*exp5/oma2
        end if

        Fx = (h1x + fxa*(h0x - h1x))*gx

        dFx_dp     = dgx_dp*(h1x + fxa*(h0x - h1x)) 
     &             + gx*dh1x_dp*(1.0d0 - fxa)
        dFx_dalpha = gx*(dh1x_dalpha + dfxa_dalpha*(h0x - h1x) 
     &             - fxa*dh1x_dalpha)

        dFx_dn   = dFx_dalpha*dalpha_dn   + dFx_dp*dp_dn
        dFx_dagr = dFx_dalpha*dalpha_dagr + dFx_dp*dp_dagr
        dFx_dtau = dFx_dalpha*dalpha_dtau

        ex0  = Cx*n_13
        nex0 = n*ex0

        eupx      = ex0*Fx
        fnupx     = nex0*dFx_dn     + frthrd*eupx
        fdnupx    = nex0*dFx_dagr
        fdtauupx  = nex0*dFx_dtau

*       ***** DOWN *****
        n   = 2.0d0*ndn
        agr = 2.0d0*agrdn
        tau = 2.0d0*taudn

        n_13  = n**thrd
        n_53  = n_13*n_13*n
        n_83  = n_53*n
        agr2  = agr*agr

        p       =  agr2/(4.0d0*P23*n_83)
        p_14    =  dsqrt(dsqrt(p))
        dp_dn   = -etthrd*p/n
        dp_dagr =  2.0d0*p/agr
c       dp_dtau =  0.0d0

        tauW  = 0.125d0*agr2/n
        tauU  = 0.3d0*P23*n_53
        z    = tauW/tau
        dz_dn   = -z/n
        dz_dagr =  2.0d0*z/agr
        dz_dtau = -z/tau

        z2 = z*z

        alpha       = 0.0d0
        dalpha_dn   = 0.0d0
        dalpha_dagr = 0.0d0
        dalpha_dtau = 0.0d0

        thr0 = agr2/(8.0d0*tau)
        dthr0_dagr = 2.0d0*thr0/agr
        dthr0_dtau = -thr0/tau
        delta2 = delta*delta
        diff = n - thr0
        diff2 = diff*diff
        diff3 = diff2*diff
        dod = diff/delta
        dod2 = dod*dod
        tauUn = tauU*n

c       if (n .le. thr0) then
c         alpha       = 0.0d0
c         dalpha_dn   = 0.0d0
c         dalpha_dagr = 0.0d0
c         dalpha_dtau = 0.0d0
c       else if ((n .gt. thr0) .and. (n .lt. (thr0 + delta))) then
c         falpha = tau*(-diff3/delta2 + 2.0d0*diff2/delta)
c         dfalpha_dn =  tau*(-3.0d0*dod2 + 4.0d0*dod)
c         dfalpha_dagr = tau*(3.0d0*dod2 - 4.0d0*dod)*dthr0_dagr
c         dfalpha_dtau = falpha/tau
c    &                 + tau*(3.0d0*dod2 - 4.0d0*dod)*dthr0_dtau
c         alpha = falpha/tauUn
c         dalpha_dn = -etthrd*falpha/(tauUn*n) + dfalpha_dn/tauUn
c         dalpha_dagr = dfalpha_dagr/tauUn
c         dalpha_dtau = dfalpha_dtau/tauUn
c       else if (n .ge. (thr0 + delta)) then
c         falpha = tau*n - agr2/8.0d0
c         dfalpha_dn = tau
c         dfalpha_dagr = -agr/4.0d0
c         dfalpha_dtau = n
c         alpha = falpha/tauUn
c         dalpha_dn = -etthrd*falpha/(tauUn*n) + dfalpha_dn/tauUn
c         dalpha_dagr = dfalpha_dagr/tauUn
c         dalpha_dtau = dfalpha_dtau/tauUn
c       end if

        alpha = (tau - tauW)/tauU

        if (alpha .lt. 0.0d0) then
          alpha       = 0.0d0
          dalpha_dn   = 0.0d0
          dalpha_dagr = 0.0d0
          dalpha_dtau = 0.0d0
        else
          falpha = tau*n - agr2/8.0d0
          dfalpha_dn = tau
          dfalpha_dagr = -agr/4.0d0
          dfalpha_dtau = n
          alpha = falpha/tauUn
          dalpha_dn = -etthrd*falpha/(tauUn*n) + dfalpha_dn/tauUn
          dalpha_dagr = dfalpha_dagr/tauUn
          dalpha_dtau = dfalpha_dtau/tauUn
        end if
        oma  = 1.0d0 - alpha
        oma2 = oma*oma

        exp1 = dexp(-b4*p/muAK)
        exp2 = dexp(-b3*oma2)

        x1 = muAK*p*(1.0d0 + (b4*p/muAK)*exp1)
        x2 = b1*p + b2*oma*exp2
        x = x1 + x2*x2

        denh1x = K1 + x
        numh1x = denh1x + K1*x
        h1x    = numh1x/denh1x

        dx1_dp    = muAK + b4*p*exp1*(2.0d0 - p*b4/muAK)
        dx2_dp    = b1
        dx_dp     = dx1_dp + 2.0d0*x2*dx2_dp
        dx_dalpha = 2.0d0*b2*exp2*x2*(2.0d0*b3*oma2 - 1.0d0)

        dh1x_dx     = (K1/denh1x)**2.0d0
        dh1x_dp     = dh1x_dx*dx_dp
        dh1x_dalpha = dh1x_dx*dx_dalpha

        gx     = 1.0d0
        dgx_dp =  0.0d0
        if (p_14 .lt. 0.002d0) then
           exp3   = 0.0d0
           gx     = 1.0d0
           dgx_dp = 0.0d0
        else
           exp3   = dexp(-a1/p_14)
           gx     = 1.0d0 - exp3
           dgx_dp = -0.25d0*a1*exp3/(p*p_14)
        endif

        fxa = 0.d0
        dfxa_dalpha = 0.0d0

        if (alpha .le. thr1) then
           exp4 = dexp(-c1x*alpha/oma)
           exp5 = 0.0d0
           fxa  = exp4
           dfxa_dalpha = -c1x*exp4/oma2
        else if ((alpha .gt. thr1) .and. (alpha .lt. thr2)) then
           exp4 = 0.0d0
           exp5 = 0.0d0
           fxa = 0.0d0
           dfxa_dalpha = 0.0d0
        else if (alpha .ge. thr2) then
           exp4 = 0.0d0
           exp5 = dexp(c2x/oma)
           fxa = -dx*exp5
           dfxa_dalpha = -dx*c2x*exp5/oma2
        end if

        Fx = (h1x + fxa*(h0x - h1x))*gx

        dFx_dp     = dgx_dp*(h1x + fxa*(h0x - h1x)) 
     &             + gx*dh1x_dp*(1.0d0 - fxa)
        dFx_dalpha = gx*(dh1x_dalpha + dfxa_dalpha*(h0x - h1x) 
     &             - fxa*dh1x_dalpha)

        dFx_dn   = dFx_dalpha*dalpha_dn   + dFx_dp*dp_dn
        dFx_dagr = dFx_dalpha*dalpha_dagr + dFx_dp*dp_dagr
        dFx_dtau = dFx_dalpha*dalpha_dtau

        ex0  = Cx*n_13
        nex0 = n*ex0

        ednx      = ex0*Fx
        fndnx     = nex0*dFx_dn     + frthrd*ednx
        fdndnx    = nex0*dFx_dagr
        fdtaudnx  = nex0*dFx_dtau

        n  = nup + ndn

        ex = (eupx*nup + ednx*ndn)/n

*       ***** SCAN Correlation *****
        agr       = agr_in(i,3) + ETA

        tau  = tauup + taudn
        n_13 = n**thrd
        n_53 = n*n_13*n_13
        n_83 = n*n_53
        n2   = n*n
        agr2 = agr*agr

        rs     =  P13/n_13
        drs_dn = -thrd*rs/n
        rs_12  =  dsqrt(rs)

        zeta = (nup - ndn)/n

        if (dabs(nup-ndn) .lt. tol) then
           zeta       = 0.0d0
           dzeta_dnup = 0.0d0
           dzeta_dndn = 0.0d0
        else 
           dzeta_dnup =  2.0d0*ndn/n2
           dzeta_dndn = -2.0d0*nup/n2 
        end if

        p       =  agr2/(4.0d0*P23*n_83)
        dp_dn   = -etthrd*p/n
        dp_dagr =  2.0d0*p/agr

        opz    = 1.0d0 + zeta
        omz    = 1.0d0 - zeta
        opz_23 = opz**twthrd
        omz_23 = omz**twthrd

        phi  = 0.5d0*(opz_23 + omz_23)
        phi2 = phi*phi
        phi3 = phi2*phi

        if (omz .lt. tol) then
           dphi_dzeta =  0.5d0*twthrd*(opz_23/opz)
        else if (opz .lt. tol) then
           dphi_dzeta = -0.5d0*twthrd*(omz_23/omz)
        else
           dphi_dzeta =  0.5d0*twthrd*(opz_23/opz - omz_23/omz)
        end if

        beta      = beta0*(1.0d0 + beta1*rs)/(1.0d0 + beta2*rs)
        dbeta_drs = beta0*(beta1 - beta2)/((1.0d0 + beta2*rs)**2.0d0)

        call gen_PW91_c_rz(tol,rs,zeta,ecLDA1,decLDA1_drs,decLDA1_dzeta)

        w1fac     =  ecLDA1/(gamma*phi3)
        expw1     =  dexp(-w1fac)
        w1        =  expw1 - 1.0d0
        dw1_drs   = -expw1*decLDA1_drs/(gamma*phi3)
        dw1_dzeta =  (3.0d0*w1fac*dphi_dzeta/phi 
     &            -  decLDA1_dzeta/(gamma*phi3))*expw1

        A        =  beta/(gamma*w1)
        dA_drs   =  dbeta_drs/(gamma*w1) - A*dw1_drs/w1
        dA_dzeta = -A*dw1_dzeta/w1

        t2        =  P23t*p/(phi2*rs)
        dt2_drs   = -t2/rs
        dt2_dzeta = -2.0d0*t2*dphi_dzeta/phi
        dt2_dp    =  P23t/(phi2*rs)

        At2     = A*t2
        dengAt2 = 1.0d0 + 4.0d0*At2
        gAt2    = 1.0d0/(dsqrt(dsqrt(dengAt2)))

        tmp1        = 1.0d0 + w1*(1.0d0 - gAt2)
        dtmp1_drs   = dw1_drs*(1.0d0 - gAt2) +
     &                gAt2*w1*(t2*dA_drs + A*dt2_drs)/dengAt2
        dtmp1_dzeta = dw1_dzeta*(1.0d0 - gAt2) 
     &              + gAt2*w1*(t2*dA_dzeta + A*dt2_dzeta)/dengAt2
        dtmp1_dp    = gAt2*w1*A*dt2_dp/dengAt2


        H1        = gamma*phi3*dlog(tmp1)
        dH1_drs   = gamma*phi3*dtmp1_drs/tmp1
        dH1_dzeta = 3.0d0*H1*dphi_dzeta/phi 
     &            + gamma*phi3*dtmp1_dzeta/tmp1
        dH1_dp    = gamma*phi3*dtmp1_dp/tmp1

        ec1        = ecLDA1 + H1
        dec1_drs   = decLDA1_drs + dH1_drs
        dec1_dzeta = decLDA1_dzeta + dH1_dzeta
        dec1_dp    = dH1_dp

        ecLDA0      = -b1c/(1.0d0 + b2c*rs_12 + b3c*rs)
        decLDA0_drs =  b1c*(b3c + 0.5d0*b2c/rs_12)/
     &                 ((1.0d0 + b2c*rs_12 + b3c*rs)**2.0d0)

        dxz        =  0.5d0*(opz**frthrd + omz**frthrd)
        ddxz_dzeta =  0.5d0*frthrd*(opz**thrd - omz**thrd)
        zeta12     =  zeta**12.0d0
        omz12      =  1.0d0 - zeta12
        zeta11     =  zeta**11.0d0
        gc         =  (1.0d0 - dxc*(dxz - 1.0d0))*omz12
        dgc_dzeta  = -dxc*ddxz_dzeta*omz12 
     &             -  12.0d0*zeta11*(1.0d0 - dxc*(dxz - 1.0d0))

        w0fac   =  ecLDA0/b1c
        expw0   =  dexp(-w0fac)
        w0      =  expw0 - 1.0d0
        dw0_drs = -decLDA0_drs*expw0/b1c

        ginf     =  1.0d0/(dsqrt(dsqrt(1.0d0 + 4.0d0*xi*p)))
        dginf_dp = -xi*ginf/(1.0d0 + 4.0d0*xi*p)

        tmp0      =  1.0d0 + w0*(1.0d0 - ginf)
        dtmp0_drs =  dw0_drs*(1.0d0 - ginf)
        dtmp0_dp  = -w0*dginf_dp

        H0      = b1c*dlog(tmp0)
        dH0_drs = b1c*dtmp0_drs/tmp0
        dH0_dp  = b1c*dtmp0_dp/tmp0

        ec0        = (ecLDA0 + H0)*gc
        dec0_drs   = gc*(decLDA0_drs + dH0_drs)
        dec0_dzeta = dgc_dzeta*(ecLDA0 + H0)
        dec0_dp    = gc*dH0_dp

        ds        = 0.5d0*(opz**fvthrd + omz**fvthrd)
        dds_dzeta = 0.5d0*fvthrd*(opz**twthrd - omz**twthrd)

        tauW  = 0.125d0*agr2/n
        tauU  = 0.3d0*P23*ds*n_53
        z    = tauW/tau
        dz_dn   = -z/n
        dz_dagr =  2.0d0*z/agr
        dz_dtau = -z/tau

        z2 = z*z

        alpha         = 0.0d0
        dalpha_dnup   = 0.0d0
        dalpha_dndn   = 0.0d0
        dalpha_dagr   = 0.0d0
        dalpha_dtauup = 0.0d0
        dalpha_dtaudn = 0.0d0

        thr0 = agr2/(8.0d0*tau)
        dthr0_dagr = 2.0d0*thr0/agr
        dthr0_dtau = -thr0/tau
        delta2 = delta*delta
        diff = n - thr0
        diff2 = diff*diff
        diff3 = diff2*diff
        dod = diff/delta
        dod2 = dod*dod
        tauUn = tauU*n

        if (n .le. thr0) then
          alpha       = 0.0d0
          dalpha_dn   = 0.0d0
          dalpha_dagr = 0.0d0
          dalpha_dtau = 0.0d0
        else if ((n .gt. thr0) .and. (n .lt. (thr0 + delta))) then
          falpha       = tau*(-diff3/delta2 + 2.0d0*diff2/delta)
          dfalpha_dn   = tau*(-3.0d0*dod2 + 4.0d0*dod)
          dfalpha_dagr = tau*(3.0d0*dod2 - 4.0d0*dod)*dthr0_dagr
          dfalpha_dtau = falpha/tau
     &                 + tau*(3.0d0*dod2 - 4.0d0*dod)*dthr0_dtau
          alpha = falpha/tauUn
          dalpha_dnup = -etthrd*falpha/(tauUn*n) + dfalpha_dn/tauUn
     &              - (alpha/ds)*dds_dzeta*dzeta_dnup
          dalpha_dndn = -etthrd*falpha/(tauUn*n) + dfalpha_dn/tauUn
     &              - (alpha/ds)*dds_dzeta*dzeta_dndn
          dalpha_dagr = dfalpha_dagr/tauUn
          dalpha_dtauup = dfalpha_dtau/tauUn
          dalpha_dtaudn = dalpha_dtauup
        else if (n .ge. (thr0 + delta)) then
          falpha = tau*n - agr2/8.0d0
          dfalpha_dn = tau
          dfalpha_dagr = -agr/4.0d0
          dfalpha_dtau = n
          alpha = falpha/tauUn
          dalpha_dnup = -etthrd*falpha/(tauUn*n) + dfalpha_dn/tauUn
     &                - (alpha/ds)*dds_dzeta*dzeta_dnup 
          dalpha_dndn = -etthrd*falpha/(tauUn*n) + dfalpha_dn/tauUn
     &                - (alpha/ds)*dds_dzeta*dzeta_dndn
          dalpha_dagr = dfalpha_dagr/tauUn
          dalpha_dtauup = dfalpha_dtau/tauUn
          dalpha_dtaudn = dalpha_dtauup 
        end if

        oma  = 1.0d0 - alpha
        oma2 = oma*oma

        fca = 0.d0
        dfca_dalpha = 0.0d0

        if (alpha .le. thr1) then
           exp6 = dexp(-c1c*alpha/oma)
           exp7 = 0.0d0
           fca  = exp6
           dfca_dalpha = -c1c*exp6/oma2
        else if ((alpha .gt. thr1) .and. (alpha .lt. thr2)) then
           exp6 = 0.0d0
           exp7 = 0.0d0
           fca  = 0.0d0
           dfca_dalpha = 0.0d0
        else if (alpha .ge. thr2) then
           exp6 = 0.0d0
           exp7 = dexp(c2c/oma)
           fca = -dc*exp7
           dfca_dalpha = -dc*c2c*exp7/oma2
        end if
 
        dec1_dnup = dec1_drs*drs_dn + dec1_dzeta*dzeta_dnup 
     &            + dec1_dp*dp_dn
        dec1_dndn = dec1_drs*drs_dn + dec1_dzeta*dzeta_dndn
     &            + dec1_dp*dp_dn
        dec1_dagr = dec1_dp*dp_dagr
        
        dec0_dnup = dec0_drs*drs_dn + dec0_dzeta*dzeta_dnup
     &            + dec0_dp*dp_dn
        dec0_dndn = dec0_drs*drs_dn + dec0_dzeta*dzeta_dndn
     &            + dec0_dp*dp_dn
        dec0_dagr = dec0_dp*dp_dagr

        dfca_dnup   = dfca_dalpha*dalpha_dnup
        dfca_dndn   = dfca_dalpha*dalpha_dndn
        dfca_dagr   = dfca_dalpha*dalpha_dagr
        dfca_dtauup = dfca_dalpha*dalpha_dtauup
        dfca_dtaudn = dfca_dalpha*dalpha_dtaudn
          
        ec       = ec1 + fca*(ec0 - ec1)
        fnupc    = n*(dec1_dnup + dfca_dnup*(ec0 - ec1) 
     &           + fca*(dec0_dnup - dec1_dnup)) + ec
        fndnc    = n*(dec1_dndn + dfca_dndn*(ec0 - ec1) 
     &           + fca*(dec0_dndn - dec1_dndn)) + ec
        fdnupc   = 0.0d0
        fdndnc   = 0.0d0
        fdnc     = n*(dec1_dagr + dfca_dagr*(ec0 - ec1) 
     &           + fca*(dec0_dagr - dec1_dagr)) 
        fdtauupc = n*dfca_dtauup*(ec0 - ec1)
        fdtaudnc = n*dfca_dtaudn*(ec0 - ec1)
       
 
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
