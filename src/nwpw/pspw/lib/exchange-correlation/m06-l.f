*     ************************************************
*     *                                              *
*     *              gen_M06L_restricted              *
*     *                                              *
*     ************************************************

*    This function returns the M06L exchange-correlation
*  energy density, xce, and its derivatives with respect
*  to n, |grad n|, tau.

*
*     Entry - n2ft3d   : number of grid points
*             rho_in(*) :  density (nup+ndn)
*             agr_in(*): |grad rho_in|
*             tau_in(*): tau
*             x_parameter: scale parameter for exchange
*             c_parameter: scale parameter for correlation
*
*     Exit  - xce(n2ft3d) : M06L exchange correlation energy density
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
      real*8 n,agr,tau
      real*8 Cx,P23
      real*8 nup,agrup,tauup
      real*8 inv_nup,nup_13,nup_43,nup_53,nup_83,agrup2
      real*8 xup,dxup_dnup,dxup_dagrup,zup,dzup_dnup,dzup_dtauup
      real*8 ex1,fnx1,fdnx1,fdtaux1
      real*8 ex2,fnx2,fdnx2,fdtaux2
      real*8 ex,fnx,fdnx,fdtaux
      real*8 ec1,fnc1,fdnc1,fdtauc1
      real*8 rs,drs_dn
      real*8 ec2,fnc2,fdnc2,fdtauc2
      real*8 ec,fnc,fdnc,fdtauc
*     ***** constants *****
      real*8 pi,thrd,twthrd,frthrd,fvthrd,etthrd
      parameter (pi     = 3.14159265358979311599d0)
      parameter (thrd   = 1.0d0/3.0d0)
      parameter (twthrd = 2.0d0/3.0d0)
      parameter (frthrd = 4.0d0/3.0d0)
      parameter (fvthrd = 5.0d0/3.0d0)
      parameter (etthrd = 8.0d0/3.0d0)
*     ***** density cutoff parameters *****
      real*8 tol,ETA
      parameter (tol = 1.0d-18)
      parameter (ETA            =      1.0d-20)
*     ***** VS98 constants *****
      real*8 aax,bbx,ccx,ddx,eex,ffx
      real*8 aass,bbss,ccss,ddss,eess,ffss
      real*8 aaopp,bbopp,ccopp,ddopp,eeopp,ffopp
      parameter (aax   = 6.012244d-1*(-0.9305257363491))
      parameter (bbx   = 4.748822d-3*(-0.9305257363491))
      parameter (ccx   = -8.635108d-3*(-0.9305257363491))
      parameter (ddx   = -9.308062d-6*(-0.9305257363491))
      parameter (eex   = 4.482811d-5*(-0.9305257363491))
      parameter (ffx   = 0.000000d0)
      parameter (aass  = 4.650534d-1)
      parameter (bbss  = 1.617589d-1)
      parameter (ccss  = 1.833657d-1)
      parameter (ddss  = 4.692100d-4)
      parameter (eess  = -4.990573d-3)
      parameter (ffss  = 0.000000d0)
      parameter (aaopp = 3.957626d-1)
      parameter (bbopp = -5.614546d-1)
      parameter (ccopp = 1.403963d-2)
      parameter (ddopp = 9.831442d-4)
      parameter (eeopp = -3.577176d-3)
      parameter (ffopp = 0.000000d0)
*     ***** M06 constants *****
      real*8 cf,ma0,ma1,ma2,ma3,ma4,ma5,ma6,ma7,ma8,ma9,ma10,ma11
      real*8 C1,C2
      real*8 css,ccss0,ccss1,ccss2,ccss3,ccss4
      real*8 copp,ccopp0,ccopp1,ccopp2,ccopp3,ccopp4
      parameter (cf     = 9.115599720d0)         
c     cf = 3.0d0/5.0d0*(6.0d0*pi*pi)**(2.0d0/3.0d0))
      parameter (ma0    =  3.987756d-1)
      parameter (ma1    = 2.548219d-1)
      parameter (ma2    = 3.923994d-1)
      parameter (ma3    = -2.103655d0)
      parameter (ma4    = -6.302147d0)
      parameter (ma5    = 1.097615d1)
      parameter (ma6    = 3.097273d1)
      parameter (ma7    = -2.318489d1)
      parameter (ma8    = -5.673480d1)
      parameter (ma9    = 2.160364d1)
      parameter (ma10   = 3.421814d1)
      parameter (ma11   = -9.049762d0)
      parameter (C1     =  3.36116d-3)
      parameter (C2     =  4.49267d-3)
      parameter (css    = 0.06d0)
      parameter (ccss0  = 5.349466d-1)
      parameter (ccss1  = 5.396620d-1)
      parameter (ccss2  = -3.161217d1)
      parameter (ccss3  =  5.149592d1)
      parameter (ccss4  =  -2.919613d1)
      parameter (copp   = 0.0031d0)
      parameter (ccopp0 = 6.042374d-1)
      parameter (ccopp1 = 1.776783d2)
      parameter (ccopp2 = -2.513252d2)
      parameter (ccopp3 = 7.635173d1)
      parameter (ccopp4 = -1.255699d1)
      

      Cx = -1.50d0*(0.75d0/pi)**thrd
      P23 = 0.60d0*(6.0d0*pi*pi)**twthrd

      do i=1,n2ft3d
        n       =       rho_in(i) + ETA
        agr     =       agr_in(i) + ETA
        tau     = 2.0d0*tau_in(i) + ETA

        nup   = 0.50d0*n
        agrup = 0.50d0*agr
        tauup = 0.50d0*tau

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

        rs     =  (0.75d0/(pi*n))**thrd
        drs_dn = -thrd*rs/n

*       ***** VS98 Exchange *****
        call nwpw_vs98_x(aax,bbx,ccx,ddx,eex,ffx,
     >                   frthrd,nup_13,nup_43,
     >                   xup,zup,
     >                   dxup_dnup,dzup_dnup,dxup_dagrup,dzup_dtauup,
     >                   ex1,fnx1,fdnx1,fdtaux1)

*       ***** M06 Exchange *****
        call nwpw_m06_x(pi,thrd,fvthrd,etthrd,
     >                  ma0,ma1,ma2,ma3,ma4,ma5,
     >                  ma6,ma7,ma8,ma9,ma10,ma11,
     >                  C1,C2,Cx,P23,
     >                  nup,agrup,tauup,
     >                  ex2,fnx2,fdnx2,fdtaux2)

        ex     = ex1     + ex2
        fnx    = fnx1    + fnx2
        fdnx   = fdnx1   + fdnx2
        fdtaux = fdtaux1 + fdtaux2

*       ***** VS98 Correlation *****
        call nwpw_vs98_c_restricted(aass,bbss,ccss,ddss,eess,ffss,
     >                              aaopp,bbopp,ccopp,ddopp,eeopp,ffopp,
     >                              n,nup,agrup,tauup,
     >                              xup,dxup_dnup,dxup_dagrup,
     >                              zup,dzup_dnup,dzup_dtauup,
     >                              rs,drs_dn,
     >                              ec1,fnc1,fdnc1,fdtauc1)

*       ***** M06 Correlation *****
        call nwpw_m06_c_restricted(css,ccss0,ccss1,ccss2,ccss3,ccss4,
     >                          copp,ccopp0,ccopp1,ccopp2,ccopp3,ccopp4,
     >                             n,nup,
     >                             xup,dxup_dnup,dxup_dagrup,
     >                             zup,dzup_dnup,dzup_dtauup,
     >                             rs,drs_dn,
     >                             ec2,fnc2,fdnc2,fdtauc2)

        ec     = ec1     + ec2
        fnc    = fnc1    + fnc2
        fdnc   = fdnc1   + fdnc2
        fdtauc = fdtauc1 + fdtauc2

        xce(i)   = x_parameter*ex     + c_parameter*ec
        fn(i)    = x_parameter*fnx    + c_parameter*fnc
        fdn(i)   = x_parameter*fdnx   + c_parameter*fdnc
        fdtau(i) = x_parameter*fdtaux + c_parameter*fdtauc

      end do

      return
      end

*     ************************************************
*     *                                              *
*     *             gen_M06L_unrestricted            *
*     *                                              *
*     ************************************************

*    This function returns the M06L exchange-correlation
*  energy density, xce, and its derivatives with respect
*  to nup, ndn, |grad nup|, |grad ndn|, tauup, taudn.

*
*     Entry - n2ft3d   : number of grid points
*             rho_in(*,2) :  density (nup and ndn)
*             agr_in(*,3): |grad rho_in| (nup, ndn and n)
*             tau_in(*,2): tau (nup and ndn)
*             x_parameter: scale parameter for exchange
*             c_parameter: scale parameter for correlation
*
*     Exit  - xce(n2ft3d) : M06L exchange correlation energy density
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
      real*8 n,agr,tau
      real*8 nup,agrup,tauup
      real*8 ndn,agrdn,taudn
      real*8 Cx,P23
      real*8 inv_nup,nup_13,nup_43,nup_53,nup_83,agrup2
      real*8 inv_ndn,ndn_13,ndn_43,ndn_53,ndn_83,agrdn2
      real*8 xup,dxup_dnup,dxup_dagrup,zup,dzup_dnup,dzup_dtauup
      real*8 xdn,dxdn_dndn,dxdn_dagrdn,zdn,dzdn_dndn,dzdn_dtaudn
      real*8 ex1,eupx1,fnupx1,fdnupx1,fdtauupx1
      real*8 ednx1,fndnx1,fdndnx1,fdtaudnx1
      real*8 ex2,eupx2,fnupx2,fdnupx2,fdtauupx2
      real*8 ednx2,fndnx2,fdndnx2,fdtaudnx2
      real*8 ex,eupx,fnupx,fdnupx,fdtauupx,ednx,fndnx,fdndnx,fdtaudnx
      real*8 ec1,fnupc1,fndnc1,fdnupc1,fdndnc1,fdtauupc1,fdtaudnc1
      real*8 rs,drs_dn,zeta,dzeta_dnup,dzeta_dndn
      real*8 ec2,fnupc2,fndnc2,fdnupc2,fdndnc2,fdtauupc2,fdtaudnc2
      real*8 ec,fnupc,fndnc,fdnupc,fdndnc,fdtauupc,fdtaudnc
*     ***** constants *****
      real*8 pi,thrd,twthrd,frthrd,fvthrd,etthrd
      parameter (pi     = 3.14159265358979311599d0)
      parameter (thrd   = 1.0d0/3.0d0)
      parameter (twthrd = 2.0d0/3.0d0)
      parameter (frthrd = 4.0d0/3.0d0)
      parameter (fvthrd = 5.0d0/3.0d0)
      parameter (etthrd = 8.0d0/3.0d0)
*     ***** density cutoff parameters *****
      real*8 tol,ETA
      parameter (tol = 1.0d-18)
      parameter (ETA            =      1.0d-20)
*     ***** VS98 constants *****
      real*8 aax,bbx,ccx,ddx,eex,ffx
      real*8 aass,bbss,ccss,ddss,eess,ffss
      real*8 aaopp,bbopp,ccopp,ddopp,eeopp,ffopp
      parameter (aax   = 6.012244d-1*(-0.9305257363491))
      parameter (bbx   = 4.748822d-3*(-0.9305257363491))
      parameter (ccx   = -8.635108d-3*(-0.9305257363491))
      parameter (ddx   = -9.308062d-6*(-0.9305257363491))
      parameter (eex   = 4.482811d-5*(-0.9305257363491))
      parameter (ffx   = 0.000000d0)
      parameter (aass  = 4.650534d-1)
      parameter (bbss  = 1.617589d-1)
      parameter (ccss  = 1.833657d-1)
      parameter (ddss  = 4.692100d-4)
      parameter (eess  = -4.990573d-3)
      parameter (ffss  = 0.000000d0)
      parameter (aaopp = 3.957626d-1)
      parameter (bbopp = -5.614546d-1)
      parameter (ccopp = 1.403963d-2)
      parameter (ddopp = 9.831442d-4)
      parameter (eeopp = -3.577176d-3)
      parameter (ffopp = 0.000000d0)
*     ***** M06 constants *****
      real*8 cf,ma0,ma1,ma2,ma3,ma4,ma5,ma6,ma7,ma8,ma9,ma10,ma11
      real*8 C1,C2
      real*8 css,ccss0,ccss1,ccss2,ccss3,ccss4
      real*8 copp,ccopp0,ccopp1,ccopp2,ccopp3,ccopp4
      parameter (cf     = 9.115599720d0)         
c     cf = 3.0d0/5.0d0*(6.0d0*pi*pi)**(2.0d0/3.0d0))
      parameter (ma0    =  3.987756d-1)
      parameter (ma1    = 2.548219d-1)
      parameter (ma2    = 3.923994d-1)
      parameter (ma3    = -2.103655d0)
      parameter (ma4    = -6.302147d0)
      parameter (ma5    = 1.097615d1)
      parameter (ma6    = 3.097273d1)
      parameter (ma7    = -2.318489d1)
      parameter (ma8    = -5.673480d1)
      parameter (ma9    = 2.160364d1)
      parameter (ma10   = 3.421814d1)
      parameter (ma11   = -9.049762d0)
      parameter (C1     =  3.36116d-3)
      parameter (C2     =  4.49267d-3)
      parameter (css    = 0.06d0)
      parameter (ccss0  = 5.349466d-1)
      parameter (ccss1  = 5.396620d-1)
      parameter (ccss2  = -3.161217d1)
      parameter (ccss3  =  5.149592d1)
      parameter (ccss4  =  -2.919613d1)
      parameter (copp   = 0.0031d0)
      parameter (ccopp0 = 6.042374d-1)
      parameter (ccopp1 = 1.776783d2)
      parameter (ccopp2 = -2.513252d2)
      parameter (ccopp3 = 7.635173d1)
      parameter (ccopp4 = -1.255699d1)

      Cx = -1.50d0*(0.75d0/pi)**thrd
      P23 = 0.60d0*(6.0d0*pi*pi)**twthrd

      do i=1,n2ft3d
        nup        = rho_in(i,1) + ETA
        ndn        = rho_in(i,2) + ETA
        agrup      = agr_in(i,1) + ETA
        agrdn      = agr_in(i,2) + ETA
        agr        = agr_in(i,3) + ETA
        tauup      = tau_in(i,1) + ETA
        taudn      = tau_in(i,2) + ETA

        n  = nup + ndn

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

        rs         =  (0.75d0/(pi*n))**thrd
        drs_dn     = -thrd*rs/n
        zeta       =  (nup - ndn)/n
        dzeta_dnup =  ( 1.0d0 - zeta)/n
        dzeta_dndn =  (-1.0d0 - zeta)/n

*       ***** UP *****
*       ***** VS98 Exchange *****
        call nwpw_vs98_x(aax,bbx,ccx,ddx,eex,ffx,
     >                   frthrd,nup_13,nup_43,xup,zup,
     >                   dxup_dnup,dzup_dnup,dxup_dagrup,dzup_dtauup,
     >                   eupx1,fnupx1,fdnupx1,fdtauupx1)

*       ***** M06 Exchange *****
        call nwpw_m06_x(pi,thrd,fvthrd,etthrd,
     >                  ma0,ma1,ma2,ma3,ma4,ma5,
     >                  ma6,ma7,ma8,ma9,ma10,ma11,
     >                  C1,C2,Cx,P23,
     >                  nup,agrup,tauup,
     >                  eupx2,fnupx2,fdnupx2,fdtauupx2)

*       ***** DOWN *****
*       ***** VS98 Exchange *****
        call nwpw_vs98_x(aax,bbx,ccx,ddx,eex,ffx,
     >                   frthrd,ndn_13,ndn_43,xdn,zdn,
     >                   dxdn_dndn,dzdn_dndn,dxdn_dagrdn,dzdn_dtaudn,
     >                   ednx1,fndnx1,fdndnx1,fdtaudnx1)

*       ***** M06 Exchange *****
        call nwpw_m06_x(pi,thrd,fvthrd,etthrd,
     >                  ma0,ma1,ma2,ma3,ma4,ma5,
     >                  ma6,ma7,ma8,ma9,ma10,ma11,
     >                  C1,C2,Cx,P23,
     >                  ndn,agrdn,taudn,
     >                  ednx2,fndnx2,fdndnx2,fdtaudnx2)

        ex1 = (eupx1*nup + ednx1*ndn)/n
        ex2 = (eupx2*nup + ednx2*ndn)/n

        ex       = ex1       + ex2
        fnupx    = fnupx1    + fnupx2
        fndnx    = fndnx1    + fndnx2
        fdnupx   = fdnupx1   + fdnupx2
        fdndnx   = fdndnx1   + fdndnx2
        fdtauupx = fdtauupx1 + fdtauupx2
        fdtaudnx = fdtaudnx1 + fdtaudnx2
        
*       ***** VS98 Correlation *****
        call nwpw_vs98_c_unrestricted(aass,bbss,ccss,ddss,eess,ffss,
     >              aaopp,bbopp,ccopp,ddopp,eeopp,ffopp,
     >              n,nup,agrup,tauup,ndn,agrdn,taudn,
     >              xup,dxup_dnup,dxup_dagrup,zup,dzup_dnup,dzup_dtauup,
     >              xdn,dxdn_dndn,dxdn_dagrdn,zdn,dzdn_dndn,dzdn_dtaudn,
     >              rs,drs_dn,zeta,dzeta_dnup,dzeta_dndn,
     >            ec1,fnupc1,fndnc1,fdnupc1,fdndnc1,fdtauupc1,fdtaudnc1)

*       ***** M06 Correlation *****
        call nwpw_m06_c_unrestricted(css,ccss0,ccss1,ccss2,ccss3,ccss4,
     >              copp,ccopp0,ccopp1,ccopp2,ccopp3,ccopp4,
     >              n,nup,agrup,tauup,ndn,agrdn,taudn,
     >              xup,dxup_dnup,dxup_dagrup,zup,dzup_dnup,dzup_dtauup,
     >              xdn,dxdn_dndn,dxdn_dagrdn,zdn,dzdn_dndn,dzdn_dtaudn,
     >              rs,drs_dn,zeta,dzeta_dnup,dzeta_dndn,
     >            ec2,fnupc2,fndnc2,fdnupc2,fdndnc2,fdtauupc2,fdtaudnc2)

        ec       = ec1       + ec2
        fnupc    = fnupc1    + fnupc2
        fndnc    = fndnc1    + fndnc2
        fdnupc   = fdnupc1   + fdnupc2
        fdndnc   = fdndnc1   + fdndnc2
        fdtauupc = fdtauupc1 + fdtauupc2
        fdtaudnc = fdtaudnc1 + fdtaudnc2

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
