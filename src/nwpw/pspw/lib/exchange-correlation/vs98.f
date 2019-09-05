*     ************************************************
*     *                                              *
*     *                   nwpw_GVT4                  *
*     *                                              *
*     ************************************************
      subroutine nwpw_GVT4(ag,bg,cg,dg,eg,fg,alphag,alphagt,
     >                     xg,zg,gammag,G,dGdx,dGdz)
      implicit none
      real*8 ag,bg,cg,dg,eg,fg,alphag,alphagt,xg,zg,gammag
      real*8 G,dGdx,dGdz
      real*8 xg2,zg2,gammag2

      xg2     = xg*xg
      zg2     = zg*zg
      gammag2 = gammag*gammag

      G = (ag + bg*xg + cg*zg 
     &     + dg*xg2 + eg*xg*zg + fg*zg2)/gammag

      dGdx = ( -ag*alphag
     &        + bg*(1.0d0 - 2.0d0*alphag*xg)
     &        - 2.0d0*cg*alphag*zg
     &        + dg*(2.0d0*xg - 3.0d0*alphag*xg2)
     &        + eg*(zg - 3.0d0*alphag*zg*xg)
     &        - 3.0d0*fg*alphag*zg2)/gammag2
      
      dGdz = ( -ag*alphagt
     &        - 2.0d0*bg*alphagt*xg
     &        + cg*(1.0d0 - 2.0d0*alphagt*zg)
     &        - 3.0d0*dg*alphagt*xg2
     &        + eg*(xg - 3.0d0*alphagt*xg*zg)
     &        + fg*(2.0d0*zg - 3.0d0*alphagt*zg2))/gammag2

      return
      end
*     ************************************************
*     *                                              *
*     *                   nwpw_SS                    *
*     *                                              *
*     ************************************************
      subroutine nwpw_ss(as,bs,cs,ds,es,fs,alphas,n,agr,tau,e0c,f0nc,
     >                   xs,zs,gammas,xsr,zsr,
     >                   dxsdn,dxsdagr,dzsdn,dzsdtau,
     >                   ce,dfnc,dfdnc,dfdtauc)
      implicit none
      real*8 as,bs,cs,ds,es,fs,alphas,n,agr,tau,e0c,f0nc            !*** input ***
      real*8 xs,zs,gammas,xsr,zsr
      real*8 dxsdn,dxsdagr,dzsdn,dzsdtau
      real*8 ce,dfnc,dfdnc,dfdtauc                                  !*** output ***
********* local declarations******
      real*8 D,fD,thresD,dDdx,dDdz,dfDdD,dDdn,dDdagr,dDdtau
      real*8 dxdn,dxdagr,dzdn,dzdtau
      real*8 GS,dGdx,dGdz,dGdn,dGdagr,dGdtau
********* constants ******
      real*8 fvthrd,etthrd
      parameter (fvthrd = 5.0d0/3.0d0)
      parameter (etthrd = 8.0d0/3.0d0)
********* VS98 constants ******
      real*8 cf
      parameter (cf = 9.115599720d0)         !***cf = 3.0d0/5.0d0*(3.0d0*pi*pi)**(2.0d0/3.0d0))***
      
      D      = 1.0d0 - xs/(4.0d0*(zs + cf))      
      dDdx   = -1.0d0/(4.0d0*(zs + cf))
      dDdz   = xs/(4.0d0*(zs + cf)*(zs + cf))
      thresD = 0.05d0
      fD     = 0.0d0
      dDdn   = 0.0d0
      dDdagr = 0.0d0
      dDdtau = 0.0d0

      if (D .ge. thresD) then

        fD     = D
        dDdn   = dDdx*dxsdn + dDdz*dzsdn
        dDdagr = dDdx*dxsdagr
        dDdtau = dDdz*dzsdtau

      else if (D .lt. 0.0d0) then

        fD     = 0.0d0
        dDdn   = 0.0d0
        dDdagr = 0.0d0
        dDdtau = 0.0d0

      else if (D .ge. 0.0d0 .and. D .lt. thresD) then

        fD     = 2.0d0*D*D/thresD - D*D*D/(thresD*thresD)
        dfDdD  = (4.0d0*D/thresD - 3.0d0*D*D/(thresD*thresD))
        dDdn   = dfDdD*(dDdx*dxsdn + dDdz*dzsdn)
        dDdagr = dfDdD*(dDdx*dxsdagr)
        dDdtau = dfDdD*(dDdz*dzsdtau)

      end if

      call nwpw_GVT4(as,bs,cs,ds,es,fs,alphas,alphas,
     >               xsr,zsr,gammas,GS,dGdx,dGdz)


      dGdn   = dGdx*dxsdn + dGdz*dzsdn
      dGdagr = dGdx*dxsdagr
      dGdtau = dGdz*dzsdtau

      ce      = GS*fD*e0c
      dfnc    = (dGdn*fD*n + GS*dDdn*n + GS*fD)*e0c + n*GS*fD*f0nc
      dfdnc   = n*e0c*(dGdagr*fD + GS*dDdagr)
      dfdtauc = n*e0c*(dGdtau*fD + GS*dDdtau)
      
        
      return
      end


*     ************************************************
*     *                                              *
*     *              gen_VS98_restricted             *
*     *                                              *
*     ************************************************

*    This function returns the VS98 exchange-correlation
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
*     Exit  - xce(n2ft3d) : VS98 exchange correlation energy density
*             fn(n2ft3d)  : d(n*xce)/dn
*             fdn(n2ft3d) : d(n*xce)/d|grad n|
*             fdtau(n2ft3d) : d(n*xce)/dtau
*

      subroutine gen_VS98_restricted(n2ft3d,rho_in,agr_in,tau_in,
     >                               x_parameter,c_parameter,
     >                               xce,fn,fdn,fdtau)
      implicit none
      integer n2ft3d
      real*8 rho_in(*),agr_in(*),tau_in(*)     !*** input  ****
      real*8 x_parameter,c_parameter
      real*8 xce(*),fn(*),fdn(*),fdtau(*)      !*** output ****

********* local declarations******
      integer i
      real*8 n,agr,tau
      real*8 inv_n,n_13,n_43,n_53,n_83,agr2
      real*8 x,dxdn,dxdagr,z,dzdn,dzdtau
      real*8 gamma,xr,zr,GG,dGdx,dGdz
      real*8 ec0(n2ft3d),fnc0(n2ft3d),ec_LDA,fnc_LDA
      real*8 ex,fnx,fdnx,fdtaux
      real*8 ec,fnc,fdnc,fdtauc
********* constants ***************
      real*8 thrd,frthrd,fvthrd,etthrd
      parameter (thrd   = 1.0d0/3.0d0)
      parameter (frthrd = 4.0d0/3.0d0)
      parameter (fvthrd = 5.0d0/3.0d0)
      parameter (etthrd = 8.0d0/3.0d0)
********* density cutoff parameters******
      real*8 tol,ETA
      parameter (tol = 1.0d-18)
      parameter (ETA            =      1.0d-20)

********* VS98 constants ******
      real*8 cf,alpha,aa,bb,cc,dd,ee,ff
      parameter (cf = 9.115599720d0)         !***cf = 3.0d0/5.0d0*(6.0d0*pi*pi)**(2.0d0/3.0d0))***

      call gen_PW92_BW_restricted(n2ft3d,rho_in,ec0,fnc0)

      do i=1,n2ft3d
        n       = 0.50d0*rho_in(i) + ETA
        agr     = 0.50d0*agr_in(i)
        tau     = 0.50d0*tau_in(i)
        ec_LDA  = 0.50d0*ec0(i)
        fnc_LDA = 0.50d0*fnc0(i)

       !if (n .gt. tol)  then

        agr2  = agr*agr
        inv_n = 1.0d0/n
        n_13  = n**thrd
        n_43  = n_13*n
        n_53  = n_43*n_13
        n_83  = n_53*n

        x      = agr2/n_83
        dxdn   = -etthrd*x*inv_n
        dxdagr = 2.0d0*agr/n_83
        z      = tau/n_53 - cf
        dzdn   = -fvthrd*tau/n_83
        dzdtau = 1.0d0/n_53
        
********* VS98 Exchange *********
        aa    = -9.800683d-1
        bb    = -3.556788d-3
        cc    =  6.250326d-3
        dd    = -2.354518d-5
        ee    = -1.282732d-4
        ff    =  3.574822d-4
        alpha =  0.00186726d0

        gamma =  1.0d0 + alpha*(x + z)
        xr    = x/gamma
        zr    = z/gamma

        call nwpw_GVT4(aa,bb,cc,dd,ee,ff,alpha,alpha,
     >              xr,zr,gamma,GG,dGdx,dGdz)

        ex     = 2.0d0*n_13*GG
        fnx    = 2.0d0*(frthrd*n_13*GG + n_43*(dGdx*dxdn + dGdz*dzdn))
        fdnx   = 2.0d0*n_43*(dGdx*dxdagr)
        fdtaux = 2.0d0*n_43*(dGdz*dzdtau) 

********* VS98 Correlation *********
********* Same-spin *********
        aa    =  3.270912d-1
        bb    = -3.228915d-2
        cc    = -2.942406d-2
        dd    =  2.134222d-3
        ee    = -5.451559d-3
        ff    =  1.577575d-2
        alpha =  0.00515088d0

        gamma =  1.0d0 + alpha*(x + z)
        xr    = x/gamma
        zr    = z/gamma
      
        call nwpw_ss(aa,bb,cc,dd,ee,ff,alpha,n,agr,tau,ec_LDA,fnc_LDA,
     >               x,z,gamma,xr,zr,
     >               dxdn,dxdagr,dzdn,dzdtau,
     >               ec,fnc,fdnc,fdtauc)

        ec     = 2.0d0*ec
        fnc    = 2.0d0*fnc
        fdnc   = 2.0d0*fdnc
        fdtauc = 2.0d0*fdtauc

c      else
c         ex = 0.0d0
c         fnx = 0.0d0
c         fdnx = 0.0d0
c         fdtaux = 0.0d0
c         ec = 0.0d0
c         fnc = 0.0d0
c         fdnc = 0.0d0
c         fdtauc = 0.0d0
c      end if

        xce(i)   = x_parameter*ex     + c_parameter*ec
        fn(i)    = x_parameter*fnx    + c_parameter*fnc
        fdn(i)   = x_parameter*fdnx   + c_parameter*fdnc
        fdtau(i) = x_parameter*fdtaux + c_parameter*fdtauc


      end do
      return
      end

*     ************************************************
*     *                                              *
*     *              gen_VS98_unrestricted             *
*     *                                              *
*     ************************************************

*    This function returns the VS98 exchange-correlation
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
*     Exit  - xce(n2ft3d) : VS98 exchange correlation energy density
*             fn(n2ft3d,2)  : d(n*xce)/dnup, d(n*xce)/dndn
*             fdn(n2ft3d,3) : d(n*xce)/d|grad nup|, d(n*xce)/d|grad ndn|, d(n*xce)/d|grad n|
*             fdtau(n2ft3d,2) : d(n*xce)/dtauup, d(n*xce)/dtaudn
*
      subroutine gen_VS98_unrestricted(n2ft3d,rho_in,agr_in,tau_in,
     >                               x_parameter,c_parameter,
     >                               xce,fn,fdn,fdtau)
      implicit none

      integer n2ft3d
      real*8 rho_in(n2ft3d,2),agr_in(n2ft3d,3),tau_in(n2ft3d,2)      !*** input  ****
      real*8 x_parameter,c_parameter
      real*8 xce(n2ft3d),fn(n2ft3d,2),fdn(n2ft3d,3),fdtau(n2ft3d,2)  !*** output ****

********* local declarations******
      integer i
      real*8 rhouu_in(n2ft3d,2),euu0(n2ft3d),fuu0n(n2ft3d,2)
      real*8 rhodd_in(n2ft3d,2),edd0(n2ft3d),fdd0n(n2ft3d,2)
      real*8 eud0(n2ft3d),fud0n(n2ft3d,2)
      real*8 nup,agrup,tauup
      real*8 ndn,agrdn,taudn
      real*8 inv_nup,nup_13,nup_43,nup_53,nup_83,agrup2
      real*8 inv_ndn,ndn_13,ndn_43,ndn_53,ndn_83,agrdn2
      real*8 xup,dxup_dnup,dxup_dagrup,zup,dzup_dnup,dzup_dtauup
      real*8 xdn,dxdn_dndn,dxdn_dagrdn,zdn,dzdn_dndn,dzdn_dtaudn
      real*8 xud,zud
      real*8 gammaup,xrup,zrup,GGup,dGup_dxup,dGup_dzup
      real*8 gammadn,xrdn,zrdn,GGdn,dGdn_dxdn,dGdn_dzdn
      real*8 gammaud,xrud,zrud,GGud,dGud_dxud,dGud_dzud
      real*8 dGud_dnup,dGud_dagrup,dGud_dtauup
      real*8 dGud_dndn,dGud_dagrdn,dGud_dtaudn
      real*8 eupx,fnupx,fdnupx,fdtauupx
      real*8 ednx,fndnx,fdndnx,fdtaudnx
      real*8 euu0c,fuu0_nupc,euuc
      real*8 fuu_nupc,fuu_dnupc,fuu_dtauupc
      real*8 edd0c,fdd0_ndnc,eddc
      real*8 fdd_ndnc,fdd_dndnc,fdd_dtaudnc
      real*8 eud0c,fud0_nupc,fud0_ndnc,eudc
      real*8 fud_nupc,fud_dnupc,fud_dtauupc
      real*8 fud_ndnc,fud_dndnc,fud_dtaudnc
      real*8 eupc,fnupc,fdnupc,fdtauupc
      real*8 ednc,fndnc,fdndnc,fdtaudnc
      real*8 ex,ec,n
********* constants ***************
      real*8 thrd,frthrd,fvthrd,etthrd
      parameter (thrd   = 1.0d0/3.0d0)
      parameter (frthrd = 4.0d0/3.0d0)
      parameter (fvthrd = 5.0d0/3.0d0)
      parameter (etthrd = 8.0d0/3.0d0)
********* density cutoff parameters******
      real*8 tol,ETA
      parameter (tol = 1.0d-10)
      parameter (ETA            =      1.0d-20)

********* VS98 constants ******
      real*8 cf,alpha,aa,bb,cc,dd,ee,ff
      parameter (cf = 9.115599720d0)         !***cf = 3.0d0/5.0d0*(6.0d0*pi*pi)**(2.0d0/3.0d0))***

      do i=1,n2ft3d
         rhouu_in(i,1) = rho_in(i,1)
         rhouu_in(i,2) = 0.0d0
         rhodd_in(i,1) = 0.0d0 
         rhodd_in(i,2) = rho_in(i,2)
      end do
      call gen_PW92_BW_unrestricted(n2ft3d,rhouu_in,euu0,fuu0n) 
      call gen_PW92_BW_unrestricted(n2ft3d,rhodd_in,edd0,fdd0n) 
      call gen_PW92_BW_unrestricted(n2ft3d,rho_in,eud0,fud0n) 
     
      do i=1,n2ft3d

        nup        = rho_in(i,1) + ETA
        agrup      = agr_in(i,1)
        tauup      = tau_in(i,1)
        euu0c      = euu0(i)
        fuu0_nupc  = fuu0n(i,1)

        ndn        = rho_in(i,2) + ETA
        agrdn      = agr_in(i,2)
        taudn      = tau_in(i,2)
        edd0c      = edd0(i)
        fdd0_ndnc  = fdd0n(i,2)

        eud0c      = eud0(i)
        fud0_nupc  = fud0n(i,1)
        fud0_ndnc  = fud0n(i,2)

        n          = nup + ndn
     
        agrup2  = agrup*agrup
        inv_nup = 1.0d0/nup
        nup_13  = nup**thrd
        nup_43  = nup_13*nup
        nup_53  = nup_43*nup_13
        nup_83  = nup_53*nup

        xup         = agrup2/nup_83
        dxup_dnup   = -etthrd*xup*inv_nup
        dxup_dagrup = 2.0d0*agrup/nup_83
        zup         = tauup/nup_53 - cf
        dzup_dnup   = -fvthrd*tauup/nup_83
        dzup_dtauup = 1.0d0/nup_53
        
        agrdn2  = agrdn*agrdn
        inv_ndn = 1.0d0/ndn
        ndn_13  = ndn**thrd
        ndn_43  = ndn_13*ndn
        ndn_53  = ndn_43*ndn_13
        ndn_83  = ndn_53*ndn

        xdn         = agrdn2/ndn_83
        dxdn_dndn   = -etthrd*xdn*inv_ndn
        dxdn_dagrdn = 2.0d0*agrdn/ndn_83
        zdn         = taudn/ndn_53 - cf
        dzdn_dndn   = -fvthrd*taudn/ndn_83
        dzdn_dtaudn = 1.0d0/ndn_53
********* VS98 Exchange *********
        aa    = -9.800683d-1
        bb    = -3.556788d-3
        cc    =  6.250326d-3
        dd    = -2.354518d-5
        ee    = -1.282732d-4
        ff    =  3.574822d-4
        alpha =  0.00186726d0

********* Spin-up *********
        gammaup = 1.0d0 + alpha*(xup + zup)
        xrup    = xup/gammaup
        zrup    = zup/gammaup

        call nwpw_GVT4(aa,bb,cc,dd,ee,ff,alpha,alpha,
     >                 xrup,zrup,gammaup,GGup,dGup_dxup,dGup_dzup)

        eupx     = nup_13*GGup
        fnupx    = frthrd*nup_13*GGup 
     &           + nup_43*(dGup_dxup*dxup_dnup + dGup_dzup*dzup_dnup)
        fdnupx   = nup_43*(dGup_dxup*dxup_dagrup)
        fdtauupx = nup_43*(dGup_dzup*dzup_dtauup) 
********* Spin-down *********
        gammadn = 1.0d0 + alpha*(xdn + zdn)
        xrdn    = xdn/gammadn
        zrdn    = zdn/gammadn

        call nwpw_GVT4(aa,bb,cc,dd,ee,ff,alpha,alpha,
     >                 xrdn,zrdn,gammadn,GGdn,dGdn_dxdn,dGdn_dzdn)

        ednx     = nup_13*GGdn
        fndnx    = frthrd*ndn_13*GGdn 
     &           + ndn_43*(dGdn_dxdn*dxdn_dndn + dGdn_dzdn*dzdn_dndn)
        fdndnx   = ndn_43*(dGdn_dxdn*dxdn_dagrdn)
        fdtaudnx = ndn_43*(dGdn_dzdn*dzdn_dtaudn) 
 
        ex = eupx + ednx
********* VS98 Correlation *********
********* Same-spin *********
        aa    =  3.270912d-1
        bb    = -3.228915d-2
        cc    = -2.942406d-2
        dd    =  2.134222d-3
        ee    = -5.451559d-3
        ff    =  1.577575d-2
        alpha =  0.00515088d0

********* alpha-alpha *********
        gammaup = 1.0d0 + alpha*(xup + zup)
        xrup    = xup/gammaup
        zrup    = zup/gammaup

        call nwpw_ss(aa,bb,cc,dd,ee,ff,alpha,nup,agrup,tauup,
     >               euu0c,fuu0_nupc,xup,zup,gammaup,xrup,zrup,
     >               dxup_dnup,dxup_dagrup,dzup_dnup,dzup_dtauup,
     >               euuc,fuu_nupc,fuu_dnupc,fuu_dtauupc) 
********* beta-beta *********
        gammadn = 1.0d0 + alpha*(xdn + zdn)
        xrdn    = xdn/gammadn
        zrdn    = zdn/gammadn

        call nwpw_ss(aa,bb,cc,dd,ee,ff,alpha,ndn,agrdn,taudn,
     >               edd0c,fdd0_ndnc,xdn,zdn,gammadn,xrdn,zrdn,
     >               dxdn_dndn,dxdn_dagrdn,dzdn_dndn,dzdn_dtaudn,
     >               eddc,fdd_ndnc,fdd_dndnc,fdd_dtaudnc) 
*********Opposite-spin *********
        aa    =  7.035010d-1
        bb    =  7.694574d-3
        cc    =  5.152765d-2
        dd    =  3.394308d-5
        ee    = -1.269420d-3
        ff    =  1.296118d-3
        alpha =  0.00304966d0

        xud     = xup + xdn
        zud     = zup + zdn
        gammaud = 1.0d0 + alpha*(xud + zud)
        xrud    = xud/gammaud
        zrud    = zud/gammaud

        call nwpw_GVT4(aa,bb,cc,dd,ee,ff,alpha,alpha,
     >                 xrud,zrud,gammaud,GGud,dGud_dxud,dGud_dzud)

        eudc        = GGud*eud0c
        dGud_dnup   = dGud_dxud*dxup_dnup + dGud_dzud*dzup_dnup
        dGud_dagrup = dGud_dxud*dxup_dagrup 
        dGud_dtauup = dGud_dzud*dzup_dtauup
        dGud_dndn   = dGud_dxud*dxdn_dndn + dGud_dzud*dzdn_dndn
        dGud_dagrdn = dGud_dxud*dxdn_dagrdn 
        dGud_dtaudn = dGud_dzud*dzdn_dtaudn

        !fudc = n*GGud*eud0c
        fud_nupc    = GGud*eud0c 
     &                + nup*(dGud_dnup*eud0c + GGud*fud0_nupc)
        fud_dnupc   = nup*eud0c*dGud_dagrup
        fud_dtauupc = nup*eud0c*dGud_dtauup
        fud_ndnc    = GGud*eud0c 
     &                + ndn*(dGud_dndn*eud0c + GGud*fud0_ndnc)
        fud_dndnc   = ndn*eud0c*dGud_dagrdn
        fud_dtaudnc = ndn*eud0c*dGud_dtaudn

        ec       = euuc + eddc + eudc
        fnupc    = fuu_nupc + fud_nupc
        fdnupc   = fuu_dnupc + fud_dnupc
        fdtauupc = fuu_dtauupc + fud_dtauupc 
        fndnc    = fdd_ndnc + fud_ndnc
        fdndnc   = fdd_dndnc + fud_dndnc
        fdtaudnc = fdd_dtaudnc + fud_dtaudnc 
        
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

