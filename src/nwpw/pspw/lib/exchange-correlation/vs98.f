*     ************************************************
*     *                                              *
*     *                   nwpw_GVT4                  *
*     *                                              *
*     ************************************************
      subroutine nwpw_GVT4(ag,bg,cg,dg,eg,fg,alphag,alphagt,
     >                     xg,zg,gammag,G,dG_dx,dG_dz)
      implicit none
*     ***** input *****
      real*8 ag,bg,cg,dg,eg,fg,alphag,alphagt,xg,zg,gammag
*     ***** output *****
      real*8 G,dG_dx,dG_dz
*     ***** local declarations *****
      real*8 xg2,zg2,gammag2

      xg2     = xg*xg
      zg2     = zg*zg
      gammag2 = gammag*gammag

      G = (ag + bg*xg + cg*zg 
     &  + dg*xg2 + eg*xg*zg + fg*zg2)/gammag

      dG_dx = (-ag*alphag
     &      + bg*(1.0d0 - 2.0d0*alphag*xg)
     &      - 2.0d0*cg*alphag*zg
     &      + dg*(2.0d0*xg - 3.0d0*alphag*xg2)
     &      + eg*(zg - 3.0d0*alphag*zg*xg)
     &      - 3.0d0*fg*alphag*zg2)/gammag2
      
      dG_dz = (-ag*alphagt
     &      - 2.0d0*bg*alphagt*xg
     &      + cg*(1.0d0 - 2.0d0*alphagt*zg)
     &      - 3.0d0*dg*alphagt*xg2
     &      + eg*(xg - 3.0d0*alphagt*xg*zg)
     &      + fg*(2.0d0*zg - 3.0d0*alphagt*zg2))/gammag2

      return
      end

*     ************************************************
*     *                                              *
*     *                nwpw_vs98c_ss                 *
*     *                                              *
*     ************************************************
      subroutine nwpw_vs98c_ss(as,bs,cs,ds,es,fs,alphas,n,agr,tau,
     >                   xs,zs,gammas,xsr,zsr,
     >                   dxs_dn,dxs_dagr,dzs_dn,dzs_dtau,
     >                   e0c,de0c_dn,
     >                   ce,dfnc,dfdnc,dfdtauc)
      implicit none
*     ***** input *****
      real*8 as,bs,cs,ds,es,fs,alphas,n,agr,tau
      real*8 xs,zs,gammas,xsr,zsr
      real*8 dxs_dn,dxs_dagr,dzs_dn,dzs_dtau
*     ***** output *****
      real*8 e0c,de0c_dn
      real*8 ce,dfnc,dfdnc,dfdtauc            
*     ***** local declarations *****
      real*8 rs,drs_dn,ne0c,pwc,dpwc_drs,dummy
      real*8 D,dD_dx,dD_dz,dD_dn,dD_dagr,dD_dtau
      real*8 GS,dG_dx,dG_dz,dG_dn,dG_dagr,dG_dtau
*     ***** constants *****
      real*8 pi,thrd,fvthrd,etthrd
      parameter (pi     =  3.14159265358979311599d0)
      parameter (thrd   = 1.0d0/3.0d0)
      parameter (fvthrd = 5.0d0/3.0d0)
      parameter (etthrd = 8.0d0/3.0d0)
*     ***** density cutoff parameters *****
      real*8 tol
      parameter (tol = 1.0d-18)
*     ***** VS98 constants *****
      real*8 cf
      parameter (cf = 9.115599720d0)         
c     cf = 3.0d0/5.0d0*(6.0d0*pi*pi)**(2.0d0/3.0d0))
     
      rs     =  (0.75d0/(pi*n))**thrd
      drs_dn = -rs/(3.0d0*n)

      call gen_PW91_c_rz(tol,rs,1.0d0,pwc,dpwc_drs,dummy)

      D = 1.0d0 - xs/(4.0d0*(zs + cf))      

      if (D .lt. 0.0d0) then
        D       = 0.0d0
        dD_dn   = 0.0d0
        dD_dagr = 0.0d0
        dD_dtau = 0.0d0
      else 
        dD_dx   = -1.0d0/(4.0d0*(zs + cf))
        dD_dz   =  xs/(4.0d0*(zs + cf)*(zs + cf))
        dD_dn   =  dD_dx*dxs_dn + dD_dz*dzs_dn
        dD_dagr =  dD_dx*dxs_dagr
        dD_dtau =  dD_dz*dzs_dtau
      end if

      call nwpw_GVT4(as,bs,cs,ds,es,fs,alphas,alphas,
     >               xsr,zsr,gammas,GS,dG_dx,dG_dz)

      dG_dn   = dG_dx*dxs_dn    + dG_dz*dzs_dn
      dG_dagr = dG_dx*dxs_dagr
      dG_dtau = dG_dz*dzs_dtau

      e0c     = pwc
      de0c_dn = dpwc_drs*drs_dn
      ne0c    = n*e0c
      ce      = GS*D*e0c
      dfnc    = ne0c*(dG_dn*D   + GS*dD_dn)    
     &        + n*GS*D*de0c_dn  + ce
      dfdnc   = ne0c*(dG_dagr*D + GS*dD_dagr)
      dfdtauc = ne0c*(dG_dtau*D + GS*dD_dtau)
      
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
*     Entry - n2ft3d   : number of grid points
*             rho_in(*) :  density (nup+ndn)
*             agr_in(*): |grad rho_in|
*             tau_in(*): tau
*             x_parameter: scale parameter for exchange
*             c_parameter: scale parameter for correlation
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
*     ***** input *****
      integer n2ft3d
      real*8 rho_in(*),agr_in(*),tau_in(*)     
      real*8 x_parameter,c_parameter
*     ***** output *****
      real*8 xce(*),fn(*),fdn(*),fdtau(*)      
*     ***** local declarations *****
      integer i
      real*8 n,agr,tau
      real*8 inv_n,n_13,n_43,n_53,n_83,agr2
      real*8 x,dx_dn,dx_dagr,z,dz_dn,dz_dtau
      real*8 gamma,xr,zr,GG,dG_dx,dG_dz
      real*8 ex,fnx,fdnx,fdtaux
      real*8 nup,agrup,tauup
      real*8 inv_nup,nup_13,nup_43,nup_53,nup_83,agrup2
      real*8 xup,dxup_dnup,dxup_dagrup,zup,dzup_dnup,dzup_dtauup
      real*8 xrup,zrup,gammaup
      real*8 euu0c,deuu0c_dn
      real*8 euuc,dfuuc_dn,dfuuc_dagr,dfuuc_dtau
      real*8 rs,drs_dn
      real*8 dG_dn,dG_dagr,dG_dtau
      real*8 eud0c,deud0c_drs,dummy
      real*8 eudc,eud1c,fud1c,dfud1c_dn
      real*8 dfudc_dn,dfudc_dagr,dfudc_dtau
      real*8 ec,fnc,fdnc,fdtauc
*     ***** constants *****
      real*8 pi,thrd,frthrd,fvthrd,etthrd
      parameter (pi     = 3.14159265358979311599d0)
      parameter (thrd   = 1.0d0/3.0d0)
      parameter (frthrd = 4.0d0/3.0d0)
      parameter (fvthrd = 5.0d0/3.0d0)
      parameter (etthrd = 8.0d0/3.0d0)
*     ***** density cutoff parameters *****
      real*8 tol,ETA
      parameter (tol = 1.0d-18)
      parameter (ETA            =      1.0d-20)
*     ***** VS98 constants *****
      real*8 cf,aa,bb,cc,dd,ee,ff,alpha
      parameter (cf = 9.115599720d0)         
c     cf = 3.0d0/5.0d0*(6.0d0*pi*pi)**(2.0d0/3.0d0))

      do i=1,n2ft3d
        n       = 0.50d0*rho_in(i) + ETA
        agr     = 0.50d0*agr_in(i) + ETA
        tau     =        tau_in(i) + ETA
c       tau     = 2.0d0*tau

c       if (n .gt. tol)  then
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
     >              xr,zr,gamma,GG,dG_dx,dG_dz)

        ex     = n_13*GG
        fnx    = n_43*(dG_dx*dx_dn + dG_dz*dz_dn) + frthrd*ex
        fdnx   = n_43*(dG_dx*dx_dagr)
        fdtaux = n_43*(dG_dz*dz_dtau) 

*       ***** VS98 Correlation *****
*       ***** same-spin *****
        aa    =  3.270912d-1
        bb    = -3.228915d-2
        cc    = -2.942406d-2
        dd    =  2.134222d-3
        ee    = -5.451559d-3
        ff    =  1.577575d-2
        alpha =  0.00515088d0
*       ***** alpha-alpha (beta-beta) *****
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

        gammaup = 1.0d0 + alpha*(xup + zup)
        xrup    = xup/gammaup
        zrup    = zup/gammaup

        call nwpw_vs98c_ss(aa,bb,cc,dd,ee,ff,alpha,nup,agrup,tauup,
     >               xup,zup,gammaup,xrup,zrup,
     >               dxup_dnup,dxup_dagrup,dzup_dnup,dzup_dtauup,
     >               euu0c,deuu0c_dn,
     >               euuc,dfuuc_dn,dfuuc_dagr,dfuuc_dtau)

*       ***** opposite-spin *****
        aa    =  7.035010d-1
        bb    =  7.694574d-3
        cc    =  5.152765d-2
        dd    =  3.394308d-5
        ee    = -1.269420d-3
        ff    =  1.296118d-3
        alpha =  0.00304966d0
*       ***** alpha-beta *****
        rs     =  (0.75d0/(pi*n))**thrd
        drs_dn = -rs/(3.0d0*n)

        call gen_PW91_c_rz(tol,rs,0.0d0,eud0c,deud0c_drs,dummy)
        eud1c = eud0c - euu0c
        fud1c = n*eud1c

        x   = 2.0d0*xup
        z   = 2.0d0*zup

        gamma = 1.0d0 + alpha*(x + z)
        xr    = x/gamma
        zr    = z/gamma

        call nwpw_GVT4(aa,bb,cc,dd,ee,ff,alpha,alpha,
     >                 xr,zr,gamma,GG,dG_dx,dG_dz)

        dG_dn   = dG_dx*dxup_dnup   + dG_dz*dzup_dnup
        dG_dagr = dG_dx*dxup_dagrup
        dG_dtau = dG_dz*dzup_dtauup

        eudc   = GG*eud1c

        dfud1c_dn  = eud0c + n*deud0c_drs*drs_dn  
     &             - euu0c - nup*deuu0c_dn
        dfudc_dn   = fud1c*dG_dn    + GG*dfud1c_dn
        dfudc_dagr = fud1c*dG_dagr
        dfudc_dtau = fud1c*dG_dtau

        ec     = eudc        + euuc 
        fnc    = dfudc_dn    + dfuuc_dn
        fdnc   = dfudc_dagr  + dfuuc_dagr
        fdtauc = dfudc_dtau  + dfuuc_dtau

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
*     Entry - n2ft3d   : number of grid points
*             rho_in(*,2) :  density (nup and ndn)
*             agr_in(*,3): |grad rho_in| (nup, ndn and n)
*             tau_in(*,2): tau (nup and ndn)
*             x_parameter: scale parameter for exchange
*             c_parameter: scale parameter for correlation
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
*     ***** input *****
      integer n2ft3d
      real*8 rho_in(n2ft3d,2),agr_in(n2ft3d,3),tau_in(n2ft3d,2)      
      real*8 x_parameter,c_parameter
*     ***** output *****
      real*8 xce(n2ft3d),fn(n2ft3d,2),fdn(n2ft3d,3),fdtau(n2ft3d,2)  
*     ***** local declarations *****
      integer i
      real*8 nup,agrup,tauup
      real*8 ndn,agrdn,taudn
      real*8 n,agr,tau
      real*8 inv_n,n_13,n_43,n_53,n_83,agr2
      real*8 x,dx_dn,dx_dagr,z,dz_dn,dz_dtau
      real*8 gamma,xr,zr,GG,dG_dx,dG_dz
      real*8 inv_nup,nup_13,nup_43,nup_53,nup_83,agrup2
      real*8 inv_ndn,ndn_13,ndn_43,ndn_53,ndn_83,agrdn2
      real*8 xup,dxup_dnup,dxup_dagrup,zup,dzup_dnup,dzup_dtauup
      real*8 xdn,dxdn_dndn,dxdn_dagrdn,zdn,dzdn_dndn,dzdn_dtaudn
      real*8 xud,zud
      real*8 gammaup,xrup,zrup
      real*8 gammadn,xrdn,zrdn
      real*8 gammaud,xrud,zrud,GGud,dGud_dxud,dGud_dzud
      real*8 dGud_dnup,dGud_dagrup,dGud_dtauup
      real*8 dGud_dndn,dGud_dagrdn,dGud_dtaudn
      real*8 eupx,fnupx,fdnupx,fdtauupx
      real*8 ednx,fndnx,fdndnx,fdtaudnx
      real*8 euu0c,deuu0c_dnup
      real*8 euuc,dfuuc_dnup,dfuuc_dagrup,dfuuc_dtauup
      real*8 edd0c,dedd0c_dndn
      real*8 eddc,dfddc_dndn,dfddc_dagrdn,dfddc_dtaudn
      real*8 rs,drs_dn,zeta,dzeta_dnup,dzeta_dndn
      real*8 eud0c,deud0c_drs,deud0c_dzeta
      real*8 eudc,eud1c,fud1c,dfud1c_dnup,dfud1c_dndn
      real*8 dfudc_dnup,dfudc_dagrup,dfudc_dtauup
      real*8 dfudc_dndn,dfudc_dagrdn,dfudc_dtaudn
      real*8 fnupc,fdnupc,fdtauupc
      real*8 fndnc,fdndnc,fdtaudnc
      real*8 ex,ec
********* constants ***************
      real*8 pi,thrd,frthrd,fvthrd,etthrd
      parameter (pi     =  3.14159265358979311599d0)
      parameter (thrd   = 1.0d0/3.0d0)
      parameter (frthrd = 4.0d0/3.0d0)
      parameter (fvthrd = 5.0d0/3.0d0)
      parameter (etthrd = 8.0d0/3.0d0)
********* density cutoff parameters******
      real*8 tol,ETA
      parameter (tol = 1.0d-18)
      parameter (ETA            =      1.0d-20)

********* VS98 constants ******
      real*8 cf,aa,bb,cc,dd,ee,ff,alpha
      parameter (cf = 9.115599720d0)        
c     cf = 3.0d0/5.0d0*(6.0d0*pi*pi)**(2.0d0/3.0d0))

      do i=1,n2ft3d
        nup        = 0.50d0*rho_in(i,1) + ETA
        agrup      = 0.50d0*agr_in(i,1) + ETA
        tauup      = 0.50d0*tau_in(i,1) + ETA
        ndn        = 0.50d0*rho_in(i,2) + ETA
        agrdn      = 0.50d0*agr_in(i,2) + ETA
        taudn      = 0.50d0*tau_in(i,2) + ETA

*       ***** VS98 Exchange *****
        aa    = -9.800683d-1
        bb    = -3.556788d-3
        cc    =  6.250326d-3
        dd    = -2.354518d-5
        ee    = -1.282732d-4
        ff    =  3.574822d-4
        alpha =  0.00186726d0
*       ***** UP *****
        n     = 2.0d0*nup
        agr   = 2.0d0*agrup
        tau   = 2.0d0*tauup

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
        
        gamma = 1.0d0 + alpha*(x + z)
        xr    = x/gamma
        zr    = z/gamma

        call nwpw_GVT4(aa,bb,cc,dd,ee,ff,alpha,alpha,
     >              xr,zr,gamma,GG,dG_dx,dG_dz)

        eupx     = n_13*GG
        fnupx    = frthrd*eupx + n_43*(dG_dx*dx_dn + dG_dz*dz_dn)
        fdnupx   = n_43*(dG_dx*dx_dagr)
        fdtauupx = n_43*(dG_dz*dz_dtau) 

*       ***** DOWN *****
        n     = 2.0d0*ndn
        agr   = 2.0d0*agrdn
        tau   = 2.0d0*taudn

        agr2  = agr*agr
        inv_n = 1.0d0/n
        n_13  = n**thrd
        n_43  = n_13*n
        n_53  = n_43*n_13
        n_83  = n_53*n

        x       =  agr2/n_83
        dx_dn   = -etthrd*x*inv_n
        dx_dagr =  2.0d0*agr/n_83

        z      =  tau/n_53 - cf
        dz_dn   = -fvthrd*tau/n_83
        dz_dtau =  1.0d0/n_53
        
        gamma = 1.0d0 + alpha*(x + z)
        xr    = x/gamma
        zr    = z/gamma

        call nwpw_GVT4(aa,bb,cc,dd,ee,ff,alpha,alpha,
     >              xr,zr,gamma,GG,dG_dx,dG_dz)

        ednx     = n_13*GG
        fndnx    = n_43*(dG_dx*dx_dn    + dG_dz*dz_dn) + frthrd*ednx
        fdndnx   = n_43*(dG_dx*dx_dagr)
        fdtaudnx = n_43*(dG_dz*dz_dtau) 

        n = nup + ndn

        ex = (eupx*nup + ednx*ndn)/n
*       ***** VS98 Correlation *****
        agr        = 0.50d0*agr_in(i,3) + ETA

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
*       ***** same-spin *****
        aa    =  3.270912d-1
        bb    = -3.228915d-2
        cc    = -2.942406d-2
        dd    =  2.134222d-3
        ee    = -5.451559d-3
        ff    =  1.577575d-2
        alpha =  0.00515088d0
*       ***** alpha-alpha *****
        gammaup = 1.0d0 + alpha*(xup + zup)
        xrup    = xup/gammaup
        zrup    = zup/gammaup

        call nwpw_vs98c_ss(aa,bb,cc,dd,ee,ff,alpha,nup,agrup,tauup,
     >               xup,zup,gammaup,xrup,zrup,
     >               dxup_dnup,dxup_dagrup,dzup_dnup,dzup_dtauup,
     >               euu0c,deuu0c_dnup,
     >               euuc,dfuuc_dnup,dfuuc_dagrup,dfuuc_dtauup)

*       ***** beta-beta *****
        gammadn = 1.0d0 + alpha*(xdn + zdn)
        xrdn    = xdn/gammadn
        zrdn    = zdn/gammadn

        call nwpw_vs98c_ss(aa,bb,cc,dd,ee,ff,alpha,ndn,agrdn,taudn,
     >               xdn,zdn,gammadn,xrdn,zrdn,
     >               dxdn_dndn,dxdn_dagrdn,dzdn_dndn,dzdn_dtaudn,
     >               edd0c,dedd0c_dndn,
     >               eddc,dfddc_dndn,dfddc_dagrdn,dfddc_dtaudn)

*       ***** opposite-spin *****
        aa    =  7.035010d-1
        bb    =  7.694574d-3
        cc    =  5.152765d-2
        dd    =  3.394308d-5
        ee    = -1.269420d-3
        ff    =  1.296118d-3
        alpha =  0.00304966d0
*       ***** alpha-beta *****
        rs     =  (0.75d0/(pi*n))**thrd
        drs_dn = -rs/(3.0d0*n)

        zeta       = (nup - ndn)/n
        dzeta_dnup = ( 1.0d0 - zeta)/n
        dzeta_dndn = (-1.0d0 - zeta)/n

        call gen_PW91_c_rz(tol,rs,zeta,eud0c,deud0c_drs,deud0c_dzeta)
        eud1c = eud0c - (euu0c*nup + edd0c*ndn)/n
        fud1c = n*eud1c

        xud     = xup + xdn
        zud     = zup + zdn

        gammaud = 1.0d0 + alpha*(xud + zud)
        xrud    = xud/gammaud
        zrud    = zud/gammaud

        call nwpw_GVT4(aa,bb,cc,dd,ee,ff,alpha,alpha,
     >                 xrud,zrud,gammaud,GGud,dGud_dxud,dGud_dzud)

        dGud_dnup   = dGud_dxud*dxup_dnup   + dGud_dzud*dzup_dnup
        dGud_dndn   = dGud_dxud*dxdn_dndn   + dGud_dzud*dzdn_dndn
        dGud_dagrup = dGud_dxud*dxup_dagrup
        dGud_dagrdn = dGud_dxud*dxdn_dagrdn
        dGud_dtauup = dGud_dzud*dzup_dtauup
        dGud_dtaudn = dGud_dzud*dzdn_dtaudn

        eudc   = GGud*eud1c

        dfud1c_dnup = eud0c + n*(deud0c_drs*drs_dn 
     &              + deud0c_dzeta*dzeta_dnup) - euu0c - nup*deuu0c_dnup
        dfud1c_dndn = eud0c + n*(deud0c_drs*drs_dn 
     &              + deud0c_dzeta*dzeta_dndn) - edd0c - ndn*dedd0c_dndn

        dfudc_dnup   = fud1c*dGud_dnup    + GGud*dfud1c_dnup
        dfudc_dndn   = fud1c*dGud_dndn    + GGud*dfud1c_dndn
        dfudc_dagrup = fud1c*dGud_dagrup
        dfudc_dagrdn = fud1c*dGud_dagrdn
        dfudc_dtauup = fud1c*dGud_dtauup
        dfudc_dtaudn = fud1c*dGud_dtaudn

        ec       = eudc          + (euuc*nup  + eddc*ndn)/n 
        fnupc    = dfudc_dnup    + dfuuc_dnup
        fndnc    = dfudc_dndn    + dfddc_dndn
        fdnupc   = dfudc_dagrup  + dfuuc_dagrup
        fdndnc   = dfudc_dagrdn  + dfddc_dagrdn
        fdtauupc = dfudc_dtauup  + dfuuc_dtauup
        fdtaudnc = dfudc_dtaudn  + dfddc_dtaudn

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

