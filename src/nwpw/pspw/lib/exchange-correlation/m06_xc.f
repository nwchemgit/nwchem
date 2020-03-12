*     ************************************************
*     *                                              *
*     *                 nwpw_m06_x                   *
*     *                                              *
*     ************************************************
      subroutine nwpw_m06_x(pi,thrd,fvthrd,etthrd,
     >                      a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,
     >                      C1,C2,Cx,P23,
     >                      n,agr,tau,
     >                      xe,dfdnx,dfdagrx,dfdtaux)
      implicit none
*     ***** input *****
      real*8 pi,thrd,fvthrd,etthrd
      real*8 a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11
      real*8 C1,C2,Cx,P23
      real*8 n,agr,tau
*     ***** output *****
      real*8 xe,dfdnx,dfdagrx,dfdtaux
*     ***** local declarations *****
      real*8 n_13,n_23,n_53,n_83,agr2
      real*8 tauU,t,dt_dn,dt_dtau
      real*8 w,dw_dt,w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11
      real*8 fw,dfw_dw,dfw_dn,dfw_dtau
      real*8 x,dx_dn,dx_dagr
      real*8 enum,eden,eden2,etmp,detmp_dx
      real*8 Fpbe,dFpbe_dn,dFpbe_dagr

      n_13  = n**thrd
      n_23  = n_13*n_13
      n_53  = n_23*n
      n_83  = n_53*n
      agr2  = agr*agr

      tauU    = P23*n_53 
      t       = tauU/tau
      dt_dn   = fvthrd*t/n
      dt_dtau = -t/tau

      w     = (t - 1.0d0)/(t + 1.0d0)
      dw_dt = 2.0d0/((1.0d0 + t)**2.0d0)

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

      fw       = a0    + a1*w1 + a2*w2 + a3*w3 + a4*w4   + a5*w5 
     &         + a6*w6 + a7*w7 + a8*w8 + a9*w9 + a10*w10 + a11*w11 
      dfw_dw   = a1            + 2.0d0*a2*w1     + 3.0d0*a3*w2     
     &         + 4.0d0*a4*w3   + 5.0d0*a5*w4     + 6.0d0*a6*w5 
     &         + 7.0d0*a7*w6   + 8.0d0*a8*w7     + 9.0d0*a9*w8 
     &         + 10.0d0*a10*w9 + 11.0d0*a11*w10
      dfw_dn   = dfw_dw*dw_dt*dt_dn
      dfw_dtau = dfw_dw*dw_dt*dt_dtau

      x       =  agr2/n_83
      dx_dn   = -etthrd*x/n
      dx_dagr =  2.0d0*x/agr

      enum     =  C1*x
      eden     =  1.0d0 + C2*x
      eden2    =  eden*eden
      etmp     =  -enum/eden
      detmp_dx =  -(C1*eden - C2*enum)/eden2

      Fpbe       = n_13*(Cx + etmp)
      dFpbe_dn   = n_13*detmp_dx*dx_dn   + thrd*Fpbe/n
      dFpbe_dagr = n_13*detmp_dx*dx_dagr

      xe      = Fpbe*fw
      dfdnx   = n*(dFpbe_dn*fw + Fpbe*dfw_dn) + xe
      dfdagrx = n*dFpbe_dagr*fw
      dfdtaux = n*Fpbe*dfw_dtau

      return
      end
*     ************************************************
*     *                                              *
*     *                nwpw_m06_c_ss                 *
*     *                                              *
*     ************************************************
      subroutine nwpw_m06_c_ss(cgss,css0,css1,css2,css3,css4,
     >                         n,xs,dxs_dn,dxs_dagr,zs,dzs_dn,dzs_dtau,
     >                         ess0c,dess0c_dn,
     >                         cess,dfss_dn,dfss_dagr,dfss_dtau)
      implicit none
*     ***** input *****
      real*8 cgss,css0,css1,css2,css3,css4
      real*8 n,xs,dxs_dn,dxs_dagr,zs,dzs_dn,dzs_dtau
*     ***** output *****
      real*8 ess0c,dess0c_dn
      real*8 cess,dfss_dn,dfss_dagr,dfss_dtau
*     ***** local declarations *****
      real*8 D,dD_dx,dD_dz,dD_dn,dD_dagr,dD_dtau
      real*8 U,dU_dx,U2,U3,U4,gss,dgss_dU
      real*8 rs,drs_dn
      real*8 pwc,dpwc_drs,dummy,ness0c
*     ***** constants *****
      real*8 pi,thrd
      parameter (pi     = 3.14159265358979311599d0)
      parameter (thrd   = 1.0d0/3.0d0)
*     ***** density cutoff parameters *****
      real*8 tol
      parameter (tol = 1.0d-18)
*     ***** M06 constants *****
      real*8 cf
      parameter (cf = 9.115599720d0)

      rs = (0.75d0/(pi*n))**thrd
      drs_dn = -thrd*rs/n

      call gen_PW91_c_rz(tol,rs,1.0d0,pwc,dpwc_drs,dummy)

      ess0c = pwc
      ness0c = n*ess0c
      dess0c_dn = dpwc_drs*drs_dn

      U = cgss*xs/(1.0d0 + cgss*xs)
      dU_dx = cgss/((1.0d0 + cgss*xs)**2.0d0)

      U2 = U*U
      U3 = U2*U
      U4 = U3*U

      gss = css0 + css1*U + css2*U2 + css3*U3 + css4*U4
      dgss_dU = css1 + 2.0d0*css2*U + 3.0d0*css3*U2 + 4.0d0*css4*U3

      D = 1.0d0 - 0.25d0*xs/(zs+cf) 
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

      cess = ess0c*gss*D
      dfss_dn = cess + n*dess0c_dn*gss*D 
     &        + ness0c*dgss_dU*dU_dx*dxs_dn*D + ness0c*gss*dD_dn
      dfss_dagr = ness0c*(dgss_dU*dU_dx*dxs_dagr*D + gss*dD_dagr)
      dfss_dtau = ness0c*gss*dD_dtau

      return
      end

*     ************************************************
*     *                                              *
*     *              nwpw_m06_c_restricted           *
*     *                                              *
*     ************************************************
      subroutine nwpw_m06_c_restricted(cgss,css0,css1,css2,css3,css4,
     >                             cgopp,copp0,copp1,copp2,copp3,copp4,
     >                             n,nup,xup,dxup_dnup,dxup_dagrup,
     >                             zup,dzup_dnup,dzup_dtauup,rs,drs_dn,
     >                             ce,fnc,fdnc,fdtauc)
      implicit none
*     ***** input *****
      real*8 cgss,css0,css1,css2,css3,css4
      real*8 cgopp,copp0,copp1,copp2,copp3,copp4
      real*8 n,nup,xup,dxup_dnup,dxup_dagrup
      real*8 zup,dzup_dnup,dzup_dtauup
      real*8 rs,drs_dn
*     ***** output *****
      real*8 ce,fnc,fdnc,fdtauc
*     ***** local declarations *****
      real*8 euu0c,deuu0c_dn
      real*8 euuc,dfuuc_dn,dfuuc_dagr,dfuuc_dtau
      real*8 eud0c,deud0c_drs,dummy,eud1c,fud1c,dfud1c_dn
      real*8 x,U,dU_dx,U2,U3,U4,gopp,dgopp_dU
      real*8 eudc,dfudc_dn,dfudc_dagr,dfudc_dtau
*     ***** density cutoff parameters *****
      real*8 tol
      parameter (tol = 1.0d-18)
      
*     ***** same spin (alpha-alpha or beta-beta) *****
      call nwpw_m06_c_ss(cgss,css0,css1,css2,css3,css4,
     >                   nup,xup,dxup_dnup,dxup_dagrup,
     >                   zup,dzup_dnup,dzup_dtauup,
     >                   euu0c,deuu0c_dn,
     >                   euuc,dfuuc_dn,dfuuc_dagr,dfuuc_dtau)

*     ***** opposite spin (alpha-beta) *****
      call gen_PW91_c_rz(tol,rs,0.0d0,eud0c,deud0c_drs,dummy)
      eud1c = eud0c - euu0c
      fud1c = n*eud1c
      dfud1c_dn  = eud0c + n*deud0c_drs*drs_dn
     &           - euu0c - nup*deuu0c_dn

      x   = 2.0d0*xup
      U = cgopp*x/(1.0d0 + cgopp*x)
      dU_dx = cgopp/((1.0d0 + cgopp*x)**2.0d0)

      U2 = U*U
      U3 = U2*U
      U4 = U3*U

      gopp = copp0 + copp1*U + copp2*U2 + copp3*U3 + copp4*U4
      dgopp_dU = copp1          + 2.0d0*copp2*U 
     &         + 3.0d0*copp3*U2 + 4.0d0*copp4*U3

      eudc = eud1c*gopp
      dfudc_dn = dfud1c_dn*gopp + fud1c*dgopp_dU*dU_dx*dxup_dnup
      dfudc_dagr = fud1c*dgopp_dU*dU_dx*dxup_dagrup
      dfudc_dtau = 0.0d0

      ce     = eudc        + euuc
      fnc    = dfudc_dn    + dfuuc_dn
      fdnc   = dfudc_dagr  + dfuuc_dagr
      fdtauc =               dfuuc_dtau
      
      return
      end

*     ************************************************
*     *                                              *
*     *            nwpw_m06_c_unrestricted           *
*     *                                              *
*     ************************************************
      subroutine nwpw_m06_c_unrestricted(cgss,css0,css1,css2,css3,css4,
     >                             cgopp,copp0,copp1,copp2,copp3,copp4,
     >                             n,nup,agrup,tauup,ndn,agrdn,taudn,
     >              xup,dxup_dnup,dxup_dagrup,zup,dzup_dnup,dzup_dtauup,
     >              xdn,dxdn_dndn,dxdn_dagrdn,zdn,dzdn_dndn,dzdn_dtaudn,
     >                   rs,drs_dn,zeta,dzeta_dnup,dzeta_dndn,
     >                   ce,fnupc,fndnc,fdnupc,fdndnc,fdtauupc,fdtaudnc)
      implicit none
*     ***** input *****
      real*8 cgss,css0,css1,css2,css3,css4
      real*8 cgopp,copp0,copp1,copp2,copp3,copp4
      real*8 n,nup,agrup,tauup,ndn,agrdn,taudn
      real*8 xup,dxup_dnup,dxup_dagrup,zup,dzup_dnup,dzup_dtauup
      real*8 xdn,dxdn_dndn,dxdn_dagrdn,zdn,dzdn_dndn,dzdn_dtaudn
      real*8 rs,drs_dn,zeta,dzeta_dnup,dzeta_dndn
*     ***** output *****
      real*8 ce,fnupc,fndnc,fdnupc,fdndnc,fdtauupc,fdtaudnc
*     ***** local declarations *****
      real*8 euu0c,deuu0c_dnup
      real*8 euuc,dfuuc_dnup,dfuuc_dagrup,dfuuc_dtauup
      real*8 edd0c,dedd0c_dndn
      real*8 eddc,dfddc_dndn,dfddc_dagrdn,dfddc_dtaudn
      real*8 eud0c,deud0c_drs,deud0c_dzeta
      real*8 x,U,dU_dx,U2,U3,U4,gopp,dgopp_dU
      real*8 eudc,eud1c,fud1c,dfud1c_dnup,dfud1c_dndn
      real*8 dfudc_dnup,dfudc_dagrup,dfudc_dtauup
      real*8 dfudc_dndn,dfudc_dagrdn,dfudc_dtaudn
*     ***** density cutoff parameters *****
      real*8 tol
      parameter (tol = 1.0d-18)

*     ***** same spin *****
*     ***** alpha-alpha *****
      call nwpw_m06_c_ss(cgss,css0,css1,css2,css3,css4,
     >                   nup,xup,dxup_dnup,dxup_dagrup,
     >                   zup,dzup_dnup,dzup_dtauup,
     >                   euu0c,deuu0c_dnup,
     >                   euuc,dfuuc_dnup,dfuuc_dagrup,dfuuc_dtauup)
*     ***** beta-beta *****
      call nwpw_m06_c_ss(cgss,css0,css1,css2,css3,css4,
     >                   ndn,xdn,dxdn_dndn,dxdn_dagrdn,
     >                   zdn,dzdn_dndn,dzdn_dtaudn,
     >                   edd0c,dedd0c_dndn,
     >                   eddc,dfddc_dndn,dfddc_dagrdn,dfddc_dtaudn)

*     ***** opposite spin (alpha-beta) *****
      call gen_PW91_c_rz(tol,rs,zeta,eud0c,deud0c_drs,deud0c_dzeta)
      eud1c = eud0c - (euu0c*nup + edd0c*ndn)/n
      fud1c = n*eud1c

      x   = xup + xdn
      U = cgopp*x/(1.0d0 + cgopp*x)
      dU_dx = cgopp/((1.0d0 + cgopp*x)**2.0d0)

      U2 = U*U
      U3 = U2*U
      U4 = U3*U

      gopp = copp0 + copp1*U + copp2*U2 + copp3*U3 + copp4*U4
      dgopp_dU = copp1          + 2.0d0*copp2*U 
     &         + 3.0d0*copp3*U2 + 4.0d0*copp4*U3

      eudc = eud1c*gopp
      dfud1c_dnup = eud0c + n*(deud0c_drs*drs_dn
     &            + deud0c_dzeta*dzeta_dnup) - euu0c 
     &            - nup*deuu0c_dnup 
      dfud1c_dndn = eud0c + n*(deud0c_drs*drs_dn
     &            + deud0c_dzeta*dzeta_dndn) - edd0c 
     &            - ndn*dedd0c_dndn

      dfudc_dnup = dfud1c_dnup*gopp + fud1c*dgopp_dU*dU_dx*dxup_dnup
      dfudc_dndn = dfud1c_dndn*gopp + fud1c*dgopp_dU*dU_dx*dxdn_dndn
      dfudc_dagrup = fud1c*dgopp_dU*dU_dx*dxup_dagrup
      dfudc_dagrdn = fud1c*dgopp_dU*dU_dx*dxdn_dagrdn
      dfudc_dtauup = 0.0d0
      dfudc_dtaudn = 0.0d0

      ce       = eudc          + (euuc*nup  + eddc*ndn)/n
      fnupc    = dfudc_dnup    + dfuuc_dnup
      fndnc    = dfudc_dndn    + dfddc_dndn
      fdnupc   = dfudc_dagrup  + dfuuc_dagrup
      fdndnc   = dfudc_dagrdn  + dfddc_dagrdn
      fdtauupc =                 dfuuc_dtauup
      fdtaudnc =                 dfddc_dtaudn

      return
      end
