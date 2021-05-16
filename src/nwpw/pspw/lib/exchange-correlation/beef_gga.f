
*    ************************************
*    *					*
*    *	    gen_BEEF_BW_unrestricted	*
*    *					*
*    ************************************
*
*    This function returns the BEEF exchange-correlation
*  energy density, xce, and its derivatives with respect
*  to nup, ndn, |grad nup|, |grad ndn|, and |grad n|.
*
*   Entry - n2ft3d     : number of grid points
*           dn_in(*,2) : spin densites nup and ndn
*           agr_in(*,3): |grad nup|, |grad ndn|, and |grad n|
*           x_parameter: scale parameter for exchange
*           c_parameter: scale parameter for correlation
*
*   Exit - xce(*)  : PBE96 energy density
*        - fn(*,2) : d(n*xce)/dnup, d(n*xce)/dndn
*        - fdn(*,3): d(n*xce)/d|grad nup|, d(n*xce)/d|grad ndn|
*                    d(n*xce)/d|grad n|

      subroutine gen_BEEF_BW_unrestricted(n2ft3d,
     >                           dn_in,agr_in,
     >                           x_parameter,c_parameter,alphac,
     >                           xce,fn,fdn)
      implicit none
      
      integer n2ft3d
      real*8 dn_in(n2ft3d,2)
      real*8 agr_in(n2ft3d,3)
      real*8 x_parameter,c_parameter,alphac
      real*8 xce(n2ft3d)
      real*8 fn(n2ft3d,2)
      real*8 fdn(n2ft3d,3)
      
*     **** Density cutoff parameter ****
      real*8 DNS_CUT,ETA,ETA2,alpha_zeta,alpha_zeta2
      parameter (DNS_CUT = 1.0d-20)
      parameter (ETA=1.0d-20)
      parameter (ETA2=1.0d-14)
      parameter (alpha_zeta=(1.0d0-ETA2))
      parameter (alpha_zeta2=(1.0d0-ETA2))

c     ***** PBE96 GGA exchange constants ******
      real*8 MU,KAPPA
      parameter (MU    = 0.2195149727645171d0)
      parameter (KAPPA = 0.8040000000000000d0)
 
c     ****** PBE96 GGA correlation constants ******
      real*8 GAMMA,BETA,BOG
      parameter (GAMMA	= 0.031090690869655d0)
      parameter (BETA	= 0.066724550603149d0)
      !parameter (BETA	= 0.066725d0)
      parameter (BOG    = BETA/GAMMA)


c     ****** Perdew-Wang92 LDA correlation coefficients *******
      real*8 GAM,iGAM,FZZ,iFZZ
      parameter (GAM  	= 0.519842099789746329d0)
      parameter (iGAM  	= 1.0d0/GAM)
      parameter (FZZ    = (8.0d0/(9.0d0*GAM)) )
      parameter (iFZZ    = 0.125d0*9.0d0*GAM)

      real*8 A_1,A1_1,B1_1,B2_1,B3_1,B4_1     
      parameter (A_1  = 0.0310907d0)
      !parameter (A_1  = 0.031091d0)
      parameter (A1_1 =	0.2137000d0)
      parameter (B1_1 =	7.5957000d0)
      parameter (B2_1 =	3.5876000d0)
      parameter (B3_1 =	1.6382000d0)
      parameter (B4_1 =	0.4929400d0)

      real*8 A_2,A1_2,B1_2,B2_2,B3_2,B4_2     
      parameter (A_2  =  0.01554535d0)
      !parameter (A_2  =  0.015545d0)
      parameter (A1_2 =	 0.20548000d0)
      parameter (B1_2 =	14.11890000d0)
      parameter (B2_2 =	 6.19770000d0)
      parameter (B3_2 =	 3.36620000d0)
      parameter (B4_2 =	 0.62517000d0)
      
      real*8 A_3,A1_3,B1_3,B2_3,B3_3,B4_3     
      parameter (A_3  =  0.0168869d0)
      !parameter (A_3  =  0.016887d0)
      parameter (A1_3 =	 0.1112500d0)
      parameter (B1_3 =	10.3570000d0)
      parameter (B2_3 =	 3.6231000d0)
      parameter (B3_3 =	 0.8802600d0)
      parameter (B4_3 =	 0.4967100d0)

c     **** other constants ****
      real*8 onethird,fourthird,fivethird,onesixthm
      real*8 twothird,sevensixthm
      real*8 onethirdm
      parameter (onethird=1.0d0/3.0d0)
      parameter (onethirdm=-1.0d0/3.0d0)
      parameter (twothird=2.0d0/3.0d0)
      parameter (fourthird=4.0d0/3.0d0)
      parameter (fivethird=5.0d0/3.0d0)
      parameter (onesixthm=-1.0d0/6.0d0)
      parameter (sevensixthm=-7.0d0/6.0d0)

*---- Beef expansion parameters given by Wellendorff et al-----------------*
      real*8 malphac
c      real*8 alphac,malphac
c      parameter (alphac=0.6001664769d0)
c      parameter (malphac=1.0d0-alphac)
      real*8 am(0:29)
      data am(0:29)/1.516501714d0,  4.413532099d-1,-9.182135241d-2,
     >            -2.352754331d-2, 3.418828455d-2, 2.411870076d-3,
     >            -1.416381352d-2, 6.975895581d-4, 9.859205137d-3,
     >            -6.737855051d-3,-1.573330824d-3, 5.036146253d-3,
     >            -2.569472453d-3,-9.874953976d-4, 2.033722895d-3,
     >            -8.018718848d-4,-6.688078723d-4, 1.030936331d-3,
     >            -3.673838660d-4,-4.213635394d-4,5.761607992d-4,
     >            -8.346503735d-5,-4.458447585d-4,4.601290092d-4,
     >            -5.231775398d-6,-4.239570471d-4,3.750190679d-4,
     >             2.114938125d-5, -1.904911565d-4,7.384362421d-5/

c     **** local variables ****
      integer i,j
      real*8 n,agr
      real*8 nup,agrup
      real*8 ndn,agrdn
      real*8 kf,ks,s,n_onethird,pi,rs_scale
      real*8 rs         ! Wigner radius
      real*8 rss        ! rss  = sqrt(rs)
      real*8 rs_n       ! rs_n = n*drs/dn
      real*8 t,t2,t4,t6
      real*8 t_nup      ! t_nup = n*dt/dnup
      real*8 t_ndn      ! t_ndn = n*dt/dndn
      real*8 t_agr      ! t_agr = n*dt/dagr
      real*8 zet,twoksg
      real*8 zet_nup    ! zet_nup = n*dzet/dnup
      real*8 zet_ndn    ! zet_nup = n*dzet/dnup
      real*8 zetp_1_3,zetm_1_3
      real*8 zetpm_1_3,zetmm_1_3
      real*8 phi,phi3,phi4
      real*8 phi_zet
      real*8 A,A2
      real*8 A_phi,A_ec_lda
      real*8 Q0,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8
      real*8 PON,FZ,z4
      real*8 tau
      real*8 F
      real*8 Fs                ! dF/ds
      real*8 Hpbe
      real*8 Hpbe_t            ! dHpbe/dt
      real*8 Hpbe_phi          ! dHpbe/dphi
      real*8 Hpbe_ec_lda       ! dHpbe/d(ec_lda)
      real*8 Hpbe_nup,Hpbe_ndn ! n*dHpbe/dnup, n*dHpbe/dndn
      real*8 Ipbe
      real*8 Ipbe_t,Ipbe_A     ! dIpbe/dt, dIpbe/dA

      real*8 exup,exdn,ex,ex_lda
      real*8 ecu,ecp,eca,ec,ec_lda
      real*8 ecu_rs,ecp_rs,eca_rs
      real*8 ec_lda_rs,ec_lda_zet  ! d(ec_lda)/drs, d(ec_lda)/dzet
      real*8 ec_lda_nup,ec_lda_ndn ! n*d(ec_lda)/dnup, n*d(ec_lda)/dndn
      real*8 fnxup,fdnxup          ! d(n*ex)/dnup, d(n*ex)/dndn
      real*8 fnxdn,fdnxdn          ! d(n*ex)/d|grad nup|, d(n*ex)/d|grad ndn|
      real*8 fncup,fncdn           ! d(n*ec)/dnup, d(n*ec)/dndn
      real*8 fdnx_const

      real*8 p0,p1,p2,dp,s2,oneovers2,Ft,dp2,sgn

      malphac=1.0d0-alphac

      pi = 4.0d0*datan(1.0d0)
      rs_scale = (0.75d0/pi)**onethird
      fdnx_const = -3.0d0/(8.0d0*pi)

      
!$OMP DO
      do i=1,n2ft3d
         nup     = dn_in(i,1)+ETA
         agrup   = agr_in(i,1)
 
         ndn     = dn_in(i,2)+ETA
         agrdn   = agr_in(i,2)
 
c        ****************************************************************
c        ***** calculate polarized Exchange energies and potentials *****
c        ****************************************************************

c        ************
c        **** up ****
c        ************
         n     = 2.0d0*nup
         agr   = 2.0d0*agrup

         n_onethird = (3.0d0*n/pi)**onethird
         ex_lda     = -0.75d0*n_onethird

         kf = (3.0d0*pi*pi*n)**onethird
         s  = agr/(2.0d0*kf*n)

         t  = 2.0d0*s*s/(4.0d0+s*s) - 1.0d0
         if (dabs(t).gt.1.0d0) t = dsign(1.0d0,t)
         F  = 0.0d0
         Ft = 0.0d0
         
         s2   = t*t - 1.0d0
         if (dabs(s2).lt.1.0d-12) then
            if (t.gt.0.0d0) then
               do j=0,29
                  F  = F  + am(j)
                  Ft = Ft + am(j)*0.5d0*dble(j*(j+1))
               end do
            else
               sgn = 1.0d0
               do j=0,29
                  F  = F  + sgn*am(j)
                  Ft = Ft + sgn*am(j)*0.5d0*dble(j*(j+1))
                  sgn = -sgn
               end do
            end if
         else
            oneovers2 = 1.0d0/s2
            p0 = 1.0d0
            dp = 0
            F  = F  + am(0)*p0
            Ft = Ft + am(0)*dp
            p1 = t
            dp = 1.0d0
            F  = F  + am(1)*t
            Ft = Ft + am(1)*dp
            do j=2,29
               p2    = (1.0d0/dble(j+1))*((2*j+1)*t*p1 - dble(j)*p0)
               dp2   = dble(j)*oneovers2*(t*p2-p1)
               F  = F + am(j)*p2
               Ft = Ft + am(j)*dp2
               p0 = p1
               p1 = p2
            end do
         end if
         Fs  = (16.0d0*s/(4.0d0+s*s)**2)*Ft
        
         exup = ex_lda*F
         fnxup = fourthird*(exup - ex_lda*Fs*s)
         fdnxup = fdnx_const*Fs

c        **************
c        **** down ****
c        **************
         n     = 2.0d0*ndn
         agr   = 2.0d0*agrdn

         n_onethird = (3.0d0*n/pi)**onethird
         ex_lda     = -0.75d0*n_onethird

         kf = (3.0d0*pi*pi*n)**onethird
         s  = agr/(2.0d0*kf*n)

         t  = 2.0d0*s*s/(4.0d0+s*s) - 1.0d0
         if (dabs(t).gt.1.0d0) t = dsign(1.0d0,t)
         F  = 0.0d0
         Ft = 0.0d0
         
         s2   = t*t - 1.0d0
         if (dabs(s2).lt.1.0d-12) then
            if (t.gt.0.0d0) then
               do j=0,29
                  F  = F  + am(j)
                  Ft = Ft + am(j)*0.5d0*dble(j*(j+1))
               end do
            else
               sgn = 1.0d0
               do j=0,29
                  F  = F  + sgn*am(j)
                  Ft = Ft + sgn*am(j)*0.5d0*dble(j*(j+1))
                  sgn = -sgn
               end do
            end if
         else
            oneovers2 = 1.0d0/s2
            p0 = 1.0d0
            dp = 0
            F  = F  + am(0)*p0
            Ft = Ft + am(0)*dp
            p1 = t
            dp = 1.0d0
            F  = F  + am(1)*t
            Ft = Ft + am(1)*dp
            do j=2,29
               p2    = (1.0d0/dble(j+1))*((2*j+1)*t*p1 - dble(j)*p0)
               dp2   = dble(j)*oneovers2*(t*p2-p1)
               F  = F + am(j)*p2
               Ft = Ft + am(j)*dp2
               p0 = p1
               p1 = p2
            end do
         end if
         Fs  = (16.0d0*s/(4.0d0+s*s)**2)*Ft

         exdn   = ex_lda*F
         fnxdn  = fourthird*(exdn - ex_lda*Fs*s)
         fdnxdn = fdnx_const*Fs

         n = nup+ndn

         ex = (exup*nup+ exdn*ndn)/ n
                  
c        *******************************************************************
c        ***** calculate polarized correlation energies and potentials ***** 
c        *******************************************************************
         agr   = agr_in(i,3)

         zet = (nup-ndn)/n
c         if (zet.gt.0.0d0) zet = zet - ETA2
c         if (zet.lt.0.0d0) zet = zet + ETA2
c        if (dabs(dn_in(i,2)).gt.DNS_CUT) zet_nup =  2*ndn/n**2
c        if (dabs(dn_in(i,1)).gt.DNS_CUT) zet_ndn = -2*nup/n**2
c        if (dabs(dn_in(i,2)).gt.DNS_CUT) zet_nup =  2*ndn/n
c        zet_nup =  2*ndn/n
c        zet_ndn = -2*nup/n
         zet_nup = -(zet - 1.0d0)
         zet_ndn = -(zet + 1.0d0)
         zetpm_1_3 = (1.0d0+zet*alpha_zeta)**onethirdm
         zetmm_1_3 = (1.0d0-zet*alpha_zeta)**onethirdm
         zetp_1_3  = (1.0d0+zet*alpha_zeta)*zetpm_1_3**2
         zetm_1_3  = (1.0d0-zet*alpha_zeta)*zetmm_1_3**2


         phi = 0.5d0*( zetp_1_3**2 + zetm_1_3**2)
         phi_zet = alpha_zeta*( zetpm_1_3 - zetmm_1_3)/3.0d0
         F =(  (1.0d0+zet*alpha_zeta)*zetp_1_3
     >       + (1.0d0-zet*alpha_zeta)*zetm_1_3
     >       - 2.0d0)*iGAM

         FZ = (zetp_1_3 - zetm_1_3)*(alpha_zeta*fourthird*iGAM)



*        **** calculate Wigner radius ****
         rs    = rs_scale/(n**onethird)
         rss   = dsqrt(rs)

*        **** calculate n*drs/dn ****
c        rs_n = onethirdm*rs/n
         rs_n = onethirdm*rs



c        **** calculate t ****
         kf = (3.0d0*pi*pi*n)**onethird
         ks = dsqrt(4.0d0*kf/pi)
        
         twoksg = 2.0d0*ks*phi
       
         t  = agr/(twoksg*n)

*        *** calculate n*dt/dnup, n*dt/dndn, n*dt/d|grad n| ****
         t_nup = sevensixthm*t - (phi_zet)*(zet_nup)*t/phi
         t_ndn = sevensixthm*t - (phi_zet)*(zet_ndn)*t/phi
         t_agr  = 1.0d0/(twoksg)


 
 
c        **************************************************
c        ***** compute LSDA correlation energy density ****
c        **************************************************
         call LSDT(A_1,A1_1,B1_1,B2_1,B3_1,B4_1,rss,ecu,ecu_rs)
         call LSDT(A_2,A1_2,B1_2,B2_2,B3_2,B4_2,rss,ecp,ecp_rs)
         call LSDT(A_3,A1_3,B1_3,B2_3,B3_3,B4_3,rss,eca,eca_rs)
         
         z4 = zet**4
            
         ec_lda = ecu*(1.0d0-F*z4) 
     >          + ecp*F*z4 
     >          - eca*F*(1.0d0-z4)/FZZ
         
         ec_lda_rs = ecu_rs*(1.0d0-F*z4)
     >             + ecp_rs*F*z4 
     >             - eca_rs*F*(1.0d0-z4)/FZZ

         ec_lda_zet = (4.0d0*(zet**3)*F + FZ*z4)*(ecp-ecu+eca*iFZZ)
     >              - FZ*eca*iFZZ


         
     
c        ********************************************
c        **** calculate PBE96 correlation energy ****
c        ********************************************
         phi3 = phi**3
         phi4 = phi3*phi
         PON  = -ec_lda/(phi3*GAMMA)
         tau  = DEXP(PON)

         A = BOG/(tau-1.0d0+ETA)
         A2 = A*A
         t2 = t*t
         t4 = t2*t2
         t6 = t4*t2
         Q4 = 1.0d0 + A*t2
         Q5 = 1.0d0 + 2.0d0*A*t2
         Q6 = 2.0d0 + A*t2
         Q7 = 1.0d0+A*t2+A2*t4
         Q8 = Q7*Q7

         Ipbe = 1.0d0 + BOG*t2*Q4/Q7
         Hpbe = GAMMA*phi3*DLOG(Ipbe)

         Ipbe_t =  BOG*(2.0d0*t)*Q5/Q8
         Ipbe_A = -BOG*(A*t6)   *Q6/Q8

         A_ec_lda  = tau/(BETA*phi3)*A2
         A_phi     = -3.0d0*ec_lda*tau/(BETA*phi4)*A2


         Hpbe_ec_lda = (GAMMA*phi3/Ipbe)*Ipbe_A*A_ec_lda

         Hpbe_phi    = 3.0d0*Hpbe/phi 
     >               + (GAMMA*phi3/Ipbe)*Ipbe_A*A_phi
         
         Hpbe_t      = (GAMMA*phi3/Ipbe)*Ipbe_t

         ec_lda_nup = ec_lda_zet 
     >              - zet * ec_lda_zet
     >              + rs_n * ec_lda_rs
         ec_lda_ndn = -ec_lda_zet
     >              - zet  * ec_lda_zet
     >              + rs_n * ec_lda_rs



         Hpbe_nup  = ec_lda_nup   * Hpbe_ec_lda
     >          + phi_zet*zet_nup * Hpbe_phi
     >          + t_nup           * Hpbe_t

         Hpbe_ndn  = ec_lda_ndn   * Hpbe_ec_lda
     >          + phi_zet*zet_ndn * Hpbe_phi
     >          + t_ndn           * Hpbe_t



         ec = alphac*ec_lda + malphac*Hpbe

         fncup  = ec + (alphac*ec_lda_nup + malphac*Hpbe_nup)
         fncdn  = ec + (alphac*ec_lda_ndn + malphac*Hpbe_ndn)

         xce(i)   = x_parameter*ex     + c_parameter*ec
         fn(i,1)  = x_parameter*fnxup  + c_parameter*fncup
         fn(i,2)  = x_parameter*fnxdn  + c_parameter*fncdn

         fdn(i,1) = x_parameter*fdnxup 
         fdn(i,2) = x_parameter*fdnxdn 
         fdn(i,3) = c_parameter*t_agr*Hpbe_t*malphac

      end do
!$OMP END DO
      
      
      
      return
      end
      




*    ************************************
*    *                                  *
*    *        gen_BEEF_BW_restricted    *
*    *                                  *
*    ************************************

*   This routine calculates the non-Langreth terms of the BEEF-vdw exchange-correlation 
*   potential(xcp) and energy density(xce).
*
*
*   Entry - n2ft3d     : number of grid points
*           rho_in(*) :  density (nup+ndn)
*           agr_in(*): |grad rho_in|
*           x_parameter: scale parameter for exchange
*           c_parameter: scale parameter for correlation
*
*     Exit  - xce(n2ft3d) : PBE96 exchange correlation energy density
*             fn(n2ft3d)  : d(n*xce)/dn
*             fdn(n2ft3d) : d(n*xce/d|grad n|
*
      subroutine gen_BEEF_BW_restricted(n2ft3d,rho_in,agr_in,
     >                                x_parameter,c_parameter,alphac,
     >                                xce,fn,fdn)
*
      implicit none
      integer    n2ft3d
      real*8     rho_in(n2ft3d)
      real*8     agr_in(n2ft3d)
      real*8     x_parameter,c_parameter,alphac
      real*8     xce(n2ft3d)
      real*8     fn(n2ft3d)
      real*8     fdn(n2ft3d)

*     **** Density cutoff parameter ****
      real*8 DNS_CUT,ETA
      parameter (DNS_CUT = 1.0d-20)
      parameter (ETA     = 1.0d-20)

c     ***** PBE96 GGA exchange constants ******
      real*8 MU,KAPPA
      parameter (MU    = 0.2195149727645171d0)
      parameter (KAPPA = 0.8040000000000000d0)

c     ****** PBE96 GGA correlation constants ******
      real*8 GAMMA,BETA,BOG
      parameter (GAMMA  = 0.031090690869655d0)
      parameter (BETA   = 0.066724550603149d0)
      parameter (BOG    = BETA/GAMMA)

c     ****** Perdew-Wang92 LDA correlation coefficients *******
      real*8 GAM,iGAM,FZZ,iFZZ
      parameter (GAM    = 0.519842099789746329d0)
      parameter (iGAM   = 1.0d0/GAM)
      parameter (FZZ    = (8.0d0/(9.0d0*GAM)) )
      parameter (iFZZ    = 0.125d0*9.0d0*GAM)

      real*8 A_1,A1_1,B1_1,B2_1,B3_1,B4_1
      parameter (A_1  = 0.0310907d0)
      !parameter (A_1  = 0.031091d0)
      parameter (A1_1 = 0.2137000d0)
      parameter (B1_1 = 7.5957000d0)
      parameter (B2_1 = 3.5876000d0)
      parameter (B3_1 = 1.6382000d0)
      parameter (B4_1 = 0.4929400d0)

      real*8 A_2,A1_2,B1_2,B2_2,B3_2,B4_2
      parameter (A_2  =  0.01554535d0)
      !parameter (A_2  =  0.015545d0)
      parameter (A1_2 =  0.20548000d0)
      parameter (B1_2 = 14.11890000d0)
      parameter (B2_2 =  6.19770000d0)
      parameter (B3_2 =  3.36620000d0)
      parameter (B4_2 =  0.62517000d0)

      real*8 A_3,A1_3,B1_3,B2_3,B3_3,B4_3
      parameter (A_3  =  0.0168869d0)
      !parameter (A_3  =  0.016887d0)
      parameter (A1_3 =  0.1112500d0)
      parameter (B1_3 = 10.3570000d0)
      parameter (B2_3 =  3.6231000d0)
      parameter (B3_3 =  0.8802600d0)
      parameter (B4_3 =  0.4967100d0)

c     **** other constants ****
      real*8 onethird,fourthird,fivethird,onesixthm
      real*8 twothird,sevensixthm,sevensixths
      real*8 onethirdm
      parameter (onethird=1.0d0/3.0d0)
      parameter (onethirdm=-1.0d0/3.0d0)
      parameter (twothird=2.0d0/3.0d0)
      parameter (fourthird=4.0d0/3.0d0)
      parameter (fivethird=5.0d0/3.0d0)
      parameter (onesixthm=-1.0d0/6.0d0)
      parameter (sevensixthm=-7.0d0/6.0d0)
      parameter (sevensixths=7.0d0/6.0d0)


*---- Beef expansion parameters given by Wellendorff et al-----------------*
      real*8 malphac
c      real*8 alphac,malphac
c      parameter (alphac=0.6001664769d0)
c      parameter (malphac=1.0d0-alphac)
      real*8 am(0:29)
      data am(0:29)/1.516501714d0,  4.413532099d-1,-9.182135241d-2,
     >            -2.352754331d-2, 3.418828455d-2, 2.411870076d-3,
     >            -1.416381352d-2, 6.975895581d-4, 9.859205137d-3,
     >            -6.737855051d-3,-1.573330824d-3, 5.036146253d-3,
     >            -2.569472453d-3,-9.874953976d-4, 2.033722895d-3,
     >            -8.018718848d-4,-6.688078723d-4, 1.030936331d-3,
     >            -3.673838660d-4,-4.213635394d-4,5.761607992d-4,
     >            -8.346503735d-5,-4.458447585d-4,4.601290092d-4,
     >            -5.231775398d-6,-4.239570471d-4,3.750190679d-4,
     >             2.114938125d-5, -1.904911565d-4,7.384362421d-5/

c     **** local variables ****
      integer i,j
      real*8 n,agr
      real*8 kf,ks,s,n_onethird,pi,rs_scale
      real*8 fdnx_const
      real*8 rs,rss,t,t2,t4,t6
      real*8 Q0,Q1,Q2,Q3,Q4,Q5,Q8,Q9,B
      real*8 Ht
      real*8 B_ec,Hrs,H_B
      real*8 F,Fs

      real*8 ex_lda,ec_lda
      real*8 ec_lda_rs
      real*8 ex,ec,H
      real*8 fnx,fdnx,fnc,fdnc


      real*8 p0,p1,p2,dp,s2,oneovers2,Ft,dp2,sgn


      malphac=1.0d0-alphac
      pi         = 4.0d0*datan(1.0d0)
      rs_scale   = (0.75d0/pi)**onethird
      fdnx_const = -3.0d0/(8.0d0*pi)

!$OMP DO
      do i=1,n2ft3d
         n     = rho_in(i)+ETA
         agr   = agr_in(i)

c        ***** calculate unpolarized Exchange energies and potentials *****
         n_onethird = (3.0d0*n/pi)**onethird
         ex_lda     = -0.75d0*n_onethird

         kf = (3.0d0*pi*pi*n)**onethird
         s  = agr/(2.0d0*kf*n)

         t  = 2.0d0*s*s/(4.0d0+s*s) - 1.0d0
         if (dabs(t).gt.1.0d0) t = dsign(1.0d0,t)
         F  = 0.0d0
         Ft = 0.0d0
         s2   = t*t - 1.0d0
         if (dabs(s2).lt.1.0d-12) then
            if (t.gt.0.0d0) then
               do j=0,29
                  F  = F  + am(j)
                  Ft = Ft + am(j)*0.5d0*dble(j*(j+1))
               end do
            else
               sgn = 1.0d0
               do j=0,29
                  F  = F  + sgn*am(j)
                  Ft = Ft + sgn*am(j)*0.5d0*dble(j*(j+1))
                  sgn = -sgn
               end do
            end if
         else
            oneovers2 = 1.0d0/s2
            p0 = 1.0d0
            dp = 0
            F  = F  + am(0)*p0
            Ft = Ft + am(0)*dp
            p1 = t
            dp = 1.0d0
            F  = F  + am(1)*t
            Ft = Ft + am(1)*dp
            do j=2,29
               p2    = (1.0d0/dble(j+1))*((2*j+1)*t*p1 - dble(j)*p0)
               dp2   = dble(j)*oneovers2*(t*p2-p1)
               F  = F + am(j)*p2
               Ft = Ft + am(j)*dp2
               p0 = p1
               p1 = p2
            end do
         end if
         Fs  = (16.0d0*s/(4.0d0+s*s)**2)*Ft

         ex   = ex_lda*F
         fnx  = fourthird*(ex - ex_lda*Fs*s)
         fdnx = fdnx_const*Fs


*        *********************************************************************
c        ***** calculate unpolarized correlation energies and potentials *****
*        *********************************************************************

c        **** calculate rs and t ****
         rs    = rs_scale/(n**onethird)
         rss   = dsqrt(rs)

         kf = (3.0d0*pi*pi*n)**onethird
         ks = dsqrt(4.0d0*kf/pi)
         t  = agr/(2.0d0*ks*n)


c        **** unpolarized LDA correlation energy ****
c        **** ec_p = correlation energy          ****
c        ****   ec_p_rs = dec_p/drs              ****
c        ****   uc_p    = dec_p/dn               ****
         call LSDT(A_1,A1_1,B1_1,B2_1,B3_1,B4_1,rss,ec_lda,ec_lda_rs)
c        **** PBE96 correlation energy  corrections ****
         t2 = t*t
         t4 = t2*t2
         B = -ec_lda/GAMMA
         B = BOG/(dexp(B)-1.0d0+ETA)
         Q4 = 1.0d0 + B*t2
         Q5 = 1.0d0 + B*t2 + B*B*t4
         H = GAMMA*dlog(1.0d0 + BOG*Q4*t2/Q5)


c        **** PBE96 correlation fdn and fdnc derivatives ****
         t6   = t4*t2

         B_ec = (B/BETA)*(BOG+B)

         Q8  = Q5*Q5+BOG*Q4*Q5*t2
         Q9  = 1.0d0+2*B*t2
         H_B  = -BETA*B*t6*(2.0d0+B*t2)/Q8
         Hrs  = H_B*B_ec*ec_lda_rs

         Ht  = 2.0d0*BETA*Q9/Q8*t

         ec   = alphac*ec_lda + malphac*H
         fnc  = ec -  alphac*(onethird*rs*ec_lda_rs)
     >             - malphac*(onethird*rs*Hrs)
     >             - malphac*(sevensixths*t*Ht)
         fdnc = malphac*(0.5d0* Ht/ks)
      
         xce(i) = x_parameter*ex   + c_parameter*ec
         fn(i)  = x_parameter*fnx  + c_parameter*fnc
         fdn(i) = x_parameter*fdnx + c_parameter*fdnc

      end do
!$OMP END DO

      return
      end

