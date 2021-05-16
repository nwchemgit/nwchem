*    ********************************************
*    *				         	*
*    *	    gen_PBE96_c_full_unrestricted	*
*    *				        	*
*    ********************************************
*
*    This function returns the PBE96 correlation
*  energy density, ce, and its derivatives with respect
*  to nup, ndn, |grad nup|, |grad ndn|, and |grad n|.
*
*   Entry - dn1_in,dn2_in          :  spin densites nup and ndn
*           agr1_in,agr2_in,agr3_in: |grad nup|, |grad ndn|, and |grad n|
*
*   Exit - ce            : PBE96 correlation energy density
*        - fn1,fn2       : d(n*ce)/dnup, d(n*ce)/dndn
*        - fdn1,fdn2,fdn3: d(n*ce)/d|grad nup|, d(n*ce)/d|grad ndn|
*                           d(n*ce)/d|grad n|

      subroutine gen_PBE96_c_full_unrestricted(dn1_in,dn2_in,
     >                           agr1_in,agr2_in,agr3_in,
     >                           ce,fn1,fn2,fdn1,fdn2,fdn3)
      implicit none
*     ***** input *****      
      real*8 dn1_in,dn2_in,agr1_in,agr2_in,agr3_in
*     ***** output *****      
      real*8 ce,fn1,fn2,fdn1,fdn2,fdn3
  
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

c     **** local variables ****
      real*8 n,agr
      real*8 nup,agrup
      real*8 ndn,agrdn
      real*8 kf,ks,s,P0,n_onethird,pi,rs_scale
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

      pi = 4.0d0*datan(1.0d0)
      rs_scale = (0.75d0/pi)**onethird
      fdnx_const = -3.0d0/(8.0d0*pi)

      
ccc!$OMP DO
      nup     = dn1_in+ETA
      agrup   = agr1_in
 
      ndn     = dn2_in+ETA
      agrdn   = agr2_in
 
                 
c     *******************************************************************
c     ***** calculate polarized correlation energies and potentials ***** 
c     *******************************************************************
      n  =  nup + ndn

      agr   = agr3_in

      zet = (nup-ndn)/n
c      if (zet.gt.0.0d0) zet = zet - ETA2
c      if (zet.lt.0.0d0) zet = zet + ETA2
c     if (dabs(dn_in(i,2)).gt.DNS_CUT) zet_nup =  2*ndn/n**2
c     if (dabs(dn_in(i,1)).gt.DNS_CUT) zet_ndn = -2*nup/n**2
c     if (dabs(dn_in(i,2)).gt.DNS_CUT) zet_nup =  2*ndn/n
c     zet_nup =  2*ndn/n
c     zet_ndn = -2*nup/n
      zet_nup = -(zet - 1.0d0)
      zet_ndn = -(zet + 1.0d0)
      zetpm_1_3 = (1.0d0+zet*alpha_zeta)**onethirdm
      zetmm_1_3 = (1.0d0-zet*alpha_zeta)**onethirdm
      zetp_1_3  = (1.0d0+zet*alpha_zeta)*zetpm_1_3**2
      zetm_1_3  = (1.0d0-zet*alpha_zeta)*zetmm_1_3**2


      phi = 0.5d0*( zetp_1_3**2 + zetm_1_3**2)
      phi_zet = alpha_zeta*( zetpm_1_3 - zetmm_1_3)/3.0d0
      F =(  (1.0d0+zet*alpha_zeta)*zetp_1_3
     >    + (1.0d0-zet*alpha_zeta)*zetm_1_3
     >    - 2.0d0)*iGAM

      FZ = (zetp_1_3 - zetm_1_3)*(alpha_zeta*fourthird*iGAM)



*     **** calculate Wigner radius ****
      rs    = rs_scale/(n**onethird)
      rss   = dsqrt(rs)

*     **** calculate n*drs/dn ****
c     rs_n = onethirdm*rs/n
      rs_n = onethirdm*rs



c     **** calculate t ****
      kf = (3.0d0*pi*pi*n)**onethird
      ks = dsqrt(4.0d0*kf/pi)
      
      twoksg = 2.0d0*ks*phi
      
      t  = agr/(twoksg*n)

*     *** calculate n*dt/dnup, n*dt/dndn, n*dt/d|grad n| ****
      t_nup = sevensixthm*t - (phi_zet)*(zet_nup)*t/phi
      t_ndn = sevensixthm*t - (phi_zet)*(zet_ndn)*t/phi
      t_agr  = 1.0d0/(twoksg)
 
 
c     **************************************************
c     ***** compute LSDA correlation energy density ****
c     **************************************************
      call LSDT(A_1,A1_1,B1_1,B2_1,B3_1,B4_1,rss,ecu,ecu_rs)
      call LSDT(A_2,A1_2,B1_2,B2_2,B3_2,B4_2,rss,ecp,ecp_rs)
      call LSDT(A_3,A1_3,B1_3,B2_3,B3_3,B4_3,rss,eca,eca_rs)
      
      z4 = zet**4
         
      ec_lda = ecu*(1.0d0-F*z4) 
     >       + ecp*F*z4 
     >       - eca*F*(1.0d0-z4)/FZZ
      
      ec_lda_rs = ecu_rs*(1.0d0-F*z4)
     >          + ecp_rs*F*z4 
     >          - eca_rs*F*(1.0d0-z4)/FZZ

      ec_lda_zet = (4.0d0*(zet**3)*F + FZ*z4)*(ecp-ecu+eca*iFZZ)
     >           - FZ*eca*iFZZ


      
     
c     ********************************************
c     **** calculate PBE96 correlation energy ****
c     ********************************************
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
     >            + (GAMMA*phi3/Ipbe)*Ipbe_A*A_phi
      
      Hpbe_t      = (GAMMA*phi3/Ipbe)*Ipbe_t

      ec_lda_nup = ec_lda_zet 
     >           - zet * ec_lda_zet
     >           + rs_n * ec_lda_rs
      ec_lda_ndn = -ec_lda_zet
     >           - zet  * ec_lda_zet
     >           + rs_n * ec_lda_rs



      Hpbe_nup  = ec_lda_nup   * Hpbe_ec_lda
     >       + phi_zet*zet_nup * Hpbe_phi
     >       + t_nup           * Hpbe_t

      Hpbe_ndn  = ec_lda_ndn   * Hpbe_ec_lda
     >       + phi_zet*zet_ndn * Hpbe_phi
     >       + t_ndn           * Hpbe_t



      ec = ec_lda + Hpbe

      fncup  = ec + (ec_lda_nup + Hpbe_nup)
      fncdn  = ec + (ec_lda_ndn + Hpbe_ndn)

      ce   = ec
      fn1  = fncup
      fn2  = fncdn

      fdn1 = 0.0d0 
      fdn2 = 0.0d0 
      fdn3 = t_agr*Hpbe_t

ccc!$OMP END DO
      
      
      
      return
      end
      

*    ************************************
*    *					*
*    *	    gen_PBE96_c_unrestricted	*
*    *					*
*    ************************************
*
*    This function returns the PBE96 correlation
*  energy density, ce, and its derivatives with respect
*  to n_sigma, |grad n_sigma|.
*
*   Entry - dn_in :  spin densites n_sigma
*           agr_in: |grad n_sigma|
*
*   Exit - ce : PBE96 correlation energy density
*        - fn : d(n_sigma*ce)/dn_sigma,
*        - fdn: d(n_sigma*ce)/d|grad n_sigma|

      subroutine gen_PBE96_c_unrestricted(dn_in,agr_in,
     >                           ce,fn,fdn)
      implicit none
*     ***** input *****      
      real*8 dn_in,agr_in
*     ***** output *****
      real*8 ce,fn,fdn
      
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

c     **** local variables ****
      real*8 n,agr
      real*8 nup,agrup
      real*8 kf,ks,s,P0,n_onethird,pi,rs_scale
      real*8 rs         ! Wigner radius
      real*8 rss        ! rss  = sqrt(rs)
      real*8 rs_n       ! rs_n = n*drs/dn
      real*8 t,t2,t4,t6
      real*8 t_nup      ! t_nup = n*dt/dnup
      real*8 t_agr      ! t_agr = n*dt/dagr
      real*8 zet,twoksg
      real*8 zet_nup    ! zet_nup = n*dzet/dnup
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
      real*8 Hpbe_nup ! n*dHpbe/dnup, n*dHpbe/dndn
      real*8 Ipbe
      real*8 Ipbe_t,Ipbe_A     ! dIpbe/dt, dIpbe/dA

      real*8 exup,exdn,ex,ex_lda
      real*8 ecu,ecp,eca,ec,ec_lda
      real*8 ecu_rs,ecp_rs,eca_rs
      real*8 ec_lda_rs,ec_lda_zet  ! d(ec_lda)/drs, d(ec_lda)/dzet
      real*8 ec_lda_nup ! n*d(ec_lda)/dnup, n*d(ec_lda)/dndn
      real*8 fncup           ! d(n*ec)/dnup, d(n*ec)/dndn
      real*8 fdnx_const

      pi = 4.0d0*datan(1.0d0)
      rs_scale = (0.75d0/pi)**onethird
      fdnx_const = -3.0d0/(8.0d0*pi)

      
ccc!$OMP DO
      nup     = dn_in+ETA
 
                 
c     *******************************************************************
c     ***** calculate polarized correlation energies and potentials ***** 
c     *******************************************************************
      n  =  nup

      agr   = agr_in

c     zet = nup/n 
      zet = 1.0d0
c      if (zet.gt.0.0d0) zet = zet - ETA2
c      if (zet.lt.0.0d0) zet = zet + ETA2
c     if (dabs(dn_in(i,2)).gt.DNS_CUT) zet_nup =  2*ndn/n**2
c     if (dabs(dn_in(i,1)).gt.DNS_CUT) zet_ndn = -2*nup/n**2
c     if (dabs(dn_in(i,2)).gt.DNS_CUT) zet_nup =  2*ndn/n
c     zet_nup =  2*ndn/n
c     zet_ndn = -2*nup/n

c     zet_nup = -(zet - 1.0d0)
      zet_nup   = 0.0d0
      zetpm_1_3 = (1.0d0+zet*alpha_zeta)**onethirdm
      zetmm_1_3 = (1.0d0-zet*alpha_zeta)**onethirdm
      zetp_1_3  = (1.0d0+zet*alpha_zeta)*zetpm_1_3**2
      zetm_1_3  = (1.0d0-zet*alpha_zeta)*zetmm_1_3**2


      phi = 0.5d0*( zetp_1_3**2 + zetm_1_3**2)
      phi_zet = alpha_zeta*( zetpm_1_3 - zetmm_1_3)/3.0d0
      F =(  (1.0d0+zet*alpha_zeta)*zetp_1_3
     >    + (1.0d0-zet*alpha_zeta)*zetm_1_3
     >    - 2.0d0)*iGAM

      FZ = (zetp_1_3 - zetm_1_3)*(alpha_zeta*fourthird*iGAM)



*     **** calculate Wigner radius ****
      rs    = rs_scale/(n**onethird)
      rss   = dsqrt(rs)

*     **** calculate n*drs/dn ****
c     rs_n = onethirdm*rs/n
      rs_n = onethirdm*rs



c     **** calculate t ****
      kf = (3.0d0*pi*pi*n)**onethird
      ks = dsqrt(4.0d0*kf/pi)
      
      twoksg = 2.0d0*ks*phi
      
      t  = agr/(twoksg*n)

*     *** calculate n*dt/dnup, n*dt/dndn, n*dt/d|grad n| ****
      t_nup = sevensixthm*t - (phi_zet)*(zet_nup)*t/phi
      t_agr  = 1.0d0/(twoksg)
 
 
c     **************************************************
c     ***** compute LSDA correlation energy density ****
c     **************************************************
      call LSDT(A_1,A1_1,B1_1,B2_1,B3_1,B4_1,rss,ecu,ecu_rs)
      call LSDT(A_2,A1_2,B1_2,B2_2,B3_2,B4_2,rss,ecp,ecp_rs)
      call LSDT(A_3,A1_3,B1_3,B2_3,B3_3,B4_3,rss,eca,eca_rs)
      
      z4 = zet**4
         
      ec_lda = ecu*(1.0d0-F*z4) 
     >       + ecp*F*z4 
     >       - eca*F*(1.0d0-z4)/FZZ
      
      ec_lda_rs = ecu_rs*(1.0d0-F*z4)
     >          + ecp_rs*F*z4 
     >          - eca_rs*F*(1.0d0-z4)/FZZ

      ec_lda_zet = (4.0d0*(zet**3)*F + FZ*z4)*(ecp-ecu+eca*iFZZ)
     >           - FZ*eca*iFZZ


      
     
c     ********************************************
c     **** calculate PBE96 correlation energy ****
c     ********************************************
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
     >            + (GAMMA*phi3/Ipbe)*Ipbe_A*A_phi
      
      Hpbe_t      = (GAMMA*phi3/Ipbe)*Ipbe_t

      ec_lda_nup = ec_lda_zet 
     >           - zet * ec_lda_zet
     >           + rs_n * ec_lda_rs



      Hpbe_nup  = ec_lda_nup   * Hpbe_ec_lda
     >       + phi_zet*zet_nup * Hpbe_phi
     >       + t_nup           * Hpbe_t




      ec = ec_lda + Hpbe

      fncup  = ec + (ec_lda_nup + Hpbe_nup)

      ce   = ec
      fn  = fncup

      fdn = t_agr*Hpbe_t

CCC!$OMP END DO
      
      
      
      return
      end
      




*    ************************************
*    *					*
*    *	    gen_PBE96_c_restricted	*
*    *					*
*    ************************************
*
*   This routine calculates the PBE96 correlation 
*   potential(cp) and energy density(ce).
*
*
*    Entry - rho_in:  density (nup+ndn)
*            agr_in: |grad rho_in|
*
*    Exit  - ce : PBE96 correlation energy density
*            fn : d(n*ce)/dn
*            fdn: d(n*ce/d|grad n|
*
      subroutine gen_PBE96_c_restricted(rho_in,agr_in,
     >                                ce,fn,fdn)
      implicit none
*     ****** input ******
      real*8     rho_in,agr_in
*     ****** output ******
      real*8     ce,fn,fdn

      
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
      parameter (GAMMA	= 0.031090690869655d0)
      parameter (BETA	= 0.066724550603149d0)
      parameter (BOG    = BETA/GAMMA)


c     ****** Perdew-Wang92 LDA correlation coefficients *******
      real*8 A_1,A1_1,B1_1,B2_1,B3_1,B4_1     
      parameter (A_1  = 0.0310907d0)
      parameter (A1_1 =	0.2137000d0)
      parameter (B1_1 =	7.5957000d0)
      parameter (B2_1 =	3.5876000d0)
      parameter (B3_1 =	1.6382000d0)
      parameter (B4_1 =	0.4929400d0)

      real*8 A_2,A1_2,B1_2,B2_2,B3_2,B4_2     
      parameter (A_2  =  0.01554535d0)
      parameter (A1_2 =	 0.20548000d0)
      parameter (B1_2 =	14.11890000d0)
      parameter (B2_2 =	 6.19770000d0)
      parameter (B3_2 =	 3.36620000d0)
      parameter (B4_2 =	 0.62517000d0)
      
      real*8 A_3,A1_3,B1_3,B2_3,B3_3,B4_3     
      parameter (A_3  =  0.0168869d0)
      parameter (A1_3 =	 0.1112500d0)
      parameter (B1_3 =	10.3570000d0)
      parameter (B2_3 =	 3.6231000d0)
      parameter (B3_3 =	 0.8802600d0)
      parameter (B4_3 =	 0.4967100d0)

c     **** other constants ****
      real*8 onethird,fourthird,sevensixths
      parameter (onethird=1.0d0/3.0d0)
      parameter (fourthird=4.0d0/3.0d0)
      parameter (sevensixths=7.0d0/6.0d0)

c     **** local variables ****
      integer i
      real*8 n,agr
      real*8 kf,ks,s,P0,n_onethird,pi,rs_scale
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


      pi         = 4.0d0*datan(1.0d0)
      rs_scale   = (0.75d0/pi)**onethird
      fdnx_const = -3.0d0/(8.0d0*pi)
      
ccc!$OMP DO
      n     = rho_in+ETA
      agr   = agr_in
      

*     *********************************************************************
c     ***** calculate unpolarized correlation energies and potentials *****
*     *********************************************************************

c     **** calculate rs and t ****
      rs    = rs_scale/(n**onethird)
      rss   = dsqrt(rs)

      kf = (3.0d0*pi*pi*n)**onethird
      ks = dsqrt(4.0d0*kf/pi)
      t  = agr/(2.0d0*ks*n)


c     **** unpolarized LDA correlation energy ****
c     **** ec_p = correlation energy          ****
c     ****   ec_p_rs = dec_p/drs              ****
c     ****   uc_p    = dec_p/dn               ****
      call LSDT(A_1,A1_1,B1_1,B2_1,B3_1,B4_1,rss,ec_lda,ec_lda_rs)
c     **** PBE96 correlation energy  corrections ****
      t2 = t*t
      t4 = t2*t2
      B = -ec_lda/GAMMA
      B = BOG/(exp(B)-1.0d0+ETA)
      Q4 = 1.0d0 + B*t2
      Q5 = 1.0d0 + B*t2 + B*B*t4
      H = GAMMA*dlog(1.0d0 + BOG*Q4*t2/Q5)


c     **** PBE96 correlation fdn and fdnc derivatives ****
      t6   = t4*t2

      B_ec = (B/BETA)*(BOG+B)

      Q8  = Q5*Q5+BOG*Q4*Q5*t2
      Q9  = 1.0d0+2*B*t2
      H_B  = -BETA*B*t6*(2.0d0+B*t2)/Q8 
      Hrs  = H_B*B_ec*ec_lda_rs

      Ht  = 2.0d0*BETA*Q9/Q8*t

      ec   = ec_lda + H
      fnc = ec  - (onethird*rs*ec_lda_rs)
     >          - (onethird*rs*Hrs)
     >          - (sevensixths*t*Ht)
      fdnc = 0.5d0* Ht/ks

      ce = ec
      fn  = fnc
      fdn = fdnc
      

c       write(*,*) "pbe96:",i,ec,fnc,fdnc


ccc!$OMP END DO

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      
 

 
