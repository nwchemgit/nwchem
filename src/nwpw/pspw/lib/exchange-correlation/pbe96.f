
*    ************************************
*    *					*
*    *	    gen_PBE96_unrestricted	*
*    *					*
*    ************************************
*
*    This function returns the PBE96 exchange-correlation
*  energy density, xce, and its derivatives with respect
*  to nup, ndn, |grad nup|, |grad ndn|, and |grad n|.
*
*   Entry - n2ft3d     : number of grid points
*           dn_in(*,2) : spin densites nup and ndn
*           agr_in(*,3): |grad nup|, |grad ndn|, and |grad n|
*
*   Exit - xce(*)  : PBE96 energy density
*        - fn(*,2) : d(n*xce)/dnup, d(n*xce)/dndn
*        - fdn(*,3): d(n*xce)/d|grad nup|, d(n*xce)/d|grad ndn|
*                    d(n*xce)/d|grad n|

      subroutine gen_PBE96_BW_unrestricted(n2ft3d,
     >                           dn_in,agr_in,
     >                           xce,fn,fdn)
      implicit none
      
      integer n2ft3d
      real*8 dn_in(n2ft3d,2)
      real*8 agr_in(n2ft3d,3)
      real*8 xce(n2ft3d)
      real*8 fn(n2ft3d,2)
      real*8 fdn(n2ft3d,3)
      
*     **** Density cutoff parameter ****
      real*8 DNS_CUT,ETA,ETA2
      parameter (DNS_CUT = 1.0d-20)
      parameter (ETA=1.0d-20)
      parameter (ETA2=1.0d-14)

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
      real*8 GAM,FZZ
      parameter (GAM  	= 0.519842099789746329d0)
      parameter (FZZ    = (8.0d0/(9.0d0*GAM)) )

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
      real*8 onethird,fourthird,onesixthm
      real*8 twothird,sevensixthm
      real*8 onethirdm
      parameter (onethird=1.0d0/3.0d0)
      parameter (onethirdm=-1.0d0/3.0d0)
      parameter (twothird=2.0d0/3.0d0)
      parameter (fourthird=4.0d0/3.0d0)
      parameter (onesixthm=-1.0d0/6.0d0)
      parameter (sevensixthm=-7.0d0/6.0d0)

c     **** local variables ****
      integer i
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
         P0 = 1.0d0 + (MU/KAPPA)*s*s

         F   = (1.0d0 + KAPPA - KAPPA/P0)
         Fs  = 2.0d0*MU/(P0*P0)*s

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
         P0 = 1.0d0 + (MU/KAPPA)*s*s

         F   = (1.0d0 + KAPPA - KAPPA/P0)
         Fs  = 2.0d0*MU/(P0*P0)*s

         exdn   = ex_lda*F
         fnxdn  = fourthird*(exdn - ex_lda*Fs*s)
         fdnxdn = fdnx_const*Fs

			
         ex = (exup*nup+ exdn*ndn)/ (nup+ndn)
         
         
c        *******************************************************************
c        ***** calculate polarized correlation energies and potentials ***** 
c        *******************************************************************
         n     = dn_in(i,1) + dn_in(i,2) + ETA
         agr   = agr_in(i,3)

         zet = (nup-ndn)/n
         if (zet.gt.0.0d0) zet = zet - ETA2
         if (zet.lt.0.0d0) zet = zet + ETA2
c        if (dabs(dn_in(i,2)).gt.DNS_CUT) zet_nup =  2*ndn/n**2
c        if (dabs(dn_in(i,1)).gt.DNS_CUT) zet_ndn = -2*nup/n**2
c        if (dabs(dn_in(i,2)).gt.DNS_CUT) zet_nup =  2*ndn/n
c        zet_nup =  2*ndn/n
c        zet_ndn = -2*nup/n
         zet_nup = -(zet - 1.0d0)
         zet_ndn = -(zet + 1.0d0)

         phi = ( (1.0d0+zet)**twothird + (1.0d0-zet)**twothird)/2.0d0

         phi_zet = ( (1.0d0+zet)**onethirdm 
     >             - (1.0d0-zet)**onethirdm)/3.0d0


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
c        t  = 0.25d0*(pi/3.0)**onesixth*agr*(n**sevensixthm)/phi

*        *** calculate n*dt/dnup, n*dt/dndn, n*dt/d|grad n| ****
c        t_nup = sevensixthm*t/n - (phi_zet)*(zet_nup)*t/phi
c        t_ndn = sevensixthm*t/n - (phi_zet)*(zet_ndn)*t/phi
c        t_agr  = 1.0d0/(twoksg*n)
         t_nup = sevensixthm*t - (phi_zet)*(zet_nup)*t/phi
         t_ndn = sevensixthm*t - (phi_zet)*(zet_ndn)*t/phi
         t_agr  = 1.0d0/(twoksg)


 
 
c        **************************************************
c        ***** compute LSDA correlation energy density ****
c        **************************************************
         Q0 = -2.0d0*A_1*(1.0d0+A1_1*rs)
         Q1 =  2.0d0*A_1*rss*(B1_1+rss*(B2_1+rss*(B3_1+B4_1*rss)))
         Q2 = dlog(1.0d0+1.0d0/Q1)
         Q3 = A_1*(B1_1/rss + 2.0d0*B2_1 
     >                      + rss*(3.0d0*B3_1+4.0d0*B4_1*rss))

         ecu    = Q0*Q2
         ecu_rs = -2.0d0*A_1*A1_1*Q2-Q0*Q3/(Q1*(1.0d0+Q1))
      

      
         Q0 = -2.0d0*A_2*(1.0d0+A1_2*rs)
         Q1 =  2.0d0*A_2*rss*(B1_2+rss*(B2_2+rss*(B3_2+B4_2*rss)))
         Q2 = dlog(1.0d0+1.0d0/Q1)
         Q3 = A_2*(B1_2/rss + 2.0d0*B2_2 
     >                      + rss*(3.0d0*B3_2+4.0d0*B4_2*rss))

         ecp    = Q0*Q2
         ecp_rs = -2.0d0*A_2*A1_2*Q2-Q0*Q3/(Q1*(1+Q1))
            
                  
         Q0 = -2.0d0*A_3*(1.0d0+A1_3*rs)
         Q1 =  2.0d0*A_3*rss*(B1_3+rss*(B2_3+rss*(B3_3+B4_3*rss)))
         Q2 = dlog(1.0d0+1.0d0/Q1)
         Q3 = A_3*(B1_3/rss + 2.0d0*B2_3 
     >                      + rss*(3.0d0*B3_3+4.0d0*B4_3*rss))

         eca    = Q0*Q2
         eca_rs = -2.0d0*A_3*A1_3*Q2-Q0*Q3/(Q1*(1+Q1))
         
         z4 = zet**4
         F =(  (1.0d0+zet)**fourthird
     >       + (1.0d0-zet)**fourthird
     >       - 2.0d0)/GAM
            
         ec_lda = ecu*(1.0d0-F*z4) + ecp*F*z4 
     >          - eca*F*(1.0d0-z4)/FZZ
         
         ec_lda_rs = ecu_rs*(1.0d0-F*z4)
     >             + ecp_rs*F*z4 
     >             - eca_rs*F*(1.0d0-z4)/FZZ
         FZ = fourthird*((1.D0+zet)**onethird-(1.D0-zet)**onethird)
     >       /GAM
         ec_lda_zet = 4.D0*(zet**3)*F*(ecp-ecu+eca/FZZ)
     >          + FZ*(z4*ecp-z4*ecu
     >          - (1.D0-z4)*eca/FZZ)

         
         
     
c        ********************************************
c        **** calculate PBE96 correlation energy ****
c        ********************************************
         phi3 = phi**3
         phi4 = phi3*phi
         PON  = -ec_lda/(phi3*GAMMA)
         tau  = DEXP(PON)

         A = BOG/(tau-1.D0+ETA)
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

        ec_lda_nup = zet_nup * ec_lda_zet
     >             + rs_n    * ec_lda_rs
        ec_lda_ndn = zet_ndn * ec_lda_zet
     >             + rs_n    * ec_lda_rs



         Hpbe_nup  = ec_lda_nup   * Hpbe_ec_lda
     >          + phi_zet*zet_nup * Hpbe_phi
     >          + t_nup           * Hpbe_t

         Hpbe_ndn  = ec_lda_ndn   * Hpbe_ec_lda
     >          + phi_zet*zet_ndn * Hpbe_phi
     >          + t_ndn           * Hpbe_t


         ec = ec_lda + Hpbe
c        fncup  = ec + n*(ec_lda_nup + Hpbe_nup)
c        fncdn  = ec + n*(ec_lda_ndn + Hpbe_ndn)
         fncup  = ec + (ec_lda_nup + Hpbe_nup)
         fncdn  = ec + (ec_lda_ndn + Hpbe_ndn)

         xce(i)   = ex + ec
         fn(i,1)  = fnxup  + fncup
         fn(i,2)  = fnxdn  + fncdn

         fdn(i,1) = fdnxup 
         fdn(i,2) = fdnxdn 
c        fdn(i,3) = n*t_agr*Hpbe_t
         fdn(i,3) =  t_agr*Hpbe_t
      end do
      
      
      
      return
      end
      


*    ************************************
*    *					*
*    *	    gen_PBE96_BW_restricted	*
*    *					*
*    ************************************
*
*   This routine calculates the PBE96 exchange-correlation 
*   potential(xcp) and energy density(xce).
*
*
*   Entry - n2ft3d     : number of grid points
*           rho_in(*) :  density (nup+ndn)
*           agr_in(*): |grad rho_in|
*
*     Exit  - xce(n2ft3d) : PBE96 exchange correlation energy density
*             fn(n2ft3d)  : d(n*xce)/dn
*             fdn(n2ft3d) : d(n*xce/d|grad n|
*
      subroutine gen_PBE96_BW_restricted(n2ft3d,rho_in,agr_in,
     >                                xce,fn,fdn)
      implicit none

      integer    n2ft3d
      real*8     rho_in(n2ft3d)
      real*8     agr_in(n2ft3d)
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
      
      do i=1,n2ft3d
         n     = rho_in(i)+ETA
         agr   = agr_in(i)
        
c        ***** calculate unpolarized Exchange energies and potentials *****
         n_onethird = (3.0d0*n/pi)**onethird
         ex_lda     = -0.75d0*n_onethird

         kf = (3.0d0*pi*pi*n)**onethird
         s  = agr/(2.0d0*kf*n)
         P0 = 1.0d0 + (MU/KAPPA)*s*s

c        if (n.gt.DNS_CUT) then
c           F   = (1.0d0 + KAPPA - KAPPA/P0)
c           Fs  = 2.0d0*MU/(P0*P0)*s
c        else
c           F   = 1.0d0
c           Fs  = 0.0d0
c        end if
         F   = (1.0d0 + KAPPA - KAPPA/P0)
         Fs  = 2.0d0*MU/(P0*P0)*s

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
         t  = agr/(2.0*ks*n)


c        **** unpolarized LDA correlation energy ****
c        **** ec_p = correlation energy          ****
c        ****   ec_p_rs = dec_p/drs              ****
c        ****   uc_p    = dec_p/dn               ****
         Q0 = -2.0d0*A_1*(1.0d0+A1_1*rs)
         Q1 =  2.0d0*A_1*rss*(B1_1+rss*(B2_1+rss*(B3_1+B4_1*rss)))
         Q2 = dlog(1.0d0+1.0d0/Q1)
         Q3 = A_1*(B1_1/rss + 2.0d0*B2_1 
     >                      + rss*(3.0d0*B3_1+4.0d0*B4_1*rss))

         ec_lda    = Q0*Q2
         ec_lda_rs = -2.0d0*A_1*A1_1*Q2-Q0*Q3/(Q1*(1+Q1))
      

c        **** PBE96 correlation energy  corrections ****
         t2 = t*t
         t4 = t2*t2
         B = -ec_lda/GAMMA
         B = BOG/(exp(B)-1.0d0+ETA)
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

         ec   = ec_lda + H
         fnc = ec  - (onethird*rs*ec_lda_rs)
     >             - (onethird*rs*Hrs)
     >             - (sevensixths*t*Ht)
         fdnc = 0.5d0* Ht/ks

         xce(i) = ex + ec
         fn(i)  = fnx  + fnc
         fdn(i) = fdnx + fdnc
         

      end do

      return
      end

 
