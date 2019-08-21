*
* $Id$
*

*    ************************************
*    *					*
*    *	    gen_PW92_unrestricted	*
*    *					*
*    ************************************
*
*    This function returns the PW92 LDA correlation
*  energy density, xce, and its derivatives with respect
*  to nup, ndn.
*
*   Entry - n2ft3d     : number of grid points
*           dn_in(*,2) : spin densites nup and ndn
*
*   Exit - xce(*)  : PBE96 energy density
*        - fn(*,2) : d(n*xce)/dnup, d(n*xce)/dndn

      subroutine gen_PW92_BW_unrestricted(n2ft3d,
     >                           dn_in,
     >                           xce,fn)
      implicit none
      
      integer n2ft3d
      real*8 dn_in(n2ft3d,2)
      real*8 xce(n2ft3d)
      real*8 fn(n2ft3d,2)
      
*     **** Density cutoff parameter ****
      real*8 DNS_CUT,ETA,ETA2,alpha_zeta,alpha_zeta2
      parameter (DNS_CUT = 1.0d-20)
      parameter (ETA=1.0d-20)
      parameter (ETA2=1.0d-14)
      parameter (alpha_zeta=(1.0d0-ETA2))
      parameter (alpha_zeta2=(1.0d0-ETA2))

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
      real*8 onethird,fourthird
      real*8 onethirdm
      parameter (onethird=1.0d0/3.0d0)
      parameter (onethirdm=-1.0d0/3.0d0)
      parameter (fourthird=4.0d0/3.0d0)

c     **** local variables ****
      integer i
      real*8 n
      real*8 nup
      real*8 ndn
      real*8 n_onethird,pi,rs_scale
      real*8 rs         ! Wigner radius
      real*8 rss        ! rss  = sqrt(rs)
      real*8 rs_n       ! rs_n = n*drs/dn
      real*8 zet
      real*8 zet_nup    ! zet_nup = n*dzet/dnup
      real*8 zet_ndn    ! zet_nup = n*dzet/dnup
      real*8 zetp_1_3,zetm_1_3
      real*8 zetpm_1_3,zetmm_1_3
      real*8 FZ,z4
      real*8 F

      real*8 ecu,ecp,eca,ec,ec_lda
      real*8 ecu_rs,ecp_rs,eca_rs
      real*8 ec_lda_rs,ec_lda_zet  ! d(ec_lda)/drs, d(ec_lda)/dzet
      real*8 ec_lda_nup,ec_lda_ndn ! n*d(ec_lda)/dnup, n*d(ec_lda)/dndn
      real*8 fncup,fncdn           ! d(n*ec)/dnup, d(n*ec)/dndn

      pi = 4.0d0*datan(1.0d0)
      rs_scale = (0.75d0/pi)**onethird

      
!$OMP DO
      do i=1,n2ft3d
         nup     = dn_in(i,1)+ETA
         ndn     = dn_in(i,2)+ETA
 
c        *******************************************************************
c        ***** calculate polarized correlation energies and potentials ***** 
c        *******************************************************************
         n = nup+ndn
         zet = (nup-ndn)/n
         zet_nup = -(zet - 1.0d0)
         zet_ndn = -(zet + 1.0d0)
         zetpm_1_3 = (1.0d0+zet*alpha_zeta)**onethirdm
         zetmm_1_3 = (1.0d0-zet*alpha_zeta)**onethirdm
         zetp_1_3  = (1.0d0+zet*alpha_zeta)*zetpm_1_3**2
         zetm_1_3  = (1.0d0-zet*alpha_zeta)*zetmm_1_3**2

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


         ec_lda_nup = ec_lda_zet 
     >              - zet * ec_lda_zet
     >              + rs_n * ec_lda_rs
         ec_lda_ndn = -ec_lda_zet
     >              - zet  * ec_lda_zet
     >              + rs_n * ec_lda_rs


         ec = ec_lda 

         fncup  = ec + (ec_lda_nup)
         fncdn  = ec + (ec_lda_ndn)

         xce(i)   = ec
         fn(i,1)  = fncup
         fn(i,2)  = fncdn
      end do
!$OMP END DO
      
      return
      end
      




*    ************************************
*    *					*
*    *	    gen_PW92_BW_restricted	*
*    *					*
*    ************************************
*
*   This routine calculates the PW92 LDA correlation 
*   potential(xcp) and energy density(xce).
*
*
*   Entry - n2ft3d     : number of grid points
*           rho_in(*) :  density (nup+ndn)
*
*     Exit  - xce(n2ft3d) : PW92 correlation energy density
*             fn(n2ft3d)  : d(n*xce)/dn
*
      subroutine gen_PW92_BW_restricted(n2ft3d,rho_in,
     >                                xce,fn)
      implicit none

      integer    n2ft3d
      real*8     rho_in(n2ft3d)
      real*8     xce(n2ft3d)
      real*8     fn(n2ft3d)

      
*     **** Density cutoff parameter ****
      real*8 DNS_CUT,ETA
      parameter (DNS_CUT = 1.0d-20)
      parameter (ETA     = 1.0d-20)

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
      real*8 onethird
      parameter (onethird=1.0d0/3.0d0)

c     **** local variables ****
      integer i
      real*8 n
      real*8 n_onethird,pi,rs_scale
      real*8 rs,rss

      real*8 ec_lda
      real*8 ec_lda_rs
      real*8 ec
      real*8 fnc


      pi         = 4.0d0*datan(1.0d0)
      rs_scale   = (0.75d0/pi)**onethird
      
!$OMP DO
      do i=1,n2ft3d
         n     = rho_in(i)+ETA
        
*        *********************************************************************
c        ***** calculate unpolarized correlation energies and potentials *****
*        *********************************************************************

c        **** calculate rs and t ****
         rs    = rs_scale/(n**onethird)
         rss   = dsqrt(rs)


c        **** unpolarized LDA correlation energy ****
c        **** ec_p = correlation energy          ****
c        ****   ec_p_rs = dec_p/drs              ****
c        ****   uc_p    = dec_p/dn               ****
         call LSDT(A_1,A1_1,B1_1,B2_1,B3_1,B4_1,rss,ec_lda,ec_lda_rs)

         ec   = ec_lda
         fnc = ec  - (onethird*rs*ec_lda_rs)

         xce(i) = ec
         fn(i)  = fnc
         

      end do
!$OMP END DO

      return
      end


 
