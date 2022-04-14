!***********************************************************************
! Module description
!***********************************************************************
!
! JWF : this subroutine computes the exchange energy of SCAN
!  and some of its regularised variants.
!
!  Three Integer flags are required to select SCAN variant:
!   IALPHA  :  iso-orbital indicator
!               0: alpha        (scan)
!               1: alpha'       (rscan)
!               2: \bar{alpha}  (r++scan, r2scan, r4scan)
!
!   IINTERP :  interpolation function
!               0: scan
!               1: rscan, r++scan, r2scan, r4scan
!
!   IDELFX  :  gradient expansion correction
!               0: scan (scan, rscan, r++scan)
!               1: 2nd order (r2scan)
!               2: 4th order (r4scan)
!
!
!  Hence, functionals are accessed by passing:
!   SCAN:       (0, 0, 0)
!   rSCAN:      (1, 1, 0)
!   r++SCAN:    (2, 1, 0)
!   r2SCAN:     (2, 1, 1)
!   r4SCAN:     (2, 1, 2)
!
! Author: James W. Furness
! eMail : jfurness@tulane.edu (james.w.furness.1@gmail.com)
! Date  : 24/06/2020
!
! This work is made available under the CC0 1.0 Universal (CC0 1.0)
! Public Domain Dedication.
! https://creativecommons.org/publicdomain/zero/1.0/
!
! The person who associated a work with this deed has dedicated the work
! to the public domain by waiving all of his or her rights to the work
! worldwide under copyright law, including all related and neighboring
! rights, to the extent allowed by law.
!
! You can copy, modify, distribute and perform the work, even for
! commercial purposes, all without asking permission. See Other Information
! below.
!
! Other Information:
!
! In no way are the patent or trademark rights of any person affected by CC0,
! nor are the rights that other persons may have in the work or in how the work
! is used, such as publicity or privacy rights.
!
! Unless expressly stated otherwise, the person who associated a work with this
! deed makes no warranties about the work, and disclaims liability for all uses
! of the work, to the fullest extent permitted by applicable law.
!
! When using or citing the work, you should not imply endorsement by the author
! or the affirmer.
!
! While we have made every effort to ensure the code's correctness, it is provided
! as is and no warranties or guarantees are given.
!***********************************************************************

!===============================================================================
!   NWChem specific wrappers, gen_r2SCAN_restricted and gen_r2SCAN_unrestricted
!   author: Aaron Kaplan (kaplan@temple.edu)
!===============================================================================

!     ************************************************
!     *                                              *
!     *              gen_r2SCAN_restricted           *
!     *                                              *
!     ************************************************

!    This function returns the r2SCAN exchange-correlation
!  energy density, xce, and its derivatives with respect
!  to n, |grad n|, tau.


!   Entry - n2ft3d   : number of grid points
!           rho_in(*) :  density (nup+ndn)
!           agr_in(*): |grad rho_in|
!           tau_in(*): tau
!           x_parameter: scale parameter for exchange
!           c_parameter: scale parameter for correlation

!     Exit  - xce(n2ft3d) : r2SCAN exchange correlation energy per electron, NOT density
!             fn(n2ft3d)  : d(n*xce)/dn
!             fdn(n2ft3d) : d(n*xce)/d|grad n|
!             fdtau(n2ft3d) : d(n*xce)/dtau


      subroutine gen_r2SCAN_restricted(n2ft3d,rho_in,agr_in,tau_in,
     &                               x_parameter,c_parameter,
     &                               xce,fn,fdn,fdtau)

      implicit none

!     ***** input *****
      integer n2ft3d
      real*8 rho_in(*),agr_in(*),tau_in(*)
      real*8 x_parameter,c_parameter

!     ***** output *****
      real*8 xce(*),fn(*),fdn(*),fdtau(*)

!     ***** local declarations *****
      integer i
      real*8 n,agr,tau, ex, dexdn,dexdg,dexdt,nh,agrh,th
      real*8 ec, decdn0, decdn1, decdg0,decdg1, decdt0,decdt1, decdg

!     ***** density cutoff parameters *****
      real*8 tol,ETA
      parameter (tol = 1.0d-18)
      parameter (ETA = 1.0d-20)

!$OMP DO
      do i=1,n2ft3d

        n = rho_in(i) + ETA
        agr = agr_in(i) + ETA
        tau = 2.0d0*tau_in(i) + ETA

        call eps_SCANx(n, agr, tau, ex, dexdn, dexdg, dexdt, 0,0,0)!2,1,1)

        nh = n/2.d0
        agrh = agr/2.d0
        th = tau/2.d0
        call vrSCANc(nh, nh, agrh, agrh, agr, th,th, ec, decdn0,decdn1,
     &     decdg0, decdg1, decdg, decdt0, decdt1, 0,0,0)!2,1,1)

        xce(i)   = (x_parameter*ex + c_parameter*ec)/n
        fn(i)    = x_parameter*dexdn + c_parameter*decdn0
        fdn(i)   = x_parameter*dexdg + c_parameter*decdg
        fdtau(i) = x_parameter*dexdt + c_parameter*decdt0

      end do
!$OMP END DO

      return
      end subroutine gen_r2SCAN_restricted

!     ************************************************
!     *                                              *
!     *            gen_r2SCAN_unrestricted           *
!     *                                              *
!     ************************************************

!    This function returns the r2SCAN exchange-correlation
!  energy density, xce, and its derivatives with respect
!  to nup, ndn, |grad nup|, |grad ndn|, |grad n|, tauup, taudn.


!   Entry - n2ft3d   : number of grid points
!           rho_in(*,2) :  density (nup and ndn)
!           agr_in(*,3): |grad rho_in| (nup, ndn and n)
!           tau_in(*,2): tau (nup and ndn)
!           x_parameter: scale parameter for exchange
!           c_parameter: scale parameter for correlation

!     Exit  - xce(n2ft3d) : r2SCAN exchange correlation energy per electron, NOT density
!             fn(n2ft3d,2)  : d(n*xce)/dnup, d(n*xce)/dndn
!             fdn(n2ft3d,3) : d(n*xce)/d|grad nup|, d(n*xce)/d|grad ndn|, d(n*xce)/d|grad n|
!             fdtau(n2ft3d,2) : d(n*xce)/dtauup, d(n*xce)/dtaudn


      subroutine gen_r2SCAN_unrestricted(n2ft3d,rho_in,agr_in,tau_in,
     &                               x_parameter,c_parameter,
     &                               xce,fn,fdn,fdtau)
      implicit none
!     ***** input *****
      integer n2ft3d
      real*8 rho_in(n2ft3d,2),agr_in(n2ft3d,3),tau_in(n2ft3d,2)
      real*8 x_parameter,c_parameter
!     ***** output *****
      real*8 xce(n2ft3d),fn(n2ft3d,2),fdn(n2ft3d,3),fdtau(n2ft3d,2)
!     ***** local declarations *****
      integer i
      real*8 agr,n0,agr0,t0,n1,agr1,t1
      real*8 ex0, dexdn0,dexdg0,dexdt0
      real*8 ex1, dexdn1,dexdg1,dexdt1
      real*8 ec, decdn0, decdn1, decdg0,decdg1,decdt0,decdt1, decdg

!     ***** density cutoff parameters *****
      real*8 tol,ETA
      parameter (tol = 1.0d-18)
      parameter (ETA = 1.0d-20)


!$OMP DO
      do i=1,n2ft3d

        n0 = rho_in(i,1) + ETA
        agr0 = agr_in(i,1) + ETA
        t0 = tau_in(i,1) + ETA

        n1 = rho_in(i,2) + ETA
        agr1 = agr_in(i,2) + ETA
        t1 = tau_in(i,2) + ETA

        agr = agr_in(i,3) + ETA

        call eps_SCANx(2.d0*n0, 2.d0*agr0, 2.d0*t0, ex0, dexdn0,dexdg0,
     &     dexdt0, 0,0,0)!2,1,1)
        call eps_SCANx(2.d0*n1, 2.d0*agr1, 2.d0*t1, ex1, dexdn1,dexdg1,
     &     dexdt1, 0,0,0)!2,1,1)

        call vrSCANc(n0, n1, agr0, agr1, agr, t0, t1, ec, decdn0,
     &     decdn1, decdg0, decdg1, decdg, decdt0, decdt1, 0,0,0)!2,1,1)

        xce(i) = (x_parameter*(ex0 + ex1)/2.d0 
     &         + c_parameter*ec)/(n0 + n1)

        fn(i,1) = x_parameter*dexdn0 + c_parameter*decdn0
        fn(i,2) = x_parameter*dexdn1 + c_parameter*decdn1

        fdn(i,1) = x_parameter*dexdg0
        fdn(i,2) = x_parameter*dexdg1

        fdn(i,3) = c_parameter*decdg

        fdtau(i,1) = x_parameter*dexdt0 + c_parameter*decdt0
        fdtau(i,2) = x_parameter*dexdt1 + c_parameter*decdt1

      end do
!$OMP END DO

      return
      end subroutine gen_r2SCAN_unrestricted


!================================================================
!   SCAN family exchange
!================================================================

      subroutine eps_SCANx(den, grd, tau, eps_x, dedd, dedg, dedt, 
     &                 IALPHA, IINTERP, IDELFX)
          IMPLICIT NONE
      !       Note that this routine expects spin scaling relation to be
      !       used for total exchange energy.
      !       e.g.
      !           ex_total = (eps_SCANx(2*d0, 2*g0, 2*t0) + eps_SCANx(2*d1, 2*g1, 2*t1))/2.0
      !               where d0 is density for spin 0 etc.
      !
      !       Inputs:
      !           den : electron density (n)
      !           grd : |grad n|. Note this is not squared!
      !           tau : Kinetic energy density
      !       Outputs:
      !           eps_x : energy density
      !           dedd  : energy density derivative w.r.t. density
      !           dedg  : energy density derivative w.r.t. gradient
      !           dedt  : energy density derivative w.r.t. tau
      !       Parameters: (see top of file)
      !           IALPHA  : selects alpha regularisation
      !           IINTERP : selects interpolation function
      !           IDELFX  : selects gradient expansion corrections


          real(8), intent(in) :: den, grd, tau
          real(8), intent(out) :: eps_x, dedd, dedg, dedt
          integer, intent(in) :: IALPHA, IINTERP, IDELFX

          real(8) :: cp, den83, den53
          real(8) :: exlda, dexldadd, fx, dfxdp, dfxda
          real(8) :: p, dpdd, dpdg
          real(8) :: alpha, dadd, dadg, dadt, reg, dreg
          real(8) :: tueg, dtuegdd, tueg_con
          real(8) :: tauw, dtauwdd, dtauwdg

          real(8), parameter :: AX = -0.7385587663820224058842300326808360d0
          real(8), parameter :: PI = 3.1415926535897932384626433832795d0
          real(8), parameter :: PI2 = PI*PI
          real(8), parameter :: ETA = 1.0d-3
          real(8), parameter :: TAU_R = 1.0d-4
          real(8), parameter :: A_REG = 1.0d-3

      !       Reduced density gradient
          den83 = den**(8.d0/3.d0)
          cp = 4.d0*(3.d0*PI2)**(2.d0/3.d0)
          p = grd**2/(cp*den83)
          dpdd = -8.d0/3.d0*p/den
          dpdg = 2.d0*p/grd

      !       Regularised Alpha
          tueg_con = 3.d0/10.d0*(3.d0*PI2)**(2.d0/3.d0)
          den53 = den**(5.d0/3.d0)
          if (IALPHA .eq. 1) then ! regularised tau_ueg for rSCAN \tilde{alpha}
              tueg = tueg_con*den53 + TAU_R
          else                    ! Unregularised tau_ueg
              tueg = tueg_con*den53
          end if
          dtuegdd = 5.d0/3.d0*tueg_con*den53/den

          tauw = grd**2/8.d0/den
          dtauwdd = -grd**2/(8.d0*den**2)
          dtauwdg = 2.d0*grd/(8.d0*den)

          if (IALPHA .eq. 0) then ! SCAN unregularised alpha
              alpha = (tau - tauw)/tueg
              dadd = -(tau - tauw)*dtuegdd/tueg**2 - dtauwdd/tueg
              dadg = -dtauwdg/tueg
              dadt = 1.d0/tueg

          else if (IALPHA .eq. 1) then ! rSCAN alpha'
              alpha = (tau - tauw)/tueg
              dadd = ((tauw - tau)*dtuegdd - tueg*dtauwdd)/tueg**2
              dadg = -grd/(4.d0*den*tueg)
              dadt = 1.d0/tueg

              if (A_REG .gt. 0.d0) then
                  reg = alpha**3/(alpha**2 + A_REG)
                  dreg = (alpha**4 + 3.d0*alpha**2*A_REG)
     &                 /(alpha**2 + A_REG)**2
                  dadd = dadd*dreg
                  dadg = dadg*dreg
                  dadt = dadt*dreg
                  alpha = reg
              endif

          else if (IALPHA .eq. 2) then ! \bar{alpha} for r2scan and r4scan
              alpha = (tau - tauw)/(tueg + ETA*tauw)
              dadd = -dtauwdd/(tueg + ETA*tauw) 
     &       - (tau - tauw)*(dtuegdd + ETA*dtauwdd)/(tueg + ETA*tauw)**2
              dadg = -ETA*(tau - tauw)*dtauwdg/(tueg + ETA*tauw)**2 
     &       - dtauwdg/(tueg + ETA*tauw)
              dadt = 1.d0/(tueg + ETA*tauw)

          else
              write(*,*) 'ERROR: Unknown IDELFX in SCAN'
              stop
          end if

      !       UEG exchange density [FD]
          exlda = AX*den**(4.d0/3.d0)
          dexldadd = 4.d0*AX*den**(1.d0/3.d0)/3.d0

          call exchange_enhancement(p, alpha, fx, dfxdp, dfxda, ETA, 
     &                              IINTERP, IDELFX)

      !       Exchange Energy[FD]
          eps_x = exlda*fx
          dedd = dexldadd*fx + exlda*(dfxdp*dpdd + dfxda*dadd)
          dedg = exlda*(dfxdp*dpdg + dfxda*dadg)
          dedt = exlda*dfxda*dadt

      end subroutine

      subroutine exchange_enhancement(p, alpha, Fx, dfxdp, dfxda, ETA, 
     &                                IINTERP, IDELFX)
          IMPLICIT NONE

          real(8), intent(in) :: p, alpha
          real(8), intent(out) :: Fx, dfxdp, dfxda
          real(8), intent(in) :: ETA
          integer, intent(in) :: IINTERP, IDELFX


          real(8) :: ief, diefda, oma
          real(8) :: h0x
          real(8) :: h1x, dh1xdp, dh1xda
          real(8) :: gx, dgxdp
          real(8) :: del_f2, C2, ALPHA_GE, damp, ddampdp
          real(8) :: del_fx, ddel_fxdp, ddel_fxda
          real(8) :: wfac, dwfacdp, vfac, dvfacdp, dvfacda
          real(8) :: yfac, dyfacdp, dyfacda

          real(8), parameter :: A1 = 4.9479d0
          real(8), parameter :: K0 = 0.174d0
          real(8), parameter :: MU = 10.d0/81.d0
          real(8), parameter :: K1 = 0.065d0

          real(8), parameter :: cfx1 = 0.667d0
          real(8), parameter :: cfx2 = 0.800d0
          real(8), parameter :: cfdx1 = 1.24d0

          !real(8), parameter :: B1 = 0.156632d0
          !real(8), parameter :: B2 = 0.12083d0
          real(8), parameter :: B3 = 0.5d0
          !real(8), parameter :: B4 = MU*MU/K1 - 0.112654d0

          real(8), parameter :: b2 = (5913.0d0/405000.0d0)**(0.5d0)
          real(8), parameter :: b1 = (511.0d0/13500.0d0)/(2.0d0*b2)
          real(8), parameter :: b4 = mu*mu/K1 - 1606.0d0/18225.0d0 
     &                             -b1*b1


          real(8), parameter, dimension(8) :: 
     &             PARAMS = (/-0.023185843322d0,0.234528941479d0,
     &             -0.887998041597d0,1.451297044490d0,
     &            -0.663086601049d0,-0.4445555d0,-0.667d0,1.d0/)
          integer, dimension(8), parameter :: f_x_e = (/7,6,5,4,
     &                                                  3,2,1,0/)
          real(8), parameter :: D_DAMP2 = 0.361d0

          ALPHA_GE = 20.d0/27.d0 + ETA*5.d0/3.d0

          ief = 0.d0
          diefda = 0.d0
          oma = 1.d0 - alpha
          if (IINTERP .eq. 0) then  ! scan interpolation function
              if (alpha .lt. 1.d0) then
                  ief = exp(-cfx1*alpha/oma)
                  diefda = -cfx1*exp(-cfx1*alpha/oma)/oma**2
              else if (alpha .gt. 1.d0) then
                  ief = -cfdx1*exp(cfx2/oma)
                  diefda = -cfx2*cfdx1*exp(cfx2/oma)/oma**2
              endif

          else if (IINTERP .eq. 1) then ! rscan, r2scan, r4scan interpolation function
              if (alpha .lt. 1.0d-13) then
                  ief = exp(-cfx1*alpha/oma)
                  diefda = -cfx1*exp(-cfx1*alpha/oma)/oma**2
              else if( alpha .lt. 2.5d0) then
                  ief = dot_product(alpha**f_x_e, PARAMS)
                  diefda = dot_product(alpha**f_x_e(2:), 
     &                                 f_x_e(:7)*PARAMS(:7))
              else if (alpha .ge. 2.5d0) then
                  ief = -cfdx1*exp(cfx2/oma)
                  diefda = -cfx2*cfdx1*exp(cfx2/oma)/oma**2
              endif

          else
              write(*,*) 'ERROR: Unknown IINTERP in SCAN'
              stop
          end if

      !       Single orbital enhancement
          h0x = 1.d0 + K0

      !       Slowly varying enhancement
          if (IDELFX .eq. 0) then     ! scan, rscan
              wfac = B4*p**2*exp(-B4*p/MU)
              dwfacdp = B4*p*exp(-B4*p/MU)*(2.d0 - B4*p/MU)

              vfac = B1*p + B2*oma*exp(-B3*oma**2)
              yfac = MU*p + wfac + vfac**2
              h1x = 1.d0 + K1 - K1/(1.d0 + yfac/K1)

              dvfacdp = B1
              dyfacdp = MU + dwfacdp + 2.d0*vfac*dvfacdp
              dh1xdp = dyfacdp/(1.d0 + yfac/K1)**2

              dvfacda = -B2*(1.d0 - 2.d0*B3*oma**2)*exp(-B3*oma**2)
              dyfacda = 2.d0*vfac*dvfacda
              dh1xda = dyfacda/(1.d0 + yfac/K1)**2

          else if (IDELFX .eq. 1 .or. IDELFX .eq. 2) then  ! Second order corrections
              del_f2 = dot_product(f_x_e(:7), PARAMS(:7))
              C2 = -del_f2*(1.d0 - h0x)

      !       Damping
              damp = exp(-p**2/D_DAMP2**4)
              ddampdp = -2.d0*damp*p/D_DAMP2**4

      !       Slowly varying contribution [FD]
              h1x = 1.d0 + K1 - K1/(1.d0 + p*(MU + ALPHA_GE*C2*damp)/K1)
              dh1xdp = K1**2*(MU + ALPHA_GE*C2*(damp + p*ddampdp)) 
     &               /(K1 + MU*p + ALPHA_GE*C2*p*damp)**2
              dh1xda = 0.d0

          else
              write(*,*) 'ERROR: Unknown IDELFX in SCAN'
              stop
          end if

          gx = 1.d0 - exp(-A1/p**(1.d0/4.d0))
          dgxdp = -A1*exp(-A1/p**(1.d0/4.d0))/(4.d0*p**(5.d0/4.d0))

          if (IDELFX .eq. 2) then  ! 4th order corrections for r4scan
              call get_del_fx(p, alpha, del_fx, ddel_fxdp, ddel_fxda, 
     &                        K0, K1, C2, ETA, ALPHA_GE, 
     &                        PARAMS, f_x_e)
          else
              del_fx = 0.d0
              ddel_fxdp = 0.d0
              ddel_fxda = 0.d0
          end if

          fx = (h1x + ief*(h0x - h1x) + del_fx)*gx
          dfxdp = (del_fx + h1x + (h0x - h1x)*ief)*dgxdp 
     &       + gx*(ddel_fxdp + dh1xdp - ief*dh1xdp)
          dfxda = gx*((h0x - h1x)*diefda + ddel_fxda + dh1xda 
     &          - ief*dh1xda)

      end subroutine

      subroutine get_del_fx(p, alpha, del_fx, ddel_fxdp, ddel_fxda, 
     &                        K0, K1, C2, ETA, ALPHA_GE, PARAMS, f_x_e)
          IMPLICIT NONE

          real(8), intent(in) :: p, alpha, K0, K1, C2, ETA, ALPHA_GE
          real(8), dimension(8), intent(in) :: PARAMS
          integer, dimension(8), intent(in) :: f_x_e
          real(8), intent(out) :: del_fx, ddel_fxdp, ddel_fxda
          real(8) :: order_1, dorder_1dp, dorder_1da, C_pa, C_aa, C_pp
          real(8) :: oma, damp, t1, dt1dp, dt1da, ddampdp, ddampda

          call get_dx_terms(C_aa, C_pa, C_pp, 
     &              K0, K1, C2, ETA, PARAMS, f_x_e)

          oma = 1.d0 - alpha

          order_1 = C2*(oma - ALPHA_GE*p)
          dorder_1dp = -C2*ALPHA_GE
          dorder_1da = -C2

      !       Correcting contribution
          t1 = order_1 + C_aa*oma**2 + C_pa*p*oma + C_pp*p**2
          dt1dp = dorder_1dp + 2*C_pp*p + C_pa*oma
          dt1da = dorder_1da - 2*C_aa*oma - C_pa*p

          call get_fourth_order_damp(p, alpha, damp, ddampdp, ddampda)

          del_fx = t1*damp
          ddel_fxdp = damp*dt1dp + t1*ddampdp
          ddel_fxda = t1*ddampda + dt1da*damp

      end subroutine get_del_fx

      subroutine get_dx_terms(C_aa, C_pa, C_pp, K0, K1, C2, ETA, 
     &                  PARAMS, f_x_e)
          IMPLICIT NONE

          real(8), intent(in) :: K0, K1, C2, ETA
          real(8), dimension(8), intent(in) :: PARAMS
          integer, dimension(8), intent(in) :: f_x_e
          real(8), intent(out) :: C_aa, C_pa, C_pp
          real(8) :: h0x, eta_term, del_f2, del_f4
          real(8) :: ALPHA_GE

          real(8), parameter :: MU = 10.d0/81.d0

          ALPHA_GE = 20.d0/27.d0 + ETA*5.d0/3.d0

          eta_term = ETA*3.d0/4.d0 + 2.d0/3.d0
          h0x = 1.d0 + K0

          del_f2 = dot_product(f_x_e(:7), PARAMS(:7))
          del_f4 = dot_product(f_x_e(:7)*(f_x_e(:7) - 1.d0), PARAMS(:7))

          C_aa = 73.d0/5000.d0 - 0.5d0*del_f4*(h0x - 1.d0)

          C_pa = 511.d0/13500.d0 - 73.d0/1500.d0*ETA 
     &            - del_f2*(ALPHA_GE*C2 + MU)

          C_pp = 146.d0/2025.d0*eta_term**2 - 73.d0/405.d0*eta_term 
     &            + (ALPHA_GE*C2 + MU)**2/K1

      end subroutine

      subroutine get_fourth_order_damp(p, alpha, damp4, 
     &            ddamp4dp, ddamp4da)
          implicit NONE

          real(8), intent(in) :: p, alpha
          real(8), intent(out) :: damp4, ddamp4dp, ddamp4da
          real(8) :: t1, t2, dt1da, dt2dp, dt2da, oma
          real(8), parameter :: DX_DAMP4_P = 0.802d0
          real(8), parameter :: DX_DAMP4_A = 0.178d0

          oma = 1.d0 - alpha

          t1 = 2.d0*alpha**2/(1.d0 + alpha**4)
          dt1da = -4.d0*alpha*(alpha**4 - 1.d0)/(1.d0 + alpha**4)**2

          t2 = exp(-oma**2/DX_DAMP4_A**2 - p**2/DX_DAMP4_P**4)
          dt2dp = -2.d0*t2*p/DX_DAMP4_P**4
          dt2da = 2*oma*t2/DX_DAMP4_A**2

          damp4 = t1*t2
          ddamp4dp = t1*dt2dp
          ddamp4da = t1*dt2da + t2*dt1da

      end subroutine

!================================================================
!   SCAN family correlation
!================================================================

      subroutine vrSCANc(d0, d1, g0, g1, gt, t0, t1, eps_c, 
     &                   dedd0, dedd1, dedg0, dedg1, dedg, dedt0, 
     &                   dedt1, IALPHA, IINTERP, IDELEC)
          IMPLICIT NONE

      !       Inputs:
      !           d0 : electron density for spin 0
      !           d1 : electron density for spin 1
      !           g0 : |grad d0|. Not this is not squared!
      !           g1 : |grad d1|. Not this is not squared!
      !           gt : |grad d|. Total density gradient. Note: not g0 + g1
      !           t0 : Kinetic energy density spin 0
      !           t1 : Kinetic energy density spin 1
      !       Outputs:
      !           eps_c : energy density
      !           dedd0  : energy density derivative w.r.t. d0
      !           dedd1  : energy density derivative w.r.t. d1
      !           dedg0  : energy density derivative w.r.t. g0
      !           dedg1  : energy density derivative w.r.t. g1
      !           dedg   : energy density derivative w.r.t. gt
      !           dedt0  : energy density derivative w.r.t. t0
      !           dedt1  : energy density derivative w.r.t. t1
      !       Parameters: (see top of file)
      !           IALPHA  : selects alpha regularisation
      !           IINTERP : selects interpolation function
      !           IDELEC  : selects gradient expansion corrections

          real(8), intent(in) :: d0, d1, g0, g1, gt, t0, t1
          real(8), intent(out) :: eps_c, dedd0, dedg0, dedt0, dedd1, 
     &                            dedg1
          real(8), intent(out) :: dedt1, dedg
          integer, intent(in) :: IALPHA, IINTERP, IDELEC

          real(8) :: dt, tt, dthrd, rs, drsdd0, drsdd1
          real(8) :: zeta, dzetadd0, dzetadd1, ds_z, dds_zdd0, dds_zdd1
          real(8) :: s, dsdd0, dsdg0, dsdd1, dsdg1, dsdg, gngn0, gngn1
          real(8) :: y0, y1, y, yc, dycdg0, dycdg1
          real(8) :: alpha, dadd0, dadg0, dadt0, dadd1, dadg1, dadt1, 
     &               oma
          real(8) :: reg, dreg, dadg
          real(8) :: tueg, tueg_con, dtuegdd0, dtuegdd1
          real(8) :: tauw, dtauwdd0, dtauwdg0, dtauwdd1, dtauwdg1, 
     &               dtauwdg
          real(8) :: ief, diefda, diefdd0, diefdd1, diefdg0, diefdg1, 
     &               diefdg
          real(8) :: diefdt0, diefdt1
          real(8) :: ec0, dec0drs, dec0ds, dec0dz
          real(8) :: dec0dd0, dec0dd1, dec0dg0, dec0dg1
          real(8) :: ec1, dec1drs, dec1ds, dec1dz!, dec1da
          real(8) :: dec1dd0, dec1dd1, dec1dg0, dec1dg1
          real(8) :: dec0dg, dec1dg, s_den

          real(8), parameter :: CFDC1 = 0.7d0
          real(8), parameter :: CFC1 = 0.640d0
          real(8), parameter :: CFC2 = 1.5d0

          real(8), parameter, dimension(7) :: IE_PARAMS_C = 
     &    (/-0.64d0, -0.4352d0, -1.535685604549d0, 3.061560252175d0, 
     &        -1.915710236206d0, 0.516884468372d0, -0.051848879792d0/)

          real(8), parameter :: ETA = 1.0d-3
          real(8), parameter :: TAU_R = 1.0d-4
          real(8), parameter :: A_REG = 1.0d-3

          integer :: i

          real(8), parameter :: BETA_MB = 0.06672455060314922d0!0.066725d0
          real(8), parameter :: GAMMA = 0.031090690869655d0
          real(8), parameter :: AFACTOR = 0.1d0
          real(8), parameter :: BFACTOR = 0.1778d0
          real(8), parameter :: PI = 3.1415926535897932384626433832795d0
          real(8), parameter :: PI2 = PI*PI

          dt = d0 + d1
          y0 = g0**2
          y1 = g1**2
          y = gt**2
          yc = (y - y0 - y1)/2.d0
          tt = t0 + t1

          dycdg0 = yc/g0
          dycdg1 = yc/g1

          !        Zeta
          !zeta = min(max((d0 - d1)/dt, -0.99999999999990d0), 0.99999999999990d0)
          zeta = min(max((d0 - d1)/dt, -1.d0), 1.d0)

          dzetadd0 = 2.d0*d1/dt**2
          dzetadd1 = -2.d0*d0/dt**2

          !        Wigner-Seitz radius
          dthrd = dt**(1.d0/3.d0)
          rs = (0.75d0/PI)**(1.d0/3.d0)/dthrd
          drsdd0 = -(0.75d0/PI)**(1.d0/3.d0)/(3.d0*dthrd*dt)
          drsdd1 = drsdd0

          !        Reduced density gradient
          s_den = 1.d0/(2.d0*(3.d0*PI2)**(1.d0/3.d0)*dt**(4.d0/3.d0))
          s = gt*s_den
          dsdd0 = -4.d0/3.d0*s/dt
          dsdd1 = dsdd0
          gngn0 = 1.d0/(2.d0*gt)*(2.d0*g0+2.d0*yc/g0)
          gngn1 = 1.d0/(2.d0*gt)*(2.d0*g1+2.d0*yc/g1)
          dsdg0 = s/gt*gngn0
          dsdg1 = s/gt*gngn1
          dsdg = s_den

          !        ds_zeta
          ds_z = ((1.d0 + zeta)**(5.d0/3.d0) 
     &         + (1.d0 - zeta)**(5.d0/3.d0))/2.d0
          dds_zdd0 = 5.d0/3.d0*((1.d0 + zeta)**(2.d0/3.d0) 
     &             - (1.d0 - zeta)**(2.d0/3.d0))*dzetadd0/2.d0
          dds_zdd1 = 5.d0/3.d0*((1.d0 + zeta)**(2.d0/3.d0) 
     &             - (1.d0 - zeta)**(2.d0/3.d0))*dzetadd1/2.d0

          !       alpha
          tueg_con = 3.d0/10.d0*(3.d0*PI2)**(2.d0/3.d0)

          if (IALPHA .eq. 1) then
              tueg = (tueg_con*dt**(5.d0/3.d0) + TAU_R)*ds_z
              dtuegdd0 = 5.d0/3.d0*tueg_con*dt**(2.d0/3.d0)*ds_z 
     &           + tueg*dds_zdd0/ds_z
              dtuegdd1 = 5.d0/3.d0*tueg_con*dt**(2.d0/3.d0)*ds_z 
     &           + tueg*dds_zdd1/ds_z
          else
              tueg = tueg_con*dt**(5.d0/3.d0)*ds_z
              dtuegdd0 = 5.d0/3.d0*tueg/dt + tueg*dds_zdd0/ds_z
              dtuegdd1 = 5.d0/3.d0*tueg/dt + tueg*dds_zdd1/ds_z
          end if

          !tauw = min(y/(8.d0*dt), tt)
          tauw = y/(8.d0*dt)
          dtauwdd0 = -tauw/dt
          dtauwdd1 = dtauwdd0
          dtauwdg0 = (g0 + dycdg0)/(4.d0*dt)
          dtauwdg1 = (g1 + dycdg1)/(4.d0*dt)
          dtauwdg = gt/(4.d0*dt)

          if (IALPHA .eq. 0 .or. IALPHA .eq. 1) then
              ! Conventional alpha
              alpha = (tt - tauw)/tueg
              dadd0 = -(dtauwdd0 + alpha*dtuegdd0)/tueg
              dadd1 = -(dtauwdd1 + alpha*dtuegdd1)/tueg
              dadg0 = -dtauwdg0/tueg
              dadg1 = -dtauwdg1/tueg
              dadg = -dtauwdg/tueg
              dadt0 = 1.d0/tueg
              dadt1 = 1.d0/tueg

          else if (IALPHA .eq. 2) then
              ! regularised alpha of r++scan, r2scan, r4scan
              alpha = (tt - tauw)/(tueg + ETA*tauw)
              dadd0 = -dtauwdd0/(tueg + ETA*tauw) 
     &    - (tt - tauw)*(dtuegdd0 + ETA*dtauwdd0)/(tueg + ETA*tauw)**2
              dadd1 = -dtauwdd1/(tueg + ETA*tauw) 
     &    - (tt - tauw)*(dtuegdd1 + ETA*dtauwdd1)/(tueg + ETA*tauw)**2
              dadg0 = -ETA*(tt - tauw)*dtauwdg0/(tueg + ETA*tauw)**2 
     &    - dtauwdg0/(tueg + ETA*tauw)
              dadg1 = -ETA*(tt - tauw)*dtauwdg1/(tueg + ETA*tauw)**2 
     &    - dtauwdg1/(tueg + ETA*tauw)
              dadg = -ETA*(tt - tauw)*dtauwdg/(tueg + ETA*tauw)**2 
     &    - dtauwdg/(tueg + ETA*tauw)
              dadt0 = 1.d0/(tueg + ETA*tauw)
              dadt1 = dadt0
          else
              write(*,*) 'ERROR: Unknown IALPHA in SCAN'
              stop
          end if

          if (IALPHA .eq. 1) then
              ! alpha regularisation for rscan
              reg = alpha**3/(alpha**2 + A_REG)
              dreg = (alpha**4 + 3.d0*alpha**2*A_REG)
     &             /(alpha**2 + A_REG)**2
              dadd0 = dadd0*dreg
              dadd1 = dadd1*dreg
              dadg0 = dadg0*dreg
              dadg1 = dadg1*dreg
              dadg = dadg*dreg
              dadt0 = dadt0*dreg
              dadt1 = dadt1*dreg
              alpha = reg
          end if

          ief = 0.d0
          diefda = 0.d0
          oma = 1.d0 - alpha

          if (IINTERP .eq. 0) then
              ! SCAN interpolation function
              if (alpha .lt. 1.d0) then
                  ief = exp(-CFC1*alpha/oma)
                  diefda = -CFC1*ief/oma**2
              else if (alpha .gt. 1.0d0) then
                  ief = -CFDC1*exp(CFC2/oma)
                  diefda = CFC2*ief/oma**2
              endif

          else if (IINTERP .eq. 1) then
              ! rSCAN interpolation function
              if (alpha .lt. 1.0d-13) then
                  ief = exp(-CFC1*alpha/oma)
                  diefda = -CFC1*ief/oma**2
              else if (alpha .lt. 2.5d0) then
                  ief = 1.d0
                  do i = 1, 7
                      ief = ief + IE_PARAMS_C(i)*alpha**(i)
                      diefda = diefda + i*IE_PARAMS_C(i)*alpha**(i - 1)
                  end do
              else if (alpha .ge. 2.5d0) then
                  ief = -CFDC1*exp(CFC2/oma)
                  diefda = CFC2*ief/oma**2
              endif
          else
              write(*,*) 'ERROR: Unknown IINTERP in SCAN'
              stop
          end if

          diefdd0 = diefda*dadd0
          diefdd1 = diefda*dadd1
          diefdg0 = diefda*dadg0
          diefdg1 = diefda*dadg1
          diefdg = diefda*dadg
          diefdt0 = diefda*dadt0
          diefdt1 = diefda*dadt1

          !        Single Orbital Correlation
          ec0 = 0.d0
          dec0drs = 0.d0
          dec0ds = 0.d0
          dec0dz = 0.d0
          call scan_ec0(rs, s, zeta, ec0, dec0drs, dec0ds, dec0dz)
          dec0dd0 = dec0drs*drsdd0 + dec0dz*dzetadd0 + dec0ds*dsdd0
          dec0dd1 = dec0drs*drsdd1 + dec0dz*dzetadd1 + dec0ds*dsdd1
          dec0dg0 = dec0ds*dsdg0
          dec0dg1 = dec0ds*dsdg1
          dec0dg = dec0ds*dsdg

          !        Slowly Varying Correlation
          ec1 = 0.d0
          dec1drs = 0.d0
          dec1ds = 0.d0
          dec1dz = 0.d0
          call scan_ec1(rs, s, zeta, ec1, dec1drs, dec1dz, dec1ds, 
     &                  IE_PARAMS_C, ETA, IDELEC)
          dec1dd0 = dec1drs*drsdd0 + dec1dz*dzetadd0 + dec1ds*dsdd0
          dec1dd1 = dec1drs*drsdd1 + dec1dz*dzetadd1 + dec1ds*dsdd1
          dec1dg0 = dec1ds*dsdg0 !+ dec1da*dadg0
          dec1dg1 = dec1ds*dsdg1 !+ dec1da*dadg1
          dec1dg = dec1ds*dsdg

          !        Full correlation functional
          eps_c = (ec1 + ief*(ec0 - ec1))*dt
          dedd0 = ec1 + (ec0 - ec1)*ief 
     &    + dt*(ief*(dec0dd0 - dec1dd0) + dec1dd0 + (ec0 - ec1)*diefdd0)
          dedd1 = ec1 + (ec0 - ec1)*ief 
     &    + dt*(ief*(dec0dd1 - dec1dd1) + dec1dd1 + (ec0 - ec1)*diefdd1)
          dedg0 = dt* 
     &      (ief*(dec0dg0 - dec1dg0) + dec1dg0 + (ec0 - ec1)*diefdg0)
          dedg1 = dt* 
     &      (ief*(dec0dg1 - dec1dg1) + dec1dg1 + (ec0 - ec1)*diefdg1)
          dedg = dt* 
     &      (ief*(dec0dg - dec1dg) + dec1dg + (ec0 - ec1)*diefdg)
          dedt0 = dt*(ec0 - ec1)*diefdt0
          dedt1 = dt*(ec0 - ec1)*diefdt1

      end subroutine

      subroutine scan_ec0(rs, s, zeta, ec0, dec0drs, dec0ds, dec0dz)
          IMPLICIT NONE

          real(8), intent(in) :: rs, s, zeta
          real(8), intent(out) :: ec0, dec0drs, dec0ds, dec0dz

          real(8) :: eclda, dldadrs, dldadrsrs
          real(8) :: dx_z, ddx_zdz, gc_z, dgc_zdz
          real(8) :: w0, dw0drs, ginf, dginfds
          real(8) :: h0, dh0drs, dh0ds

          real(8), parameter :: B1C = 0.0285764d0
          real(8), parameter :: CHI_LD = 0.12802585262625815d0

          call lda_0(rs, eclda, dldadrs, dldadrsrs, B1C)

          dx_z = ((1.d0 + zeta)**(4.d0/3.d0) 
     &         + (1.d0 - zeta)**(4.d0/3.d0))/2.d0
          ddx_zdz = -2.d0*((1.d0 - zeta)**(1.d0/3.d0) 
     &            - (1.d0 + zeta)**(1.d0/3.d0))/3.d0

          gc_z = (1.d0 - 2.3631d0*(dx_z - 1.d0))*(1.d0 - zeta**12)
          dgc_zdz = -(1.d0 - 2.3631d0*(dx_z - 1.d0))*12.d0*zeta**11
          dgc_zdz = dgc_zdz - 2.3631d0*ddx_zdz*(1.d0 - zeta**12)

          w0 = exp(-eclda/B1C) - 1.d0
          dw0drs = -(w0 + 1.d0)*dldadrs/B1C

          ginf = 1.d0/(1.d0 + 4.d0*CHI_LD*s*s)**(1.d0/4.d0)
          dginfds = -2.d0*CHI_LD*s/(1.d0 + 4.d0*CHI_LD*s*s)**(5.d0/4.d0)

          h0 = B1C*log(1.d0 + w0*(1.d0 - ginf))
          dh0drs = B1C*(1.d0 - ginf)*dw0drs/(1.d0 + (1.d0 - ginf)*w0)
          dh0ds = -B1C*w0*dginfds/(1.d0 + (1.d0 - ginf)*w0)

          ec0 = (eclda + h0)*gc_z
          dec0drs = (dldadrs + dh0drs)*gc_z
          dec0dz = (h0 + eclda)*dgc_zdz
          dec0ds = dh0ds*gc_z

      end subroutine

      subroutine scan_ec1(rs, s, zeta, ec1, dec1drs, dec1dz, dec1ds, 
     &                    IE_PARAMS_C, ETA, IDELEC)
          IMPLICIT NONE

          real(8), intent(in) :: rs, s, zeta
          real(8), intent(out) :: ec1, dec1drs, dec1dz, dec1ds
          real(8), dimension(7), intent(in) :: IE_PARAMS_C
          real(8), intent(in) :: ETA
          integer, intent(in) :: IDELEC

          real(8) :: sqrt_rs, dx_z, ddx_zdz, gc_z, dgc_zdz, phi, dphidz
          real(8) :: phi3, dphi3dz
          real(8) :: eclda0, declda0drs, declda0drsrs
          real(8) :: eclsda1, declsda1drs, declsda1dz, declsda1drsrs, 
     &               declsda1drsz
          real(8) :: t, dtdrs, dtdz, dtds
          real(8) :: w1, dw1drs, dw1dz, y, dydrs, dydz, dyds
          real(8) :: del_y, ddel_ydrs, ddel_ydz, ddel_yds
          real(8) :: g_y, dg_ydrs, dg_ydz, dg_yds
          real(8) :: h1, dh1drs, dh1dz, dh1ds

          real(8), parameter :: B1C = 0.0285764d0
          real(8), parameter :: PI = 3.1415926535897932384626433832795d0
          real(8), parameter :: GAMMA = 0.0310906908696d0
          real(8), parameter :: BETA_MB = 0.06672455060314922d0!0.066724550603149220d0
          real(8), parameter :: AFACTOR = 0.1d0
          real(8), parameter :: BFACTOR = 0.1778d0
          real(8) :: AFIX_T
          AFIX_T = sqrt(PI/4.d0)*(9.d0*PI/4.d0)**(1.d0/6.d0)

          dx_z = ((1.d0 + zeta)**(4.d0/3.d0) 
     &         + (1.d0 - zeta)**(4.d0/3.d0))/2.d0
          ddx_zdz = -2.d0*((1.d0 - zeta)**(1.d0/3.d0) 
     &            - (1.d0 + zeta)**(1.d0/3.d0))/3.d0

          gc_z = (1.d0 - 2.3631d0*(dx_z - 1.d0))*(1.d0 - zeta**12)
          dgc_zdz = -(1.d0 - 2.3631d0*(dx_z - 1.d0))*12.d0*zeta**11
          dgc_zdz = dgc_zdz - 2.3631d0*ddx_zdz*(1.d0 - zeta**12)

          phi = (exp((2.d0/3.d0)*log(1.d0 + zeta)) 
     &        + exp((2.d0/3.d0)*log(1.d0 - zeta)))/2.d0
          !dphidz = (1.d0/3.d0)*((1.d0 + zeta)**(-1.d0/3.d0)-(1.d0 - zeta)**(-1.d0/3.d0))

          ! AK: this computation of d phi / d z is wrong, but this is what's done
          ! in the SCAN subroutines.
          ! Basically the problem is that phi'(z) diverges when |z| to 1.
          ! In VASP, we regularize |zeta| <= 0.99999999999990
          ! FHI-AIMS and NWChem use this alternative regularization
          if (1.d0 - zeta < 1.d-18) then
              dphidz = (1.d0/3.d0)*(1.d0 + zeta)**(-1.d0/3.d0)
          elseif (1.d0 + zeta < 1.d-18) then
              dphidz = - (1.d0/3.d0)*(1.d0 - zeta)**(-1.d0/3.d0)
          else
              dphidz = (1.d0/3.d0)*((1.d0 + zeta)**(-1.d0/3.d0)
     &               -(1.d0 - zeta)**(-1.d0/3.d0))
          end if

          phi3 = phi**3
          dphi3dz = 3.d0*phi**2*dphidz

          call lda_0(rs, eclda0, declda0drs, declda0drsrs, B1C)
          call lsda_1(rs, zeta, eclsda1, declsda1drs, declsda1dz, 
     &                declsda1drsrs, declsda1drsz)

          sqrt_rs = sqrt(rs)

          t = AFIX_T*s/(sqrt_rs*phi)
          dtdrs = -AFIX_T*s/(2.d0*phi*rs**(3.0/2.0))
          dtdz = -dphidz*AFIX_T*s/(sqrt_rs*phi**2)
          dtds = AFIX_T/(sqrt_rs*phi)

          w1 = exp(-eclsda1/(GAMMA*phi3)) - 1.d0
          dw1drs = -(w1 + 1.d0)*declsda1drs/(GAMMA*phi3)
          dw1dz = -(w1 + 1.d0)/(GAMMA*phi3)*(declsda1dz 
     &          - 3.d0*eclsda1*dphidz/phi)

          call get_y(rs, t, dtdrs, dtdz, dtds, 
     &    w1, dw1drs, dw1dz, GAMMA, 
     &    y, dydrs, dydz, dyds)

          if (IDELEC .eq. 0) then
              del_y = 0.d0
              ddel_ydrs = 0.d0
              ddel_ydz = 0.d0
              ddel_yds = 0.d0

          else if (IDELEC .eq. 1 .or. IDELEC .eq. 2) then
              call get_del_y(rs, s, zeta, 
     &                       eclda0, declda0drs, declda0drsrs, 
     &                       eclsda1, declsda1drs, declsda1dz, 
     &                       declsda1drsrs, declsda1drsz, 
     &                       gc_z, dgc_zdz, 
     &                       phi3, dphi3dz, w1, dw1drs, dw1dz, GAMMA, 
     &                       del_y, ddel_ydrs, ddel_ydz, ddel_yds, 
     &                       IE_PARAMS_C, ETA)

          else
              write(*,*) 'ERROR: Unknown IDELEC in SCAN'
              stop
          end if

          g_y = 1.d0/(1.d0 + 4.d0*(y - del_y))**(1.0/4.0)
          dg_ydrs = -(dydrs - ddel_ydrs) 
     &              /(1.d0 + 4.d0*(y - del_y))**(5.d0/4.d0)
          dg_ydz = -(dydz - ddel_ydz) 
     &              /(1.d0 + 4.d0*(y - del_y))**(5.d0/4.d0)
          dg_yds = -(dyds - ddel_yds) 
     &              /(1.d0 + 4.d0*(y - del_y))**(5.d0/4.d0)

          h1 = GAMMA*phi3*log(1.d0 + w1*(1.d0 - g_y))
          dh1drs = GAMMA*phi3*((1.d0 - g_y)*dw1drs - w1*dg_ydrs) 
     &          /(1.d0 + (1.d0 - g_y)*w1)
          dh1dz = GAMMA*log(1.d0 + (1.d0 - g_y)*w1)*dphi3dz 
     &        + GAMMA*phi3*((1.d0 - g_y)*dw1dz - w1*dg_ydz) 
     &          /(1.d0 + (1.d0 - g_y)*w1)
          dh1ds = -GAMMA*phi3*w1*dg_yds/(1.d0 + (1.d0 - g_y)*w1)

          ec1 = eclsda1 + h1
          dec1drs = declsda1drs + dh1drs
          dec1dz = declsda1dz + dh1dz
          dec1ds = dh1ds
      end subroutine

      subroutine get_y(rs, t, dtdrs, dtdz, dtds, 
     &                 w1, dw1drs, dw1dz, GAMMA, 
     &                 y, dydrs, dydz, dyds)
          implicit NONE

          real(8), intent(in) :: rs, t, dtdrs, dtdz, dtds
          real(8), intent(in) :: w1, dw1drs, dw1dz, GAMMA
          real(8), intent(out) :: y, dydrs, dydz, dyds
          real(8) :: beta, dbetadrs

          real(8), parameter :: BETA_MB = 0.06672455060314922d0!0.066725d0
          real(8), parameter :: AFACTOR = 0.1d0
          real(8), parameter :: BFACTOR = 0.1778d0

          beta = BETA_MB*(1.d0 + AFACTOR*rs)/(1.d0 + BFACTOR*rs)
          dbetadrs = BETA_MB*(AFACTOR - BFACTOR)/(1.d0 + BFACTOR*rs)**2

          y = beta/(GAMMA*w1)*t**2
          dydrs = t**2*dbetadrs/(GAMMA*w1) 
     &          - beta*t**2*dw1drs/(GAMMA*w1**2) 
     &          + 2.d0*beta*t*dtdrs/(GAMMA*w1)
          dydz = 2.d0*beta*t*dtdz/(GAMMA*w1) 
     &         - beta*t**2*dw1dz/(GAMMA*w1**2)
          dyds = 2.d0*beta*t*dtds/(GAMMA*w1)
      end subroutine

      subroutine get_del_y(rs, s, zeta, lda0, dlda0drs, dlda0drsrs, 
     &                     lsda1, dlsda1drs, dlsda1dz, dlsda1drsrs, 
     &                     dlsda1drsz, gc, dgcdz, 
     &                     phi3, dphi3dz, w1, dw1drs, dw1dz, GAMMA, 
     &                     del_y, ddel_ydrs, ddel_ydz, ddel_yds, 
     &                     IE_PARAMS_C, ETA)
          IMPLICIT NONE

          real(8), intent(in) :: rs, s, zeta, lda0, dlda0drs, dlda0drsrs
          real(8), intent(in) :: lsda1, dlsda1drs, dlsda1dz, 
     &                           dlsda1drsrs, dlsda1drsz
          real(8), intent(in) :: gc, dgcdz
          real(8), intent(in) :: phi3, dphi3dz, w1, dw1drs, dw1dz, GAMMA
          real(8), intent(out) :: del_y, ddel_ydrs, ddel_ydz, ddel_yds
          real(8), dimension(7), intent(in) :: IE_PARAMS_C
          real(8), intent(in) :: ETA

          real(8) :: del_f2
          real(8) :: lsda0, dlsda0drs, dlsda0dz, dlsda0drsrs, dlsda0drsz
          real(8) :: K, dKdrs, dKdz
          real(8) :: t1, dt1drs, dt1dz, t2, dt2drs, dt2dz, t3, dt3drs, 
     &               dt3dz
          real(8) :: damp, ddampds
          real(8) :: p, ds_z, dds_zdz
          integer :: i

          real(8), parameter :: D_DAMP2 = 0.361

          p = s*s
          ds_z = ((1.d0 + zeta)**(5.d0/3.d0) 
     &         + (1.d0 - zeta)**(5.d0/3.d0))/2.d0
          dds_zdz = -5.d0/6.d0*((1.d0 - zeta)**(2.d0/3.d0) 
     &                         - (1.d0 + zeta)**(2.d0/3.d0))

          lsda0 = lda0*gc
          dlsda0drs = dlda0drs*gc
          dlsda0dz = lda0*dgcdz
          dlsda0drsrs = dlda0drsrs*gc
          dlsda0drsz = dlda0drs*dgcdz

          del_f2 = 0.d0
          do i = 1, 7
              del_f2 = del_f2 + i*IE_PARAMS_C(i)
          end do

          t1 = del_f2/(27.d0*GAMMA*ds_z*phi3*w1)
          dt1drs = -t1*dw1drs/w1
          dt1dz = -del_f2*(w1*(phi3*dds_zdz + ds_z*dphi3dz) 
     &          + ds_z*phi3*dw1dz)/(27.d0*GAMMA*(ds_z*phi3*w1)**2)

          t2 = 20.d0*rs*(dlsda0drs - dlsda1drs)
          dt2drs = 20.d0*(dlsda0drs - dlsda1drs 
     &           + rs*(dlsda0drsrs - dlsda1drsrs))
          dt2dz = 20.d0*rs*(dlsda0drsz - dlsda1drsz)

          t3 = 45.d0*ETA*(lsda0 - lsda1)
          dt3drs = 45.d0*ETA*(dlsda0drs - dlsda1drs)
          dt3dz = 45.d0*ETA*(dlsda0dz - dlsda1dz)

          K = t1*(t2 - t3)
          dKdrs = dt1drs*(t2 - t3) + t1*(dt2drs - dt3drs)
          dKdz = dt1dz*(t2 - t3) + t1*(dt2dz - dt3dz)

          damp = exp(-p**2/D_DAMP2**4)
          ddampds = -4.d0*damp*s**3/D_DAMP2**4

          del_y = K*p*damp
          ddel_ydrs = p*damp*dKdrs
          ddel_ydz = p*damp*dKdz
          ddel_yds = K*s*(2.d0*damp + s*ddampds)
      end subroutine get_del_y

      subroutine lda_0(rs, elda_0, dldadrs, dldadrsrs, B1C)
          IMPLICIT NONE

          real(8), intent(in) :: rs
          real(8), intent(out) :: elda_0, dldadrs, dldadrsrs
          real(8), intent(in) :: B1C

          real(8) :: sqrtrs

          real(8), parameter :: B2C = 0.0889d0
          real(8), parameter :: B3C = 0.125541d0

          sqrtrs = sqrt(rs)
          elda_0 = -B1C/(1.d0 + B2C*sqrtrs + B3C*rs)
          dldadrs = (B3C + B2C/(2.d0*sqrtrs))*elda_0**2/B1C
          dldadrsrs = -B1C*(B2C + 3.d0*B2C**2*sqrtrs + 9.d0*B2C*B3C*rs 
     &              + 8.d0*B3C**2*rs*sqrtrs) 
     &              /(4.d0*rs*sqrtrs*(1.d0 + B2C*sqrtrs + B3C*rs)**3)
      end subroutine

      subroutine lsda_1(rs, zeta, eclda1, declda1drs, declda1dz, 
     &                  declda1drsrs, declda1drsz)
          IMPLICIT NONE

          real(8), intent(in) :: rs, zeta
          real(8), intent(out) :: eclda1, declda1drs, declda1dz
          real(8), intent(out) :: declda1drsrs, declda1drsz
          real(8) :: z3, z4
          real(8) :: eu, deudrs, deudrsrs
          real(8) :: ep, depdrs, depdrsrs
          real(8) :: alfm, dalfmdrs, dalfmdrsrs
          real(8) :: F, dFdz

          real(8), parameter :: GAM = 0.51984209978974632953442121455650d0
          real(8), parameter :: FZZ = 8.d0/(9.d0*GAM)

          call grcor2(0.03109070d0, 0.213700d0, 7.59570d0, 3.58760d0, 
     &          1.63820d0, 0.492940d0, rs, eu, deudrs, deudrsrs)
          call grcor2(0.015545350d0, 0.205480d0, 14.11890d0, 6.19770d0,
     &          3.36620d0, 0.625170d0, rs, ep, depdrs, depdrsrs)
          call grcor2(0.01688690d0, 0.111250d0, 10.3570d0, 3.62310d0, 
     &          0.880260d0, 0.496710d0, rs, alfm, dalfmdrs, dalfmdrsrs)

          z3 = zeta**3
          z4 = zeta*z3

          F = ((1.d0 + zeta)**(4.d0/3.d0) + (1.d0 - zeta)**(4.d0/3.d0) 
     &      - 2.d0)/GAM
          dFdz = -4.d0*((1.d0 - zeta)**(1.d0/3.d0) 
     &           - (1.d0 + zeta)**(1.d0/3.d0))/(3.d0*GAM)

          eclda1 = EU*(1.d0 - F*z4) + EP*F*z4 - ALFM*F*(1.d0 - z4)/FZZ
          declda1drs = (1.d0 - z4*F)*deudrs + z4*F*depdrs 
     &              - (1.d0 - z4)*F*dalfmdrs/FZZ
          declda1dz = EU*(-4.d0*z3*F - z4*dFdz) + z4*EP*dFdz 
     &              + 4.d0*z3*EP*F 
     &              + 4.d0*z3*alfm*F/FZZ - (1.d0 - z4)*alfm*dFdz/FZZ

          ! Some second derivatives are required for r2scan/r4scan correlation.
          declda1drsrs = z4*F*depdrsrs + (1.d0 - z4*F)*deudrsrs 
     &                 - (1.d0 - z4)*F*dalfmdrsrs/FZZ
          declda1drsz = (4.d0*z3*F*(dalfmdrs + FZZ*(depdrs - deudrs)) 
     &                + ((z4 - 1.d0)*dalfmdrs 
     &                + FZZ*z4*(depdrs - deudrs))*dFdz)/FZZ
      end subroutine

      subroutine grcor2(A, A1, B1, B2, B3, B4, rs, GG, GGRS, GGRSRS)
          IMPLICIT NONE

          real(8), intent(in) :: A, A1, B1, B2, B3, B4, rs
          real(8), intent(out) :: GG, GGRS, GGRSRS
          real(8) :: rtrs, Q0, Q1, Q2
          real(8) :: Q0RS, Q1RS, Q1RSRS, Q2RS, Q2RSRS

          rtrs = sqrt(rs)

          Q0 = -2.d0*A*(1.d0 + A1*rs)
          Q0RS = -2.d0*A*A1

          Q1 = 2.d0*A*rtrs*(B1 + rtrs*(B2 + rtrs*(B3 + B4*rtrs)))
          Q1RS = A*(2.d0*B2 + B1/rtrs + 3.d0*B3*rtrs + 4.d0*B4*rs)
          Q1RSRS = A*(4.d0*B4 - B1/(2.d0*rs*rtrs) + 3.d0*B3/(2.d0*rtrs))

          Q2 = log(1.d0 + 1.d0/Q1)
          Q2RS = -Q1RS/((1.d0 + 1.d0/Q1)*Q1**2)
          Q2RSRS = ((1.d0 + 2.d0*Q1)*Q1RS**2 - Q1*(1.d0 + Q1)*Q1RSRS)
     &           /(Q1**2*(1.d0 + Q1)**2)

          GG = Q0*Q2
          GGRS = Q0*Q2RS + Q2*Q0RS
          GGRSRS = 2.d0*Q0RS*Q2RS + Q0*Q2RSRS
      end subroutine
