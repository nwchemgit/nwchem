c
c     Wrapper routine for evaluating XC functional and derivatives,
c     and combining them with quadrature weights.
c     Copied from xc_quadv0_a.
c
c     BGJ (8/98)
c
      Subroutine xc_eval_fnl(rho, delrho, Amat, Amat2, Cmat, Cmat2,
     &                       nq, Ex, Ec, qwght, GRAD, ldew, func,
     &                       do_2nd)
c
c     $Id: xc_eval_fnl.F,v 1.26 2001-11-21 01:57:13 edo Exp $
c
      implicit none
c
#include "cdft.fh"
#include "dft2drv.fh"
#ifdef BETTIS
#include "screened_op.fh"
#endif
c
      integer nq
      double precision rho(*), delrho(*)
      double precision Amat(*), Cmat(*), Amat2(*), Cmat2(*)
      double precision Ex, Ec
      double precision qwght(nq), func(nq)
      logical GRAD, ldew, do_2nd
      logical dohcth
      character*4 whichf
c
      double precision eps
      integer nx,nc
      parameter (eps=1.e-8)
c
c     Initialize the XC potential and energy sampling matrices.
c
      call dfill(ipol*nq, 0.d0, Amat, 1)
      call dfill(3*nq*ipol, 0.d0, Cmat, 1)
      if (do_2nd) then
         call dfill(nq*NCOL_AMAT2, 0.d0, Amat2, 1)
         if (GRAD) call dfill(nq*NCOL_CMAT2, 0.d0, Cmat2, 1)
      endif
      if (ldew) call dfill(nq, 0.d0, func, 1)
#ifdef BETTIS
c     Update quantities that help in determining screening constant
      call density_average(rho, qwght, nq, ipol, 1)
#endif
c
c     warning!! xc_dirac has to be called before all the other
c     XC routines
c     
      if (abs(xfac(2)).gt.eps)then
         if (.not. do_2nd) then
            call xc_dirac(tol_rho, xfac(2), lxfac(2), nlxfac(2), rho, 
     &           Amat, nq, ipol, Ex, qwght, 
     &           ldew, func)
         else
            call xc_dirac_d2(tol_rho, xfac(2), lxfac(2), nlxfac(2), rho, 
     &           Amat, Amat2, nq, ipol, Ex, qwght, 
     &           .false., func)
         endif
#ifdef BETTIS
         if (screened_exchange()) then
c        Add screened LDA exchange contribution
         if (.not. do_2nd) then
            call xc_dirac_screened(tol_rho, -xfac(2), lxfac(2),
     &           nlxfac(2), rho, Amat, nq, ipol, Ex,
     &           screened_operator_Ks(), qwght, ldew, func)
         else
            call xc_dirac_screened_d2(tol_rho, -xfac(2), lxfac(2),
     &           nlxfac(2), rho, Amat, Amat2, nq, ipol, Ex,
     &           screened_operator_Ks(), qwght, ldew, func)
         endif
         endif
#endif
      endif
      if (abs(xfac(3)).gt.eps)then
         if (.not. do_2nd) then
            call xc_becke88(tol_rho, xfac(3), lxfac(3), nlxfac(3), 
     &           rho, delrho, Amat, Cmat, nq, ipol, 
     &           Ex, qwght,ldew,func)
         else
            call xc_becke88_d2(tol_rho, xfac(3), lxfac(3), nlxfac(3), 
     &           rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol, 
     &           Ex, qwght,ldew,func)
         endif
      endif
      if (abs(xfac(16)).gt.eps)then
         call xc_optx(tol_rho, xfac(16), 2,
     &        rho, delrho, Amat, Cmat, nq, ipol, 
     &        Ex, qwght,ldew,func)
      endif
C
C     hcth , becke97s functionals
C
      nx = 4       ! take care of compiler warnings
      nc = 13

      dohcth=.false.
      if (abs(xfac(4)).gt.eps.or.abs(cfac(13)).gt.eps)then
         whichf='hcth'
         dohcth=.true.
         nx=4
         nc=13
      elseif (abs(xfac(10)).gt.eps.or.abs(cfac(16)).gt.eps)then
         whichf='h120'
         dohcth=.true.
         nx=10
         nc=16
      elseif (abs(xfac(11)).gt.eps.or.abs(cfac(17)).gt.eps)then
         whichf='h147'
         dohcth=.true.
         nx=11
         nc=17
      elseif (abs(xfac(5)).gt.eps.or.abs(cfac(14)).gt.eps)then
         whichf='b970'
         dohcth=.true.
         nx=5
         nc=14
      elseif (abs(xfac(6)).gt.eps.or.abs(cfac(15)).gt.eps)then
         whichf='b971'
         dohcth=.true.
         nx=6
         nc=15
      elseif (abs(xfac(12)).gt.eps.or.abs(cfac(18)).gt.eps)then
         whichf='b980'
         dohcth=.true.
         nx=12
         nc=18
      elseif (abs(xfac(13)).gt.eps.or.abs(cfac(19)).gt.eps)then
         whichf='b97g'
         dohcth=.true.
         nx=13
         nc=19
      elseif (abs(xfac(14)).gt.eps.or.abs(cfac(20)).gt.eps)then
         whichf='h407'
         dohcth=.true.
         nx=14
         nc=20
      elseif (abs(xfac(15)).gt.eps.or.abs(cfac(21)).gt.eps)then
         whichf='hp14'
         dohcth=.true.
         nx=17
         nc=23
      elseif (abs(xfac(17)).gt.eps.or.abs(cfac(23)).gt.eps)then
         whichf='b972'
         dohcth=.true.
         nx=17
         nc=23
      endif
       if(dohcth) then  
         if (.not. do_2nd) then
            call xc_hcth(tol_rho, xfac(nx), lxfac(nx), nlxfac(nx), 
     ,           cfac(nc), lcfac(nc), nlxfac(nc), rho, 
     &           delrho, Amat, Cmat, nq, ipol, Ex, Ec, qwght, 
     &           ldew, func,whichf)
         else
            call errquit('xc_eval_fnl: 2nd derivs not coded yet',0)
         endif
      endif
c     
c     compute partial derivatives of the correlation energy functional.
c     
      if (abs(cfac(1)).gt.eps)then
         if (.not. do_2nd) then
            call xc_vwn_5(tol_rho, cfac(1), rho, 
     &           Amat, nq, ipol, Ec, qwght, 
     &           ldew, func)
         else
            call xc_vwn_5_d2(tol_rho, cfac(1),  rho, 
     &           Amat, Amat2, nq, ipol, Ec, qwght, 
     &           ldew, func)
         endif
      endif
c     
      if (abs(cfac(7)).gt.eps)then
         if (.not. do_2nd) then
            call xc_vwn_1_rpa(tol_rho, cfac(7),  
     &           rho, Amat, nq, ipol, Ec, 
     &           qwght, ldew, func)
         else
            call xc_vwn_1_rpa_d2(tol_rho, cfac(7), 
     &           rho, Amat, Amat2, nq, ipol, Ec, 
     &           qwght, ldew, func)
         endif
      endif
c     
      if (abs(cfac(8)).gt.eps)then
         if (.not. do_2nd) then
            call xc_vwn_1(tol_rho, cfac(8), rho, 
     &           Amat, nq, ipol, Ec, qwght, 
     &           ldew, func)
         else
            call xc_vwn_1_d2(tol_rho, cfac(8),  rho, 
     &           Amat, Amat2, nq, ipol, Ec, qwght, 
     &           ldew, func)
         endif
      endif
c     
      if (abs(cfac(9)).gt.eps)then
         if (.not. do_2nd) then
            call xc_vwn_2(tol_rho, cfac(9),  rho, 
     &           Amat, nq, ipol, Ec, qwght, 
     &           ldew, func)
         else
            call xc_vwn_2_d2(tol_rho, cfac(9),  rho, 
     &           Amat, Amat2, nq, ipol, Ec, qwght, 
     &           ldew, func)
         endif
      endif
c     
      if (abs(cfac(10)).gt.eps)then
         if (.not. do_2nd) then
            call xc_vwn_3(tol_rho, cfac(10), 
     &           rho, Amat, nq, ipol, Ec, qwght, 
     &           ldew, func)
         else
            call xc_vwn_3_d2(tol_rho, cfac(10), 
     &           rho, Amat, Amat2, nq, ipol, Ec, qwght, 
     &           ldew, func)
         endif
      endif
c     
      if (abs(cfac(11)).gt.eps)then
         if (.not. do_2nd) then
            call xc_vwn_4(tol_rho, cfac(11), 
     &           rho, Amat, nq, ipol, Ec, qwght, 
     &           ldew, func)
         else
            call xc_vwn_4_d2(tol_rho, cfac(11), 
     &           rho, Amat, Amat2, nq, ipol, Ec, qwght, 
     &           ldew, func)
         endif
      endif
c     
      if (abs(cfac(6)).gt.eps)then
#ifdef BETTIS
         call errquit('pw91lda on in input',0)
#endif
         if (.not. do_2nd) then
            call xc_pw91lda(tol_rho, cfac(6), lcfac(6), nlcfac(6), 
     &           rho, Amat, nq, ipol, Ec, 
     &           qwght, ldew, func)
         else
            call xc_pw91lda_d2(tol_rho, cfac(6), lcfac(6), nlcfac(6),
     &           rho, Amat, Amat2, nq, ipol, Ec, qwght, 
     &           ldew, func)
         endif
      endif
c     
      if (abs(cfac(2)).gt.eps)then
         if (.not. do_2nd) then
            call xc_lyp88(tol_rho, cfac(2), 
     &           rho, delrho, Amat, Cmat, nq, ipol, Ec, 
     &           qwght, ldew, func)
         else
            call xc_lyp88_d2(tol_rho, cfac(2), 
     &           rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol, Ec, 
     &           qwght, ldew, func)
         endif
      endif
c     
      if (abs(cfac(3)).gt.eps)then
         if (.not. do_2nd) then
            call xc_p81(tol_rho, cfac(3), lcfac(3), nlcfac(3), rho, 
     &           Amat, nq, ipol, Ec, qwght, ldew, func)
         else
            call xc_p81_d2(tol_rho, cfac(3), lcfac(3), nlcfac(3), rho, 
     &           Amat, Amat2, nq, ipol, Ec, qwght, ldew, func)
         endif
      endif
c     
      if (abs(cfac(4)).gt.eps)then
         if (.not. do_2nd) then
            call xc_perdew86(tol_rho, cfac(4), lcfac(4), nlcfac(4), 
     &           rho, delrho, Amat, Cmat, nq, ipol, 
     &           Ec, qwght, ldew, func)
         else
            call xc_perdew86_d2(tol_rho, cfac(4), lcfac(4), nlcfac(4), 
     &           rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol, 
     &           Ec, qwght, ldew, func)
         endif
      endif
c     
c     PW91 is special in that the GGA part is dependent on
c     the E(LDA) ... so more info has to be passed in.
c     
      if (abs(cfac(5)).gt.eps)then
         if (.not. do_2nd) then
            call xc_perdew91(tol_rho, cfac, lcfac, nlcfac, rho, 
     &           delrho, Amat, Cmat, nq, ipol, 
     &           Ec, qwght, ldew, func)
         else
            call errquit('xc_eval_fnl: 2nd derivs not coded yet',0)
         endif
      endif
c     
c     PBE96 is special in that the GGA part is dependent on
c     the E(LDA) ... so more info has to be passed in.
c     
      if (abs(cfac(12)).gt.eps)then
         if (.not. do_2nd) then
            call xc_cpbe96(tol_rho, cfac, lcfac, nlcfac, rho, 
     &           delrho, Amat, Cmat, nq, ipol, 
     &           Ec, qwght, ldew, func)
         else
            call errquit('xc_eval_fnl: 2nd derivs not coded yet',0)
         endif
      endif
      if (abs(xfac(7)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xpbe96(tol_rho, xfac(7), lxfac(7), nlxfac(7), 
     &           rho, delrho, Amat, Cmat, nq, ipol, 
     &           Ex, qwght,ldew,func)
         else
            call xc_xpbe96_d2(tol_rho, xfac(7), lxfac(7), nlxfac(7), 
     &           rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol, 
     &           Ex, qwght,ldew,func)
         endif
      endif
      if (abs(xfac(8)).gt.eps)then
         if (.not. do_2nd) then
            call xc_gill96(tol_rho, xfac(8), lxfac(8), nlxfac(8), 
     &           rho, delrho, Amat, Cmat, nq, ipol, 
     &           Ex, qwght,ldew,func)
         else
            call xc_gill96_d2(tol_rho, xfac(8), lxfac(8), nlxfac(8), 
     &           rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol, 
     &           Ex, qwght,ldew,func)
         endif
      endif
      if (abs(xfac(9)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xpw91(tol_rho, xfac(9), lxfac(9), nlxfac(9), 
     &           rho, delrho, Amat, Cmat, nq, ipol, 
     &           Ex, qwght,ldew,func)
         else
            call xc_xpw91_d2(tol_rho, xfac(9), lxfac(9), nlxfac(9), 
     &           rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol, 
     &           Ex, qwght,ldew,func)
         endif
      endif
      if (abs(cfac(22)).gt.eps)then
         call xc_optc(rho, delrho, 
     &                      Amat, Cmat, nq, Ec, qwght,ldew,func)
      endif
c     
c     Combine with quadrature weights
c
      if (.not. do_2nd) then
         call setACmat(delrho, Amat, Cmat, qwght, ipol, nq, GRAD,
     &                 .not. do_2nd)
      else
         call setACmat_d2(delrho, Amat, Amat2, Cmat, Cmat2, qwght, ipol,
     &                    nq, GRAD, .not. do_2nd)
      endif
c
      return
      end
