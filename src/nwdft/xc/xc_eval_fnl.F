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
c     $Id: xc_eval_fnl.F,v 1.1 1998-09-09 13:32:36 bjohnson Exp $
c
      implicit none
c
#include "cdft.fh"
#include "dft2drv.fh"
c
      integer nq
      double precision rho(*), delrho(*)
      double precision Amat(*), Cmat(*), Amat2(*), Cmat2(*)
      double precision Ex, Ec
      double precision qwght(nq), func(nq)
      logical GRAD, ldew, do_2nd
c
      double precision eps
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
      endif
      if (abs(xfac(3)).gt.eps)then
         call xc_becke88(tol_rho, xfac(3), lxfac(3), nlxfac(3), 
     &        rho, delrho, Amat, Cmat, nq, ipol, 
     &        Ex, qwght)
      endif
c     
c     compute partial derivatives of the correlation energy functional.
c     
      if (abs(cfac(1)).gt.eps)then
         call xc_vwn_5(tol_rho, cfac(1), lcfac(1), nlcfac(1), rho, 
     &        Amat, nq, ipol, Ec, qwght, 
     &        ldew, func)
      endif
c     
      if (abs(cfac(7)).gt.eps)then
         call xc_vwn_1_rpa(tol_rho, cfac(7), lcfac(7), nlcfac(7), 
     &        rho, Amat, nq, ipol, Ec, 
     &        qwght, ldew, func)
      endif
c     
      if (abs(cfac(8)).gt.eps)then
         call xc_vwn_1(tol_rho, cfac(8), lcfac(8), nlcfac(8), rho, 
     &                       Amat, nq, ipol, Ec, qwght, 
     &        ldew, func)
      endif
c     
      if (abs(cfac(9)).gt.eps)then
         call xc_vwn_2(tol_rho, cfac(9), lcfac(9), nlcfac(9), rho, 
     &        Amat, nq, ipol, Ec, qwght, 
     &        ldew, func)
      endif
c     
      if (abs(cfac(10)).gt.eps)then
         call xc_vwn_3(tol_rho, cfac(10), lcfac(10), nlcfac(10), 
     &        rho, Amat, nq, ipol, Ec, qwght, 
     &        ldew, func)
      endif
c     
      if (abs(cfac(11)).gt.eps)then
         call xc_vwn_4(tol_rho, cfac(11), lcfac(11), nlcfac(11),
     &        rho, Amat, nq, ipol, Ec, qwght, 
     &        ldew, func)
      endif
c     
      if (abs(cfac(6)).gt.eps)then
         call xc_pw91lda(tol_rho, cfac(6), lcfac(6), nlcfac(6), 
     &        rho, Amat, nq, ipol, Ec, 
     &        qwght)
      endif
c     
      if (abs(cfac(2)).gt.eps)then
         call xc_lyp88(tol_rho, cfac(2), lcfac(2), nlcfac(2), 
     &        rho, delrho, Amat, Cmat, nq, ipol, Ec, 
     &        qwght)
      endif
c     
      if (abs(cfac(3)).gt.eps)then
         call xc_p81(tol_rho, cfac(3), lcfac(3), nlcfac(3), rho, 
     &        Amat, nq, ipol, Ec, qwght)
      endif
c     
      if (abs(cfac(4)).gt.eps)then
         call xc_perdew86(tol_rho, cfac(4), lcfac(4), nlcfac(4), 
     &        rho, delrho, Amat, Cmat, nq, ipol, 
     &        Ec, qwght)
      endif
c     
c     PW91 is special in that the GGA part is dependent on
c     the E(LDA) ... so more info has to be passed in.
c     
      if (abs(cfac(5)).gt.eps)then
         call xc_perdew91(tol_rho, cfac, lcfac, nlcfac, rho, 
     &        delrho, Amat, Cmat, nq, ipol, 
     &        Ec, qwght)
      endif
c     
c     Combine with quadrature weights
c
      if (.not. do_2nd) then
         call setACmat(Amat, Cmat, qwght, ipol, nq, GRAD)
      else
         call setACmat_d2(Amat, Amat2, Cmat, Cmat2, qwght, ipol, nq,
     &                    GRAD)
      endif
c
      return
      end
