c
c     Wrapper routine for evaluating XC functional and derivatives,
c     and combining them with quadrature weights.
c     Copied from xc_quadv0_a.
c
c     BGJ (8/98)
c
      Subroutine xc_eval_fnl(rho, delrho, Amat, Amat2, Cmat, Cmat2,
     &                       nq, Ex, Ec, qwght, grad, ldew, func,
     &                       do_2nd, tau, ttau, kske, Mmat, Mmat2,
     &                       StericEnergy)
c
c     $Id: xc_eval_fnl.F,v 1.65 2009-02-18 02:34:15 niri Exp $
c
      implicit none
c
#include "cdft.fh"
#include "dft2drv.fh"
#include "stdio.fh"
#include "steric.fh"
c
      integer nq
      double precision rho(*), delrho(*)
c
      double precision tau(*), Mmat(*), ttau(*)
c
      double precision Amat(*), Cmat(*), Amat2(*), Cmat2(*), Mmat2(*)
      double precision Ex, Ec
      double precision StericEnergy
      double precision qwght(nq), func(nq)
      logical grad, ldew, do_2nd, kske
      logical dohcth
      character*4 whichf
c
      double precision eps,dumd
      integer nx,nc,dumi
      parameter (eps=1.e-8)
c
c     Initialize the XC potential and energy sampling matrices.
c
      call dfill(ipol*nq, 0.d0, Amat, 1)
      call dfill(3*nq*ipol, 0.d0, Cmat, 1)
      if(kske) call dfill(nq*ipol, 0.d0, Mmat, 1)
      if (do_2nd) then
         call dfill(nq*NCOL_AMAT2, 0.d0, Amat2, 1)
         if (grad) call dfill(nq*NCOL_CMAT2, 0.d0, Cmat2, 1)
         if (kske) call dfill(nq*NCOL_MMAT2, 0.d0, Mmat2, 1)
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
c
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
c
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
         nx=15
         nc=21
      elseif (abs(xfac(17)).gt.eps.or.abs(cfac(23)).gt.eps)then
         whichf='b972'
         dohcth=.true.
         nx=17
         nc=23
      elseif (abs(xfac(20)).gt.eps.or.abs(cfac(26)).gt.eps)then
         whichf='407p'
         dohcth=.true.
         nx=20
         nc=26
      elseif (abs(xfac(22)).gt.eps.or.abs(cfac(28)).gt.eps)then
         whichf='b973'
         dohcth=.true.
         nx=22
         nc=28
      elseif (abs(xfac(39)).gt.eps.or.abs(cfac(41)).gt.eps)then
         whichf='b97d'
         dohcth=.true.
         nx=39
         nc=41
      endif
       if(dohcth) then  
         if (.not. do_2nd) then
            call xc_hcth(tol_rho, xfac(nx), lxfac(nx), nlxfac(nx), 
     ,           cfac(nc), lcfac(nc), nlxfac(nc), rho, 
     &           delrho, Amat, Cmat, nq, ipol, Ex, Ec, qwght, 
     &           ldew, func,whichf)
         else
            call xc_hcth_d2(tol_rho, xfac(nx), lxfac(nx), nlxfac(nx),
     ,           cfac(nc), lcfac(nc), nlxfac(nc), rho,
     &           delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol, Ex, Ec,
     &           qwght, ldew, func,whichf)
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
            call xc_perdew91_d2(tol_rho, cfac, lcfac, nlcfac, rho, 
     &           delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol, 
     &           Ec, qwght, ldew, func)
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
            call xc_cpbe96_d2(tol_rho, cfac, lcfac, nlcfac, rho, 
     &           delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol, 
     &           Ec, qwght, ldew, func)
         endif
      endif
      if (abs(xfac(7)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xpbe96('orig',
     T           tol_rho, xfac(7), lxfac(7), nlxfac(7), 
     &           rho, delrho, Amat, Cmat, nq, ipol, 
     &           Ex, qwght,ldew,func)
         else
            call xc_xpbe96_d2('orig',
     T           tol_rho, xfac(7), lxfac(7), nlxfac(7), 
     &           rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol, 
     &           Ex, qwght,ldew,func)
         endif
      endif
      if (abs(xfac(30)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xpbe96('rpbe',
     T           tol_rho, xfac(30), lxfac(30), nlxfac(30), 
     &           rho, delrho, Amat, Cmat, nq, ipol, 
     &           Ex, qwght,ldew,func)
         else
            call xc_xpbe96_d2('rpbe',
     T           tol_rho, xfac(30), lxfac(30), nlxfac(30), 
     &           rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol, 
     &           Ex, qwght,ldew,func)
         endif
      endif
      if (abs(xfac(31)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xpbe96('revp',
     T           tol_rho, xfac(31), lxfac(31), nlxfac(31), 
     &           rho, delrho, Amat, Cmat, nq, ipol, 
     &           Ex, qwght,ldew,func)
         else
            call xc_xpbe96_d2('revp',
     T           tol_rho, xfac(31), lxfac(31), nlxfac(31), 
     &           rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol, 
     &           Ex, qwght,ldew,func)
         endif
      endif
c
c     SSB-D
c     consists of three parts, SSB-1 (correction to PBEx),
c     portion of KT1 gradient correction, and sPBEc
c     (it also includes a portion of Grimme's dispersion correction)
c     see: Swart, Sola, Bickelhaupt  JCP 2009, 131, 094103
c
c     sPBEc is special in that the GGA part is dependent on
c     the E(LDA) ... so more info has to be passed in.
c
      if (abs(cfac(46)).gt.eps)then
         if (.not. do_2nd) then
            call xc_spbe96(tol_rho, cfac, lcfac, nlcfac,
     &           rho, delrho, Amat, Cmat, nq, ipol, 
     &           Ec, qwght, ldew, func)
         else    
            call xc_spbe96_d2(tol_rho, cfac, lcfac, nlcfac,
     &           rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol,
     &           Ec, qwght, ldew, func)
         endif
      endif
      if (abs(xfac(46)).gt.eps)then
         if (.not. do_2nd) then
c
c           first the part that depends on s (correction to PBEx)
c
            call xc_ssbD_1(tol_rho, xfac(46), lxfac(46), nlxfac(46),
     &           rho, delrho, Amat, Cmat, nq, ipol,
     &           Ex, qwght, ldew, func)
         else
c
c           first the part that depends on s (correction to PBEx)
c
            call xc_ssbD_1_d2(tol_rho, xfac(46), lxfac(46), nlxfac(46),
     &           rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol,
     &           Ex, qwght, ldew, func)
         endif
      endif
c
c kt1
c
      if (abs(xfac(47)).gt.eps)then
         if (.not. do_2nd) then
            call xc_kt1(tol_rho, xfac(47), rho, delrho,
     &                     Amat, Cmat, nq, ipol, Ex, qwght, ldew, func)
         
         else
            call xc_kt1_d2(tol_rho, xfac(47), rho, delrho,
     &                     Amat, Amat2, Cmat, Cmat2, nq, ipol,
     &                     Ex, qwght, ldew, func)
         endif 
      endif
c
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
      if (abs(xfac(19)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xmpw91(tol_rho,xfac(19),lxfac(19),nlxfac(19), 
     &           rho, delrho, Amat, Cmat, nq, ipol, 
     &           Ex, qwght,ldew,func)
         else
            call xc_xmpw91_d2(tol_rho,xfac(19),lxfac(19),nlxfac(19), 
     &           rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol, 
     &           Ex, qwght,ldew,func)
         endif
      endif
c
      if (abs(xfac(25)).gt.eps.or.abs(cfac(24)).gt.eps)then
         if (.not. do_2nd) then
            call xc_ft97(tol_rho,xfac(25),lxfac(25),nlxfac(25), 
     .           cfac(24),lcfac(24),nlcfac(24), 
     &           rho, delrho, Amat, Cmat, nq, ipol, 
     &           Ex, Ec, qwght,ldew,func)
         else
            call errquit('2nd derivative not available 
     &for this xc functional',0,0)
!           call xc_ft97_d2()
         endif
      endif
c
      if ((abs(cfac(36)).gt.eps).or.(abs(cfac(37)).gt.eps))then
         if(abs(cfac(36)).gt.eps) then
            nc=36
            whichf='be88'
         endif
         if(abs(cfac(37)).gt.eps) then
            nc=37
            whichf='pb96'
         endif
         if (.not. do_2nd) then
            call xc_op(tol_rho,whichf,
     .           cfac(nc),lcfac(nc),nlcfac(nc), 
     &           rho, delrho, Amat, Cmat, nq, ipol, 
     &           Ec, qwght,ldew,func)
         else
            call errquit('2nd derivative not available 
     &for this xc functional',0,0)
!           call xc_op_d2()
         endif
      endif
c
c     meta GGA
c
      if (abs(xfac(18)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xpkzb99(tol_rho, xfac(18), lxfac(18), nlxfac(18), 
     &           rho, delrho, Amat, Cmat, nq, ipol, 
     &           Ex, qwght,ldew,func,tau, Mmat)
         else
            call xc_xpkzb99_d2()
         endif
      endif
c
c     LB94 or CS00 correction is added to xc potential only
c     (xc functional and functional 2nds are unchanged)
c
      if (cs00) then
         call xc_cs00(tol_rho, xfac(1), rho, delrho, Amat, nq, ipol, 
     &   delta_ac, e_homo)
      else if (lb94) then
         call xc_lb94(tol_rho, xfac(1), rho, delrho, Amat, nq, ipol)
      endif
c     
c     PKZB99-COR is special in that the GGA part is
c     defined to be  PBE COR GGA  and also is dependent on
c     the E(LDA) ... 
c     the decision has been made to use the PW91-LDA as the
c     LDA-correlation.  at present, this LDA  cannot be
c     set by the user
c     
      if (abs(cfac(25)).gt.eps)then
         if (.not. do_2nd) then
            call xc_cpkzb99(tol_rho, cfac(25), lcfac(25), nlcfac(25),  
     &           rho, delrho,  nq, ipol, 
     &           Ec, qwght, ldew, func, tau,Amat,Cmat,Mmat)
         else
            call xc_cpkzb99_d2()
         endif
      endif
c
c   TPSS  meta GGA
c
      if (abs(xfac(21)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xtpss03(tol_rho, xfac(21),  
     &           rho, delrho, Amat, Cmat, nq, ipol, 
     &           Ex, qwght,ldew,func,tau,Mmat)
         else
            call xc_xtpss03_d2()
         endif
      endif
c     
c     TPSS03-COR is special in that the GGA part is
c     defined to be  PBE COR GGA  and also is dependent on
c     the E(LDA) ... 
c     the decision has been made to use the PW91-LDA as the
c     LDA-correlation.  at present, this LDA  cannot be
c     set by the user
c     
      if (abs(cfac(27)).gt.eps)then
         if (.not. do_2nd) then
            call xc_ctpss03(tol_rho, cfac(27), lcfac(27), nlcfac(27),  
     &           rho, delrho,  nq, ipol, 
     &           Ec, qwght, ldew, func, tau,Amat,Cmat,Mmat)

         else
            call xc_ctpss03_d2()
         endif
      endif
c     
c     Bc95-COR is special in that the GGA part is
c     defined to be dependent on
c     the E(LDA) ... 
c     the decision has been made to use the PW91-LDA as the
c     LDA-correlation.  at present, this LDA  cannot be
c     set by the user
c     
c     Note that bc95, cpw6b95, cpwb6k use the same subroutine xc_bc95()
c
      if (abs(cfac(31)).gt.eps)then
         if (.not. do_2nd) then
            call xc_bc95(tol_rho, cfac(31), lcfac(31), nlcfac(31),
     &           rho, delrho,  nq, ipol,
     &           Ec, qwght, ldew,func, tau,Amat,Cmat,Mmat,0)
      
         else
            call xc_bc95_d2()
         endif
      endif
c   
c   PW6B95 Exchange part
c
      if (abs(xfac(26)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xpw6(tol_rho,xfac(26),lxfac(26),nlxfac(26),
     &           rho, delrho, Amat, Cmat, nq, ipol,
     &           Ex, qwght,ldew,func,1)
         else
            call xc_xpw6_d2(tol_rho,xfac(26),lxfac(26),nlxfac(26),
     &           rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol,
     &           Ex, qwght,ldew,func,1)
         endif
      endif
c
c   PWB6K Exchange part
c
      if (abs(xfac(27)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xpw6(tol_rho,xfac(27),lxfac(27),nlxfac(27),
     &           rho, delrho, Amat, Cmat, nq, ipol,
     &           Ex, qwght,ldew,func,2)
         else
            call xc_xpw6_d2(tol_rho,xfac(27),lxfac(27),nlxfac(27),
     &           rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol,
     &           Ex, qwght,ldew,func,2)
         endif
      endif
c
c M05   meta GGA Exchange
c
      if (abs(xfac(28)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xm05(tol_rho, xfac(28), lxfac(28), nlxfac(28),
     &           rho, delrho, Amat, Cmat, nq, ipol,
     &           Ex, qwght,ldew,func,tau,Mmat,1)
         else
            call xc_xm05_d2()
         endif
      endif
c
c M05-2X   meta GGA Exchange
c      
      if (abs(xfac(29)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xm05(tol_rho, xfac(29), lxfac(29), nlxfac(29),
     &           rho, delrho, Amat, Cmat, nq, ipol,
     &           Ex, qwght,ldew, func,tau,Mmat,2)
         else
            call xc_xm05_d2()
         endif
      endif
c
c VSXC   meta GGA Exchange
c
      if (abs(xfac(33)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xvs98(tol_rho, xfac(33), lxfac(33), nlxfac(33),
     &           rho, delrho, Amat, Cmat, nq, ipol,
     &           Ex, qwght,ldew, func,tau,Mmat,1)
         else
            call xc_xvs98_d2()
         endif
      endif
c
c M06-L   meta GGA Exchange
c
      if (abs(xfac(34)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xm06(tol_rho, xfac(34), lxfac(34), nlxfac(34),
     &           rho, delrho, Amat, Cmat, nq, ipol,
     &           Ex, qwght,ldew, func,tau,Mmat,1)
         else
            call xc_xvs98_d2()
         endif
      endif
c
c M06-HF   meta GGA Exchange
c
      if (abs(xfac(35)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xm06(tol_rho, xfac(35), lxfac(35), nlxfac(35),
     &           rho, delrho, Amat, Cmat, nq, ipol,
     &           Ex, qwght,ldew, func,tau,Mmat,2)
         else
            call xc_xm06_d2()
         endif
      endif
c
c M06   meta GGA Exchange
c
      if (abs(xfac(36)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xm06(tol_rho, xfac(36), lxfac(36), nlxfac(36),
     &           rho, delrho, Amat, Cmat, nq, ipol,
     &           Ex, qwght,ldew, func,tau,Mmat,3)
         else
            call xc_xm06_d2()
         endif
      endif
c
c M06-2X  meta GGA Exchange
c
      if (abs(xfac(37)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xm06(tol_rho, xfac(37), lxfac(37), nlxfac(37),
     &           rho, delrho, Amat, Cmat, nq, ipol,
     &           Ex, qwght,ldew, func,tau,Mmat,4)
         else
            call xc_xm06_d2()
         endif
      endif
c
c LC-BNL 2007 Exchange 
c
      if (abs(xfac(38)).gt.eps)then
       if (.not. do_2nd) then
        call xc_rodaes_erf(tol_rho, xfac(38), lxfac(38), nlxfac(38), 
     &    rho, Amat, nq, ipol, Ex, qwght, ldew, func)
       else
        call xc_rodaes_erf_d2(tol_rho, xfac(38), lxfac(38), nlxfac(38), 
     &    rho, Amat, Amat2, nq, ipol, Ex, qwght, ldew, func)
       endif
      endif
c
c CAM-B88 Exchange
c
      if (abs(xfac(40)).gt.eps)then
       if (.not. do_2nd) then
        call xc_camb88(tol_rho, xfac(40), lxfac(40), nlxfac(40),
     &    rho, delrho, Amat, Cmat, nq, ipol,
     &    Ex, qwght,ldew,func)
       else
        call xc_camb88_d2(tol_rho, xfac(40), lxfac(40), nlxfac(40),
     &    rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol,
     &    Ex, qwght,ldew,func)
       endif
      endif
c
c CAM-PBE96 Exchange
c
      if (abs(xfac(41)).gt.eps)then
       if (.not. do_2nd) then
        call xc_camxpbe96('orig',
     T           tol_rho, xfac(41), lxfac(41), nlxfac(41), 
     &           rho, delrho, Amat, Cmat, nq, ipol, 
     &           Ex, qwght,ldew,func)
       else
        call xc_camxpbe96_d2('orig',
     T           tol_rho, xfac(41), lxfac(41), nlxfac(41), 
     &           rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol, 
     &           Ex, qwght,ldew,func)
       endif
      endif
c
c CAM-LSD Exchange
c
      if (abs(xfac(42)).gt.eps)then
         if (.not. do_2nd) then
           call xc_camxlsd(tol_rho, xfac(42), lxfac(42), nlxfac(42), 
     &           rho, Amat, nq, ipol, Ex, qwght, ldew, func)
         else
           call xc_camxlsd_d2(tol_rho, xfac(42), lxfac(42), nlxfac(42), 
     &          rho, Amat, Amat2, nq, ipol, Ex, qwght, .false., func)
         endif
      endif
c
c     cpw6b95 is special in that the GGA part is
c     defined to be dependent on
c     the E(LDA) ...
c     the decision has been made to use the PW91-LDA as the
c     LDA-correlation.  at present, this LDA  cannot be
c     set by the user
c     Note that bc95, cpw6b95, cpwb6k use the same subroutine xc_bc95()
c
      if (abs(cfac(32)).gt.eps)then
         if (.not. do_2nd) then
            call xc_bc95(tol_rho, cfac(32), lcfac(32), nlcfac(32),
     &           rho, delrho,  nq, ipol,
     &           Ec, qwght, ldew,func, tau,Amat,Cmat,Mmat,1)

         else
            call xc_bc95_d2()
         endif
      endif
c
c     cpwb6k is special in that the GGA part is
c     defined to be dependent on
c     the E(LDA) ...
c     the decision has been made to use the PW91-LDA as the
c     LDA-correlation.  at present, this LDA  cannot be
c     set by the user
c
      if (abs(cfac(33)).gt.eps)then
         if (.not. do_2nd) then
            
            call xc_bc95(tol_rho, cfac(33), lcfac(33), nlcfac(33),
     &           rho, delrho,  nq, ipol,
     &           Ec, qwght, ldew,func,tau,Amat,Cmat,Mmat,2)

         else
            call xc_bc95_d2()
         endif
      endif
c
c     cm05 is special in that the GGA part is
c     defined to be dependent on
c     the E(LDA) ...
c     the decision has been made to use the PW91-LDA as the
c     LDA-correlation.  at present, this LDA  cannot be
c     set by the user
c     
      if (abs(cfac(34)).gt.eps)then
         if (.not. do_2nd) then
            call xc_cm05(tol_rho, cfac(34), lcfac(34), nlcfac(34),
     &           rho, delrho,  nq, ipol,
     &           Ec, qwght, ldew, func,tau,Amat,Cmat,Mmat,1)

         else
            call xc_cm05_d2()
         endif
      endif
c
c     cm05-2x is special in that the GGA part is
c     defined to be dependent on
c     the E(LDA) ...
c     the decision has been made to use the PW91-LDA as the
c     LDA-correlation.  at present, this LDA  cannot be
c     set by the user
c
      if (abs(cfac(35)).gt.eps)then
         if (.not. do_2nd) then
            call xc_cm05(tol_rho, cfac(35), lcfac(35), nlcfac(35),
     &           rho, delrho,  nq, ipol,
     &           Ec, qwght, ldew,func,tau,Amat,Cmat,Mmat,2)

         else
            call xc_cm05_d2()
         endif
      endif
c
c     cvs98 is special in that the GGA part is
c     defined to be dependent on
c     the E(LDA) ...
c     the decision has been made to use the PW91-LDA as the
c     LDA-correlation.  at present, this LDA  cannot be
c     set by the user
c
      if (abs(cfac(29)).gt.eps)then
         if (.not. do_2nd) then
            call xc_cvs98(tol_rho, cfac(29), lcfac(29), nlcfac(29),
     &           rho, delrho,  nq, ipol,
     &           Ec, qwght, ldew,func,tau,Amat,Cmat,Mmat,1)

         else
            call xc_cvs98_d2()
         endif
      endif
c
c     cm06-L is special in that the GGA part is
c     defined to be dependent on
c     the E(LDA) ...
c     the decision has been made to use the PW91-LDA as the
c     LDA-correlation.  at present, this LDA  cannot be
c     set by the user
c
      if (abs(cfac(30)).gt.eps)then
         if (.not. do_2nd) then
            call xc_cm06(tol_rho, cfac(30), lcfac(30), nlcfac(30),
     &           rho, delrho,  nq, ipol,
     &           Ec, qwght, ldew,func,tau,Amat,Cmat,Mmat,1)

         else
            call xc_cm06_d2()
         endif
      endif
c
c     cm06-hf is special in that the GGA part is
c     defined to be dependent on
c     the E(LDA) ...
c     the decision has been made to use the PW91-LDA as the
c     LDA-correlation.  at present, this LDA  cannot be
c     set by the user
c
      if (abs(cfac(38)).gt.eps)then
         if (.not. do_2nd) then
            call xc_cm06(tol_rho, cfac(38), lcfac(38), nlcfac(38),
     &           rho, delrho,  nq, ipol,
     &           Ec, qwght, ldew,func,tau,Amat,Cmat,Mmat,2)

         else
            call xc_cm06_d2()
         endif
      endif
c
c     cm06 is special in that the GGA part is
c     defined to be dependent on
c     the E(LDA) ...
c     the decision has been made to use the PW91-LDA as the
c     LDA-correlation.  at present, this LDA  cannot be
c     set by the user
c
      if (abs(cfac(39)).gt.eps)then
         if (.not. do_2nd) then
            call xc_cm06(tol_rho, cfac(39), lcfac(39), nlcfac(39),
     &           rho, delrho,  nq, ipol,
     &           Ec, qwght, ldew,func,tau,Amat,Cmat,Mmat,3)

         else
            call xc_cm06_d2()
         endif
      endif
c
c     cm06-2x is special in that the GGA part is
c     defined to be dependent on
c     the E(LDA) ...
c     the decision has been made to use the PW91-LDA as the
c     LDA-correlation.  at present, this LDA  cannot be
c     set by the user
c
      if (abs(cfac(40)).gt.eps)then
         if (.not. do_2nd) then
            call xc_cm06(tol_rho, cfac(40), lcfac(40), nlcfac(40),
     &           rho, delrho,  nq, ipol,
     &           Ec, qwght, ldew,func,tau,Amat,Cmat,Mmat,4)

         else
            call xc_cm06_d2()
         endif
      endif
c
c     Calculate the steric energy
c
      if (lsteric) then
        StericEnergy = 0.d0
        call steric_energy(tol_rho,xfac(1),rho,delrho,nq,
     &   qwght,ipol,StericEnergy)
      endif
c
c     Combine with quadrature weights
c
      if (.not. do_2nd) then
         call setACmat(delrho, Amat, Cmat, qwght, ipol, nq, grad,
     &         (.not. do_2nd), kske, Mmat)
      else
         call setACmat_d2(delrho, Amat, Amat2, Cmat, Cmat2, qwght, ipol,
     &         nq, grad, (.not. do_2nd), kske, Mmat2)
      endif
c
      return
      end
