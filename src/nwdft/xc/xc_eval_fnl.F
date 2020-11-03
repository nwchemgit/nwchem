c
c     Wrapper routine for evaluating XC functional and derivatives,
c     and combining them with quadrature weights.
c     Copied from xc_quadv0_a.
c
c     BGJ (8/98)
c
C> \ingroup nwdft_xc
C> @{
C>
C> \brief the functional evaluation routine
C>
C> This routine iterates over the components of the current density
C> functional and evaluates and accumulates all density dependent 
C> terms.
C> At present it will evaluate
C>
C> - the exchange-correlation energy
C>
C> - the exchange-correlation 1st order partial derivatives
C>
C> - the exchange-correlation 2nd order partial derivatives
C>
C> as desired.
C>
      Subroutine xc_eval_fnl(rho, delrho, Amat, Amat2, Cmat, Cmat2,
     &                       nq, Ex, Ec, qwght, grad, ldew, func,
     &                       do_2nd, ttau, kske, Mmat, Mmat2,
     &                       laprho, kslap, Lmat,
     &                       StericEnergy,
     &                       do_3rd, Amat3, Cmat3, ldew2)
c
      implicit none
c
#include "cdft.fh"
c xc-second derivative header file
#include "dft2drv.fh"
c xc-third derivative header file
#include "dft3drv.fh"
#include "stdio.fh"
#include "steric.fh"
cHvD
#include "errquit.fh"
#include "mafdecls.fh"
#include "util.fh"
      integer l_f, k_f
      integer l_g, k_g
      integer ii,iq
cHvD
c
c List of functionals with implemented third derivatives.
c
c    Exchange:    Slater/Dirac
c                 PBE 96 (original, RPBE, and revPBE)
c                 Becke 88
c                 CAM-Becke 88
c                 CAM-PBE 96
c                 BNL 
c                 LC-wPBE
c    Correlation: VWN (1-5 and 1 RPA)
c                 LYP
c                 PBE 96 
c                 PW91 LDA (needed for the default local part of PBE)
c                 Perdew 86 (restricted only at the moment)
c                 Perdew 81 LDA (needed for the default local part of Perdew 86)
c
      integer nq, is12x
      double precision rho(*)    !< [Input] The electron density
                                 !< \f$\rho^\alpha\f$ and
                                 !< \f$\rho^\beta\f$.
      double precision delrho(*) !< [Input] The electron density gradient
                                 !< \f$\partial\rho^\alpha/\partial r\f$
                                 !< and
                                 !< \f$\partial\rho^\beta/\partial r\f$
      double precision laprho(*)
c
      double precision Amat(*), Cmat(*)    
      double precision Amat2(*),Cmat2(*)
      double precision Amat3(*),Cmat3(*)
      double precision ttau(*) !< [Input] The kinetic energy density
      double precision Mmat(*), Mmat2(*), Mmat3(2)
      double precision Lmat(*)
c
      double precision Ex !< [Output] The exchange energy
      double precision Ec !< [Output] The correlation energy
      double precision StericEnergy
      double precision qwght(nq), func(nq)
      logical grad, kske, dohcth, use_nwxc, kslap
c
      character*4 whichf

      logical do_2nd, do_3rd
      logical ldew, ldew2, ldew3
cc AJL/Unused
c      logical out1
cc AJL/End
c
      double precision eps,dumd
      integer nx,nc,dumi
      parameter (eps=1.e-8)
c
c     Initialize the XC potential and energy sampling matrices.
c
      call dcopy(ipol*nq, 0.0d0, 0, Amat, 1)
      if(grad) call dcopy(3*nq*ipol, 0.0d0, 0, Cmat, 1)
      if(kske) call dcopy(nq*ipol, 0.0d0, 0, Mmat, 1)
      if(kslap) call dcopy(nq*ipol, 0.0d0, 0, Lmat, 1)
      if (do_2nd) then
         call dcopy(nq*NCOL_AMAT2, 0.0d0, 0, Amat2, 1)
         if (grad) call dcopy(nq*NCOL_CMAT2, 0.0d0, 0, Cmat2, 1)
         if (kske) call dcopy(nq*NCOL_MMAT2, 0.0d0, 0, Mmat2, 1)
      endif
c
c     Initialize the 3rd derivatives.  Note that we may need 
c     the 2nd derivative matrices in some cases as well, so
c     those are initialized also.
      if (do_3rd) then
         call dcopy(nq*NCOL_AMAT2, 0.0d0, 0, Amat2, 1)
         call dcopy(nq*NCOL_AMAT3, 0.0d0, 0, Amat3, 1)
         if (grad) call dcopy(nq*NCOL_CMAT2, 0.0d0, 0, Cmat2, 1)
         if (grad) call dcopy(nq*NCOL_CMAT3, 0.0d0, 0, Cmat3, 1)
         if (kske) call dcopy(nq*NCOL_MMAT2, 0.0d0, 0, Mmat2, 1)
         if (kske) call errquit("xc_eval_fnc: do_3rd and kske",0,
     +                          CAPMIS_ERR)
      endif      
c
c     This prevents the XC-kernel from being multiplied
c     with the quadrature weights in TDDFT gradients. 
c
      ldew3 = (do_3rd.and.ldew)
      if (ldew.or.ldew2.or.ldew3) call dcopy(nq, 0.d0, 0, func, 1)
c
      use_nwxc = util_module_avail("nwxc")
      if (use_nwxc) then
         call nwxc_getvals("nwxc_is_on",use_nwxc)
      endif
      if (use_nwxc) then
c
c       Use the nwxc functional library
c
        if (.not.ma_push_get(mt_dbl,nq,"ffunc",l_f,k_f))
     +    call errquit("xc_eval_fnc: no memory",nq,MA_ERR)
        if (.not.ma_push_get(mt_dbl,3*nq,"rgamma",l_g,k_g))
     +    call errquit("xc_eval_fnc: no memory",3*nq,MA_ERR)
        if (grad) call calc_rgamma(ipol,nq,delrho,dbl_mb(k_g))
        if (ipol.eq.1) then
          if (do_3rd) then
            call nwxc_eval_df3(ipol,nq,rho(1),dbl_mb(k_g),ttau,
     +                         dbl_mb(k_f),
     +                         Amat,Amat2,Amat3,Cmat,Cmat2,Cmat3,
     +                         Mmat,Mmat2,Mmat3)
          else if (do_2nd) then
            call nwxc_eval_df2(ipol,nq,rho(1),dbl_mb(k_g),ttau,
     +                         dbl_mb(k_f),
     +                         Amat,Amat2,Cmat,Cmat2,Mmat,Mmat2)
          else
            call nwxc_eval_df(ipol,nq,rho(1),dbl_mb(k_g),ttau,
     +                        dbl_mb(k_f),
     +                        Amat,Cmat,Mmat)
          endif
        else
          if (do_3rd) then
            call nwxc_eval_df3(ipol,nq,rho(nq+1),dbl_mb(k_g),ttau,
     +                         dbl_mb(k_f),
     +                         Amat,Amat2,Amat3,Cmat,Cmat2,Cmat3,
     +                         Mmat,Mmat2,Mmat3)
          else if (do_2nd) then
            call nwxc_eval_df2(ipol,nq,rho(nq+1),dbl_mb(k_g),ttau,
     +                         dbl_mb(k_f),
     +                         Amat,Amat2,Cmat,Cmat2,Mmat,Mmat2)
          else
            call nwxc_eval_df(ipol,nq,rho(nq+1),dbl_mb(k_g),ttau,
     +                        dbl_mb(k_f),
     +                        Amat,Cmat,Mmat)
          endif
        endif
c
c       The NWXC module calculates derivatives wrt the proper 
c       kinetic energy density. It seems that the original subroutines
c       calculate derivatives wrt. twice the kinetic energy density
c       (in many cases people choose to absorb the factor half in the
c       kinetic energy density in the functional expression). As a
c       result there is factor half between what NWXC calculates and
c       what NWChem expects. Hence we need to scale Mmat here.
c
        if (kske) then
          call dscal(nq*ipol,0.5d0,mmat,1)
        endif
c
        Ex = Ex + ddot(nq,dbl_mb(k_f),1,qwght,1)
        if (ldew) then
          call dcopy(nq,dbl_mb(k_f),1,func,1)
        endif
        if (.not.ma_pop_stack(l_g))
     +    call errquit("xc_eval_fnc: no free",3*nq,MA_ERR)
        if (.not.ma_pop_stack(l_f))
     +    call errquit("xc_eval_fnc: no free",nq,MA_ERR)
c
c       Combine with quadrature weights
c
c Daniel (1-11-13): Added XC-third derivative stuff.  This currently
c doesn't include any functionality for meta-GGAs.
c Daniel (1-17-13): Because we need derivatives of the quadrature
c weights for TDDFT gradients, we don't want to multiply the A 
c and C matrices with the quadrature weights here (which would be
c double counting).  This prevents us from multiplying the
c quadrature weights into the matrices and then dividing them out
c later, which is prone to numerical problems (and slow). 
c Daniel (2-4-13): ldew3 is for the dfxc*(X+Y)*(X+Y) term in the TDDFT 
c gradients.  ldew2 is for the dVxc*P term in the TDDFT gradients.
        if (.not.(ldew2.or.ldew3)) then

          if (.not. do_2nd .and. .not. do_3rd) then
             call setACmat(delrho, Amat, Cmat, qwght, ipol, nq, grad,
     &             (.not. do_2nd), kske, Mmat, kslap, Lmat)
          else if (.not. do_3rd) then
             call setACmat_d2(delrho, Amat, Amat2, Cmat, Cmat2, qwght, 
     &             ipol, nq, grad, (.not. do_2nd), kske, Mmat, Mmat2,
     &             .false.)
          else
             call setACmat_d3(delrho, Amat, Amat2, Amat3, Cmat, Cmat2,
     &             Cmat3, qwght, ipol, nq, grad, (.not. do_3rd),
     &             .false.,.false.)
          endif
        endif
        return
      endif
c      if (ldew.or.ldew3) call dfill(nq, 0.d0, func, 1)
c      if (ldew) call dfill(nq, 0.d0, func, 1)
c
c     warning!! xc_dirac has to be called before all the other
c     XC routines
c     
      if (abs(xfac(2)).gt.eps)then
         if (.not. do_2nd .and. .not. do_3rd) then
            call xc_dirac(tol_rho, xfac(2), lxfac(2), nlxfac(2), rho,
     &           Amat, nq, ipol, Ex, qwght,
     &           ldew, func)
         else if (.not. do_3rd) then
            call xc_dirac_d2(tol_rho, xfac(2), lxfac(2), nlxfac(2), rho,
     &           Amat, Amat2, nq, ipol, Ex, qwght,
     &           .false., func)
         else
            call xc_dirac_d3(tol_rho, xfac(2), lxfac(2), nlxfac(2), rho,
     &           Amat, Amat2, Amat3, nq, ipol, Ex, qwght,
     &           .false., func)
         endif
      endif
c
c     LB94 potential and energy (not variational!)
c
      if (abs(xfac(70)).gt.eps)then
         if (.not. do_2nd) then
            call xc_lb94_e(tol_rho, xfac(70), rho, delrho, Amat, nq,
     I          ipol, Ex, qwght,ldew,func,ttau)
         else
        call xc_lb94_e_d2(tol_rho, xfac(70), rho, delrho, Amat, 
     I           Amat2, nq, ipol, Ex, qwght,.false.,func,ttau)
c            call errquit( ' lb94 2nds not ready yet',0,0)
         endif
      endif
c
      if (abs(xfac(3)).gt.eps)then
         if (.not. do_2nd .and. .not. do_3rd) then
            call xc_becke88(tol_rho, xfac(3), lxfac(3), nlxfac(3), 
     &           rho, delrho, Amat, Cmat, nq, ipol, 
     &           Ex, qwght,ldew,func)
         else if (.not. do_3rd) then
            call xc_becke88_d2(tol_rho, xfac(3), lxfac(3), nlxfac(3), 
     &           rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol, 
     &           Ex, qwght,ldew,func)
         else
            call xc_becke88_d3(tol_rho, xfac(3), lxfac(3), nlxfac(3), 
     &           rho, delrho, Amat, Amat2, Amat3, Cmat, Cmat2, Cmat3,
     &           nq, ipol, Ex, qwght,ldew,func)
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
c
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
      elseif (abs(cfac(45)).gt.eps)then
         whichf='n120'
         dohcth=.true.
         nx=45
         nc=45
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
         endif  ! do_2nd
      endif  ! dohctch
c     
c     compute partial derivatives of the correlation energy functional.
c     
      if (abs(cfac(1)).gt.eps)then
         if (.not. do_2nd .and..not. do_3rd) then
            call xc_vwn_5(tol_rho, cfac(1), rho,
     &           Amat, nq, ipol, Ec, qwght,
     &           ldew, func)
         else if (.not. do_3rd) then
            call xc_vwn_5_d2(tol_rho, cfac(1),  rho,
     &           Amat, Amat2, nq, ipol, Ec, qwght,
     &           ldew, func)
         else
            call xc_vwn_5_d3(tol_rho, cfac(1),  rho,
     &           Amat, Amat2, Amat3, nq, ipol, Ec, qwght,
     &           ldew, func)
         endif
      endif
c     
      if (abs(cfac(7)).gt.eps)then
         if (.not. do_2nd .and..not. do_3rd) then
            call xc_vwn_1_rpa(tol_rho, cfac(7),
     &           rho, Amat, nq, ipol, Ec,
     &           qwght, ldew, func)
         else if (.not. do_3rd) then
            call xc_vwn_1_rpa_d2(tol_rho, cfac(7),
     &           rho, Amat, Amat2, nq, ipol, Ec,
     &           qwght, ldew, func)
         else
            call xc_vwn_1_rpa_d3(tol_rho, cfac(7),
     &           rho, Amat, Amat2, Amat3, nq, ipol, Ec,
     &           qwght, ldew, func)
         endif
      endif
c     
      if (abs(cfac(8)).gt.eps)then
         if (.not. do_2nd .and..not. do_3rd) then
            call xc_vwn_1(tol_rho, cfac(8), rho,
     &           Amat, nq, ipol, Ec, qwght,
     &           ldew, func)
         else if (.not. do_3rd) then
            call xc_vwn_1_d2(tol_rho, cfac(8),  rho,
     &           Amat, Amat2, nq, ipol, Ec, qwght,
     &           ldew, func)
         else
            call xc_vwn_1_d3(tol_rho, cfac(8),  rho,
     &           Amat, Amat2, Amat3, nq, ipol, Ec, qwght,
     &           ldew, func)
         endif
      endif
c     
      if (abs(cfac(9)).gt.eps)then
         if (.not. do_2nd .and..not. do_3rd) then
            call xc_vwn_2(tol_rho, cfac(9),  rho,
     &           Amat, nq, ipol, Ec, qwght,
     &           ldew, func)
         else if (.not. do_3rd) then
            call xc_vwn_2_d2(tol_rho, cfac(9),  rho,
     &           Amat, Amat2, nq, ipol, Ec, qwght,
     &           ldew, func)
         else
            call xc_vwn_2_d3(tol_rho, cfac(9),  rho,
     &           Amat, Amat2, Amat3, nq, ipol, Ec, qwght,
     &           ldew, func)
         endif
      endif
c     
      if (abs(cfac(10)).gt.eps)then
         if (.not. do_2nd .and..not. do_3rd) then
            call xc_vwn_3(tol_rho, cfac(10),
     &           rho, Amat, nq, ipol, Ec, qwght,
     &           ldew, func)
         else if (.not. do_3rd) then
            call xc_vwn_3_d2(tol_rho, cfac(10),
     &           rho, Amat, Amat2, nq, ipol, Ec, qwght,
     &           ldew, func)
         else
            call xc_vwn_3_d3(tol_rho, cfac(10),
     &           rho, Amat, Amat2, Amat3, nq, ipol, Ec, qwght,
     &           ldew, func)
         endif
      endif
c     
      if (abs(cfac(11)).gt.eps)then
         if (.not. do_2nd .and..not. do_3rd) then
            call xc_vwn_4(tol_rho, cfac(11),
     &           rho, Amat, nq, ipol, Ec, qwght,
     &           ldew, func)
         else if (.not. do_3rd) then
            call xc_vwn_4_d2(tol_rho, cfac(11),
     &           rho, Amat, Amat2, nq, ipol, Ec, qwght,
     &           ldew, func)
         else
            call xc_vwn_4_d3(tol_rho, cfac(11),
     &           rho, Amat, Amat2, Amat3, nq, ipol, Ec, qwght,
     &           ldew, func)
         endif
      endif
c     
      if (abs(cfac(6)).gt.eps)then
         if (.not. do_2nd .and..not. do_3rd) then
            call xc_pw91lda(tol_rho, cfac(6), lcfac(6), nlcfac(6), 
     &           rho, Amat, nq, ipol, Ec, 
     &           qwght, ldew, func)
         else if (.not. do_3rd) then
            call xc_pw91lda_d2(tol_rho, cfac(6), lcfac(6), nlcfac(6),
     &           rho, Amat, Amat2, nq, ipol, Ec, qwght, 
     &           ldew, func)
         else
            call xc_pw91lda_d3(tol_rho, cfac(6), lcfac(6), nlcfac(6),
     &           rho, Amat, Amat2, Amat3, nq, ipol, Ec, qwght, 
     &           ldew, func)
         endif
      endif
c     
      if (abs(cfac(2)).gt.eps)then
         if (.not. do_2nd .and. .not. do_3rd) then
            call xc_lyp88(tol_rho, cfac(2), 
     &           rho, delrho, Amat, Cmat, nq, ipol, Ec, 
     &           qwght, ldew, func)
         else if (.not. do_3rd) then
            call xc_lyp88_d2(tol_rho, cfac(2), 
     &           rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol, Ec, 
     &           qwght, ldew, func)
         else
            call xc_lyp88_d3(tol_rho, cfac(2), 
     &           rho, delrho, Amat, Amat2, Amat3, Cmat, Cmat2, Cmat3,
     &           nq, ipol, Ec, qwght, ldew, func)
         endif
      endif
c     
      if (abs(cfac(3)).gt.eps)then
         if (.not. do_2nd .and..not. do_3rd) then
            call xc_p81(tol_rho, cfac(3), lcfac(3), nlcfac(3), rho, 
     &           Amat, nq, ipol, Ec, qwght, ldew, func)
         else if (.not. do_3rd) then
            call xc_p81_d2(tol_rho, cfac(3), lcfac(3), nlcfac(3), rho, 
     &           Amat, Amat2, nq, ipol, Ec, qwght, ldew, func)
         else
            call xc_p81_d3(tol_rho, cfac(3), lcfac(3), nlcfac(3), rho, 
     &           Amat, Amat2, Amat3, nq, ipol, Ec, qwght, ldew, func)
         endif
      endif
c     
      if (abs(cfac(4)).gt.eps)then
         if (.not. do_2nd .and..not. do_3rd) then
            call xc_perdew86(tol_rho, cfac(4), lcfac(4), nlcfac(4), 
     &           rho, delrho, Amat, Cmat, nq, ipol, 
     &           Ec, qwght, ldew, func)
         else if (.not. do_3rd) then
            call xc_perdew86_d2(tol_rho, cfac(4), lcfac(4), nlcfac(4), 
     &           rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol, 
     &           Ec, qwght, ldew, func)
         else 
            call xc_perdew86_d3(tol_rho, cfac(4), lcfac(4), nlcfac(4),
     &           rho, delrho, Amat, Amat2, Amat3, Cmat, Cmat2, Cmat3, 
     &           nq, ipol, Ec, qwght, ldew, func)
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
c     same is true for revTPSS-variant of it
c
      if (abs(cfac(12)).gt.eps.or.abs(cfac(64)).gt.eps)then
         if (.not. do_2nd .and..not. do_3rd) then
            call xc_cpbe96(tol_rho, cfac, lcfac, nlcfac, rho, 
     &           delrho, Amat, Cmat, nq, ipol, 
     &           Ec, qwght, ldew, func)
         else if (.not. do_3rd) then
            call xc_cpbe96_d2(tol_rho, cfac, lcfac, nlcfac, rho, 
     &           delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol, 
     &           Ec, qwght, ldew, func)
         else
            call xc_cpbe96_d3(tol_rho, cfac, lcfac, nlcfac, rho, 
     &           delrho, Amat, Amat2, Amat3, Cmat, Cmat2, Cmat3,
     &           nq, ipol, Ec, qwght, ldew, func)
         endif
      endif
      if (abs(xfac(7)).gt.eps)then
         if (.not. do_2nd .and..not. do_3rd) then
            call xc_xpbe96('orig',
     T           tol_rho, xfac(7), lxfac(7), nlxfac(7), 
     &           rho, delrho, Amat, Cmat, nq, ipol, 
     &           Ex, qwght,ldew,func)
         else if (.not. do_3rd) then
            call xc_xpbe96_d2('orig',
     T           tol_rho, xfac(7), lxfac(7), nlxfac(7), 
     &           rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol, 
     &           Ex, qwght,ldew,func)
         else
            call xc_xpbe96_d3('orig',
     T           tol_rho, xfac(7), lxfac(7), nlxfac(7), 
     &           rho, delrho, Amat, Amat2, Amat3, Cmat, Cmat2, Cmat3,
     &           nq, ipol, Ex, qwght,ldew,func)
         endif
      endif
      if (abs(xfac(30)).gt.eps)then
         if (.not. do_2nd .and..not. do_3rd) then
            call xc_xpbe96('rpbe',
     T           tol_rho, xfac(30), lxfac(30), nlxfac(30), 
     &           rho, delrho, Amat, Cmat, nq, ipol, 
     &           Ex, qwght,ldew,func)
         else if (.not. do_3rd) then
            call xc_xpbe96_d2('rpbe',
     T           tol_rho, xfac(30), lxfac(30), nlxfac(30), 
     &           rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol, 
     &           Ex, qwght,ldew,func)
         else
            call xc_xpbe96_d3('rpbe',
     T           tol_rho, xfac(30), lxfac(30), nlxfac(30), 
     &           rho, delrho, Amat, Amat2, Amat3, Cmat, Cmat2, Cmat3,
     &           nq, ipol, Ex, qwght,ldew,func)
         endif
      endif
      if (abs(xfac(31)).gt.eps)then
         if (.not. do_2nd .and..not. do_3rd) then
            call xc_xpbe96('revp',
     T           tol_rho, xfac(31), lxfac(31), nlxfac(31), 
     &           rho, delrho, Amat, Cmat, nq, ipol, 
     &           Ex, qwght,ldew,func)
         else if (.not. do_3rd) then
            call xc_xpbe96_d2('revp',
     T           tol_rho, xfac(31), lxfac(31), nlxfac(31), 
     &           rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol, 
     &           Ex, qwght,ldew,func)
         else
            call xc_xpbe96_d3('revp',
     T           tol_rho, xfac(31), lxfac(31), nlxfac(31), 
     &           rho, delrho, Amat, Amat2, Amat3, Cmat, Cmat2, Cmat3, 
     &           nq, ipol, Ex, qwght,ldew,func)
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
c s12g
c
      if (abs(xfac(60)).gt.eps) then
         is12x = 1
         if (.not. do_2nd) then
            call xc_s12x(tol_rho, xfac(60), lxfac(60), nlxfac(60),
     &           rho, delrho, Amat, Cmat, nq, ipol,
     &           Ex, qwght, ldew, func, is12x)
         else
            call xc_s12x_d2(tol_rho, xfac(60), lxfac(60), nlxfac(60),
     &         rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol,
     &           Ex, qwght, ldew, func, is12x)
         endif
      endif
c
c s12h
c
      if (abs(xfac(61)).gt.eps) then
         is12x = 2
         if (.not. do_2nd) then
            call xc_s12x(tol_rho, xfac(61), lxfac(61), nlxfac(61),
     &           rho, delrho, Amat, Cmat, nq, ipol,
     &           Ex, qwght, ldew, func, is12x)
         else
            call xc_s12x_d2(tol_rho, xfac(61), lxfac(61), nlxfac(61),
     &         rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol,
     &           Ex, qwght, ldew, func, is12x)
         endif
      endif
c
c cam-s12g
c
      if (abs(xfac(62)).gt.eps) then
         is12x = 1
         if (.not. do_2nd) then
            call xc_cams12x(tol_rho, xfac(62), lxfac(62), nlxfac(62),
     &           rho, delrho, Amat, Cmat, nq, ipol,
     &           Ex, qwght, ldew, func, is12x)
         else
            call xc_cams12x_d2(tol_rho, xfac(62), lxfac(62), nlxfac(62),
     &         rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol,
     &           Ex, qwght, ldew, func, is12x)
         endif
      endif
c
c cam-s12h
c
      if (abs(xfac(63)).gt.eps) then
         is12x = 2
         if (.not. do_2nd) then
            call xc_cams12x(tol_rho, xfac(63), lxfac(63), nlxfac(63),
     &           rho, delrho, Amat, Cmat, nq, ipol,
     &           Ex, qwght, ldew, func, is12x)
         else
            call xc_cams12x_d2(tol_rho, xfac(63), lxfac(63), nlxfac(63),
     &         rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol,
     &           Ex, qwght, ldew, func, is12x)
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
     &           Ex, qwght,ldew,func,ttau, Mmat)
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
      else if (ncap) then
         if (delta_ac.eq.1.0d99) then
           Amat(1:nq) = Amat(1:nq) - 
     &     0.053805222d0*(1d0 + dsqrt(1d0 - 2d0*e_homo/0.053805222d0))
           if (ipol.eq.2) then
             Amat(nq+1:2*nq) = Amat(nq+1:2*nq) -
     &       0.053805222d0*(1d0 + dsqrt(1d0 - 2d0*e_homo/0.053805222d0))
           endif
         else
           Amat(1:nq) = Amat(1:nq) - delta_ac
           if (ipol.eq.2) Amat(nq+1:2*nq) = Amat(nq+1:2*nq) - delta_ac
         endif
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
     &           Ec, qwght, ldew, func, ttau,Amat,Cmat,Mmat)
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
     &           Ex, qwght,ldew,func,ttau,Mmat)
         else
            call xc_xtpss03_d2()
         endif
      endif
c
c   MVS  meta GGA
c
      if (abs(xfac(64)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xmvs15(tol_rho, xfac(64),
     &           rho, delrho, Amat, Cmat, nq, ipol,
     &           Ex, qwght,ldew,func,ttau,Mmat)
         else
            call xc_xmvs15_d2()
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
     &           Ec, qwght, ldew, func, ttau,Amat,Cmat,Mmat)

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
     &           Ec, qwght, ldew,func, ttau,Amat,Cmat,Mmat,0)
      
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
     &           Ex, qwght,ldew,func,ttau,Mmat,1)
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
     &           Ex, qwght,ldew, func,ttau,Mmat,2)
         else
            call xc_xm05_d2()
         endif
      endif
c
c dlDF    meta GGA Exchange
c      
      if (abs(xfac(32)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xdldf(tol_rho, xfac(32), lxfac(32), nlxfac(32),
     &           rho, delrho, Amat, Cmat, nq, ipol,
     &           Ex, qwght,ldew, func,ttau,Mmat)
         else
            call xc_xdldf_d2()
         endif
      endif
c
c VSXC   meta GGA Exchange
c
      if (abs(xfac(33)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xvs98(tol_rho, xfac(33), lxfac(33), nlxfac(33),
     &           rho, delrho, Amat, Cmat, nq, ipol,
     &           Ex, qwght,ldew, func,ttau,Mmat,1)
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
     &           Ex, qwght,ldew, func,ttau,Mmat,1)
         else
            call xc_xm06_d2()
         endif
      endif
c
c revM06   meta GGA Exchange
c
      if (abs(xfac(68)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xm06(tol_rho, xfac(68), lxfac(68), nlxfac(68),
     &           rho, delrho, Amat, Cmat, nq, ipol,
     &           Ex, qwght,ldew, func,ttau,Mmat,6)
         else
            call xc_xm06_d2()
         endif
      endif
c
c revM06-L   meta GGA Exchange
c
      if (abs(xfac(69)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xm06(tol_rho, xfac(69), lxfac(69), nlxfac(69),
     &           rho, delrho, Amat, Cmat, nq, ipol,
     &           Ex, qwght,ldew, func,ttau,Mmat,5)
         else
            call xc_xm06_d2()
         endif
      endif
c
c M06-HF   meta GGA Exchange
c
      if (abs(xfac(35)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xm06(tol_rho, xfac(35), lxfac(35), nlxfac(35),
     &           rho, delrho, Amat, Cmat, nq, ipol,
     &           Ex, qwght,ldew, func,ttau,Mmat,2)
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
     &           Ex, qwght,ldew, func,ttau,Mmat,3)
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
     &           Ex, qwght,ldew, func,ttau,Mmat,4)
         else
            call xc_xm06_d2()
         endif
      endif
c
c M08-HX   meta GGA Exchange
c
      if (abs(xfac(48)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xm11(tol_rho, xfac(48), lxfac(48), nlxfac(48),
     &           rho, delrho, Amat, Cmat, nq, ipol,
     &           Ex, qwght,ldew, func,ttau,Mmat,1)
         else
            call xc_xm11_d2()
         endif
      endif
c
c M08-SO   meta GGA Exchange
c
      if (abs(xfac(49)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xm11(tol_rho, xfac(49), lxfac(49), nlxfac(49),
     &           rho, delrho, Amat, Cmat, nq, ipol,
     &           Ex, qwght,ldew, func,ttau,Mmat,2)
         else
            call xc_xm11_d2()
         endif
      endif
c
c M11   meta GGA Exchange
c
      if (abs(xfac(50)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xm11(tol_rho, xfac(50), lxfac(50), nlxfac(50),
     &           rho, delrho, Amat, Cmat, nq, ipol,
     &           Ex, qwght,ldew, func,ttau,Mmat,3)
         else
            call xc_xm11_d2()
         endif
      endif
c
c M11-L   meta GGA Exchange
c
      if (abs(xfac(51)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xm11(tol_rho, xfac(51), lxfac(51), nlxfac(51),
     &           rho, delrho, Amat, Cmat, nq, ipol,
     &           Ex, qwght,ldew, func,ttau,Mmat,4)
         else
            call xc_xm11_d2()
         endif
      endif
c
c SOGGA GGA Exchange
c
      if (abs(xfac(52)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xsogga(tol_rho, xfac(48), lxfac(48), nlxfac(48),
     &           rho, delrho, Amat, Cmat, nq, ipol,
     &           Ex, qwght,ldew, func,1)
         else
            call xc_xsogga_d2()
         endif
      endif
c
c SOGGA11 GGA Exchange
c
      if (abs(xfac(53)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xsogga(tol_rho, xfac(49), lxfac(49), nlxfac(49),
     &           rho, delrho, Amat, Cmat, nq, ipol,
     &           Ex, qwght,ldew, func,2)
         else
            call xc_xsogga_d2()
         endif
      endif
c
c SOGGA11-X GGA Exchange
c
      if (abs(xfac(54)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xsogga(tol_rho, xfac(50), lxfac(50), nlxfac(50),
     &           rho, delrho, Amat, Cmat, nq, ipol,
     &           Ex, qwght,ldew, func,3)
         else
            call xc_xsogga_d2()
         endif
      endif
c
      if (abs(xfac(55)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xb86b(tol_rho, xfac(55), lxfac(55), nlxfac(55), 
     &           rho, delrho, Amat, Cmat, nq, ipol, 
     &           Ex, qwght,ldew,func)
         else
            call xc_xb86b_d2(tol_rho, xfac(7), lxfac(7), nlxfac(7), 
     &           rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol, 
     &           Ex, qwght,ldew,func)
         endif
      endif
      if (abs(xfac(56)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xpw86(tol_rho, xfac(56), lxfac(56), nlxfac(56), 
     &           rho, delrho, Amat, Cmat, nq, ipol, 
     &           Ex, qwght,ldew,func)
         else
            call xc_xpw86(tol_rho, xfac(56), lxfac(56), nlxfac(56), 
     &           rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol, 
     &           Ex, qwght,ldew,func)
         endif
      endif
c
c     cm08-hx is special in that the GGA part is
c     defined to be dependent on
c     the E(LDA) ...
c     the decision has been made to use the PW91-LDA as the
c     LDA-correlation.  at present, this LDA  cannot be
c     set by the user
c
      if (abs(cfac(48)).gt.eps)then
         if (.not. do_2nd) then
            call xc_cm11(tol_rho, cfac(48), lcfac(48), nlcfac(48),
     &           rho, delrho,  nq, ipol,
     &           Ec, qwght, ldew,func,ttau,Amat,Cmat,Mmat,1)

         else
            call xc_cm11_d2()
         endif
      endif
c
c     cm08-so is special in that the GGA part is
c     defined to be dependent on
c     the E(LDA) ...
c     the decision has been made to use the PW91-LDA as the
c     LDA-correlation.  at present, this LDA  cannot be
c     set by the user
c
      if (abs(cfac(49)).gt.eps)then
         if (.not. do_2nd) then
            call xc_cm11(tol_rho, cfac(49), lcfac(49), nlcfac(49),
     &           rho, delrho,  nq, ipol,
     &           Ec, qwght, ldew,func,ttau,Amat,Cmat,Mmat,2)

         else
            call xc_cm11_d2()
         endif
      endif
c
c     cm11 is special in that the GGA part is
c     defined to be dependent on
c     the E(LDA) ...
c     the decision has been made to use the PW91-LDA as the
c     LDA-correlation.  at present, this LDA  cannot be
c     set by the user
c
      if (abs(cfac(50)).gt.eps)then
         if (.not. do_2nd) then
            call xc_cm11(tol_rho, cfac(50), lcfac(50), nlcfac(50),
     &           rho, delrho,  nq, ipol,
     &           Ec, qwght, ldew,func,ttau,Amat,Cmat,Mmat,3)

         else
            call xc_cm11_d2()
         endif
      endif
c
c     cm11-l is special in that the GGA part is
c     defined to be dependent on
c     the E(LDA) ...
c     the decision has been made to use the PW91-LDA as the
c     LDA-correlation.  at present, this LDA  cannot be
c     set by the user
c
      if (abs(cfac(51)).gt.eps)then
         if (.not. do_2nd) then
            call xc_cm11(tol_rho, cfac(51), lcfac(51), nlcfac(51),
     &           rho, delrho,  nq, ipol,
     &           Ec, qwght, ldew,func,ttau,Amat,Cmat,Mmat,4)

         else
            call xc_cm11_d2()
         endif
      endif
c
c     csogga is special in that the GGA part is
c     defined to be dependent on
c     the E(LDA) ...
c     the decision has been made to use the PW91-LDA as the
c     LDA-correlation.  at present, this LDA  cannot be
c     set by the user
c
      if (abs(cfac(52)).gt.eps)then
         if (.not. do_2nd) then
            call xc_cpbe96(tol_rho, cfac, lcfac, nlcfac, rho,
     &           delrho, Amat, Cmat, nq, ipol,
     &           Ec, qwght, ldew, func)
         else
            call xc_csogga_d2()
         endif
      endif
c
c     csogga11 is special in that the GGA part is
c     defined to be dependent on
c     the E(LDA) ...
c     the decision has been made to use the PW91-LDA as the
c     LDA-correlation.  at present, this LDA  cannot be
c     set by the user
c
      if (abs(cfac(53)).gt.eps)then
         if (.not. do_2nd) then
            call xc_csogga(tol_rho, cfac(49), lcfac(49), nlcfac(49),
     &           rho, delrho,  nq, ipol,
     &           Ec, qwght, ldew,func,Amat,Cmat,1)

         else
            call xc_csogga_d2()
         endif
      endif
c
c     csogga11-x is special in that the GGA part is
c     defined to be dependent on
c     the E(LDA) ...
c     the decision has been made to use the PW91-LDA as the
c     LDA-correlation.  at present, this LDA  cannot be
c     set by the user
c
      if (abs(cfac(54)).gt.eps)then
         if (.not. do_2nd) then
            call xc_csogga(tol_rho, cfac(50), lcfac(50), nlcfac(50),
     &           rho, delrho,  nq, ipol,
     &           Ec, qwght, ldew,func,Amat,Cmat,2)

         else
            call xc_csogga_d2()
         endif
      endif
c
c LC-BNL 2007 Exchange 
c
      if (abs(xfac(38)).gt.eps)then
        if (.not. do_2nd .and..not. do_3rd) then
          call xc_bnl(tol_rho, xfac(38), lxfac(38), nlxfac(38), 
     &      rho, Amat, nq, ipol, Ex, qwght, ldew, func)
        else if (.not. do_3rd) then
          call xc_bnl_d2(tol_rho, xfac(38), lxfac(38), nlxfac(38), 
     &      rho, Amat, Amat2, nq, ipol, Ex, qwght, ldew, func)
        else
          call xc_bnl_d3(tol_rho, xfac(38), lxfac(38), nlxfac(38), 
     &      rho, Amat, Amat2, Amat3, nq, ipol, Ex, qwght, ldew, func)
        endif
      endif
c
c CAM-B88 Exchange
c
      if (abs(xfac(40)).gt.eps)then
        if (.not. do_2nd .and..not. do_3rd) then
          call xc_camb88(tol_rho, xfac(40), lxfac(40), nlxfac(40),
     &      rho, delrho, Amat, Cmat, nq, ipol,
     &      Ex, qwght,ldew,func)
        else if (.not. do_3rd) then
          call xc_camb88_d2(tol_rho, xfac(40), lxfac(40), nlxfac(40),
     &      rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol,
     &      Ex, qwght,ldew,func)
        else
          call xc_camb88_d3(tol_rho, xfac(40), lxfac(40), nlxfac(40),
     &      rho, delrho, Amat, Amat2, Amat3, Cmat, Cmat2, Cmat3, 
     &      nq, ipol, Ex, qwght,ldew,func)
        endif
      endif
c
c CAM-PBE96 Exchange
c
      if (abs(xfac(41)).gt.eps)then
        if (.not. do_2nd .and..not. do_3rd) then
          call xc_camxpbe96('orig',
     T             tol_rho, xfac(41), lxfac(41), nlxfac(41), 
     &             rho, delrho, Amat, Cmat, nq, ipol, 
     &             Ex, qwght,ldew,func)
        else if (.not. do_3rd) then
          call xc_camxpbe96_d2('orig',
     T             tol_rho, xfac(41), lxfac(41), nlxfac(41), 
     &             rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol, 
     &             Ex, qwght,ldew,func)
        else
          call xc_camxpbe96_d3('orig',
     T             tol_rho, xfac(41), lxfac(41), nlxfac(41), 
     &             rho, delrho, Amat, Amat2, Amat3, Cmat, Cmat2, Cmat3,
     &             nq, ipol, Ex, qwght,ldew,func)
        endif
      endif
c
c CAM-LSD Exchange
c
      if (abs(xfac(42)).gt.eps)then
        if (.not. do_2nd .and..not. do_3rd) then
          call xc_camxlsd(tol_rho, xfac(42), lxfac(42), nlxfac(42), 
     &        rho, Amat, nq, ipol, Ex, qwght, ldew, func)
        else if (.not. do_3rd) then
          call xc_camxlsd_d2(tol_rho, xfac(42), lxfac(42), nlxfac(42), 
     &        rho, Amat, Amat2, nq, ipol, Ex, qwght, .false., func)
        else
          call xc_camxlsd_d3(tol_rho, xfac(42), lxfac(42), nlxfac(42), 
     &        rho, Amat, Amat2, Amat3, nq, ipol, Ex, qwght, .false., 
     &        func)
        endif
      endif
c
c xwpbe exchange: HSE screened exchange
c 
      if (abs(xfac(43)).gt.eps)then
         if (.not. do_2nd .and..not. do_3rd) then
           call xc_xwpbe(tol_rho, xfac(43), lxfac(43), nlxfac(43), 
     &         rho, delrho, Amat, Cmat, nq, ipol, Ex, qwght,ldew,func)
         else if (.not. do_3rd) then
           call xc_xwpbe_d2(tol_rho, xfac(43), lxfac(43), nlxfac(43), 
     &         rho, delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol, Ex,
     &                        qwght,ldew,func)
         else
           call xc_xwpbe_d3(tol_rho, xfac(43), lxfac(43), nlxfac(43), 
     &         rho, delrho, Amat, Amat2, Amat3, Cmat, Cmat2, Cmat3, 
     &         nq, ipol, Ex, qwght, ldew, func)
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
     &           Ec, qwght, ldew,func, ttau,Amat,Cmat,Mmat,1)

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
     &           Ec, qwght, ldew,func,ttau,Amat,Cmat,Mmat,2)

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
     &           Ec, qwght, ldew, func,ttau,Amat,Cmat,Mmat,1)

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
     &           Ec, qwght, ldew,func,ttau,Amat,Cmat,Mmat,2)

         else
            call xc_cm05_d2()
         endif
      endif
c
c     dlDF Correlation
c
c     cdldf is special in that the GGA part is
c     defined to be dependent on
c     the E(LDA) ...
c     the decision has been made to use the PW91-LDA as the
c     LDA-correlation.  at present, this LDA  cannot be
c     set by the user
c
      if (abs(cfac(42)).gt.eps)then
         if (.not. do_2nd) then
            call xc_cdldf(tol_rho, cfac(42), lcfac(42), nlcfac(42),
     &           rho, delrho,  nq, ipol,
     &           Ec, qwght, ldew,func,ttau,Amat,Cmat,Mmat)

         else
            call xc_cdldf_d2()
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
     &           Ec, qwght, ldew,func,ttau,Amat,Cmat,Mmat,1)

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
     &           Ec, qwght, ldew,func,ttau,Amat,Cmat,Mmat,1)

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
     &           Ec, qwght, ldew,func,ttau,Amat,Cmat,Mmat,2)

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
     &           Ec, qwght, ldew,func,ttau,Amat,Cmat,Mmat,3)

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
     &           Ec, qwght, ldew,func,ttau,Amat,Cmat,Mmat,4)

         else
            call xc_cm06_d2()
         endif
      endif
c
c     N12x
c
      if (abs(xfac(45)).gt.eps)then
         if (.not. do_2nd) then
            call xc_xn12(tol_rho, xfac(45), lxfac(45), nlxfac(45),
     &           rho, delrho,  Amat, Cmat, nq, ipol,
     &           Ex, qwght, ldew,func,1)

         else
            call xc_xn12_d2()
         endif
      endif
c
c     revM06 correlation
c
      if (abs(cfac(68)).gt.eps)then
         if (.not. do_2nd) then
            call xc_cm06(tol_rho, cfac(68), lcfac(68), nlcfac(68),
     &           rho, delrho,  nq, ipol,
     &           Ec, qwght, ldew,func,ttau,Amat,Cmat,Mmat,6)

         else
            call xc_cm06_d2()
         endif
      endif
c
c     revM06-L correlation
c
      if (abs(cfac(69)).gt.eps)then
         if (.not. do_2nd) then
            call xc_cm06(tol_rho, cfac(69), lcfac(69), nlcfac(69),
     &           rho, delrho,  nq, ipol,
     &           Ec, qwght, ldew,func,ttau,Amat,Cmat,Mmat,5)

         else
            call xc_cm06_d2()
         endif
      endif
c
c     SCAN meta GGA
c
      if (abs(xfac(66)).gt.eps) then
         if(.not. do_2nd) then
            call xc_xscan('orig', tol_rho, xfac(66), rho, delrho, Amat, 
     &                    Cmat, nq, ipol, Ex, qwght, ldew, func, ttau, 
     &                    Mmat)
         else
            call xc_xscan_d2()
         endif
      endif
c
      if (abs(cfac(66)).gt.eps) then
         if(.not. do_2nd) then
            call xc_cscan('orig', tol_rho, cfac(66), rho, delrho, Amat, 
     &                    Cmat, nq, ipol, Ec, qwght, ldew, func, ttau, 
     &                    Mmat)
         else
            call xc_cscan_d2()
         endif
      endif
c
c     regularized SCAN
c
      if (abs(xfac(71)).gt.eps) then
         if(.not. do_2nd) then
            call xc_xscan('regu', tol_rho, xfac(71), rho, delrho, Amat,
     &                    Cmat, nq, ipol, Ex, qwght, ldew, func, ttau, 
     &                    Mmat)
         else
            call xc_xscan_d2()
         endif
      endif
c
      if (abs(cfac(71)).gt.eps) then
         if(.not. do_2nd) then
            call xc_cscan('regu', tol_rho, cfac(71), rho, delrho, Amat, 
     &                    Cmat, nq, ipol, Ec, qwght, ldew, func, ttau, 
     &                    Mmat)
         else
            call xc_cscan_d2()
         endif
      endif
c
      if (abs(xfac(72)).gt.eps)then
         if (.not. do_2nd .and..not. do_3rd) then
            call xc_xncap(tol_rho, xfac(72), rho,
     &           delrho, Amat, Cmat, nq, ipol,
     &           Ex, qwght, ldew, func)
         else if (.not. do_3rd) then
            call xc_xncap_d2(tol_rho, xfac(72), rho,
     &           delrho, Amat, Amat2, Cmat, Cmat2, nq, ipol,
     &           Ex, qwght, ldew, func)
         else
            call xc_xncap_d3(tol_rho, xfac(72), rho,
     &           delrho, Amat, Amat2, Amat3, Cmat, Cmat2, Cmat3,
     &           nq, ipol, Ex, qwght, ldew, func)
         endif
      endif
c
c     SCAN-L Laplacian meta GGA
c
      if (abs(xfac(67)).gt.eps) then
         if(.not. do_2nd) then
            call xc_xscanl(tol_rho, xfac(67), rho, delrho, laprho,
     &                     Amat, Cmat, Lmat, nq, ipol, Ex, qwght, ldew, 
     &                     func)
         else
            call xc_xscanl_d2()
         endif
      endif
c
      if (abs(cfac(67)).gt.eps) then
         if(.not. do_2nd) then
            call xc_cscanl(tol_rho, cfac(67), rho, delrho, laprho, 
     &                     Amat, Cmat, Lmat, nq, ipol, Ec, qwght, ldew, 
     &                     func)
         else
            call xc_cscanl_d2()
         endif
      endif
c
c     r^2SCAN meta GGA
c
      if (abs(xfac(73)).gt.eps) then
         if(.not. do_2nd) then
            call xc_xr2scan(tol_rho, xfac(73), rho, delrho, Amat, 
     &                    Cmat, nq, ipol, Ex, qwght, ldew, func, ttau, 
     &                    Mmat)
         else
            call xc_xr2scan_d2()
         endif
      endif
c
      if (abs(cfac(73)).gt.eps) then
         if(.not. do_2nd) then
            call xc_cr2scan(tol_rho, cfac(73), rho, delrho, Amat, 
     &                    Cmat, nq, ipol, Ec, qwght, ldew, func, ttau, 
     &                    Mmat)
         else
            call xc_cr2scan_d2()
         endif
      endif
c
c     r^2SCAN-L meta GGA
c
      if (abs(xfac(74)).gt.eps) then
         if(.not. do_2nd) then
            call xc_xr2scanl(tol_rho, xfac(74), rho, delrho, laprho, 
     &                    Amat, Cmat, Lmat, nq, ipol, Ex, qwght, ldew, 
     &                    func)
         else
            call xc_xr2scanl_d2()
         endif
      endif
c
      if (abs(cfac(74)).gt.eps) then
         if(.not. do_2nd) then
            call xc_cr2scanl(tol_rho, cfac(74), rho, delrho, laprho, 
     &                    Amat, Cmat, Lmat, nq, ipol, Ec, qwght, ldew, 
     &                    func)
         else
            call xc_cr2scanl_d2()
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
c Daniel (1-11-13): Added XC-third derivative stuff.  This currently
c doesn't include any functionality for meta-GGAs.
c Daniel (1-17-13): Because we need derivatives of the quadrature
c weights for TDDFT gradients, we don't want to multiply the A 
c and C matrices with the quadrature weights here (which would be
c double counting).  This prevents us from multiplying the
c quadrature weights into the matrices and then dividing them out
c later, which is prone to numerical problems (and slow). 
c Daniel (2-4-13): ldew3 is for the dfxc*(X+Y)*(X+Y) term in the TDDFT 
c gradients.  ldew2 is for the dVxc*P term in the TDDFT gradients.
      if (.not.(ldew2.or.ldew3)) then

        if (.not. do_2nd .and. .not. do_3rd) then
           call setACmat(delrho, Amat, Cmat, qwght, ipol, nq, grad,
     &           (.not. do_2nd), kske, Mmat, kslap, Lmat)
        else if (.not. do_3rd) then
           call setACmat_d2(delrho, Amat, Amat2, Cmat, Cmat2, qwght, 
     &           ipol, nq, grad, (.not. do_2nd), kske, Mmat, Mmat2,
     &           .false.)
        else
           call setACmat_d3(delrho, Amat, Amat2, Amat3, Cmat, Cmat2,
     &           Cmat3, qwght, ipol, nq, grad, (.not. do_3rd), .false.,
     &           .false.)
        endif
      endif
c
      return
      end

cc AJL/Begin/FDE
c
      Subroutine xc_eval_fnl_fde(rho, delrho, Amat, Amat2, Cmat, Cmat2,
     &                       nq, Ex, Ec, qwght, grad, ldew, func,
     &                       do_2nd, ttau, kske, Mmat, Mmat2,
     &                       StericEnergy, do_3rd, Amat3, Cmat3, ldew2)
c
      implicit none
c
#include "cdft.fh"
c xc-second derivative header file
c#include "dft2drv.fh"
c xc-third derivative header file
c#include "dft3drv.fh"
c#include "stdio.fh"
#include "steric.fh"
c
      integer nq
      double precision rho(*)    !< [Input] The electron density
                                 !< \f$\rho^\alpha\f$ and
                                 !< \f$\rho^\beta\f$.
      double precision delrho(*) !< [Input] The electron density gradient
                                 !< \f$\partial\rho^\alpha/\partial r\f$
                                 !< and
                                 !< \f$\partial\rho^\beta/\partial r\f$
c
      double precision Amat(*), Cmat(*)
      double precision Amat2(*),Cmat2(*)
      double precision Amat3(*),Cmat3(*)
      double precision ttau(*) !< [Input] The kinetic energy density
      double precision Mmat(*), Mmat2(*)
c
      double precision Ex !< [Output] The exchange energy
      double precision Ec !< [Output] The correlation energy
      double precision StericEnergy
      double precision qwght(nq), func(nq)
      logical grad, kske
c
      logical do_2nd, do_3rd
      logical ldew, ldew2
c Local
      double precision cfac_temp(numfunc), xfac_temp(numfunc)
      logical lcfac_temp(numfunc), nlcfac_temp(numfunc)
      logical lxfac_temp(numfunc), nlxfac_temp(numfunc)
      integer i
c
      do i=1,numfunc
        cfac_temp(i)   = cfac(i)
        xfac_temp(i)   = xfac(i)
        lcfac_temp(i)  = lcfac(i)
        nlcfac_temp(i) = nlcfac(i)
        lxfac_temp(i)  = lxfac(i)
        nlxfac_temp(i) = nlxfac(i)

        cfac(i)   = cfac_fde(i)
        xfac(i)   = xfac_fde(i)
        lcfac(i)  = lcfac_fde(i)
        nlcfac(i) = nlcfac_fde(i)
        lxfac(i)  = lxfac_fde(i)
        nlxfac(i) = nlxfac_fde(i)
      enddo

      call xc_eval_fnl(rho, delrho, Amat, Amat2, Cmat, Cmat2,
     &                       nq, Ex, Ec, qwght, grad, ldew, func,
     &                       do_2nd, ttau, kske, Mmat, Mmat2,
     &                       0d0, .false., 0d0,
     &                       StericEnergy, do_3rd, Amat3, Cmat3, ldew2)

      do i=1,numfunc
        cfac(i)   = cfac_temp(i)
        xfac(i)   = xfac_temp(i)
        lcfac(i)  = lcfac_temp(i)
        nlcfac(i) = nlcfac_temp(i)
        lxfac(i)  = lxfac_temp(i)
        nlxfac(i) = nlxfac_temp(i)
      enddo

      return
      end
c
cc AJL/End

      subroutine calc_rgamma(ipol,nq,delrho,rgamma)
      implicit none
cinclude "nwxc_param.fh"
#define G_TT 1
#define G_AA 1
#define G_AB 2
#define G_BB 3
      integer ipol !< [Input] The number of spin channels
      integer nq   !< [Input] The number of grid points
      double precision delrho(nq,3,ipol) !< [Input] The density gradient
      double precision rgamma(nq,3)      !< [Output] The density gradient norm
c
      integer iq
c
      if (ipol.eq.1) then
        do iq = 1, nq
          rgamma(iq,G_TT) = delrho(iq,1,1)*delrho(iq,1,1)
     &                    + delrho(iq,2,1)*delrho(iq,2,1)
     &                    + delrho(iq,3,1)*delrho(iq,3,1)
        enddo
      else
        do iq = 1, nq
          rgamma(iq,G_AA) = delrho(iq,1,1)*delrho(iq,1,1)
     &                    + delrho(iq,2,1)*delrho(iq,2,1)
     &                    + delrho(iq,3,1)*delrho(iq,3,1)
          rgamma(iq,G_BB) = delrho(iq,1,2)*delrho(iq,1,2)
     &                    + delrho(iq,2,2)*delrho(iq,2,2)
     &                    + delrho(iq,3,2)*delrho(iq,3,2)
          rgamma(iq,G_AB) = delrho(iq,1,1)*delrho(iq,1,2)
     &                    + delrho(iq,2,1)*delrho(iq,2,2)
     &                    + delrho(iq,3,1)*delrho(iq,3,2)
        enddo
      endif
c
      end
C>
C> @}
c $Id$
