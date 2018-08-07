cc AJL/Begin/FDE
c
c     Wrapper routine for evaluating Ts functional and derivatives,
c     and combining them with quadrature weights.
c     Copied from xc_eval_fnl, with minor adaptation
c
c     AJL (10/2014)
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
C> - the kinetic energy energy
C>
C> - the kinetic energy 1st order partial derivatives wrt density
C>
C> as desired.
C>
      Subroutine ts_eval_fnl(rho, delrho, Amat, Amat2, Cmat, Cmat2,
     &                       nq, Ts, qwght, grad, ldew, func,
     &                       do_2nd, ttau, kske, Mmat, Mmat2,
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
c#include "steric.fh"
c
cc AJL/Unused
c      integer nq, is12x
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
      double precision Ts !< [Output] The exchange energy
      double precision qwght(nq), func(nq)
cc AJL/Unused
c      logical grad, kske, dohcth
      logical grad, kske
c
c      character*4 whichf

      logical do_2nd, do_3rd
      logical ldew, ldew2, ldew3
c
      double precision eps
c AJL/Unused
c      integer nx,nc,dumi
      parameter (eps=1.e-8)
c
c     Initialize the Ts potential and energy sampling matrices.
c
      call dfill(ipol*nq, 0.d0, Amat, 1)
      if(grad) call dfill(3*nq*ipol, 0.d0, Cmat, 1)
      if(kske) call dfill(nq*ipol, 0.d0, Mmat, 1)
      if (do_2nd) then
         call dfill(nq*NCOL_AMAT2, 0.d0, Amat2, 1)
         if (grad) call dfill(nq*NCOL_CMAT2, 0.d0, Cmat2, 1)
         if (kske) call dfill(nq*NCOL_MMAT2, 0.d0, Mmat2, 1)
      endif
c
c     Initialize the 3rd derivatives.  Note that we may need 
c     the 2nd derivative matrices in some cases as well, so
c     those are initialized also.
      if (do_3rd) then
         call dfill(nq*NCOL_AMAT2, 0.0d0, Amat2, 1)
         call dfill(nq*NCOL_AMAT3, 0.0d0, Amat3, 1)
         if (grad) call dfill(nq*NCOL_CMAT2, 0.0d0, Cmat2, 1)         
         if (grad) call dfill(nq*NCOL_CMAT3, 0.0d0, Cmat3, 1)
      endif      
c
c     This prevents the Ts-kernel from being multiplied
c     with the quadrature weights in TDDFT gradients. 
c
      ldew3 = (do_3rd.and.ldew)
      if (ldew.or.ldew2.or.ldew3) call dfill(nq, 0.d0, func, 1)
c      if (ldew.or.ldew3) call dfill(nq, 0.d0, func, 1)
c      if (ldew) call dfill(nq, 0.d0, func, 1)
c
      if (tsfac(1).gt.eps) then
         if (.not. do_2nd .and. .not. do_3rd) then
            call ts_tf(tol_rho, tsfac(1), 0.0d0, 0.0d0, rho,
     &           Amat, nq, ipol, Ts, qwght,
     &           ldew, func)
c         else if (.not. do_3rd) then
c            call ts_tf_d2(tol_rho, 1.d0, 0.d0, 0.d0, rho,
c     &           Amat, Amat2, nq, ipol, Ts, qwght,
c     &           .false., func)
c         else
c            call ts_tf_d3(tol_rho, 1.d0, 0.d0, 0.d0, rho,
c     &           Amat, Amat2, Amat3, nq, ipol, Ts, qwght,
c     &           .false., func)
         endif
      endif
 
      if(tsfac(2).gt.eps) then 
         if (.not. do_2nd .and. .not. do_3rd) then
            call ts_vw(tol_rho, tsfac(2), 0.0d0, 0.0d0, rho,
     &           delrho, Amat, Cmat, nq, ipol, 
     &           Ts, qwght, ldew, func)
!         else if (.not. do_3rd) then
!            call ts_vw_d2(tol_rho, 1.d0, 0.d0, 0.d0, rho,
!     &           Amat, Amat2, nq, ipol, Ts, qwght,
!     &           .false., func)
!         else
!            call ts_vw_d3(tol_rho, 1.d0, 0.d0, 0.d0, rho,
!     &           Amat, Amat2, Amat3, nq, ipol, Ts, qwght,
!     &           .false., func)
         endif
      endif
c
c     Calculate the steric energy
c
c      if (lsteric) then
c        StericEnergy = 0.d0
c        call steric_energy(tol_rho,xfac(1),rho,delrho,nq,
c     &   qwght,ipol,StericEnergy)
c      endif
c
c     Combine with quadrature weights
c
      if (.not.(ldew2.or.ldew3)) then

        if (.not. do_2nd .and. .not. do_3rd) then
           call setACmat(delrho, Amat, Cmat, qwght, ipol, nq, grad,
     &           (.not. do_2nd), kske, Mmat, .false., 0d0)
c        else if (.not. do_3rd) then
c           call setACmat_d2(delrho, Amat, Amat2, Cmat, Cmat2, qwght, 
c     &           ipol, nq, grad, (.not. do_2nd), kske, Mmat, Mmat2, 
c     &           .false.)
c        else
c           call setACmat_d3(delrho, Amat, Amat2, Amat3, Cmat, Cmat2,
c     &           Cmat3, qwght, ipol, nq, grad, (.not. do_3rd), .false., 
c     &.          false.)
        endif
      endif
c
      return
      end 

cC>
cC> @}
cc AJL/End
