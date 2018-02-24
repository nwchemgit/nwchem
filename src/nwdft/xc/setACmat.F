c
c     Modified to handle second derivatives while reusing code
c
c     BGJ - 8/98
c
C> \ingroup nwdft_xc
C> @{
C>
C> \brief Do clever stuff to the partial derivatives
C>
C> The density functionals return plain partial derivatives of the
C> of the functionals in addition to the energy. This routine modifies
C> those quantities to get them in a form that is more "convenient"
C> in the construction of various matrices. 
C>
C> In all cases this modification includes multiplying the partial 
C> derivatives for each grid point with the corresponding quadrature
C> weight.
C>
C> If `xform_Cmat` is set to `.TRUE.` then the 1st derivatives with
C> respect to the norm of the electron density gradient are replaced
C> in the closed shell case with
C> \f{eqnarray*}{
C>   Cmat(i,1,1) &=& \frac{\partial\rho(r)}{\partial x} \left(
C>       \frac{\partial f}{\partial\gamma^{\alpha\alpha}} +
C>       \frac{1}{2}\frac{\partial f}{\partial\gamma^{\alpha\beta}}
C>       \right) w_i\\\\
C>   Cmat(i,2,1) &=& \frac{\partial\rho(r)}{\partial y} \left(
C>       \frac{\partial f}{\partial\gamma^{\alpha\alpha}} +
C>       \frac{1}{2}\frac{\partial f}{\partial\gamma^{\alpha\beta}}
C>       \right) w_i\\\\
C>   Cmat(i,3,1) &=& \frac{\partial\rho(r)}{\partial z} \left(
C>       \frac{\partial f}{\partial\gamma^{\alpha\alpha}} +
C>       \frac{1}{2}\frac{\partial f}{\partial\gamma^{\alpha\beta}}
C>       \right) w_i\\\\
C> \f}
C> and in the open shell case with
C> \f{eqnarray*}{
C>   Cmat(i,1,1) &=& \left(2\frac{\partial\rho^\alpha(r)}{\partial x}
C>       \frac{\partial f}{\partial\gamma^{\alpha\alpha}} +
C>       \frac{\partial\rho^\beta(r)}{\partial x}
C>       \frac{\partial f}{\partial\gamma^{\alpha\beta}}\right)w_i \\\\
C>   Cmat(i,2,1) &=& \left(2\frac{\partial\rho^\alpha(r)}{\partial y}
C>       \frac{\partial f}{\partial\gamma^{\alpha\alpha}} +
C>       \frac{\partial\rho^\beta(r)}{\partial y}
C>       \frac{\partial f}{\partial\gamma^{\alpha\beta}}\right)w_i \\\\
C>   Cmat(i,3,1) &=& \left(2\frac{\partial\rho^\alpha(r)}{\partial z} 
C>       \frac{\partial f}{\partial\gamma^{\alpha\alpha}} +
C>       \frac{\partial\rho^\beta(r)}{\partial z}
C>       \frac{\partial f}{\partial\gamma^{\alpha\beta}}\right)w_i \\\\
C>   Cmat(i,1,2) &=& \left(2\frac{\partial\rho^\beta(r)}{\partial x}
C>       \frac{\partial f}{\partial\gamma^{\beta\beta}} +
C>       \frac{\partial\rho^\alpha(r)}{\partial x}
C>       \frac{\partial f}{\partial\gamma^{\alpha\beta}}\right)w_i \\\\
C>   Cmat(i,2,2) &=& \left(2\frac{\partial\rho^\beta(r)}{\partial y}
C>       \frac{\partial f}{\partial\gamma^{\beta\beta}} +
C>       \frac{\partial\rho^\alpha(r)}{\partial y}
C>       \frac{\partial f}{\partial\gamma^{\alpha\beta}}\right)w_i \\\\
C>   Cmat(i,3,2) &=& \left(2\frac{\partial\rho^\beta(r)}{\partial z} 
C>       \frac{\partial f}{\partial\gamma^{\beta\beta}} +
C>       \frac{\partial\rho^\alpha(r)}{\partial z}
C>       \frac{\partial f}{\partial\gamma^{\alpha\beta}}\right)w_i \\\\
C> \f}
C> The grid point is indexed with `i`, the density functional is
C> referred to with \f$f\f$, and \f$w\f$ is the quadrature weight.
C>
#if !defined SECOND_DERIV && !defined THIRD_DERIV
      Subroutine setACmat(delrho, Amat, Cmat, qwght, ipol, nq, GRAD,
     &            xform_Cmat, kske, Mmat, kslap, Lmat)
#elif defined(SECOND_DERIV) && !defined THIRD_DERIV
c Second derivatives are not yet implemented for the meta-GGAs
      Subroutine setACmat_d2(delrho, Amat, Amat2, Cmat, Cmat2, qwght,
     &            ipol, nq, GRAD, xform_Cmat, kske, Mmat, Mmat2,
     &            kslap)
#else
c Third derivatives are not yet implemented for the meta-GGAs
      Subroutine setACmat_d3(delrho, Amat, Amat2, Amat3, Cmat, Cmat2,
     &            Cmat3, qwght, ipol, nq, GRAD, xform_Cmat, kske, kslap)
#endif
c
C$Id$
c
      implicit none
c
#include "dft2drv.fh"
#include "dft3drv.fh"
#include "stdio.fh"
c !!! BGJ test
#include "bgj.fh"
c !!! BGJ test
c
      integer ipol  !< [Input] The number of spin channels
      integer nq    !< [Input] The number of grid points
      Logical GRAD  !< [Input] If `.TRUE.` the density functional is
                    !< a GGA
      logical xform_Cmat !< [Input] If `.TRUE.` do the elaborate 
                         !< transformation
c
c     Density gradients - used for transforming fnl gamma derivatives
c
      double precision delrho(nq,3,ipol) !< [Input] The gradient of the
                                         !< density
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(nq,ipol) !< [In/Output] \f$\partial f/\partial\rho^\sigma\f$
      double precision Cmat(nq,3,ipol) !< [In/Output] \f$\partial f/\partial\gamma^{\sigma\sigma'}\f$
      double precision Mmat(nq,ipol) !< [In/Output] \f$\partial f/\partial\tau^\sigma\f$
      logical kske    !< [Input] If `.TRUE.` the density functional is
                      !< a meta-GGA
c
      logical kslap  !< [Input] If `.TRUE.` the density functional is
                     !< a laplacian-dependent meta-GGA
      double precision Lmat(nq,ipol)   !< [In/Output] \f$\partial f/\partial\lap^\sigma\f$

c
#ifdef SECOND_DERIV
      double precision Amat2(nq,NCOL_AMAT2) !< [In/Output] Similar to
      !< `Amat` but 2nd order derivatives
      double precision Cmat2(nq,NCOL_CMAT2) !< [In/Output] Similar to
      !< `Cmat` but 2nd order derivatives
      double precision Mmat2(nq,NCOL_MMAT2) !< [In/Output] Similar to
      !< `Mmat` but 2nd order derivatives
#endif
c
#ifdef THIRD_DERIV
      double precision Amat3(nq,NCOL_AMAT3)
      double precision Cmat3(nq,NCOL_CMAT3)
#endif
c
c     Quadrature Weights
c
      double precision qwght(nq) !< [Input] The quadrature weights
c
      integer ii, jj
c
c
c     Transform derivatives of functional with respect to gammas
c     to derivatives of functional with respect to density gradient
c
      if (GRAD .and. xform_Cmat) then
         call transform_Cmat(delrho, Cmat, ipol, nq)
      endif
c
c     Combine derivatives of functional with quadrature weights
c
      if (GRAD)then
         if (xform_Cmat) then
            do ii = 1, ipol
               do jj = 1, nq
                  Amat(jj,ii) = Amat(jj,ii)*qwght(jj)
                  Cmat(jj,1,ii) = Cmat(jj,1,ii)*qwght(jj)
                  Cmat(jj,2,ii) = Cmat(jj,2,ii)*qwght(jj)
                  Cmat(jj,3,ii) = Cmat(jj,3,ii)*qwght(jj)
               enddo
            enddo
         else
            if (ipol .eq. 1) then
               do jj = 1, nq
                  Amat(jj,1) = Amat(jj,1)*qwght(jj)
                  Cmat(jj,D1_GAA,1) = Cmat(jj,D1_GAA,1)*qwght(jj)
                  Cmat(jj,D1_GAB,1) = Cmat(jj,D1_GAB,1)*qwght(jj)
               enddo
            else
               do jj = 1, nq
                  Amat(jj,1) = Amat(jj,1)*qwght(jj)
                  Amat(jj,2) = Amat(jj,2)*qwght(jj)
                  Cmat(jj,D1_GAA,1) = Cmat(jj,D1_GAA,1)*qwght(jj)
                  Cmat(jj,D1_GAB,1) = Cmat(jj,D1_GAB,1)*qwght(jj)
                  Cmat(jj,D1_GBB,1) = Cmat(jj,D1_GBB,1)*qwght(jj)
               enddo
            endif
         endif
      else
         do ii = 1, ipol
            do jj = 1, nq
               Amat(jj,ii) = Amat(jj,ii)*qwght(jj)
            enddo
         enddo
      endif
      if (kske) then
         do ii = 1, ipol
            do jj = 1, nq
               Mmat(jj,ii) = Mmat(jj,ii)*qwght(jj)
            enddo
         enddo
      endif
      if (kslap) then
         do ii = 1, ipol
            do jj = 1, nq
               Lmat(jj,ii) = Lmat(jj,ii)*qwght(jj)
            enddo
         enddo
      endif
c
#if 0
      if (bgj_print() .gt. 1) then
      write(LuOut,*) ' setACmat: AMAT out' 
      call output(amat, 1, nq, 1, ipol, nq, ipol, 1)
      if (grad) then
         write(LuOut,*) ' setACmat: CMAT out ',xform_Cmat
         call output(cmat, 1, nq, 1, 3*ipol, nq, 3*ipol, 1)
      endif
      if (kske) then
         write(LuOut,*) ' setACmat: MMAT out'
         call output(mmat, 1, nq, 1, ipol, nq, ipol, 1)
      endif
      endif
#endif
c
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      do ii = 1, NCOL_AMAT2
         do jj = 1, nq
            Amat2(jj,ii) = Amat2(jj,ii)*qwght(jj)
         enddo
      enddo
      if (GRAD)then
         do ii = 1, NCOL_CMAT2
            do jj = 1, nq
               Cmat2(jj,ii) = Cmat2(jj,ii)*qwght(jj)
            enddo
         enddo
      endif
      if (kske)then
         do ii = 1, NCOL_MMAT2
            do jj = 1, nq
               Mmat2(jj,ii) = Mmat2(jj,ii)*qwght(jj)
            enddo
         enddo
      endif
#if 0
      if (bgj_print() .gt. 1) then
      write(LuOut,*) ' setACmat_d2: AMAT2 out'
      call output(amat2, 1, nq, 1, NCOL_AMAT2, nq, NCOL_AMAT2, 1)
      if (grad) then
         write(LuOut,*) ' setACmat_d2: CMAT2 out'
         call output(cmat2, 1, nq, 1, NCOL_CMAT2, nq, NCOL_CMAT2, 1)
      endif
      endif
#endif
#endif
c
#ifdef THIRD_DERIV
      do ii = 1, NCOL_AMAT3
         do jj = 1, nq
            Amat3(jj,ii) = Amat3(jj,ii)*qwght(jj)
         enddo
      enddo
      if (GRAD) then
         do ii = 1, NCOL_CMAT3
            do jj = 1, nq
               Cmat3(jj,ii) = Cmat3(jj,ii)*qwght(jj)
            enddo
         enddo
      endif
#endif
      return
      end
c
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
c
c     Transform Cmat from gamma form to delrho form
c
c     BGJ - 8/98
c
      Subroutine transform_Cmat(delrho, Cmat, ipol, nq)
c
      implicit none
c
#include "dft2drv.fh"
c
      integer ipol, nq
      double precision delrho(nq,3,ipol), Cmat(nq,3,ipol)
c
      integer n
      double precision gaa, gab, gbb
c
      if (ipol .eq. 1) then
         do n = 1, nq
c     Must account for delrho being total density gradient, not alpha
            gaa = Cmat(n,D1_GAA,1) + Cmat(n,D1_GAB,1)*0.5d0
            Cmat(n,1,1) = delrho(n,1,1)*gaa
            Cmat(n,2,1) = delrho(n,2,1)*gaa
            Cmat(n,3,1) = delrho(n,3,1)*gaa
         enddo
      else
         do n = 1, nq
            gaa = Cmat(n,D1_GAA,1)
            gab = Cmat(n,D1_GAB,1)
            gbb = Cmat(n,D1_GBB,1)
            Cmat(n,1,1) = 2d0*delrho(n,1,1)*gaa + delrho(n,1,2)*gab
            Cmat(n,2,1) = 2d0*delrho(n,2,1)*gaa + delrho(n,2,2)*gab
            Cmat(n,3,1) = 2d0*delrho(n,3,1)*gaa + delrho(n,3,2)*gab
            Cmat(n,1,2) = 2d0*delrho(n,1,2)*gbb + delrho(n,1,1)*gab
            Cmat(n,2,2) = 2d0*delrho(n,2,2)*gbb + delrho(n,2,1)*gab
            Cmat(n,3,2) = 2d0*delrho(n,3,2)*gbb + delrho(n,3,1)*gab
         enddo
      endif
c
      return
      end
#endif
c
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "setACmat.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "setACmat.F"
#endif
C> @}
