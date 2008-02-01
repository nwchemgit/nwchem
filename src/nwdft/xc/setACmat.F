c
c     Modified to handle second derivatives while reusing code
c
c     BGJ - 8/98
c
#ifndef SECOND_DERIV
      Subroutine setACmat(delrho, Amat, Cmat, qwght, ipol, nq, GRAD,
     &            xform_Cmat, kske, Mmat)
#else
      Subroutine setACmat_d2(delrho, Amat, Amat2, Cmat, Cmat2, qwght,
     &            ipol, nq, GRAD, xform_Cmat, kske, Mmat, Mmat2)
#endif
c
C$Id: setACmat.F,v 1.17 2008-02-01 22:27:55 niri Exp $
c
      implicit none
c
#include "dft2drv.fh"
#include "stdio.fh"
c !!! BGJ test
#include "bgj.fh"
c !!! BGJ test
c
      integer ipol, nq
      Logical GRAD, xform_Cmat
c
c     Density gradients - used for transforming fnl gamma derivatives
c
      double precision delrho(nq,3,ipol)
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(nq,ipol)
      double precision Cmat(nq,3,ipol)
      double precision Mmat(nq,ipol)
      logical kske    
c 
#ifdef SECOND_DERIV
      double precision Amat2(nq,NCOL_AMAT2)
      double precision Cmat2(nq,NCOL_CMAT2)
      double precision Mmat2(nq,NCOL_MMAT2)
#endif
c
c     Quadrature Weights
c
      double precision qwght(nq)
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
                  if(kske) Mmat(jj,ii) = Mmat(jj,ii)*qwght(jj)
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
c
#if 0
      if (bgj_print() .gt. 1) then
      write(LuOut,*) ' setACmat: AMAT out'
      call output(amat, 1, nq, 1, ipol, nq, ipol, 1)
      if (grad) then
         write(LuOut,*) ' setACmat: CMAT out ',xform_Cmat
         call output(cmat, 1, nq, 1, 3*ipol, nq, 3*ipol, 1)
      endif
      endif
#endif
c
#ifdef SECOND_DERIV
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
      return
      end
c
#ifndef SECOND_DERIV
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
