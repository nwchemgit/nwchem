c
c     Modified to handle second derivatives while reusing code
c
c     BGJ - 8/98
c
#ifndef SECOND_DERIV
      Subroutine setACmat(delrho, Amat, Cmat, qwght, ipol, nq, GRAD)
#else
      Subroutine setACmat_d2(delrho, Amat, Amat2, Cmat, Cmat2, qwght,
     &                       ipol, nq, GRAD)
#endif
c
C$Id: setACmat.F,v 1.10 1999-02-15 16:47:13 bjohnson Exp $
c
      implicit none
c
#include "dft2drv.fh"
c
      integer ipol, nq
      Logical GRAD
c
c     Density gradients - used for transforming fnl gamma derivatives
c
      double precision delrho(nq,3,ipol)
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(nq,ipol),Cmat(nq,3,ipol)
#ifdef SECOND_DERIV
      double precision Amat2(nq,NCOL_AMAT2),Cmat2(nq,NCOL_CMAT2)
#endif
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
      integer ii, jj
      double precision gaa, gab, gbb
c
c      write(6,*) ' setACmat: AMAT '
c      call output(amat, 1, nq, 1, ipol, nq, ipol, 1)
c      write(6,*) ' setACmat: CMAT '
c      call output(cmat, 1, 3*nq, 1, ipol, 3*nq, ipol, 1)
c      write(6,*) ' setACmat: qwght '
c      call output(qwght, 1, nq, 1, 1, nq, 1, 1)
c
c     Transform derivatives of functional with respect to gammas
c     to derivatives of functional with respect to density gradient
c
      if (GRAD) then
         if (ipol .eq. 1) then
            do jj = 1, nq
c     Must account for delrho being total density gradient, not alpha
               gaa = Cmat(jj,D1_GAA,1) + Cmat(jj,D1_GAB,1)*0.5d0
               Cmat(jj,1,1) = delrho(jj,1,1)*gaa
               Cmat(jj,2,1) = delrho(jj,2,1)*gaa
               Cmat(jj,3,1) = delrho(jj,3,1)*gaa
            enddo
         else
            do jj = 1, nq
               gaa = Cmat(jj,D1_GAA,1)
               gab = Cmat(jj,D1_GAB,1)
               gbb = Cmat(jj,D1_GBB,1)
               Cmat(jj,1,1) = 2*delrho(jj,1,1)*gaa + delrho(jj,1,2)*gab
               Cmat(jj,2,1) = 2*delrho(jj,2,1)*gaa + delrho(jj,2,2)*gab
               Cmat(jj,3,1) = 2*delrho(jj,3,1)*gaa + delrho(jj,3,2)*gab
               Cmat(jj,1,2) = 2*delrho(jj,1,2)*gbb + delrho(jj,1,1)*gab
               Cmat(jj,2,2) = 2*delrho(jj,2,2)*gbb + delrho(jj,2,1)*gab
               Cmat(jj,3,2) = 2*delrho(jj,3,2)*gbb + delrho(jj,3,1)*gab
            enddo
         endif
      endif
c
c     Combine derivatives of functional with quadrature weights
c
      if (GRAD)then
         do ii = 1, ipol
            do jj = 1, nq
               Amat(jj,ii) = Amat(jj,ii)*qwght(jj)
               Cmat(jj,1,ii) = Cmat(jj,1,ii)*qwght(jj)
               Cmat(jj,2,ii) = Cmat(jj,2,ii)*qwght(jj)
               Cmat(jj,3,ii) = Cmat(jj,3,ii)*qwght(jj)
            enddo
         enddo
      else
         do ii = 1, ipol
            do jj = 1, nq
               Amat(jj,ii) = Amat(jj,ii)*qwght(jj)
            enddo
         enddo
      endif
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
               Cmat2(jj,ii) = Cmat2(jj,ii)*qwght(jj)
               Cmat2(jj,ii) = Cmat2(jj,ii)*qwght(jj)
            enddo
         enddo
      endif
c      write(6,*) ' setACmat_d2: AMAT '
c      call output(amat, 1, nq, 1, ipol, nq, ipol, 1)
c      write(6,*) ' setACmat_d2: AMAT2 '
c      call output(amat2, 1, nq, 1, NCOL_AMAT2, nq, NCOL_AMAT2 , 1)
#endif
      return
      end
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "setACmat.F"
#endif
