c
c     Modified to handle second derivatives while reusing code
c
c     BGJ - 8/97
c
#ifndef SECOND_DERIV
      Subroutine setACmat(Amat, Cmat, qwght, ipol, nq, GRAD)
#else
#include "dft2drv.fh"
      Subroutine setACmat_d2(Amat, Amat2, Cmat, Cmat2, qwght, ipol, nq,
     &                       GRAD)
#endif
c
C$Id: setACmat.F,v 1.8 1998-09-09 13:32:35 bjohnson Exp $
c
      implicit none
      integer ipol, nq
      Logical GRAD
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
c
c      write(6,*) ' setACmat: AMAT '
c      call output(amat, 1, nq, 1, ipol, nq, ipol, 1)
c      write(6,*) ' setACmat: CMAT '
c      call output(cmat, 1, 3*nq, 1, ipol, 3*nq, ipol, 1)
c      write(6,*) ' setACmat: qwght '
c      call output(qwght, 1, nq, 1, 1, nq, 1, 1)
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
