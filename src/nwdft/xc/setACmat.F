      Subroutine setACmat(Amat, Cmat, qwght, ipol, nq, GRAD)
c
C$Id: setACmat.F,v 1.7 1998-04-23 00:10:53 d3h449 Exp $
c
      implicit none
      integer ipol, nq
      Logical GRAD
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision Amat(nq,ipol),Cmat(nq,3,ipol)
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
      return
      end


