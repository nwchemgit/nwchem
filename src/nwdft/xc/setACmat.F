      Subroutine setACmat(Amat, Cmat, qwght, ipol, nq, GRAD)
c
C$Id: setACmat.F,v 1.6 1997-03-28 18:45:45 d3h449 Exp $
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


