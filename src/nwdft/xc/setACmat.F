      Subroutine setACmat(Amat,Cmat,qwght,
     &                    ipol,nq,GRAD)

C$Id: setACmat.F,v 1.5 1996-10-18 00:53:21 d3h449 Exp $
      implicit none

      integer ipol,nq
      
      Logical GRAD

c--> Sampling Matrices for the XC Potential & Energy

      double precision Amat(ipol,nq),Cmat(ipol,3,nq)

c--> Quadrature Weights

      double precision qwght(nq)
c
      integer ii,jj
c
      if( GRAD )then
        do ii=1,ipol
          do jj=1,nq
            Amat(ii,jj)=Amat(ii,jj)*qwght(jj)
            Cmat(ii,1,jj)=Cmat(ii,1,jj)*qwght(jj)
            Cmat(ii,2,jj)=Cmat(ii,2,jj)*qwght(jj)
            Cmat(ii,3,jj)=Cmat(ii,3,jj)*qwght(jj)
          enddo
        enddo
c
      else
c
        do ii=1,ipol
          do jj=1,nq
            Amat(ii,jj)=Amat(ii,jj)*qwght(jj)
          enddo
        enddo
      end if
      return
      end


