      Subroutine setACmat(Amat,Cmat,qwght,
     &                    ipol,nq,GRAD)

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Logical GRAD

c--> Sampling Matrices for the XC Potential & Energy

      Dimension Amat(ipol,nq),Cmat(ipol,3,nq)

c--> Quadrature Weights

      Dimension qwght(nq)

      if( GRAD )then
      do ii=1,ipol

        do jj=1,nq
          
          Amat(ii,jj)=Amat(ii,jj)*qwght(jj)
          Cmat(ii,1,jj)=Cmat(ii,1,jj)*qwght(jj)
          Cmat(ii,2,jj)=Cmat(ii,2,jj)*qwght(jj)
          Cmat(ii,3,jj)=Cmat(ii,3,jj)*qwght(jj)
        
        enddo
      
      enddo
      else
      do ii=1,ipol

        do jj=1,nq
          
          Amat(ii,jj)=Amat(ii,jj)*qwght(jj)
        
        enddo
      
      enddo
      end if
      
      return
      end


