      Subroutine scat_mat(TT,T,Ni,Nj,mi,mj,ifirst,jfirst,ibf,jbf)

C$Id: scat_mat.F,v 1.4 1995-04-09 08:09:26 og845 Exp $
      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)
      dimension T(mi,mj),TT(Ni,Nj),ibf(mi),jbf(mj)

c
c     scatter into matrix
c     
      call dfill(Ni*Nj,0.d0,TT,1)
      do j = 1,mj
        jj=jbf(j)-jfirst+1
        do i =1,mi
          TT((ibf(i)-ifirst+1),jj)=T(i,j)
        enddo
      enddo
      
      return
      end
C
      Subroutine transp_mat(TT,T,Ni,Nj)

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)
      dimension T(Nj,Ni),TT(Ni,Nj)

c
c     transpose a matrix 
c     input TT
c     output T
c     
      do jj=1,Ni
        do ii=1,Nj
          T(ii,jj)=TT(jj,ii)
        enddo
      enddo
      
      return
      end
C
