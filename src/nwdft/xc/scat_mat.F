      Subroutine scat_mat(TT,T,Ni,Nj,mi,mj,ifirst,jfirst,ibf,jbf)

C$Id: scat_mat.F,v 1.3 1995-02-25 18:37:19 d3e129 Exp $
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
      Subroutine gat_mat(T,TT,Ni,Nj,mi,mj,ifirst,jfirst,ibf,jbf)

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)
      dimension T(mi,mj),TT(Ni,Nj),ibf(mi),jbf(mj)

c
c     gather into matrix 
c     
      do j = 1,mj
        jj=jbf(j)-jfirst+1
        do i =1,mi
          T(i,j)=TT((ibf(i)-ifirst+1),jj)
        enddo
      enddo
      
      return
      end
