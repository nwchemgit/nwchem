      Subroutine scat_mat(TT,T,Ni,Nj,mi,mj,ifirst,jfirst,ibf,jbf)

C$Id: scat_mat.F,v 1.6 2001-11-16 02:12:42 edo Exp $
      implicit none
      integer mi,mj,Ni,Nj
      integer ifirst,jfirst
      integer ibf(mi),jbf(mj)
      double precision T(mi,mj),TT(Ni,Nj)
c
      integer i,j,jj

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

      implicit none
      integer Ni,Nj
      double precision T(Nj,Ni),TT(Ni,Nj)
c
      integer jj,ii

c
c     transpose a matrix 
c     input TT
c     output T
c     
!      do jj=1,Ni,2
!      r1 = iand(max0(rows,0),3)
!        do ii=1,Nj,2
!          T(ii,jj)=TT(jj,ii)
!          T(ii+1,jj)=TT(jj,ii+1)
!          T(ii,jj+1)=TT(jj+1,ii)
!          T(ii+1,jj+1)=TT(jj+1,ii+1)
!        enddo
!      enddo
      do jj=1,Ni
        do ii=1,Nj
          T(ii,jj)=TT(jj,ii)
        enddo
      enddo
      
      return
      end
C
