      Subroutine scat_mat(TT,T,Ni,Nj,mi,mj,ifirst,jfirst,ibf,jbf)

C$Id$
      implicit none

#include "stdio.fh"

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

      subroutine scat_matup(TT,T,Ni,Nj,mi,mj,
     I     ifirst,jfirst,ibf,jbf)
      implicit none

      integer mi,mj,Ni,Nj
      integer ifirst,jfirst
      integer ibf(mi),jbf(mj)
      double precision T(mi,mj),TT(Ni,Nj)
c
      integer i,j,jj
c
c     scatter into matrix
c     prepare for ga_symmetrize
c     
      call dcopy(Ni*Nj,0.d0,0,TT,1)
      if(ifirst.eq.jfirst) then
         do j =1,mi
            TT((ibf(j)-jfirst+1),
     A           ibf(j)-jfirst+1)=T(j,j)
            do i = j+1,mi
               TT((ibf(i)-ifirst+1),
     A              ibf(j)-jfirst+1)=T(i,j)*2d0
            enddo
         enddo
      else
         do j = 1,mj
            jj=jbf(j)-jfirst+1
            do i =1,mi
               TT((ibf(i)-ifirst+1),jj)=T(i,j)*2d0
            enddo
         enddo
      endif
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
      Subroutine scat_mat3(n3d,i3d,
     .     TT,T,Ni,Nj,mi,mj,ifirst,jfirst,ibf,jbf)
      implicit none
      integer n3d,i3d
      integer mi,mj,Ni,Nj
      integer ifirst,jfirst
      integer ibf(mi),jbf(mj)
      double precision T(mi,mj),TT(n3d,Ni,Nj)
c
      integer i,j,jj

c
c     scatter into matrix
c     
      do j = 1,mj
        jj=jbf(j)-jfirst+1
        do i =1,mi
          TT(i3d,(ibf(i)-ifirst+1),jj)=T(i,j)
        enddo
      enddo
      
      return
      end

      Subroutine transp_mat3(n3d,TT,T,Ni,Nj)

      implicit none
      integer n3d
      integer Ni,Nj
      double precision T(n3d,Nj,Ni),TT(n3d,Ni,Nj)
c
      integer jj,ii,i3d

c
c     transpose a matrix 
c     input TT
c     output T
c     
      do jj=1,Ni
        do ii=1,Nj
           do i3d=1,n3d
              T(i3d,ii,jj)=TT(i3d,jj,ii)
           enddo
        enddo
      enddo
      
      return
      end
