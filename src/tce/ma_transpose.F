      subroutine ma_transpose(a,n)
      implicit none
      integer i,j,n
      double precision a(n,n)
      double precision swap
      do i = 1,n
        !$omp simd
        do j = i,n
          swap = a(i,j)
          a(i,j) = a(j,i)
          a(j,i) = swap
        enddo
        !$omp end simd
      enddo
      return
      end

!
! 32 works better than 16 or 64
! 22 works best so far
!
#define BLOCKING 22
!
      subroutine ma_transpose_blocked(a,n)
      implicit none
      integer i,j,n,r,ii,jj
      integer chunks,offset
      integer i_ind,j_ind
      integer i_jump,j_jump
      double precision a(n,n)
      double precision temp1(BLOCKING,BLOCKING)
      double precision temp2(BLOCKING,BLOCKING)
      double precision swap
!     blocking parameters
      r = mod(n,BLOCKING)
      chunks = (n-r)/BLOCKING
      offset = chunks*BLOCKING
!     begin blocked main
      do i = 1,chunks
        i_jump = (i-1)*BLOCKING
!       diagonal blocks
        do jj = 1,BLOCKING
          do ii = 1,BLOCKING
            temp1(ii,jj) = a(i_jump+ii,i_jump+jj)
          enddo ! jj
        enddo ! ii
        do jj = 1,BLOCKING
          do ii = 1,BLOCKING
            a(i_jump+ii,i_jump+jj) = temp1(jj,ii)
          enddo ! jj
        enddo ! ii
!       off-diagonal blocks
        do j = i+1,chunks
          j_jump = (j-1)*BLOCKING
!         copy upper triangle block
          do jj = 1,BLOCKING
            do ii = 1,BLOCKING
              temp1(ii,jj) = a(i_jump+ii,j_jump+jj)
            enddo ! jj
          enddo ! ii
!         copy lower triangle block
          do ii = 1,BLOCKING
            do jj = 1,BLOCKING
              temp2(jj,ii) = a(j_jump+jj,i_jump+ii)
            enddo ! jj
          enddo ! ii
!         read transpose of lower triangle block into upper triangle
          do jj = 1,BLOCKING
            do ii = 1,BLOCKING
              a(i_jump+ii,j_jump+jj) = temp2(jj,ii)
            enddo ! jj
          enddo ! ii
!         read transpose of upper triangle block into lower triangle
          do ii = 1,BLOCKING
            do jj = 1,BLOCKING
              a(j_jump+jj,i_jump+ii) = temp1(ii,jj)
            enddo ! jj
          enddo ! ii
!       end blocked main
        enddo ! j
      enddo ! i
!     begin blocked edge
      do i = 1,chunks
        i_jump = (i-1)*BLOCKING
        do ii = 1,BLOCKING
          do j = 1,r
            swap = a(i_jump+ii,offset+j)
            a(i_jump+ii,offset+j) = a(offset+j,i_jump+ii)
            a(offset+j,i_jump+ii) = swap
          enddo ! j
        enddo ! ii
      enddo ! i
!     end non-blocked edge
!     begin non-blocked lower corner
      do i = 1,r
        i_ind = offset+i
        do j = i,r
          j_ind = offset+j
          swap = a(i_ind,j_ind)
          a(i_ind,j_ind) = a(j_ind,i_ind)
          a(j_ind,i_ind) = swap
        enddo
      enddo
!     end non-blocked lower corner
      return
      end
