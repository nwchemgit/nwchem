c
c 32 works better than 16 or 64
c 22 works best so far
c
#define BLOCKING 22
c
c 8 works better 4
c 11 works well with 22
c
#define UNROLLING 11
c
      subroutine ma_transpose_blocked(a,n)
c
      implicit none
      integer i,j,n,r,ii,jj
      integer chunks,offset
      integer i_ind,j_ind
      integer i_jump,j_jump
      double precision a(n,n)
      double precision temp1(BLOCKING,BLOCKING)
      double precision temp2(BLOCKING,BLOCKING)
      double precision swap
c
c     blocking parameters
c
      r = mod(n,BLOCKING)
      chunks = ( n - r )/BLOCKING
      offset = chunks*BLOCKING
c
c     begin blocked main
c
!DEC$ prefetch a
!DEC$ loop count min(0), max(200), avg(32)
      do i = 1,chunks
        i_jump = (i-1)*BLOCKING
c
c       diagonal blocks
c
!DEC$ prefetch a
!DEC$ loop count min(BLOCKING), max(BLOCKING), avg(BLOCKING)
!DEC$ unroll(UNROLLING)
!DEC$ vector always
!DEC$ ivdep
        do jj = 1,BLOCKING
!DEC$ prefetch a
!DEC$ loop count min(BLOCKING), max(BLOCKING), avg(BLOCKING)
!DEC$ unroll(UNROLLING)
!DEC$ vector always
!DEC$ ivdep
          do ii = 1,BLOCKING
c
            temp1(ii,jj) = a(i_jump+ii,i_jump+jj)
c
          enddo ! jj
        enddo ! ii
!DEC$ prefetch a
!DEC$ loop count min(BLOCKING), max(BLOCKING), avg(BLOCKING)
!DEC$ unroll(UNROLLING)
!DEC$ vector always
!DEC$ ivdep
        do jj = 1,BLOCKING
!DEC$ prefetch a
!DEC$ loop count min(BLOCKING), max(BLOCKING), avg(BLOCKING)
!DEC$ unroll(UNROLLING)
!DEC$ vector always
!DEC$ ivdep
          do ii = 1,BLOCKING
c
            a(i_jump+ii,i_jump+jj) = temp1(jj,ii)
c
          enddo ! jj
        enddo ! ii
c
c       off-diagonal blocks
c
!DEC$ prefetch a
!DEC$ loop count min(0), max(200), avg(32)
        do j = i+1,chunks
          j_jump = (j-1)*BLOCKING
c
c         copy upper triangle block
c
!DEC$ prefetch a
!DEC$ loop count min(BLOCKING), max(BLOCKING), avg(BLOCKING)
!DEC$ unroll(UNROLLING)
!DEC$ vector always
!DEC$ ivdep
          do jj = 1,BLOCKING
!DEC$ prefetch a
!DEC$ loop count min(BLOCKING), max(BLOCKING), avg(BLOCKING)
!DEC$ unroll(UNROLLING)
!DEC$ vector always
!DEC$ ivdep
            do ii = 1,BLOCKING
c
              temp1(ii,jj) = a(i_jump+ii,j_jump+jj)
c
            enddo ! jj
          enddo ! ii
c
c         copy lower triangle block
c
!DEC$ prefetch a
!DEC$ prefetch temp2
!DEC$ loop count min(BLOCKING), max(BLOCKING), avg(BLOCKING)
!DEC$ unroll(UNROLLING)
!DEC$ vector always
!DEC$ ivdep
          do ii = 1,BLOCKING
!DEC$ prefetch a
!DEC$ prefetch temp2
!DEC$ loop count min(BLOCKING), max(BLOCKING), avg(BLOCKING)
!DEC$ unroll(UNROLLING)
!DEC$ vector always
!DEC$ ivdep
            do jj = 1,BLOCKING
c
              temp2(jj,ii) = a(j_jump+jj,i_jump+ii)
c
            enddo ! jj
          enddo ! ii
c
c         read transpose of lower triangle block into upper triangle
c
!DEC$ prefetch a
!DEC$ prefetch temp2
!DEC$ loop count min(BLOCKING), max(BLOCKING), avg(BLOCKING)
!DEC$ unroll(UNROLLING)
!DEC$ vector always
!DEC$ ivdep
          do jj = 1,BLOCKING
!DEC$ prefetch a
!DEC$ prefetch temp2
!DEC$ loop count min(BLOCKING), max(BLOCKING), avg(BLOCKING)
!DEC$ unroll(UNROLLING)
!DEC$ vector always
!DEC$ ivdep
            do ii = 1,BLOCKING
c
              a(i_jump+ii,j_jump+jj) = temp2(jj,ii)
c
            enddo ! jj
          enddo ! ii
c
c         read transpose of upper triangle block into lower triangle
c
!DEC$ prefetch a
!DEC$ prefetch temp1
!DEC$ loop count min(BLOCKING), max(BLOCKING), avg(BLOCKING)
!DEC$ unroll(UNROLLING)
!DEC$ vector always
!DEC$ ivdep
          do ii = 1,BLOCKING
!DEC$ prefetch a
!DEC$ prefetch temp1
!DEC$ loop count min(BLOCKING), max(BLOCKING), avg(BLOCKING)
!DEC$ unroll(UNROLLING)
!DEC$ vector always
!DEC$ ivdep
            do jj = 1,BLOCKING
c
              a(j_jump+jj,i_jump+ii) = temp1(ii,jj)
c
            enddo ! jj
          enddo ! ii
c
c       end blocked main 
c
        enddo ! j
      enddo ! i
c
c     begin blocked edge
c
!DEC$ prefetch a
!DEC$ loop count min(0), max(200), avg(32)
      do i = 1,chunks
        i_jump = (i-1)*BLOCKING
!DEC$ prefetch a
!DEC$ loop count min(BLOCKING), max(BLOCKING), avg(BLOCKING)
!DEC$ unroll(UNROLLING)
!DEC$ vector always
        do ii = 1,BLOCKING
!DEC$ prefetch a
!DEC$ loop count min(0), max(BLOCKING), avg(10)
!DEC$ unroll(UNROLLING)
!DEC$ vector always
          do j = 1,r
c
            swap = a(i_jump+ii,offset+j)
            a(i_jump+ii,offset+j) = a(offset+j,i_jump+ii)
            a(offset+j,i_jump+ii) = swap
c
          enddo ! j
        enddo ! ii
      enddo ! i
c
c     end non-blocked edge
c
c     begin non-blocked lower corner
c
!DEC$ prefetch a
!DEC$ loop count min(0), max(BLOCKING), avg(10)
!DEC$ unroll(UNROLLING)
!DEC$ vector always
      do i = 1,r
        i_ind = offset+i
!DEC$ prefetch a
!DEC$ loop count min(0), max(BLOCKING), avg(10)
!DEC$ unroll(UNROLLING)
!DEC$ vector always
        do j = i,r
          j_ind = offset+j
c
          swap = a(i_ind,j_ind)
          a(i_ind,j_ind) = a(j_ind,i_ind)
          a(j_ind,i_ind) = swap
c
        enddo
      enddo
c
c     end non-blocked lower corner
c
      return
      end
