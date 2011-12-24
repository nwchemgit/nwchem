      subroutine get_hash_block_ma(big_array,array,size,hash,key)
c
c $Id$
c
      implicit none
c
      integer size
      integer hash(*)
      integer key
      integer offset
      integer i
      double precision array(size)
      double precision big_array(*)
c
      call tce_hash(hash,key,offset)
      do i=1,size
       array(i)=big_array(offset+i)
      enddo
c
      return
      end
c
c
c
      subroutine ma_zero(array,size_ma)
      implicit none
c
      integer k_ma,size_ma,i
      double precision array(size_ma)
      do i=1,size_ma
       array(i)=0.0d0 
      enddo 
      return
      end
c
c
c
