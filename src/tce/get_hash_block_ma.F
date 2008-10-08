      subroutine get_hash_block_ma(big_array,array,size,hash,key)
c
c $Id: get_hash_block_ma.F,v 1.1 2008-10-08 20:46:07 kowalski Exp $
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "tce.fh"
#include "tce_main.fh"
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
#include "global.fh"
#include "mafdecls.fh"
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
