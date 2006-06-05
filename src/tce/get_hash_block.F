      subroutine get_hash_block(d_file,array,size,hash,key)
c
c $Id: get_hash_block.F,v 1.2 2006-06-05 18:22:57 kowalski Exp $
c
      implicit none
#include "tce.fh"
c
      integer d_file
      integer size
      integer hash(*)
      integer key
      integer offset
      integer indexc
      integer iauxo
      double precision array(size)
c
      if (.not.intorb) then
       call tce_hash(hash,key,offset)
       call get_block(d_file,array,size,offset)
      else if (intorb.and.(d_file.ne.d_v2orb)) then
       call tce_hash(hash,key,offset)
       call get_block(d_file,array,size,offset)
      else if (intorb.and.(d_file.eq.d_v2orb)) then
       iauxo=k_v2spin_offset
       call tce_hash_ind(iauxo,key,indexc)
       call get_block_ind(d_file,array,size,key,indexc)
      end if
c
c
      return
      end
