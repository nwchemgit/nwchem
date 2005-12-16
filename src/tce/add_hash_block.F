      subroutine add_hash_block(d_file,array,size,hash,key)
c
c $Id: add_hash_block.F,v 1.1 2005-12-16 22:20:36 kowalski Exp $
c
      implicit none
      integer d_file
      integer size
      integer hash(*)
      integer key
      integer offset
      double precision array(size)
c
      call tce_hash(hash,key,offset)
      call add_block(d_file,array,size,offset)
      return
      end
