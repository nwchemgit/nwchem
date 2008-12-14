      subroutine add_hash_block(d_file,array,size,hash,key)
c
c $Id: add_hash_block.F,v 1.2 2008-12-14 00:10:48 jhammond Exp $
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









      subroutine add_hash_block_nb(d_file,array,size,hash,key,nbh)
c
      implicit none
c
      double precision array(size)
c
      integer d_file
      integer size
      integer hash(*)
      integer key
      integer offset
      integer nbh
c
      call tce_hash(hash,key,offset)
      call add_block_nb(d_file,array,size,offset,nbh)
c
      return
      end
