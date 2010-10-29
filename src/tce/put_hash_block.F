      subroutine put_hash_block(d_file,array,size,hash,key)
c
c $Id$
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
      call put_block(d_file,array,size,offset)
      return
      end






      subroutine put_hash_block_sf(d_file,array,size,hash,key)
c
c $Id$
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
      call put_block_sf(d_file,array,size,offset)
      return
      end
