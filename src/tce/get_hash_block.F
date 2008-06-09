      subroutine get_hash_block(d_file,array,size,hash,key)
c
c $Id: get_hash_block.F,v 1.5 2008-06-09 21:01:24 jhammond Exp $
c
      implicit none
#include "tce.fh"
#include "tce_main.fh"
#include "tce_hetio.fh"
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
        if ((d_file.eq.d_v2orb).and.(io_v2.eq.IOSF)) then
          call tce_hash(hash,key,offset)
          call get_block_sf(d_file,array,size,offset)
        else
          call tce_hash(hash,key,offset)
          call get_block(d_file,array,size,offset)
        endif
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
c
c
      subroutine get_hash_block_i(d_file,array,size,hash,key,
     &                            g2b,g1b,g4b,g3b)
c
c $Id: get_hash_block.F,v 1.5 2008-06-09 21:01:24 jhammond Exp $
c
      implicit none
#include "tce.fh"
#include "tce_main.fh"
c
      integer d_file
      integer size
      integer hash(*)
      integer key
      integer offset
      integer indexc
      integer iauxo
      integer g2b,g1b,g4b,g3b
      integer qy2,qy1,qy4,qy3
      double precision array(size)
c
      if (.not.intorb) then
       call tce_hash(hash,key,offset)
       call get_block(d_file,array,size,offset)
      else if (intorb.and.(d_file.ne.d_v2orb)) then
       call tce_hash(hash,key,offset)
       call get_block(d_file,array,size,offset)
      else if (intorb.and.(d_file.eq.d_v2orb)) then
c       iauxo=k_v2spin_offset
c       call tce_hash_ind(iauxo,key,indexc)
c *** debug ***
c       write(6,*)'from ghbi g-s',g2b,g1b,g4b,g3b
c       call util_flush(6)
c ***************
       qy2=g2b
       qy1=g1b
       qy4=g4b
       qy3=g3b
       call get_block_ind_i(d_file,array,size,key,indexc,
     &                      qy2,qy1,qy4,qy3)
c *** debug ***
c       write(6,*)'after ghbi g-s',g2b,g1b,g4b,g3b
c       call util_flush(6)
c ***************
      end if
c
c
      return
      end
