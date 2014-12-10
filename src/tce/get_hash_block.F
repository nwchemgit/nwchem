      subroutine get_hash_block(d_file,array,size,hash,key)
c
c $Id$
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
      else if (intorb) then
        if ((d_file.eq.d_v2orb).and.(io_v2.eq.IOSF)) then
          iauxo=k_v2spin_offset
          call tce_hash_ind(iauxo,key,indexc)
          call get_block_ind_sf(d_file,array,size,key,indexc)
        else if ((d_file.eq.d_v2orb).and.(io_v2.ne.IOSF)) then
          iauxo=k_v2spin_offset
          call tce_hash_ind(iauxo,key,indexc)
          call get_block_ind(d_file,array,size,key,indexc)
        else
          call tce_hash(hash,key,offset)
          call get_block(d_file,array,size,offset)
        endif
      endif
c
c
      return
      end
c
c
      subroutine get_hash_block_i(d_file,array,size,hash,key,
     &                            g2b,g1b,g4b,g3b)
c
c $Id$
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
      integer g2b,g1b,g4b,g3b
      integer qy2,qy1,qy4,qy3
      double precision array(size)
c
      if (.not.intorb) then
       call tce_hash(hash,key,offset)
       call get_block(d_file,array,size,offset)
      else if (intorb) then
        if ((d_file.eq.d_v2orb).and.(io_v2.eq.IOGA)) then
c          iauxo=k_v2spin_offset
c          call tce_hash_ind(iauxo,key,indexc)
          qy2=g2b
          qy1=g1b
          qy4=g4b
          qy3=g3b
          call get_block_ind_i(d_file,array,size,key,indexc,
     &                         qy2,qy1,qy4,qy3)
        else if ((d_file.eq.d_v2orb).and.(io_v2.eq.IOHYBRID)) then
c          print*,'d_file  = ',d_file
c          print*,'d_v2orb = ',d_v2orb
          qy2=g2b
          qy1=g1b
          qy4=g4b
          qy3=g3b
          if (.not.((g3b.gt.noa).and.(g4b.gt.noa).and.
     1              (g1b.gt.noa).and.(g2b.gt.noa))) then
c           This is the GA in-core original
c            print*,'going to GA'
c            print*,'d_v2ga  = ',d_v2ga
            call get_block_ind_i(d_v2ga,array,size,key,indexc,
     &                           qy2,qy1,qy4,qy3)
          else
c           This is the SF out-of-core copy
c            print*,'going to SF'
c            print*,'d_v2sf  = ',d_v2sf
            call get_block_ind_sf_i(d_v2sf,array,size,key,indexc,
     &                              qy2,qy1,qy4,qy3)
          endif
        else if ((d_file.eq.d_v2orb).and.(io_v2.eq.IOSF)) then
          qy2=g2b
          qy1=g1b
          qy4=g4b
          qy3=g3b
          call get_block_ind_sf_i(d_file,array,size,key,indexc,
     &                            qy2,qy1,qy4,qy3)
        else
          call tce_hash(hash,key,offset)
          call get_block(d_file,array,size,offset)
        endif
      end if
c
c
      return
      end
