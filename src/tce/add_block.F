      subroutine add_block(d_file,array,size,offset)
c
c $Id: add_block.F,v 1.2 2002-10-23 01:38:51 sohirata Exp $
c
c     This subroutine reads a section of a file and adds
c     the content of an array to it and write it back to
c     the original place in the file.  Note that two or more
c     different processes can read/add/write to the same 
c     section, the entire read/add/write must be protected
c     from interruption.
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sf.fh"
#include "util.fh"
#include "stdio.fh"
#include "tce.fh"
      integer d_file
      integer size
      integer offset
      double precision array(size)
      integer l_temp,k_temp
      integer request
      integer i
      if (util_print('add_block',print_debug)) then
         write(LuOut,9000) ga_nodeid(),d_file,size,offset
      endif
      if (offset .eq. -1)
     1  call errquit('add_block: illegal offset',0)
      if (.not.ma_push_get(mt_dbl,size,'temporary',
     1  l_temp,k_temp))
     2  call errquit('add_block: MA problem',0)
c
c     Critical section
c
      call ga_lock(0)
      if (sf_read(d_file,dfloat(bytes*offset),
     1  dfloat(bytes*size),dbl_mb(k_temp),request).ne.0)
     2  call errquit('add_block: sf problem',0)
      if (sf_wait(request).ne.0)
     1  call errquit('add_block: sf problem',1)
      do i = 1,size
        dbl_mb(k_temp+i-1) = dbl_mb(k_temp+i-1) + array(i)
      enddo
      if (sf_write(d_file,dfloat(bytes*offset),
     1  dfloat(bytes*size),dbl_mb(k_temp),request).ne.0)
     2  call errquit('add_block: sf problem',2)
      if (sf_wait(request).ne.0)
     1  call errquit('add_block: sf problem',3)
      call ga_unlock(0)
c
c     End of critical section
c
      if (.not.ma_pop_stack(l_temp))
     1  call errquit('add_block: MA problem',0)
 9000 format(1x,'node',i3,' add_block request to file:',i10,
     1  ' size:',i10,' offset:',i10)
      return
      end
