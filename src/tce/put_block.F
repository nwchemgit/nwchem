      subroutine put_block(d_file,array,size,offset)
c
c $Id: put_block.F,v 1.4 2002-10-30 21:49:14 sohirata Exp $
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
      integer request
      integer i,j
      logical used
c
      if (util_print('put_block',print_debug)) then
         write(LuOut,9000) ga_nodeid(),d_file,size,offset
      endif
      if (offset .eq. -1)
     1  call errquit('put_block: illegal offset',0)
      if (ioalg.eq.0) then
        call ga_lock(0)
        used = .false.
        do i = 1,nfiles
          if (filehandles(i) .eq. d_file) then
            j = i
            used = .true.
          endif
        enddo
        if (.not.used) call errquit('get_block: invalid handle',d_file)
        open(d_file,file=filenames(j),access='direct',recl=bytes)
        do i = 1,size
          write(d_file,rec=offset+i) array(i)
        enddo
        call flush(d_file)
        close(d_file)
        call ga_unlock(0)
      else if (ioalg.eq.1) then
        call ga_lock(0)
        if (sf_write(d_file,dfloat(bytes*offset),
     1    dfloat(bytes*size),array,request).ne.0)
     2    call errquit('put_block: sf problem',0)
        if (sf_wait(request).ne.0)
     1    call errquit('put_block: sf problem',1)
        call ga_unlock(0)
      endif
 9000 format(1x,'node',i3,' put_block request to file:',i10,
     1  ' size:',i10,' offset:',i10)
      return
      end
