      subroutine get_block(d_file,array,size,offset)
c
c $Id: get_block.F,v 1.1 2002-10-18 23:42:06 sohirata Exp $
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
      integer dummy
      integer request
      logical nodezero
c
      nodezero = (ga_nodeid().eq.0)
      if (nodezero.and.util_print('get_block',print_debug)) then
         write(LuOut,9000) d_file,size,offset
      endif
      if (offset .eq. -1)
     1  call errquit('get_block: illegal offset',0)
      dummy=sf_read(d_file,dfloat(bytes*offset),
     1  dfloat(bytes*size),array,request)
      dummy=sf_wait(request)
 9000 format(1x,'get_block request to file:',i10,
     1  ' size:',i10,' offset:',i10)
      return
      end
