      subroutine tce_e_offset(l_e_offset,k_e_offset,size)
c
c $Id: tce_e_offset.F,v 1.1 2002-10-18 23:42:07 sohirata Exp $
c 
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "tce.fh"
      integer l_e_offset
      integer k_e_offset
      integer size
      logical nodezero
c
      nodezero = (ga_nodeid().eq.0)
      if (.not.ma_push_get(mt_int,1,'e offset',
     1  l_e_offset,k_e_offset))
     2  call errquit('tce_e_offset: MA problem',0)
      int_mb(k_e_offset) = 0
      size = 1
      if (nodezero.and.util_print('offset',print_debug)) then
        write(LuOut,9000)
        write(LuOut,9010) 1,1,0
      endif
 9000 format(/,1x,'e offsets')
 9010 format(1x,i5,' [',i10,'] offset:',i10)
      return
      end
