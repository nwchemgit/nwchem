      subroutine tce_t1_offset(l_t1_offset,k_t1_offset,size)
c
c $Id: tce_t1_offset.F,v 1.1 2002-10-18 23:42:09 sohirata Exp $
c 
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "tce.fh"
      integer l_t1_offset
      integer k_t1_offset
      integer size
      integer p1b,h2b
      integer i
      logical nodezero
c
      nodezero = (ga_nodeid().eq.0)
      if (nodezero.and.util_print('offset',print_debug))
     1   write(LuOut,9000)
      if (.not.ma_push_get(mt_int,nvab*noab,'t1 offset',
     1  l_t1_offset,k_t1_offset))
     2  call errquit('tce_t1_offset: MA problem',0)
      do i = 1,nvab*noab
        int_mb(k_t1_offset+i-1) = -1
      enddo
      size = 0
      do p1b = noab+1,noab+nvab
        do h2b = 1,noab
          if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h2b-1)) then
            if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h2b-1))
     1        .eq. 0) then
              int_mb(k_t1_offset+(p1b-noab-1)*noab+h2b-1) = size
              if (nodezero.and.util_print('offset',print_debug))
     1           write(LuOut,9010) p1b,h2b,(p1b-noab-1)*noab+h2b-1,size
              size = size + int_mb(k_range+p1b-1)
     1                    * int_mb(k_range+h2b-1)
            endif
          endif
        enddo
      enddo
 9000 format(/,1x,'t1 offsets')
 9010 format(1x,2i5,' [',i10,'] offset:',i10)
      return
      end
