      subroutine tce_t2_offset(l_t2_offset,k_t2_offset,size)
c
c $Id: tce_t2_offset.F,v 1.1 2002-10-18 23:42:09 sohirata Exp $
c 
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "tce.fh"
      integer l_t2_offset
      integer k_t2_offset
      integer size
      integer p1b,p2b,h3b,h4b
      integer i
      logical nodezero
c
      nodezero = (ga_nodeid().eq.0)
      if (nodezero.and.util_print('offset',print_debug))
     1   write(LuOut,9000)
      if (.not.ma_push_get(mt_int,nvab*nvab*noab*noab,'t2 offset',
     1  l_t2_offset,k_t2_offset))
     2  call errquit('tce_t2_offset: MA problem',0)
      do i = 1,nvab*nvab*noab*noab
        int_mb(k_t2_offset+i-1) = -1
      enddo
      size = 0
      do p1b = noab+1,noab+nvab
        do p2b = p1b,noab+nvab
          do h3b = 1,noab
            do h4b = h3b,noab
              if (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)
     1          .eq. int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) then
                if (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),
     1            ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1))))
     2            .eq. 0) then
                  int_mb(k_t2_offset+(((p1b-noab-1)*nvab
     1              +p2b-noab-1)*noab+h3b-1)*noab+h4b-1) = size
                  if (nodezero.and.util_print('offset',print_debug))
     1              write(LuOut,9010)
     2              p1b,p2b,h3b,h4b,(((p1b-noab-1)*nvab
     3              +p2b-noab-1)*noab+h3b-1)*noab+h4b-1,size
                  size = size + int_mb(k_range+p1b-1)
     1                        * int_mb(k_range+p2b-1)
     2                        * int_mb(k_range+h3b-1)
     3                        * int_mb(k_range+h4b-1)
                endif
              endif
            enddo
          enddo
        enddo
      enddo
 9000 format(/,1x,'t2 offsets')
 9010 format(1x,4i5,' [',i10,'] offset:',i10)
      return
      end
