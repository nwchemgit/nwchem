      subroutine tce_t3_offset(l_t3_offset,k_t3_offset,size)
c
c $Id: tce_t3_offset.F,v 1.2 2002-11-20 05:03:06 sohirata Exp $
c 
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "tce.fh"
      integer l_t3_offset
      integer k_t3_offset
      integer size
      integer p1b,p2b,p3b,h4b,h5b,h6b
      integer i
      logical nodezero
c
      nodezero = (ga_nodeid().eq.0)
      if (nodezero.and.util_print('offset',print_debug))
     1   write(LuOut,9000)
      if (.not.ma_push_get(mt_int,(nvab**3)*(noab**3),'t3 offset',
     1  l_t3_offset,k_t3_offset))
     2  call errquit('tce_t3_offset: MA problem',0)
      do i = 1,(nvab**3)*(noab**3)
        int_mb(k_t3_offset+i-1) = -1
      enddo
      size = 0
      do p1b = noab+1,noab+nvab
        do p2b = p1b,noab+nvab
          do p3b = p2b,noab+nvab
            do h4b = 1,noab
              do h5b = h4b,noab
                do h6b = h5b,noab
                  if (int_mb(k_spin+p1b-1)
     1               +int_mb(k_spin+p2b-1)
     2               +int_mb(k_spin+p3b-1)
     3            .eq.int_mb(k_spin+h4b-1)
     4               +int_mb(k_spin+h5b-1)
     5               +int_mb(k_spin+h6b-1)) then
                  if ((.not.restricted).or.
     1               (int_mb(k_spin+p1b-1)
     2               +int_mb(k_spin+p2b-1)
     3               +int_mb(k_spin+p3b-1)
     4               +int_mb(k_spin+h4b-1)
     5               +int_mb(k_spin+h5b-1)
     6               +int_mb(k_spin+h6b-1).ne.12)) then
                  if (ieor(int_mb(k_sym+p1b-1),
     1                ieor(int_mb(k_sym+p2b-1),
     2                ieor(int_mb(k_sym+p3b-1),
     3                ieor(int_mb(k_sym+h4b-1),
     4                ieor(int_mb(k_sym+h5b-1),
     5                     int_mb(k_sym+h6b-1)))))) .eq. 0) then
                    int_mb(k_t3_offset+(((((p1b-noab-1)*nvab
     1                +p2b-noab-1)*nvab+p3b-noab-1)*noab
     2                +h4b-1)*noab+h5b-1)*noab+h6b-1) = size
                    if (nodezero.and.util_print('offset',print_debug))
     1                write(LuOut,9010)
     2                p1b,p2b,p3b,h4b,h5b,h6b,
     3                (((((p1b-noab-1)*nvab
     4                +p2b-noab-1)*nvab+p3b-noab-1)*noab
     5                +h4b-1)*noab+h5b-1)*noab+h6b-1,size
                    size = size + int_mb(k_range+p1b-1)
     1                          * int_mb(k_range+p2b-1)
     2                          * int_mb(k_range+p3b-1)
     3                          * int_mb(k_range+h4b-1)
     4                          * int_mb(k_range+h5b-1)
     5                          * int_mb(k_range+h6b-1)
                  endif
                  endif
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
 9000 format(/,1x,'t3 offsets')
 9010 format(1x,6i5,' [',i10,'] offset:',i10)
      return
      end
