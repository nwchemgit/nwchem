      subroutine tce_mo1e_offset(l_1e_offset,k_1e_offset,size)
c
c $Id: tce_mo1e_offset.F,v 1.1 2002-10-18 23:42:08 sohirata Exp $
c 
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "tce.fh"
      integer l_1e_offset
      integer k_1e_offset
      integer Nb
      integer size
      integer g1b,g2b
      logical nodezero
c
      nodezero = (ga_nodeid().eq.0)
      if (nodezero.and.util_print('offset',print_debug))
     1   write(LuOut,9000)
      Nb = noab + nvab
      if (.not.ma_push_get(mt_int,Nb*Nb,'1e offset',
     1  l_1e_offset,k_1e_offset))
     2  call errquit('tce_mo1e_offset: MA problem',0)
      size = 0
      do g1b = 1,Nb
        do g2b = 1,Nb
          if (int_mb(k_spin+g1b-1) .eq. int_mb(k_spin+g2b-1)) then
            if (ieor(int_mb(k_sym+g1b-1),int_mb(k_sym+g2b-1))
     1        .eq. 0) then
              int_mb(k_1e_offset+(g1b-1)*Nb+g2b-1) = size
              if (nodezero.and.util_print('offset',print_debug))
     1           write(LuOut,9010) g1b,g2b,(g1b-1)*Nb+g2b-1,size
              size = size + int_mb(k_range+g1b-1)
     1                    * int_mb(k_range+g2b-1)
            endif
          endif
        enddo
      enddo
 9000 format(/,1x,'1e offsets')
 9010 format(1x,2i5,' [',i10,'] offset:',i10)
      return
      end
