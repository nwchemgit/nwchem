      subroutine tce_mod1_offset(l_d1_offset,k_d1_offset,size)
c
c $Id: tce_mod1_offset.F,v 1.1 2003-07-25 17:56:57 sohirata Exp $
c 
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "tce.fh"
      integer l_d1_offset
      integer k_d1_offset
      integer Nb
      integer size
      integer g1b,g2b
      logical nodezero
c
      nodezero = (ga_nodeid().eq.0)
      if (nodezero.and.util_print('offset',print_debug))
     1   write(LuOut,9000)
      Nb = noab + nvab
      if (.not.ma_push_get(mt_int,Nb*Nb,'d1 offset',
     1  l_d1_offset,k_d1_offset))
     2  call errquit('tce_mod1_offset: MA problem',0)
      size = 0
      do g1b = 1,Nb
        do g2b = 1,Nb
          if (int_mb(k_spin+g1b-1) .eq. int_mb(k_spin+g2b-1)) then
          if ((.not.restricted).or.(int_mb(k_spin+g1b-1)
     1      +int_mb(k_spin+g2b-1).ne.4)) then
          if (ieor(int_mb(k_sym+g1b-1),int_mb(k_sym+g2b-1))
     1      .eq. irrep_d) then
            int_mb(k_d1_offset+(g1b-1)*Nb+g2b-1) = size
            if (nodezero.and.util_print('offset',print_debug))
     1         write(LuOut,9010) g1b,g2b,(g1b-1)*Nb+g2b-1,size
            size = size + int_mb(k_range+g1b-1)
     1                  * int_mb(k_range+g2b-1)
          endif
          endif
          endif
        enddo
      enddo
 9000 format(/,1x,'d1 offsets')
 9010 format(1x,2i5,' [',i10,'] offset:',i10)
      return
      end
