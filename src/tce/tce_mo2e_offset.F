      subroutine tce_mo2e_offset(l_2e_offset,k_2e_offset,size)
c
c $Id: tce_mo2e_offset.F,v 1.2 2002-11-20 06:01:33 sohirata Exp $
c 
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "tce.fh"
      integer l_2e_offset
      integer k_2e_offset
      integer Nb
      integer size
      integer g1b,g2b,g3b,g4b
      logical nodezero
c
      nodezero = (ga_nodeid().eq.0)
      if (nodezero.and.util_print('offset',print_debug))
     1   write(LuOut,9000)
      Nb = noab + nvab
      if (.not.ma_push_get(mt_int,Nb**4,'2e offset',
     1  l_2e_offset,k_2e_offset))
     2  call errquit('tce_mo2e_offset: MA problem',0)
      size = 0
      do g1b = 1,Nb
        do g2b = g1b,Nb
          do g3b = 1,Nb
            do g4b = g3b,Nb
              if (int_mb(k_spin+g1b-1)+int_mb(k_spin+g2b-1)
     1          .eq. int_mb(k_spin+g3b-1)+int_mb(k_spin+g4b-1)) then
              if ((.not.restricted).or.
     1          (int_mb(k_spin+g1b-1)+int_mb(k_spin+g2b-1)+
     2          int_mb(k_spin+g3b-1)+int_mb(k_spin+g4b-1).ne.8)) then
              if (ieor(int_mb(k_sym+g1b-1),ieor(int_mb(k_sym+g2b-1),
     1          ieor(int_mb(k_sym+g3b-1),int_mb(k_sym+g4b-1))))
     2          .eq. 0) then
                int_mb(k_2e_offset+(((g1b-1)*Nb+g2b-1)*Nb
     1            +g3b-1)*Nb+g4b-1) = size
                if (nodezero.and.util_print('offset',print_debug))
     1            write(LuOut,9010) g1b,g2b,g3b,g4b,
     2            (((g1b-1)*Nb+g2b-1)*Nb+g3b-1)*Nb+g4b-1,size
                size = size + int_mb(k_range+g1b-1)
     1                      * int_mb(k_range+g2b-1)
     2                      * int_mb(k_range+g3b-1)
     3                      * int_mb(k_range+g4b-1)
              endif
              endif
              endif
            enddo
          enddo
        enddo
      enddo
 9000 format(/,1x,'2e offsets')
 9010 format(1x,4i5,' [',i10,'] offset:',i10)
      return
      end
