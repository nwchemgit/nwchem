      subroutine tce_x1_offset(l_x1_offset,k_x1_offset,size)
c
c $Id: tce_x1_offset.F,v 1.1 2003-07-25 17:56:58 sohirata Exp $
c 
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "tce.fh"
      integer l_x1_offset
      integer k_x1_offset
      integer size
      integer p1b,h2b
      integer i
      logical nodezero
      double precision dbl_size
c
      nodezero = (ga_nodeid().eq.0)
      if (nodezero.and.util_print('offset',print_debug))
     1   write(LuOut,9000)
      if (.not.ma_push_get(mt_int,nvab*noab,'x1 offset',
     1  l_x1_offset,k_x1_offset))
     2  call errquit('tce_x1_offset: MA problem',0)
      do i = 1,nvab*noab
        int_mb(k_x1_offset+i-1) = -1
      enddo
      size = 0
      dbl_size = 0
      do p1b = noab+1,noab+nvab
        do h2b = 1,noab
          if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h2b-1)) then
          if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
     1      +int_mb(k_spin+h2b-1).ne.4)) then
          if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h2b-1))
     1      .eq.irrep_x) then
            int_mb(k_x1_offset+(p1b-noab-1)*noab+h2b-1) = size
            if (nodezero.and.util_print('offset',print_debug))
     1         write(LuOut,9010) p1b,h2b,(p1b-noab-1)*noab+h2b-1,size
            size = size + int_mb(k_range+p1b-1)
     1                  * int_mb(k_range+h2b-1)
            dbl_size = dbl_size + dfloat(int_mb(k_range+p1b-1))
     1                          * dfloat(int_mb(k_range+h2b-1))
          endif
          endif
          endif
        enddo
      enddo
#if !defined(EXT_INT)
      if ((dbl_size .gt. 16.0d0*1024.0d0**3).and.(ioalg.eq.0)) then
        write(LuOut,*) '--------------------------------------------'
        write(LuOut,*) 'Abort: 2 GB file size limit will be exceeded'
        write(LuOut,*) 'Use GA/SF scheme or OS large file support'
        write(LuOut,*) '--------------------------------------------'
        call errquit('tce_x1_offset: abort',0)
      endif
      if (dbl_size .gt. 2.0d0*1024.0d0**3) then
        write(LuOut,*) '---------------------------------------------'
        write(LuOut,*) 'Abort: 32 bit integer overflow is anticipated'
        write(LuOut,*) 'Use a machine with a 64 bit operating system'
        write(LuOut,*) '---------------------------------------------'
        call errquit('tce_x1_offset: abort',0)
      endif
#endif
 9000 format(/,1x,'x1 offsets')
 9010 format(1x,2i5,' [',i10,'] offset:',i10)
      return
      end
