      subroutine tce_t2_offset(l_t2_offset,k_t2_offset,size)
c
c $Id: tce_t2_offset.F,v 1.3 2003-01-08 22:59:14 sohirata Exp $
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
      double precision dbl_size
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
      dbl_size = 0
      do p1b = noab+1,noab+nvab
        do p2b = p1b,noab+nvab
          do h3b = 1,noab
            do h4b = h3b,noab
              if (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)
     1          .eq. int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) then
              if ((.not.restricted).or.
     1          (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+
     2          int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) then
              if (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),
     1          ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1))))
     2          .eq. 0) then
                int_mb(k_t2_offset+(((p1b-noab-1)*nvab
     1            +p2b-noab-1)*noab+h3b-1)*noab+h4b-1) = size
                if (nodezero.and.util_print('offset',print_debug))
     1            write(LuOut,9010)
     2            p1b,p2b,h3b,h4b,(((p1b-noab-1)*nvab
     3            +p2b-noab-1)*noab+h3b-1)*noab+h4b-1,size
                size = size + int_mb(k_range+p1b-1)
     1                      * int_mb(k_range+p2b-1)
     2                      * int_mb(k_range+h3b-1)
     3                      * int_mb(k_range+h4b-1)
                dbl_size = dbl_size + dfloat(int_mb(k_range+p1b-1)
     1                      * int_mb(k_range+p2b-1)
     2                      * int_mb(k_range+h3b-1)
     3                      * int_mb(k_range+h4b-1))
              endif
              endif
              endif
            enddo
          enddo
        enddo
      enddo
#if !defined(EXT_INT)
      if ((dbl_size .gt. 16.0d0*1024.0d0**3).and.
     1   ((ioalg.eq.0).or.(ioalg.eq.1).or.(ioalg.eq.4))) then
        write(LuOut,*) '--------------------------------------------'
        write(LuOut,*) 'Abort: 2 GB file size limit will be exceeded'
        write(LuOut,*) 'Use GA/SF scheme or OS large file support'
        write(LuOut,*) '--------------------------------------------'
        call errquit('tce_t2_offset: abort',0)
      endif
      if (dbl_size .gt. 2.0d0*1024.0d0**3) then
        write(LuOut,*) '---------------------------------------------'
        write(LuOut,*) 'Abort: 32 bit integer overflow is anticipated'
        write(LuOut,*) 'Use a machine with a 64 bit operating system'
        write(LuOut,*) '---------------------------------------------'
        call errquit('tce_t2_offset: abort',0)
      endif
#endif
 9000 format(/,1x,'t2 offsets')
 9010 format(1x,4i5,' [',i10,'] offset:',i10)
      return
      end
