      subroutine tce_x3_offset(l_x3_offset,k_x3_offset,size)
c
c $Id: tce_x3_offset.F,v 1.2 2003-10-19 00:52:10 sohirata Exp $
c 
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
      integer l_x3_offset
      integer k_x3_offset
      integer size
      integer p1b,p2b,p3b,h4b,h5b,h6b
      integer i
      logical nodezero
      double precision dbl_size
c
      nodezero = (ga_nodeid().eq.0)
      if (nodezero.and.util_print('offset',print_debug))
     1   write(LuOut,9000)
      if (.not.ma_push_get(mt_int,(nvab**3)*(noab**3),'x3 offset',
     1  l_x3_offset,k_x3_offset))
     2  call errquit('tce_x3_offset: MA problem',0,MA_ERR)
      do i = 1,(nvab**3)*(noab**3)
        int_mb(k_x3_offset+i-1) = -1
      enddo
      size = 0
      dbl_size = 0
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
     5                     int_mb(k_sym+h6b-1)))))) .eq. irrep_x) then
                    int_mb(k_x3_offset+(((((p1b-noab-1)*nvab
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
                    dbl_size = dbl_size + dfloat(int_mb(k_range+p1b-1)
     1                          * int_mb(k_range+p2b-1)
     2                          * int_mb(k_range+p3b-1)
     3                          * int_mb(k_range+h4b-1)
     4                          * int_mb(k_range+h5b-1)
     5                          * int_mb(k_range+h6b-1))
                  endif
                  endif
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
#if !defined(EXT_INT)
      if ((dbl_size .gt. 16.0d0*1024.0d0**3).and.(ioalg.eq.0)) then
        write(LuOut,*) '--------------------------------------------'
        write(LuOut,*) 'Abort: 2 GB file size limit will be exceeded'
        write(LuOut,*) 'Use GA/SF scheme or OS large file support'
        write(LuOut,*) '--------------------------------------------'
        call errquit('tce_x3_offset: abort',0,UNKNOWN_ERR)
      endif
      if (dbl_size .gt. 2.0d0*1024.0d0**3) then
        write(LuOut,*) '---------------------------------------------'
        write(LuOut,*) 'Abort: 32 bit integer overflow is anticipated'
        write(LuOut,*) 'Use a machine with a 64 bit operating system'
        write(LuOut,*) '---------------------------------------------'
        call errquit('tce_x3_offset: abort',0,UNKNOWN_ERR)
      endif
#endif
 9000 format(/,1x,'x3 offsets')
 9010 format(1x,6i5,' [',i10,'] offset:',i10)
      return
      end
