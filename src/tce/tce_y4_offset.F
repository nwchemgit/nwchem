      subroutine tce_y4_offset(l_y4_offset,k_y4_offset,size)
c
c $Id: tce_y4_offset.F,v 1.2 2003-10-19 00:52:11 sohirata Exp $
c 
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
      integer l_y4_offset
      integer k_y4_offset
      integer size
      integer p1b,p2b,p3b,p4b,h5b,h6b,h7b,h8b
      integer i
      logical nodezero
      double precision dbl_size
c
      nodezero = (ga_nodeid().eq.0)
      if (nodezero.and.util_print('offset',print_debug))
     1   write(LuOut,9000)
      if (.not.ma_push_get(mt_int,(nvab**4)*(noab**4),'y4 offset',
     1  l_y4_offset,k_y4_offset))
     2  call errquit('tce_y4_offset: MA problem',0,MA_ERR)
      do i = 1,(nvab**4)*(noab**4)
        int_mb(k_y4_offset+i-1) = -1
      enddo
      size = 0
      dbl_size = 0
      do h5b = 1,noab
        do h6b = h5b,noab
          do h7b = h6b,noab
            do h8b = h7b,noab
              do p1b = noab+1,noab+nvab
                do p2b = p1b,noab+nvab
                  do p3b = p2b,noab+nvab
                    do p4b = p3b,noab+nvab
                      if (int_mb(k_spin+p1b-1)
     1                   +int_mb(k_spin+p2b-1)
     2                   +int_mb(k_spin+p3b-1)
     3                   +int_mb(k_spin+p4b-1)
     4                .eq.int_mb(k_spin+h5b-1)
     5                   +int_mb(k_spin+h6b-1)
     6                   +int_mb(k_spin+h7b-1)
     7                   +int_mb(k_spin+h8b-1)) then
                      if ((.not.restricted).or.
     1                   (int_mb(k_spin+p1b-1)
     2                   +int_mb(k_spin+p2b-1)
     3                   +int_mb(k_spin+p3b-1)
     4                   +int_mb(k_spin+p4b-1)
     5                   +int_mb(k_spin+h5b-1)
     6                   +int_mb(k_spin+h6b-1)
     7                   +int_mb(k_spin+h7b-1)
     8                   +int_mb(k_spin+h8b-1).ne.16)) then
                      if (ieor(int_mb(k_sym+p1b-1),
     1                    ieor(int_mb(k_sym+p2b-1),
     2                    ieor(int_mb(k_sym+p3b-1),
     3                    ieor(int_mb(k_sym+p4b-1),
     4                    ieor(int_mb(k_sym+h5b-1),
     5                    ieor(int_mb(k_sym+h6b-1),
     6                    ieor(int_mb(k_sym+h7b-1),
     7                    int_mb(k_sym+h8b-1)))))))) .eq. irrep_y) then
                   int_mb(k_y4_offset+(((((((h5b-1)*noab
     1             +h6b-1)*noab+h7b-1)*noab+h8b-1)*nvab+p1b-noab-1)*nvab
     3             +p2b-noab-1)*nvab+p3b-noab-1)*nvab+p4b-noab-1) = size
                        if (nodezero.and.util_print('offset',
     1                    print_debug))
     2             write(LuOut,9010) p1b,p2b,p3b,p4b,
     3             h5b,h6b,h7b,h8b,(((((((h5b-1)*noab
     4             +h6b-1)*noab+h7b-1)*noab+h8b-1)*nvab+p1b-noab-1)*nvab
     5             +p2b-noab-1)*nvab+p3b-noab-1)*nvab+p4b-noab-1,size
                        size = size + int_mb(k_range+p1b-1)
     1                              * int_mb(k_range+p2b-1)
     2                              * int_mb(k_range+p3b-1)
     3                              * int_mb(k_range+p4b-1)
     4                              * int_mb(k_range+h5b-1)
     5                              * int_mb(k_range+h6b-1)
     6                              * int_mb(k_range+h7b-1)
     7                              * int_mb(k_range+h8b-1)
         dbl_size = dbl_size + dfloat(int_mb(k_range+p1b-1)
     1                              * int_mb(k_range+p2b-1)
     2                              * int_mb(k_range+p3b-1)
     3                              * int_mb(k_range+p4b-1)
     4                              * int_mb(k_range+h5b-1)
     5                              * int_mb(k_range+h6b-1)
     6                              * int_mb(k_range+h7b-1)
     7                              * int_mb(k_range+h8b-1))
                      endif
                      endif
                      endif
                    enddo
                  enddo
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
        call errquit('tce_y4_offset: abort',0,UNKNOWN_ERR)
      endif
      if (dbl_size .gt. 2.0d0*1024.0d0**3) then
        write(LuOut,*) '---------------------------------------------'
        write(LuOut,*) 'Abort: 32 bit integer overflow is anticipated'
        write(LuOut,*) 'Use a machine with a 64 bit operating system'
        write(LuOut,*) '---------------------------------------------'
        call errquit('tce_y4_offset: abort',0,UNKNOWN_ERR)
      endif
#endif
 9000 format(/,1x,'y4 offsets')
 9010 format(1x,8i5,' [',i10,'] offset:',i10)
      return
      end
