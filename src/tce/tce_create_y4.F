      subroutine tce_create_y4(d_y4,l_y4_offset,k_y4_offset,size_y4,
     &                       irrep,title)
c
c $Id$
c
c Written by Jeff Hammond, January 2008.
c
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "errquit.fh"
c#include "tce.fh"
c#include "tce_main.fh"
c#include "tce_prop.fh"
c#include "tce_ints.fh"
c#include "tce_amps.fh"
c#include "tce_diis.fh"
      integer d_y1             ! Left 1 file
      integer d_y2             ! Left 2 file
      integer d_y3             ! Left 3 file
      integer d_y4             ! Left 4 file
      integer l_y1_offset      ! Offset for  1 file
      integer k_y1_offset      ! Offset for  1 file
      integer size_y1          ! File size in doubles
      integer l_y2_offset      ! Offset for  2 file
      integer k_y2_offset      ! Offset for  2 file
      integer size_y2          ! File size in doubles
      integer l_y3_offset      ! Offset for  3 file
      integer k_y3_offset      ! Offset for  3 file
      integer size_y3          ! File size in doubles
      integer l_y4_offset      ! Offset for  4 file
      integer k_y4_offset      ! Offset for  4 file
      integer size_y4          ! File size in doubles
      integer rtdb
      logical nodezero         ! True if node 0
      double precision cpu     ! CPU sec counter
      double precision wall    ! WALL sec counter
      integer irrep
      character*8 title
      character*255 filename
c
      nodezero=(ga_nodeid().eq.0)
      if (nodezero) write(LuOut,*) '============================='
c
      call tce_y4_offset_new(l_y4_offset,k_y4_offset,size_y4,irrep)
#if defined(DEBUG_PRINT)
        if(nodezero) write(LuOut,*) 'l_y4_offset = ',l_y4_offset
#endif
      call tce_filename(title,filename)
      call createfile(filename,d_y4,size_y4)
      call gatoeaf(d_y4)
      if (nodezero.and.util_print(title,print_default)) then
c        write(LuOut,*)
        write(LuOut,9000) title,size_y4
        write(LuOut,9010) title,filename(1:120)
#if defined(DEBUG_PRINT)
        write(LuOut,9090) title,d_y4
#endif
c       write(LuOut,9020) cpu
        call util_flush(LuOut)
      endif
c
      return
 9000 format(1x,A,' file size   = ',i16)
 9010 format(1x,A,' file name   = ',A)
 9020 format(1x,'Cpu & wall time / sec',2f15.1)
 9090 format(1x,A,' file handle = ',i10)

      end
