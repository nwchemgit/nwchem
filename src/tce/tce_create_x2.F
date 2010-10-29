      subroutine tce_create_x2(d_x2,l_x2_offset,k_x2_offset,size_x2,
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
      integer d_x1             ! Right 1 file
      integer d_x2             ! Right 2 file
      integer d_x3             ! Right 3 file
      integer d_x4             ! Right 4 file
      integer l_x1_offset      ! Offset for  1 file
      integer k_x1_offset      ! Offset for  1 file
      integer size_x1          ! File size in doubles
      integer l_x2_offset      ! Offset for  2 file
      integer k_x2_offset      ! Offset for  2 file
      integer size_x2          ! File size in doubles
      integer l_x3_offset      ! Offset for  3 file
      integer k_x3_offset      ! Offset for  3 file
      integer size_x3          ! File size in doubles
      integer l_x4_offset      ! Offset for  4 file
      integer k_x4_offset      ! Offset for  4 file
      integer size_x4          ! File size in doubles
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
      call tce_x2_offset_new(l_x2_offset,k_x2_offset,size_x2,irrep)
#if defined(DEBUG_PRINT)
        if(nodezero) write(LuOut,*) 'l_x2_offset = ',l_x2_offset
#endif
      call tce_filename(title,filename)
      call createfile(filename,d_x2,size_x2)
      call gatoeaf(d_x2)
      if (nodezero.and.util_print(title,print_default)) then
c        write(LuOut,*)
        write(LuOut,9000) title,size_x2
        write(LuOut,9010) title,filename(1:120)
#if defined(DEBUG_PRINT)
        write(LuOut,9090) title,d_x2
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
