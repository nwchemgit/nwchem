      subroutine tce_clone_y3(d_y3,size_y3,irrep,title)
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
#include "tce_main.fh"
c#include "tce_prop.fh"
c#include "tce_ints.fh"
#include "tce_amps.fh"
c#include "tce_diis.fh"
      integer rtdb
      logical nodezero         ! True if node 0
      double precision cpu     ! CPU sec counter
      double precision wall    ! WALL sec counter
      integer irrep
      character*8 title
      character*255 filename
c
      nodezero=(ga_nodeid().eq.0)
c
c      call tce_y3_offset(l_y3_offset,k_y3_offset,size_y3,irrep)
      call tce_filename(title,filename)
      call createfile(filename,d_y3,size_y3)
      call gatoeaf(d_y3)
      if (nodezero.and.util_print(title,print_default)) then
c        write(LuOut,*)
        write(LuOut,9000) title,size_y3
        write(LuOut,9010) title,filename(1:120)
#if defined(DEBUG_PRINT)
        write(LuOut,9090) title,d_y3
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
