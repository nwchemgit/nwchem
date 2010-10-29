      subroutine tce_clone_x2(d_x2,size_x2,irrep,title)
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
c      call tce_x2_offset(l_x2_offset,k_x2_offset,size_x2,irrep)
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
