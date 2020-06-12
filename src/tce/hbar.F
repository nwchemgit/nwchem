      logical function hbar_restore(ivec,jvec,hbar)
c
c $Id$
c
c     Read (if exist) hbar(ivec,jvec)
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sf.fh"
#include "msgids.fh"
#include "eaf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      character*255 filename
      integer ivec
      integer jvec
      double precision hbar
      integer d_temp
      integer n, i
      double precision dj
      double precision di
      double precision dn
c
      call tce_filename('hbar',filename)
      if (eaf_open(filename, EAF_R, d_temp).ne.0)
     1  call errquit('hbar_restore: eaf problem',0,DISK_ERR)
      if (eaf_read(d_temp,0.0d0,dn,bytes).ne.0)
     1  call errquit('hbar_restore: eaf problem',1,DISK_ERR)
      n = int(dn)
      if (n.gt.0) then
        do i = 1, n
          if (eaf_read(d_temp,dble(bytes)*dble((i-1)*3+1),
     1      di,bytes).ne.0)
     2      call errquit('hbar_restore: eaf problem',2,DISK_ERR)
          if (eaf_read(d_temp,dble(bytes)*dble((i-1)*3+2),
     1      dj,bytes).ne.0)
     2      call errquit('hbar_restore: eaf problem',3,DISK_ERR)
          if (eaf_read(d_temp,dble(bytes)*dble((i-1)*3+3),
     1      hbar,bytes).ne.0)
     2      call errquit('hbar_restore: eaf problem',4,DISK_ERR)
          if ((int(di).eq.ivec).and.(int(dj).eq.jvec)) then
            if (eaf_close(d_temp).ne.0)
     1        call errquit('hbar_restore: eaf problem',9,DISK_ERR)
c           write(LuOut,*) 'Restoring H bar(',ivec,jvec,') = ',hbar
            hbar_restore = .true.
            return
          endif
        enddo
      endif
      if (eaf_close(d_temp).ne.0)
     1  call errquit('hbar_restore: eaf problem',9,DISK_ERR)
      hbar_restore = .false.
      return
      end
      subroutine hbar_store(ivec,jvec,hbar)
c
c $Id$
c
c     Write hbar(ivec,jvec)
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sf.fh"
#include "msgids.fh"
#include "eaf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      character*255 filename
      integer ivec
      integer jvec
      double precision hbar
      integer d_temp
      integer n
      double precision dn
      double precision di
      double precision dj
c
      call tce_filename('hbar',filename)
      if (eaf_open(filename, EAF_RW, d_temp).ne.0)
     1  call errquit('hbar_store: eaf problem',0,DISK_ERR)
      if (eaf_read(d_temp,0.0d0,dn,bytes).ne.0)
     1  call errquit('hbar_store: eaf problem',1,DISK_ERR)
      n = int(dn)
      n = n + 1
      dn = dble(n)
      if (eaf_write(d_temp,0.0d0,dn,bytes).ne.0)
     1  call errquit('hbar_store: eaf problem',2,DISK_ERR)
      di = dble(ivec)
      dj = dble(jvec)
      if (eaf_write(d_temp,dble(bytes)*dble((n-1)*3+1),
     1  di,bytes).ne.0)
     2  call errquit('hbar_store: eaf problem',3,DISK_ERR)
      if (eaf_write(d_temp,dble(bytes)*dble((n-1)*3+2),
     1  dj,bytes).ne.0)
     2  call errquit('hbar_store: eaf problem',4,DISK_ERR)
      if (eaf_write(d_temp,dble(bytes)*dble((n-1)*3+3),
     1  hbar,bytes).ne.0)
     2  call errquit('hbar_store: eaf problem',5,DISK_ERR)
      if (eaf_close(d_temp).ne.0)
     1  call errquit('hbar_store: eaf problem',9,DISK_ERR)
c     write(LuOut,*) 'Storing H bar(',ivec,jvec,') = ',hbar
      return
      end
      subroutine hbar_init
c
c $Id$
c
c     Initialize hbar file
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sf.fh"
#include "msgids.fh"
#include "eaf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      character*255 filename
      integer d_temp
      double precision dn
c
      call tce_filename('hbar',filename)
      if (eaf_open(filename, EAF_RW, d_temp).ne.0)
     1  call errquit('hbar_init: eaf problem',0,DISK_ERR)
      dn = 0.0d0
      if (eaf_write(d_temp,0.0d0,dn,bytes).ne.0)
     1  call errquit('hbar_init: eaf problem',1,DISK_ERR)
      if (eaf_close(d_temp).ne.0)
     1  call errquit('hbar_init: eaf problem',2,DISK_ERR)
      return
      end
      subroutine hbar_tidy
c
c $Id$
c
c     Terminate hbar file
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sf.fh"
#include "msgids.fh"
#include "eaf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      character*255 filename
c
      call tce_filename('hbar',filename)
      if (eaf_delete(filename).ne.0)
     1  call errquit('hbar_tidy: eaf problem',0,DISK_ERR)
      return
      end
