      subroutine gatodra(d_sf)
c
c $Id$
c
c     Unload data from GA to DRA
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sf.fh"
#include "eaf.fh"
#include "dra.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer d_sf
      integer d_dra
      integer request
      logical nodezero
c
      nodezero = (ga_nodeid().eq.0)
      if (ioalg.eq.5) then
        if (.not.fileisga(d_sf))
     1    call errquit('gatodra: no GA found',0,UNKNOWN_ERR)
        if (dra_create(mt_dbl,filesizes(d_sf),1,filenames(d_sf),
     1    filenames(d_sf),DRA_RW,-1,1,d_dra).ne.0)
     2    call errquit ('gatodra: dra problem',1,DISK_ERR)
        if (dra_write(filehandles(d_sf),d_dra,request).ne.0)
     1    call errquit('gatodra: dra problem',1,DISK_ERR)
        if (dra_wait(request).ne.0)
     1    call errquit('gatodra: dra problem',2,DISK_ERR)
        if (dra_close(d_dra).ne.0)
     1    call errquit('gatodra: dra problem',3,DISK_ERR)
        if (.not.ga_destroy(filehandles(d_sf)))
     1    call errquit('gatodra: GA problem',0,GA_ERR)
        filehandles(d_sf) = d_dra
        fileisga(d_sf) = .false.
        if (nodezero.and.util_print('files',print_debug))
     1    call filesummary
      endif
      return
      end
