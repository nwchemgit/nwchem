      subroutine gatodra(d_sf)
c
c $Id: gatodra.F,v 1.2 2003-09-16 21:50:11 sohirata Exp $
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
#include "tce.fh"
      integer d_sf
      integer d_temp
      integer request
      logical nodezero
c
      nodezero = (ga_nodeid().eq.0)
      if (ioalg.eq.5) then
        if (.not.fileisga(d_sf))
     1    call errquit('gatodra: no GA found',0)
        if (dra_create(mt_dbl,filesizes(d_sf),1,filenames(d_sf),
     1    filenames(d_sf),DRA_RW,-1,1,d_temp).ne.0)
     2    call errquit ('gatodra: dra problem',1)
        if (dra_write(filehandles(d_sf),d_temp,request).ne.0)
     1    call errquit('gatodra: dra problem',1)
        if (dra_wait(request).ne.0)
     1    call errquit('gatodra: dra problem',2)
        if (dra_close(d_temp).ne.0)
     1    call errquit('gatodra: dra problem',3)
        if (.not.ga_destroy(filehandles(d_sf)))
     1    call errquit('gatodra: GA problem',0)
        filehandles(d_sf) = d_temp
        fileisga(d_sf) = .false.
        if (nodezero.and.util_print('files',print_debug))
     1    call filesummary
      endif
      return
      end
