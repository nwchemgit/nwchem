      subroutine dratoga(d_sf,size)
c
c $Id: dratoga.F,v 1.1 2003-08-01 22:59:51 sohirata Exp $
c
c     Load data from DRA to GA
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
      integer size
      integer request
      integer i
      logical nodezero
      logical used
c
      nodezero = (ga_nodeid().eq.0)
      if (ioalg.eq.5) then
        if (fileisga(d_sf))
     1    call errquit('dratoga: already in GA',0)
        if (.not.ga_create(mt_dbl,size,1,filenames(d_sf),
     1    -1,1,d_temp))
     2    call errquit('dratoga: ga problem',0)
        filehandles(d_sf) = d_temp
        fileisga(d_sf) = .true.
        if (dra_open(filenames(d_sf),DRA_R,d_temp).ne.0)
     1    call errquit('dratoga: dra problem',0)
        if (dra_read(filehandles(d_sf),d_temp,request).ne.0)
     1    call errquit('dratoga: dra problem',1)
        if (dra_wait(request).ne.0)
     1    call errquit('dratoga: dra problem',2)
        if (dra_close(d_temp).ne.0)
     1    call errquit('dratoga: dra problem',3)
        if (nodezero.and.util_print('files',print_debug))
     1    call filesummary
      endif
      return
      end



      subroutine gatodra(d_sf,size)
c
c $Id: dratoga.F,v 1.1 2003-08-01 22:59:51 sohirata Exp $
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
      integer size
      integer request
      integer i
      logical nodezero
      logical used
c
      nodezero = (ga_nodeid().eq.0)
      if (ioalg.eq.5) then
        if (.not.fileisga(d_sf))
     1    call errquit('gatodra: no GA found',0)
        if (dra_open(filenames(d_sf),DRA_W,d_temp).ne.0)
     1    call errquit('gatodra: dra problem',0)
        if (dra_write(filehandles(d_sf),d_temp,request).ne.0)
     1    call errquit('gatodra: dra problem',1)
        if (dra_wait(request).ne.0)
     1    call errquit('gatodra: dra problem',2)
        if (dra_close(d_temp).ne.0)
     1    call errquit('gatodra: dra problem',3)
        if (.not.ga_destroy(filehandles(d_sf)))
     1    call errquit('gatodra: GA problem',0)
        fileisga(d_sf) = .false.
        if (nodezero.and.util_print('files',print_debug))
     1    call filesummary
      endif
      return
      end
