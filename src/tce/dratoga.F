      subroutine dratoga(d_sf)
c
c $Id$
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
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer d_sf
      integer d_ga
      integer d_dra
      integer request
      logical nodezero
c
      nodezero = (ga_nodeid().eq.0)
      if (ioalg.eq.5) then
        if (fileisga(d_sf)) then
          write(LuOut,*) 'filename = ',filenames(d_sf)
          write(LuOut,*) 'filehandle = ',filehandles(d_sf)
          call errquit('dratoga: already in GA',0,UNKNOWN_ERR)
        endif
        if (.not.ga_create(mt_dbl,filesizes(d_sf),1,filenames(d_sf),
     1    -1,1,d_ga))
     2    call errquit('dratoga: ga problem',0,GA_ERR)
        if (dra_open(filenames(d_sf),DRA_R,d_dra).ne.0)
     1    call errquit('dratoga: dra problem',0,DISK_ERR)
        if (dra_read(d_ga,d_dra,request).ne.0)
     1    call errquit('dratoga: dra problem',1,DISK_ERR)
        if (dra_wait(request).ne.0)
     1    call errquit('dratoga: dra problem',2,DISK_ERR)
        if (dra_delete(d_dra).ne.0)
     1    call errquit('dratoga: dra problem',4,DISK_ERR)
        filehandles(d_sf) = d_ga
        fileisga(d_sf) = .true.
        if (nodezero.and.util_print('files',print_debug))
     1    call filesummary
      endif
      return
      end
