      subroutine dratoga(d_sf)
c
c $Id: dratoga.F,v 1.3 2003-09-16 21:50:11 sohirata Exp $
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
      integer request
      logical nodezero
c
      nodezero = (ga_nodeid().eq.0)
      if (ioalg.eq.5) then
        if (fileisga(d_sf)) then
          write(LuOut,*) 'filename = ',filenames(d_sf)
          write(LuOut,*) 'filehandle = ',filehandles(d_sf)
          call errquit('dratoga: already in GA',0)
        endif
        if (.not.ga_create(mt_dbl,filesizes(d_sf),1,filenames(d_sf),
     1    -1,1,d_temp))
     2    call errquit('dratoga: ga problem',0)
        if (dra_open(filenames(d_sf),DRA_R,d_temp).ne.0)
     1    call errquit('dratoga: dra problem',0)
        if (dra_read(filehandles(d_sf),d_temp,request).ne.0)
     1    call errquit('dratoga: dra problem',1)
        if (dra_wait(request).ne.0)
     1    call errquit('dratoga: dra problem',2)
        if (dra_close(filehandles(d_sf)).ne.0)
     1    call errquit('dratoga: dra problem',3)
        if (dra_delete(filehandles(d_sf)).ne.0)
     1    call errquit('dratoga: dra problem',4)
        filehandles(d_sf) = d_temp
        fileisga(d_sf) = .true.
        if (nodezero.and.util_print('files',print_debug))
     1    call filesummary
      endif
      return
      end
