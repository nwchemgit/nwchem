      logical function tce_prop_start(rtdb)
c
c $Id: tce_prop_start.F,v 1.3 2008-01-18 19:32:19 jhammond Exp $
c
c Wrapper for tce_property
c
c Written by Jeff Hammond
c
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
#include "stdio.fh"
c
      integer rtdb             ! Run-time database
      integer fragment         ! -1 if not fragment MO
      double precision energy
      integer iter, nfrags
      double precision r
      logical nodezero
      logical tce_property
      external tce_property
c
      nodezero = (ga_nodeid().eq.0)
c      if (.not.rtdb_get(rtdb,'tce:fragment',mt_int,1,fragment))
c     1  call errquit('tce_prop_start',0,RTDB_ERR)
c      if (fragment.eq.-1) then
        if(nodezero) write(LuOut,*) 'tce_prop_start'
        tce_prop_start = tce_property(rtdb,.true.)
        return
c      else
c        tce_prop_start = .false.
c        call geom_fragment(rtdb,nfrags)
c        if (nfrags.eq.1) then
c          tce_prop_start = tce_property(rtdb,.true.)
c          return
c        else
c          call errquit('tce_prop: no fragments!',0,RTDB_ERR)
c        endif
c      endif
c 9000 format(1x,'@ Fragment energies')
      end
