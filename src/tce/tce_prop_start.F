      logical function tce_prop_start(rtdb)
c
c $Id$
c
c Wrapper for tce_property - the only reason to keep this around is
c because eventually there may be code wrapped around the real property
c code to make its features work with QM/MM, etc.
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
c      integer fragment         ! -1 if not fragment MO
c      double precision energy
c      integer iter, nfrags
c      double precision r
      logical nodezero
      logical tce_property
      external tce_property
c
      nodezero = (ga_nodeid().eq.0)
c      if (.not.rtdb_get(rtdb,'tce:fragment',mt_int,1,fragment))
c     1  call errquit('tce_prop_start',0,RTDB_ERR)
c      if (fragment.eq.-1) then
        if(nodezero) write(LuOut,*) 'tce_prop_start'
        tce_prop_start = tce_property(rtdb)
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
