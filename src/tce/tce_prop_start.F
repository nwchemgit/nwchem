      logical function tce_prop_start(rtdb)
c
c $Id: tce_prop_start.F,v 1.1 2008-01-08 18:19:21 jhammond Exp $
c
c Wrapper for tce_property for fragment MO calculations.
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
      integer maxiter
      parameter (maxiter = 40)
      double precision residual
      double precision max_residual
      double precision damping
      double precision thresh
      logical nodezero
      logical tce_property
      external tce_property
c
      nodezero = (ga_nodeid().eq.0)
      if (.not.rtdb_get(rtdb,'tce:fragment',mt_int,1,fragment))
     1  call errquit('tce_prop_start',0,RTDB_ERR)
      if (fragment.eq.-1) then
        tce_prop_start = tce_property(rtdb,.true.)
        return
      else
        tce_prop_start = .false.
        call geom_fragment(rtdb,nfrags)
        if (nfrags.eq.1) then
          tce_prop_start = tce_property(rtdb,.true.)
          return
        else
          call errquit('tce_prop: no fragments!',0,RTDB_ERR)
        endif
      endif
 9000 format(1x,'@ Fragment energies')
      end
