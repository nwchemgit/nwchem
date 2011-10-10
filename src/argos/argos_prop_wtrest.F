      subroutine argos_prop_wtrest(lfn)
c
      implicit none
c
#include "argos_prop.fh"
c
      integer lfn
c
      integer i
c
      if(me.ne.0) return
c
      write(lfn,1000)
 1000 format('restart properties')
      write(lfn,1001) maxp,nsum,nsumt
 1001 format(3i7)
      write(lfn,1002) tsum,t2sum,tsumt,t2sumt
 1002 format(4e20.12)
      write(lfn,1003) (psum(i),i=1,maxp)
      write(lfn,1003) (p2sum(i),i=1,maxp)
      write(lfn,1003) (pslop(i),i=1,maxp)
      write(lfn,1003) (psumt(i),i=1,maxp)
      write(lfn,1003) (p2sumt(i),i=1,maxp)
      write(lfn,1003) (pslopt(i),i=1,maxp)
 1003 format(4e20.12)
c
      return
      end
c $Id$
