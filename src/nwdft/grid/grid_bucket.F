      subroutine grid_buckinit(g_buck,g_buckind)
c
c$Id: grid_bucket.F,v 1.6 2000-03-18 02:14:58 edo Exp $
c
      implicit none
      integer g_buck,g_buckind
c     
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "cdft.fh"
#include "cgridfile.fh"
#include "geom.fh"
c
      integer i,nctrs,dum
      integer grid_ncubes,ipts
      external grid_ncubes
c
      if(.not.geom_ncent(geom, nctrs))
     . call errquit(' gridb_buck: geom cent failed',0)
c
c     get ncubes
c
      ncubesmx=1
      do i=1,ntypes
         ncubes_type(i)=grid_ncubes(dble(rad_cutoff(1,i)),
     ,        .false.,dum)
         ncubesmx=max(ncubesmx,ncubes_type(i))
      enddo
      write(6,*) ' NCUBES ',ncubesmx
      if (.not.MA_alloc_Get(MT_int,ntypes*ncubesmx,'howmany cubes',
     ,     l_nxyz,k_nxyz))
     ,     call errquit('grid_buckinit: cannot allocate nxyz',0)
      do i=1,ntypes
         ipts=k_nxyz+(i-1)*ncubesmx
         dum=grid_ncubes(dble(rad_cutoff(1,i)),.true.,int_mb(ipts))
      enddo
      return
      end
      subroutine grid_buckend(g_buck,g_buckind)
      integer g_buck,g_buckind
c     
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "cdft.fh"
c
      if (.not. MA_free_heap(l_nxyz))
     .     call errquit(' grid_buckend:cannot popstack',0)
      return
      end
      subroutine grid_bucket()
      implicit none
#include "cgridfile.fh"
      return
      end
