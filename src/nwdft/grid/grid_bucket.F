      subroutine grid_buckinit(g_buck,g_buckind)
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
      integer ncubes,i,nctrs,dum
      integer grid_ncubes,ipts
      external grid_ncubes
c
      if(.not.geom_ncent(geom, nctrs))
     . call errquit(' gridb_buck: geom cent failed',0)
c
c     get ncubes
c
      ncubes=0
      do i=1,ntypes
         ncubes_type(i)=grid_ncubes(dble(rad_cutoff(1,i)),
     ,        .false.,dum)
         ncubes=max(ncubes,ncubes_type(i))
      enddo
      write(6,*) ' NCUBES ',ncubes
      if (.not.MA_alloc_Get(MT_int,ntypes*ncubes,'howmany cubes',
     , l_nxyz,k_nxyz))
     , call errquit('grid_buckinit: cannot allocate nxyz',0)
      do i=1,ntypes
         ipts=k_nxyz+(i-1)*ncubes
         dum=grid_ncubes(dble(rad_cutoff(1,i)),.true.,int_mb(ipts))
      enddo
      if(.not. ga_create(mt_dbl,ncenters*ncubes,mxpt_rect*4,
     ,     'grid bucket', 0,mxpt_rect*4,g_buck))
     . call errquit(' gridb_buck: ga_create failed',0)
      if(.not. ga_create(mt_int,ncenters,ncubes,
     ,     'grid bucket', 0,ncubes,g_buckind))
     . call errquit(' gridb_buck: ga_create failed',1)
      call ga_zero(g_buck)
      call ga_zero(g_buckind)
      return
      end
      subroutine grid_buckend(g_buck,g_buckind)
      integer g_buck,g_buckind
c     
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "cgridfile.fh"
c
       if (.not. MA_free_heap(l_nxyz))
     . call errquit(' grid_buckend:cannot popstack',0)

      if(.not. ga_destroy(g_buck))
     . call errquit(' gridb_buck: ga_destroy failed',0)
      if(.not. ga_destroy(g_buckind))
     . call errquit(' gridb_buck: ga_destroy failed',0)
      return
      end
      subroutine grid_bucket()
      implicit none
#include "cgridfile.fh"
      return
      end
