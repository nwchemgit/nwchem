      subroutine grid_buckinit(rtdb)
c
c$Id: grid_bucket.F,v 1.11 2000-04-05 01:41:22 edo Exp $
c
      implicit none
c     
#include "mafdecls.fh"
#include "global.fh"
#include "cdft.fh"
#include "cgridfile.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "grid_cube.fh"
c
      integer rtdb
      integer i,nctrs,dum
      integer grid_ncubes,ipts
      double precision grid_lcube
      external grid_ncubes,grid_lcube
c
      if(.not.geom_ncent(geom, nctrs))
     . call errquit(' gridb_buck: geom cent failed',0)
c
c     get cubes definition
c
      if (.not. rtdb_get(rtdb, 'dft:l_cube',
     $     mt_dbl, 1, l_cube)) l_cube = 5.5d0
      if (.not. rtdb_get(rtdb, 'dft:sub4',
     $     mt_log, 1, sub4)) sub4=.false.
      shft=0
      if(sub4) shft=64
      if(ga_nodeid().eq.0) write(LuOut,'(A,F12.6)') 
     .     ' l_cube',l_cube
      if(ga_nodeid().eq.0.and.sub4) write(LuOut,'(A,I4)') 
     .     ' central subdivision by',shft
c
c     get ncubes
c
      ncubesmx=1
      do i=1,ntypes
         lcube(i)=grid_lcube(dble(rad_cutoff(1,i)))
         ncubes_type(i)=grid_ncubes(dble(rad_cutoff(1,i)),
     ,        .false.,dum,lcube(i))
         ncubesmx=max(ncubesmx,ncubes_type(i))
      enddo
      if(ga_nodeid().eq.0) write(LuOut,'(A,I4)') ' NCUBESMX ',ncubesmx
      if (.not.MA_alloc_Get(MT_int,ntypes*ncubesmx,'howmany cubes',
     ,     l_nxyz,k_nxyz))
     ,     call errquit('grid_buckinit: cannot allocate nxyz',0)
      do i=1,ntypes
         ipts=k_nxyz+(i-1)*ncubesmx
         dum=grid_ncubes(dble(rad_cutoff(1,i)),.true.,int_mb(ipts),
     ,        lcube(i))
      enddo
      return
      end
      subroutine grid_buckend()
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
