      subroutine grid_buckinit(rtdb)
c
c$Id: grid_bucket.F,v 1.17 2000-07-21 01:46:58 edo Exp $
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
#include "util.fh"
#include "grid_cube.fh"
c
      integer rtdb
      integer i,nctrs,dum
      integer grid_ncubes,grid_nxyz,ipts
      double precision grid_lcube,max_rad
      integer nxyz
      integer grid_nsub4
      external grid_nsub4
      external grid_ncubes,grid_lcube,grid_nxyz
      integer max_xyz
      parameter (max_xyz=80000)
      logical oprint
c
      oprint = util_print('quadrature', print_high)
      if(.not.geom_ncent(geom, nctrs))
     . call errquit(' gridb_buck: geom cent failed',0)
c
c     get cubes definition
c
      if (.not. rtdb_get(rtdb, 'dft:l_cube',
     $     mt_dbl, 1, l_cube)) l_cube = 11d0
      if (.not. rtdb_get(rtdb, 'dft:sub4',
     $     mt_log, 1, sub4)) sub4=.true.
      if(ga_nodeid().eq.0.and.oprint) write(LuOut,'(A,F8.2)') 
     .     ' l_cube input',l_cube
c
c     get ncubes
c
      ncubesmx=1
      if(ga_nodeid().eq.0.and.oprint) write(LuOut,222) 
 222  format(/,
     &     1x,' atom type   l_cube   no. cubes   segmts   rad   sub'/
     &     1x,'----------- -------- ----------- -------- ----- -----')
 1    continue
c
c     check if we do not get too many cubes
c
      ncubesmx=0
      max_rad=0d0 
      do i=1,ntypes
         lcube(i)=grid_lcube(dble(rad_cutoff(1,i)))
         nxyz=grid_nxyz(dble(rad_cutoff(1,i)),lcube(i))
         ncubesmx=max(ncubesmx,nxyz**3)
         max_rad=max(max_rad,dble(rad_cutoff(1,i))) 
      enddo

      if(ntypes*ncubesmx.gt.max_xyz) then 
         l_cube=1.5d0*max_rad/max_xyz**(1d0/3d0)
         if (.not. rtdb_put(rtdb, 'dft:l_cube', 
     $     mt_dbl, 1, l_cube))  
     .        call errquit('gridbucki: rtdbput failed',0) 
      endif 
c
      do i=1,ntypes
         lcube(i)=grid_lcube(dble(rad_cutoff(1,i)))
         nxyz=grid_nxyz(dble(rad_cutoff(1,i)),lcube(i))
         shft=0
         if(sub4) then
            nsub4=grid_nsub4(nxyz)
            shft=nsub4**3
         endif
         ncubes_type(i)=grid_ncubes(dble(rad_cutoff(1,i)),
     ,        .false.,dum,lcube(i))
         ncubesmx=max(ncubesmx,ncubes_type(i))
      enddo
      if(ntypes*ncubesmx.gt.max_xyz) then
         l_cube=l_cube*2d0
         if (.not. rtdb_put(rtdb, 'dft:l_cube',
     $     mt_dbl, 1, l_cube)) 
     .        call errquit('gridbucki: rtdbput failed',0)
         goto 1
      endif
      if (.not.MA_alloc_Get(MT_int,ntypes*ncubesmx,'howmany cubes',
     ,     l_nxyz,k_nxyz))
     ,     call errquit('grid_buckinit: cannot allocate nxyz',
     ,     ntypes*ncubesmx)
      do i=1,ntypes
         ipts=k_nxyz+(i-1)*ncubesmx
         dum=grid_ncubes(dble(rad_cutoff(1,i)),.true.,int_mb(ipts),
     ,        lcube(i))
      if(ga_nodeid().eq.0.and.oprint) 
     .        write(LuOut,22) i,lcube(i),ncubes_type(i),
     ,     nxyz,rad_cutoff(1,i),shft
 22   format(5x,i4,2x,f8.2,4x,i6,7x,i3,5x,i3,2x,i4)
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
