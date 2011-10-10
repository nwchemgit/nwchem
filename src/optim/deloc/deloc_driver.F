c     ==================================
c     == Delocalized internals driver ==
c     ==================================
c
      logical function deloc_driver(irtdb)
c
      implicit none
c
c     == general includes ==
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "global.fh"
#include "stdio.fh"
c
c     == delocalized internals includes ==
#include "deloc_data.fh"
c
c     == general ==
      integer irtdb   ! input: handle to the runtime database
      integer master
      logical returnstatus
      character*32 pname
c
c     == geometry related ==
      integer igeom  ! local handle to the main geometry object
c
      integer iatom, natom, ndim
c
      integer i_coord,h_coord   ! coordinates
      integer i_grad,h_grad     ! gradients
      integer i_atsymb,h_atsymb ! atomic symbols
      integer i_covrad,h_covrad ! covalent radii
      integer i_atnum,h_atnum   ! atomic numbers
c
      double precision energy
      character*16 atsymb       ! atomic symbol of one atom
      double precision xloc(3)  ! coordinate of one atom
      double precision qloc(3)  ! charge on one atom
c
c     == external functions ==
      logical task_gradient     ! energy/gradient driver
      external task_gradient
      logical getcovalentradii  ! get the covalent radii
      external getcovalentradii
      logical getatomicnumbers  ! get the atomic numbers
      external getatomicnumbers
c
c     == algorithm related ==
      integer maxiter, igeomiter, itype, ierror
      logical converged
      logical periodic
c
c     == preliminaries ==
      pname ="deloc_driver"     ! name of this driver
      master=ga_nodeid()        ! get the master node
      if (master.eq.0) then
        write(luout,*)
        write(luout,*) 
     & "Geometry optimization using delocalized internal coordinates"
        write(luout,*)
      end if
c
c     == allocate necessary memory ==
c
c     == create a local geometry object ==
      if(.not.geom_create(igeom,'geometry'))
     & call errquit(pname//'Failed to create geometry',0, GEOM_ERR)
c
c     == load the current geometry into the geometry object ==
      if(.not.geom_rtdb_load(irtdb,igeom,'geometry'))
     & call errquit(pname//'Failed to load geometry',0, GEOM_ERR)
c
c     == get the number of coordinates in the geometry ==
      if(.not. geom_ncent(igeom,natom) )
     & call errquit(pname//'No coordinates present',0, GEOM_ERR)
c
c     == allocate the coordinates ==
      if(.not.ma_alloc_get(mt_dbl,3*natom,'coords',h_coord,i_coord))
     & call errquit( pname//'Failed to allocate memory for centers',
     & 3*natom, MA_ERR)
c
c     == allocate the gradients ==
      if(.not.ma_alloc_get(mt_dbl,3*natom,'gradients',h_grad,i_grad))
     & call errquit( pname//'Failed to allocate memory for gradients',
     & 3*natom, MA_ERR)
c
c     == allocate the covalent radii ==
      if(.not.ma_alloc_get(mt_dbl,natom,'cov. radii',h_covrad,i_covrad))
     & call errquit( pname//'Failed to allocate memory for cov. radii',
     & natom, MA_ERR)
c
c     == allocate the atomic numbers ==
      if(.not.ma_alloc_get(mt_int,natom,'at. nums',h_atnum,i_atnum))
     & call errquit( pname//'Failed to allocate memory for at. nums',
     & natom, MA_ERR)
c
c     == determine system type: molecular (itype = 0) or periodic (itype > 0) ==
      periodic = .false.
      if(.not. geom_systype_get(igeom,itype)) itype = 0  ! molecular by default
      if (itype.ne.0) periodic = .true.
c
c     == get the number of iterations ==
      if(.not.rtdb_get(irtdb,'driver:nptopt',mt_int,1,maxiter))
     &  maxiter = 20
c
c     == get the covalent radii ==
      if (.not. getcovalentradii(igeom, natom, dbl_mb(i_covrad)))
     & call errquit(pname//'Failed to get cov. radii',0, GEOM_ERR)
c
c     == get the atomic numbers ==
      if (.not. getatomicnumbers(igeom, natom, int_mb(i_atnum)))
     & call errquit(pname//'Failed to get atomic. nums',0, GEOM_ERR)
c
c     == test print out ==
c      write(luout,*) "natom:",natom
c      do iatom = 1,natom
c       write(luout,*) dbl_mb(i_coord+3*(iatom-1)),
c     &                dbl_mb(i_coord+3*(iatom-1)+1),
c     &                dbl_mb(i_coord+3*(iatom-1)+2),
c     &                dbl_mb(i_covrad+iatom-1),
c     &                dbl_mb(i_atnum+iatom-1)
c      end do
c      stop
c       
c     == start optimization loop ==
      igeomiter = 0  ! initialize iteration counter
  20  continue 
c
c     == get the coordinates from the geometry object ==
      if(.not.geom_cart_coords_get(igeom,dbl_mb(i_coord)))
     & call errquit( pname//'Failed to get coordinates',0, GEOM_ERR)
c
c     == calculate the energy/gradient ==
      if (.not. task_gradient(irtdb)) 
     &   call errquit(pname//'task_gradient failed',0, GEOM_ERR)
c
c     == get the energy & gradients from the runtime database ==
      if (.not. rtdb_get(irtdb,'task:energy', mt_dbl, 1, energy))
     &   call errquit(pname//'Could not get energy',0, RTDB_ERR)
      if (.not. rtdb_get(irtdb, 'task:gradient', mt_dbl, 3*natom,
     &  dbl_mb(i_grad))) 
     & call errquit(pname//'Could not get gradient',0,RTDB_ERR)
c
c     == call the main delocalized internals driver
      converged = .false.
      ierror = 0
      call deloc(periodic, natom, ndim, 
     &           dbl_mb(i_coord),   ! coordinates
     &           dbl_mb(i_grad),    ! gradients
     &           dbl_mb(i_covrad),  ! covalent radii
     &           int_mb(i_atnum),   ! atomic numbers
     &           converged, ierror)
c
c     == decide what to do if ierror > 0 == will handle later
c     switch to another scheme after updating the database with the 
c     last good geometry
c
c     == update the geometry ==
      if(.not.geom_cart_coords_set(igeom,dbl_mb(i_coord)))  
     & call errquit(pname//'Failed to set coordinates',0, GEOM_ERR)
c
c     == update the runtime database ==
      if(.not.geom_rtdb_store(irtdb,igeom,'geometry'))
     & call errquit(pname//'Failed to store geometry',0, GEOM_ERR)
c
c     == update geometry optimization counter ==
      igeomiter = igeomiter + 1
c
c     == convergence check ==
      if (igeomiter.le.abs(maxiter).and. (.not. converged) ) then
        goto 20  ! go back to the main optimization loop
      else
       call errquit(
     & pname//'Geometry optimization did not converge',0,RTDB_ERR)
      end if
c
c     == end optimization loop ==
c
c     == deallocate necessary memory ==
      if(.not.ma_free_heap(h_grad))
     & call errquit(pname//'Failed to deallocate stack gradients',0,
     &       MA_ERR)
c
      if(.not.ma_free_heap(h_coord))
     & call errquit(pname//'Failed to deallocate stack coords',0,
     &       MA_ERR)
c
      if(.not.ma_free_heap(h_covrad))
     & call errquit(pname//'Failed to deallocate stack cov. radii',0,
     &       MA_ERR)
c
      if(.not.ma_free_heap(h_atnum))
     & call errquit(pname//'Failed to deallocate stack atomic. nums',0,
     &       MA_ERR)
c
       if(.not.geom_destroy(igeom))
     & call errquit(pname//'Failed to destroy geometry',0, GEOM_ERR)
c
      deloc_driver = .true.
      return
      end
c $Id$
