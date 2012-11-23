c
c     write restart
c
      subroutine qmd_write_rst(rstfilename,rstfileunit,
     &  nat,iframe,esys,ekin,geom)
c
      implicit none
c
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "stdio.fh"
c
#include "qmd_common.fh"
c
      character*(*) rstfilename
      integer rstfileunit
      integer nat
      integer iframe
      double precision esys,ekin,etotal
      integer geom
c
      integer i
      double precision r(3,nat)
      double precision v(3,nat)
c
      character*32 pname
c
      pname = "qmd_write_rst: "
c
c     open file
      open(unit=rstfileunit,file=rstfilename,form='formatted',
     &        position='rewind',status='unknown')
c
c     write number of atoms, frame number and energy
      write(rstfileunit,110) nat
      etotal = esys + ekin
      write(rstfileunit,111) iframe,etotal,esys,ekin
      write(rstfileunit,112) targ_temp,final_temp
      if (.not.geom_cart_coords_get(geom,r))
     $  call errquit(pname//'geom_cart_coords_get',0,GEOM_ERR)
      if (.not.geom_vel_get(geom,v))
     $   call errquit(pname//'geom_vel_set',0,GEOM_ERR)
c
      write(rstfileunit,*) "Coordinates"
      do i=1,nat
         write(rstfileunit,113) r(1,i),r(2,i),r(3,i)
      end do
c
      write(rstfileunit,*) "Velocities"
      do i=1,nat
         write(rstfileunit,113) v(1,i),v(2,i),v(3,i)
      end do
c
      if (rstfileunit.ne.0) close(unit=rstfileunit)
c
  110 format(I5)
  111 format(I5,3f22.12)
  112 format(2f22.12)
  113 format(3f22.12)
c
      return
      end
