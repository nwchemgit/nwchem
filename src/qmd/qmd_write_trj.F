c
c     qmd_write_trj: write trajectory: xyz file in angstroms
c
      subroutine qmd_write_trj(xyzfilename,xyzfileunit,
     &           mode_append,n,iframe,etotal,geom,dipole)
c
      implicit none
c
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "stdio.fh"
#include "util_params.fh"
c
#include "qmd_common.fh"
c
      character*(*) xyzfilename
      integer xyzfileunit
      logical mode_append
      integer i, j
      integer n
      integer iframe
      double precision etotal
      integer geom
      character(16) t(n)
      double precision r(3,n)
      double precision v(3,n)
      double precision q(n)
      double precision dipole(3)
c
      character*32 pname
c
      pname = "qmd_write_traj: "
c
c     open file
      if (mode_append) then
         open(unit=xyzfileunit,file=xyzfilename,form='formatted',
     &        position='append',status='unknown')
      else
         open(unit=xyzfileunit,file=xyzfilename,form='formatted',
     &        position='rewind',status='unknown')
      end if
c
c     write number of atoms, frame number and energy
      write(xyzfileunit,110) n
      write(xyzfileunit,111) iframe,etotal,dipole(:)
c
c     get number of centers, tags, coordinates, and charges
      if (.not.geom_cart_get(geom,i,t,r,q))
     $  call errquit(pname//'geom_cart_coords_get',0,GEOM_ERR)
c     get velocities
      if (.not.geom_vel_get(geom,v))
     $   call errquit(pname//'geom_vel_set',0,GEOM_ERR)
c
      do i=1,n
       write(xyzfileunit,4) t(i),r(:,i)*cau2ang,v(:,i)*cau2ang/cau2fs
      end do
c
      if (xyzfileunit.ne.0) close(unit=xyzfileunit)
c
 4    format(1x,a12,3f15.8,2x,3f15.8)
 110  format(I5)
 111  format(I10,f20.10,3es15.6)
c
      return
      end
c $Id$
