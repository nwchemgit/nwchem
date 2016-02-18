c
c     write restart
c
      subroutine qmd_write_rst(rstfilename,rstfileunit,
     &  n,iframe,ekin,geom,state,nstates,zden,coup,do_namd)
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
      integer n
      integer iframe
      double precision ekin
      integer geom
      integer state
      integer nstates
      double precision coup(nstates,nstates)
      double complex zden(nstates,nstates)
      logical do_namd
c
      integer i, j
      double precision r(3,n)
      double precision v(3,n)
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
      write(rstfileunit,110) n
      write(rstfileunit,111) iframe,ekin
      if (.not.geom_cart_coords_get(geom,r))
     $  call errquit(pname//'geom_cart_coords_get',0,GEOM_ERR)
      if (.not.geom_vel_get(geom,v))
     $   call errquit(pname//'geom_vel_set',0,GEOM_ERR)
c
      write(rstfileunit,*) "Coordinates"
      do i=1,n
         write(rstfileunit,113) r(1,i),r(2,i),r(3,i)
      end do
c
      write(rstfileunit,*) "Velocities"
      do i=1,n
         write(rstfileunit,113) v(1,i),v(2,i),v(3,i)
      end do
c
      if (do_namd) then
        write(rstfileunit,*) "NAMD" 
c       we have input and output set up so state runs from 0 to nstates-1
c       inside the code state runs from 1 to nstates
        write(rstfileunit,114) state-1, nstates
        write(rstfileunit,*) "Density"
        do i=1,nstates
          write(rstfileunit,*) (zden(i,j),j=1,nstates)
        end do
        write(rstfileunit,*) "Coupling"
        do i=1,nstates
          write(rstfileunit,*) (coup(i,j),j=1,nstates)
        end do
      end if
c
      if (rstfileunit.ne.0) close(unit=rstfileunit)
c
  110 format(I5)
  111 format(I10,3f22.12)
  112 format(2f22.12)
  113 format(3f22.12)
  114 format(2I5)
c
      return
      end
c $Id$
