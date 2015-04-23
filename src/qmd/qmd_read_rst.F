c
c     qmd_read_rst: read rstfile
c
      subroutine qmd_read_rst(rstfilename,rstfileunit,n,v,r,ekin,iframe)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
c
#include "qmd_common.fh"
c
      character*(*) rstfilename
      integer rstfileunit
      integer ok
      integer iframe           ! last step
      integer n                ! number of atoms
      double precision v(3,n)  ! atomic velocities
      double precision r(3,n)  ! atomic positions
      double precision ekin    ! kinetic energy
c
      integer i,j
c
      character*32 pname
      pname='qmd_read_rst: '
c
c     open file
      open(unit=rstfileunit,file=rstfilename,form='formatted',
     $         status='old')
c     read through header and coordinates
c     number of atoms
      read(rstfileunit,*)
c     step number and energies
      read(rstfileunit,*) iframe,ekin
c     rst stores the last completed step
c     want to start with next step
      iframe=iframe+1
c     coordinates label
      read(rstfileunit,*)
c     read coordinates
      read(rstfileunit,*) ((r(j,i),j=1,3),i=1,n)
c     velocity label
      read(rstfileunit,*)
c     read velocities      
      read(rstfileunit,*) ((v(j,i),j=1,3),i=1,n)
c 
      if (rstfileunit.ne.0) close(unit=rstfileunit)
c
      return
      end
