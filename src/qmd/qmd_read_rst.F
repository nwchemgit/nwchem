c
c     qmd_read_rst: read rstfile
c
      subroutine qmd_read_rst(rstfilename,rstfileunit,n,v,r,ekin,iframe,
     &                        state,nstates,do_namd,zden,coup)
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
      integer state
      integer nstates
      logical do_namd
c
c     if do_namd false, these aren't allocated
      double precision coup(nstates,nstates)
      double complex zden(nstates,nstates)
c
      integer i,j
      integer nstates_verify
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
      if (do_namd) then
c       NAMD label
        read(rstfileunit,*)
c       current state and number of states
        read(rstfileunit,*) state, nstates_verify
c       we have input and output set up so state runs from 0 to nstates-1
c       inside the code state runs from 1 to nstates
        state=state+1
        if (nstates.ne.nstates_verify)
     $    call errquit(pname//"nstates changed",0,0)
c       density label
        read(rstfileunit,*)
c       read density matrix
        read(rstfileunit,*) ((zden(i,j),j=1,nstates),i=1,nstates)
c       coupling label
        read(rstfileunit,*)
c       read coupling matrix
        read(rstfileunit,*) ((coup(i,j),j=1,nstates),i=1,nstates)
      end if
c
  110 format(2es25.17)
c 
      if (rstfileunit.ne.0) close(unit=rstfileunit)
c
      return
      end
