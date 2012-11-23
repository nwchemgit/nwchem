c
c     write trajectory: xyz file in angstroms
c
      subroutine qmd_write_trj(xyzfilename,xyzfileunit,
     &           mode_append,nat,iframe,etotal,geom)
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
      character*(*) xyzfilename
      integer xyzfileunit
      logical mode_append
      integer nat
      integer iframe
      double precision etotal
      integer geom
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
      write(xyzfileunit,110) nat
      write(xyzfileunit,111) iframe,etotal
c
      if (.not.geom_print_xyzq(geom, xyzfileunit)) 
     &  call errquit(pname//'geom write failed',0,GEOM_ERR)
c
      if (xyzfileunit.ne.0) close(unit=xyzfileunit)
c
  110 format(I5)
  111 format(I5,f22.12)
c
      return
      end
