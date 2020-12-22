c
c     qmd_gradient: wrapper for gradient
c
      subroutine qmd_gradient(rtdb,nat,istep_nucl,m,g,r,energy,
     E     ext_forces)
c
      implicit none
c
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "inp.fh"
c
      integer rtdb               ! rtdb
      integer nat               ! number of atoms
      integer istep_nucl
      double precision m(nat)    ! masses
      double precision g(3,nat)  ! gradients
      double precision r(3,nat)  ! coordinates
      double precision energy    ! system energy
      character*(*) ext_forces
c     
      integer i,j
c
      character*30 pname
c
      logical task_gradient
      external task_gradient
      logical task_energy
      external task_energy
c
c     preliminaries
      pname = "qmd_gradient: "
c
c     energy/gradient calculation
      energy = 0.d0
      if (.not. task_gradient(rtdb))
     &   call errquit(pname//'task_gradient',0,CALC_ERR)
      if(.not.rtdb_get(rtdb,'task:gradient',mt_dbl,3*nat,g))
     &   call errquit(pname//'rtdb_get gradient',0,RTDB_ERR)
      if(.not.rtdb_get(rtdb,'task:energy',mt_dbl,1,energy))
     &     call errquit(pname//'rtdb_get energy',0,RTDB_ERR)

      if(ext_forces(1:inp_strlen(ext_forces)).ne.'none') then
        call qmd_ext_gradient(nat,istep_nucl,energy,m,g,r,
     E       ext_forces)
      endif
      return
      end
c $Id$
