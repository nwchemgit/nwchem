c
c     qmd_gradient: wrapper for gradient
c
      subroutine qmd_gradient(rtdb,nat,g,energy)
c
      implicit none
c
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
c
      integer rtdb               ! rtdb
      integer nat                ! number of atoms
      double precision g(3,nat)  ! gradients
      double precision energy    ! system energy
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
     &   call errquit(pname//'rtdb_get energy',0,RTDB_ERR)
c
      return
      end
c $Id$
