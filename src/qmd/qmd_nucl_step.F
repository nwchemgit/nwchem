c
c     nuclear md step
c
      subroutine qmd_nucl_step(rtdb,geom,nat,dt,m,r,v,g,esys,ekin)
c
      implicit none
c
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "geom.fh"
c
#include "qmd_common.fh"
c
      integer rtdb               ! rtdb
      integer geom               ! geom
      integer nat                ! number of atoms
      double precision dt        ! step size
      double precision m(nat)    ! mass
      double precision r(3,nat)  ! coordinates
      double precision v(3,nat)  ! velocities
      double precision g(3,nat)  ! gradients
      double precision esys      ! system energy
      double precision ekin      ! kinetic energy
c
      character*30 pname
      pname='qmd_nucl_step: '
c
c     choose the integrator
      if (do_vv) then
        call qmd_veloverlet(rtdb,geom,nat,dt,m,r,v,g,esys,ekin)
      else 
        call errquit(pname//'unknown integrator',0,CALC_ERR)
      end if
c
      return
      end
