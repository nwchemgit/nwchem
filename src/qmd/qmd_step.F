c
c     nuclear md step
c
      subroutine qmd_step(rtdb,integrator,nat,dt,m,r,v,g,esys)
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
      character*20 integrator
      integer nat                ! number of atoms
c
      double precision dt        ! step size
      double precision m(nat)    ! mass
      double precision r(3,nat)  ! coordinates
      double precision v(3,nat)  ! velocities
      double precision g(3,nat)  ! gradients
      double precision a(3,nat)  ! acceleration
c
      integer i,j
      double precision dt2
c
      double precision etot,esys,ekin
c
      character*30 pname
      pname='qmd_step: '
c
c     choose the integrator
      if (integrator.eq.'velocity_verlet') then
         call qmd_verlet(rtdb,nat,dt,m,r,v,g,esys)
      else 
         call errquit(pname//'unknown integrator',0,CALC_ERR)
      end if
c
      return
      end
