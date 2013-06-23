c
c     nuclear md step
c
      subroutine qmd_nucl_step(rtdb,geom,nat,istep,dt,m,r,v,g,esys,ekin)
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
      integer istep              ! step count
      integer nat                ! number of atoms
      double precision dt        ! step size
      double precision m(nat)    ! mass
      double precision r(3,nat)  ! coordinates
      double precision v(3,nat)  ! velocities
      double precision g(3,nat)  ! gradients
      double precision esys      ! system energy
      double precision ekin      ! kinetic energy
c
      integer i
      double precision a(3,nat)  ! acceleration
      double precision dt2
c
      character*30 pname
c
c     preliminaries
      pname='qmd_nucl_step: '
      dt2 = 0.5d0*dt
c
c     for the first step: acceleration and velocity
      if (istep.eq.1) then
        do i=1,nat
           a(1,i) = -g(1,i)/m(i)
           v(1,i) = v(1,i) + a(1,i)*dt2
c
           a(2,i) = -g(2,i)/m(i)
           v(2,i) = v(2,i) + a(2,i)*dt2
c
           a(3,i) = -g(3,i)/m(i)
           v(3,i) = v(3,i) + a(3,i)*dt2
        end do
      end if ! istep 
c
c     choose the integrator
      if (do_veloverlet) then
        call qmd_veloverlet(rtdb,geom,nat,dt,m,r,v,g,esys,ekin)
      else 
        call errquit(pname//'unknown integrator',0,CALC_ERR)
      end if
c
c     persist updated coordinates in geom
      if (.not.geom_cart_coords_set(geom,r))
     $   call errquit(pname//'geom_cart_coords_set',0,GEOM_ERR)
c
c     persist updated velocity in geom
      if (.not.geom_vel_set(geom,v))
     $   call errquit(pname//'geom_vel_set',0,GEOM_ERR)
c
      return
      end
c $Id$
