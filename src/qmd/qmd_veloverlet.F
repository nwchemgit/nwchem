c
c     velocity verlet md step
c
      subroutine qmd_veloverlet(rtdb,geom,nat,dt,m,r,v,g,esys,ekin)
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
      double precision a(3,nat)  ! acceleration
      integer i,j
      double precision dtsq,dtsq2,dt2
c
      character*30 pname
c
      logical task_gradient
      external task_gradient
c
c     preliminaries
      pname = "qmd_veloverlet: "
      dtsq = dt*dt
      dtsq2 = 0.5d0*dtsq
      dt2 = 0.5d0*dt
c
c     instantaneous kinetic energy
      call qmd_kinetic(nat,m,v,ekin)
c
c     apply the nose-hoover thermostat at half step
      if (do_nosehoover.and.ekin.gt.0.d0) 
     &    call qmd_nosehoover(nat,dt,v,ekin)
c
c     calculate acceleration, update positions (r) and velocities (v) at dt/2
      do i = 1,nat
        a(1,i) = -g(1,i)/m(i)
        r(1,i) = r(1,i) + v(1,i)*dt + a(1,i)*dtsq2
        v(1,i) = v(1,i) + a(1,i)*dt2
c
        a(2,i) = -g(2,i)/m(i)
        r(2,i) = r(2,i) + v(2,i)*dt + a(2,i)*dtsq2
        v(2,i) = v(2,i) + a(2,i)*dt2
c
        a(3,i) = -g(3,i)/m(i)
        r(3,i) = r(3,i) + v(3,i)*dt + a(3,i)*dtsq2
        v(3,i) = v(3,i) + a(3,i)*dt2
      end do
c
c     update geom with coords for energy-gradient calculation
      if (.not.geom_cart_coords_set(geom,r))
     $   call errquit(pname//'geom_cart_coords_set',0,GEOM_ERR)
c
c     store updated geometry
      if(.not.geom_rtdb_store(rtdb,geom,'geometry'))
     &   call errquit(pname//'geom_rtdb_store',0,GEOM_ERR)
c
c     energy/gradient calculation
      esys = 0.d0
      if (.not. task_gradient(rtdb))
     &   call errquit(pname//'task_gradient',0,CALC_ERR)
      if(.not.rtdb_get(rtdb,'task:gradient',mt_dbl,3*nat,g))
     &   call errquit(pname//'rtdb_get gradient',0,RTDB_ERR)
      if(.not.rtdb_get(rtdb,'task:energy',mt_dbl,1,esys))
     &   call errquit(pname//'rtdb_get energy',0,RTDB_ERR)
c
c     acceleration (a) and velocities (v) at dt
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
c
c     apply the thermostat
      if (ekin.le.0.d0) return
c
      if (do_nosehoover) then
         call qmd_nosehoover(nat,dt,v,ekin)
      else if (do_berendsen) then
         call qmd_berendsen(nat,dt,v,ekin)
      else
         call errquit(pname//'unknown thermostat',0,RTDB_ERR)
      end if
c
c     instantaneous kinetic energy
      call qmd_kinetic(nat,m,v,ekin)
c
c     instantaneous temperature
      final_temp = 2.0d0*ekin/kb/(1.d0*ndof)
c
      return
      end
