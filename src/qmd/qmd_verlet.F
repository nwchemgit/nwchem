c
c     velocity verlet md step
c
      subroutine qmd_verlet(rtdb,nat,dt,m,r,v,g,esys,ekin)
c
      implicit none
c
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
c
#include "qmd_params.fh"
c
      integer rtdb               ! rtdb
      integer nat                ! number of atoms
      double precision dt        ! step size
      double precision m(nat)    ! mass
      double precision r(3,nat)  ! coordinates
      double precision v(3,nat)  ! velocities
      double precision g(3,nat)  ! gradients
      double precision esys      ! system energy
      double precision ekin
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
      dtsq = dt*dt
      dtsq2 = 0.5d0*dtsq
      dt2 = 0.5d0*dt
c
c     instantaneous kinetic energy
      call qmd_kinetic(nat,m,v,ekin)
c
c     apply the nose-hoover thermostat
      if (do_nosehoover) 
     &   call qmd_nosehoover(nat,dt,v,ekin)
c
c     acceleration for the first step
      do i=1,nat
        do j=1,3
           a(j,i) = -g(j,i)/m(i)  ! a = F/m
        end do
      end do
c
c     update atom positions and velocities at dt/2
      do i = 1,nat
c       coordinates
        r(1,i) = r(1,i) + v(1,i)*dt + a(1,i)*dtsq2
        r(2,i) = r(2,i) + v(2,i)*dt + a(2,i)*dtsq2
        r(3,i) = r(3,i) + v(3,i)*dt + a(3,i)*dtsq2
c       velocities
        do j = 1,3
          v(j,i) = v(j,i) + a(j,i)*dt2
        end do
      end do
c
c     energy gradient calculation
      esys = 0.d0
      if (.not. task_gradient(rtdb))
     &   call errquit(pname//'task_gradient',0,CALC_ERR)
      if(.not.rtdb_get(rtdb,'task:gradient',mt_dbl,3*nat,g))
     &   call errquit(pname//'rtdb_get gradient',0,RTDB_ERR)
      if(.not.rtdb_get(rtdb,'task:energy',mt_dbl,1,esys))
     &   call errquit(pname//'rtdb_get energy',0,RTDB_ERR)
c
c     acceleration and velocities at dt
      do i=1,nat
        do j=1,3
           a(j,i) = -g(j,i)/m(i)  ! a = F/m
           v(j,i) = v(j,i) + a(j,i)*dt2
        end do
      end do
c
c     instantaneous kinetic energy
      call qmd_kinetic(nat,m,v,ekin)
c
c     apply the thermostat
      if (do_nosehoover) then
       call qmd_nosehoover(nat,dt,v,ekin)
      else if (do_berendsen) then
       call qmd_berendsen(nat,dt,v,ekin)
      else
       call errquit(pname//'unknown thermostat',0,RTDB_ERR)
      end if
c
c     instantaneous temperature
      final_temp = 2.0d0*ekin/kb/(1.d0*ndof)
c
      return
      end
