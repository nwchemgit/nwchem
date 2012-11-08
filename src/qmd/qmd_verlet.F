c
c     velocity verlet md step
      subroutine qmd_verlet(rtdb,nat,dt,m,r,v,g,esys)
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
c
      logical task_gradient
      external task_gradient
c
c     preliminaries
      dt2 = dt*dt
c
c     calculate the acceleration for the first step
      do i=1,nat
        do j=1,3
           a(j,i) = -g(j,i)/m(i)  ! a = F/m
        end do
      end do
c
c     update atom positions and velocities at dt/2
      do i = 1,nat
c       coordinates
        r(1,i) = r(1,i) + v(1,i)*dt + a(1,i)*0.5d0*dt2
        r(2,i) = r(2,i) + v(2,i)*dt + a(2,i)*0.5d0*dt2
        r(3,i) = r(3,i) + v(3,i)*dt + a(3,i)*0.5d0*dt2
c       velocities
        do j = 1,3
          v(j,i) = v(j,i) + a(j,i)*0.5d0*dt
        end do
      end do
c
c     perform energy gradient calculation
      esys = 0.d0
      if (.not. task_gradient(rtdb))
     &   call errquit(pname//'task_gradient',0,CALC_ERR)
      if(.not.rtdb_get(rtdb,'task:gradient',mt_dbl,3*nat,g))
     &   call errquit(pname//'rtdb_get gradient',0,RTDB_ERR)
      if(.not.rtdb_get(rtdb,'task:energy',mt_dbl,1,esys))
     &   call errquit(pname//'rtdb_get energy',0,RTDB_ERR)
c
c     calculate acceleration and velocities at dt
      do i=1,nat
        do j=1,3
           a(j,i) = -g(j,i)/m(i)  ! a = F/m
           v(j,i) = v(j,i) + a(j,i)*0.5d0*dt
        end do
      end do
c
      return
      end
