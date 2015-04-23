c
c     Update velocity and positions
c
      subroutine qmd_nucl_dt(nat,dt,m,r,v,g)
c
      implicit none
c
      integer nat                ! number of atoms
      double precision dt        ! step size
      double precision m(nat)    ! mass
      double precision r(3,nat)  ! coordinates
      double precision v(3,nat)  ! velocities
      double precision g(3,nat)  ! gradients
c
      integer i
      double precision a(3,nat)  ! acceleration
c
      character*30 pname
c
c     preliminaries
      pname='qmd_nucl_dt: '
c
c     calculate acceleration, update positions (r) and velocities (v) at dt/2
      do i = 1,nat
        a(1,i) = -g(1,i)/m(i)
        v(1,i) = v(1,i) + a(1,i)*0.5d0*dt
        r(1,i) = r(1,i) + v(1,i)*dt 
c
        a(2,i) = -g(2,i)/m(i)
        v(2,i) = v(2,i) + a(2,i)*0.5d0*dt
        r(2,i) = r(2,i) + v(2,i)*dt
c
        a(3,i) = -g(3,i)/m(i)
        v(3,i) = v(3,i) + a(3,i)*0.5d0*dt
        r(3,i) = r(3,i) + v(3,i)*dt
      end do
c
      return
      end
c $Id$
