c
c     qmd_nucl_dt: velocity at full dt
c
      subroutine qmd_nucl_dt(nat,dt,m,v,g)
c
      implicit none
c
      integer nat                ! number of atoms
      double precision dt        ! step size
      double precision m(nat)    ! mass
      double precision v(3,nat)  ! velocities
      double precision g(3,nat)  ! gradients
c
      double precision a(3,nat)  ! acceleration
      integer i,j
      double precision dt2
c
      character*30 pname
c
c     preliminaries
      pname = "qmd_nucl_dt: "
      dt2 = 0.5d0*dt
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
      return
      end
