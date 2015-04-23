c
c     Update positions
c
      subroutine qmd_nucl_pos(nat,dt,m,r,v)
c
      implicit none
c
      integer nat                ! number of atoms
      double precision dt        ! step size
      double precision m(nat)    ! mass
      double precision r(3,nat)  ! coordinates
      double precision v(3,nat)  ! velocities
c
      integer i
c
      character*30 pname
c
c     preliminaries
      pname='qmd_nucl_pos: '
c
c     Update positions (r)
      do i = 1,nat
        r(1,i) = r(1,i) + v(1,i)*dt 
        r(2,i) = r(2,i) + v(2,i)*dt
        r(3,i) = r(3,i) + v(3,i)*dt
      end do
c
      return
      end
c $Id: qmd_nucl_dt2.F 24348 2013-06-23 04:29:20Z d3y133 $
