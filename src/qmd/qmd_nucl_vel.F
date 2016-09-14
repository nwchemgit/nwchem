c
c     Update velocity
c
      subroutine qmd_nucl_vel(nat,dt,m,v,g)
c
      implicit none
c
#include "util_params.fh"
#include "qmd_common.fh"
c
      integer nat                ! number of atoms
      double precision dt        ! step size
      double precision m(nat)    ! mass
      double precision v(3,nat)  ! velocities
      double precision g(3,nat)  ! gradients
c
      integer i
      double precision a(3,nat)  ! acceleration
c
      character*30 pname
c
c     preliminaries
      pname='qmd_nucl_vel: '
c
c     Update velocity
      do i = 1,nat
        a(1,i) = -g(1,i)/(m(i)*camu2au)
        v(1,i) = v(1,i) + a(1,i)*dt
c
        a(2,i) = -g(2,i)/(m(i)*camu2au)
        v(2,i) = v(2,i) + a(2,i)*dt
c
        a(3,i) = -g(3,i)/(m(i)*camu2au)
        v(3,i) = v(3,i) + a(3,i)*dt
      end do
c
      return
      end
c $Id: qmd_nucl_dt2.F 24348 2013-06-23 04:29:20Z d3y133 $
