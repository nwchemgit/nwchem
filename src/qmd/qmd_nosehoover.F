c
c     nose-hoover thermostat
c
      subroutine qmd_nosehoover(nat,dt,v,ekin)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
c
#include "qmd_common.fh"
c
      integer nat
      double precision dt
      double precision v(3,nat)
      double precision ekin
c
      double precision dt2,dt4,dt8
c
      integer i,j
      double precision scalefac
c
      character*30 pname
c
c     preliminaries   
      pname='qmd_nosehoover: '
      dt2 = 0.5d0*dt
      dt4 = 0.25d0*dt
      dt8 = 0.125d0*dt
c
c     nose-hoover chain: 2 chains
      g_nh(2) = (m_nh(1)*v_nh(1)*v_nh(1)-kb*targ_temp)/m_nh(2)
      v_nh(2) = v_nh(2) + g_nh(2)*dt4
      v_nh(1) = v_nh(1) * dexp(-v_nh(2)*dt8)
c
      g_nh(1) = (2.0d0*ekin - ndof*kb*targ_temp)/m_nh(1)
      v_nh(1) = v_nh(1) + g_nh(1)*dt4
      v_nh(1) = v_nh(1) * dexp(-v_nh(2)*dt8)
c
      r_nh(1) = r_nh(1) + v_nh(1)*dt2
      r_nh(2) = r_nh(2) + v_nh(2)*dt2
c
      scalefac = dexp(-v_nh(1)*dt2)
c
c     scale kinetic energy
      ekin = ekin*scalefac*scalefac
c
c     nose-hoover chain: 2 chains
      v_nh(1) = v_nh(1) * dexp(-v_nh(2)*dt8)
      g_nh(1) = (2.0d0*ekin -ndof*kb*targ_temp)/m_nh(1)
      v_nh(1) = v_nh(1) + g_nh(1)*dt4
      v_nh(1) = v_nh(1) * dexp(-v_nh(2)*dt8)
      g_nh(2) = (m_nh(1)*v_nh(1)*v_nh(1)-kb*targ_temp)/m_nh(2)
      v_nh(2) = v_nh(2) + g_nh(2)*dt4
c
c     scale velocity
      do i = 1,nat
       do j = 1,3
         v(j,i) = scalefac*v(j,i)
       end do
      end do
c
      return
      end
