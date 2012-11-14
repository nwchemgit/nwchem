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
#include "qmd_params.fh"
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
c     nose-hoover chain
      gnh(2) = (mnh(1)*vnh(1)*vnh(1)-kb*targ_temp)/mnh(2)
      vnh(2) = vnh(2) + gnh(2)*dt4
      vnh(1) = vnh(1) * dexp(-vnh(2)*dt8)
      gnh(1) = (2.0d0*ekin-ndof*1.d0*kb*targ_temp)/mnh(1)
      vnh(1) = vnh(1) + gnh(1)*dt4
      vnh(1) = vnh(1) * dexp(-vnh(2)*dt8)
      rnh(1) = rnh(1) + vnh(1)*dt2
      rnh(2) = rnh(2) + vnh(2)*dt2
      scalefac = dexp(-vnh(1)*dt2)
c
c     scale kinetic energy
      ekin = ekin*scalefac*scalefac
c
c     nose-hoover chain
      vnh(1) = vnh(1) * dexp(-vnh(2)*dt8)
      gnh(1) = (2.0d0*ekin-ndof*1.d0*kb*targ_temp)/mnh(1)
      vnh(1) = vnh(1) + gnh(1)*dt4
      vnh(1) = vnh(1) * dexp(-vnh(2)*dt8)
      gnh(2) = (mnh(1)*vnh(1)*vnh(1)-kb*targ_temp)/mnh(2)
      vnh(2) = vnh(2) + gnh(2)*dt4
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
