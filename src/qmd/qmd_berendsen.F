c
c     berendsen thermostat
c
      subroutine qmd_berendsen(dt,tau,inst_temp,targ_temp,nat,v)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
c
      double precision dt
      double precision tau
      double precision inst_temp
      double precision targ_temp
      integer nat
      double precision v(3,nat)
c
      integer i,j
      double precision ratio_temp
      double precision velscale
c
      character*30 pname
      pname='qmd_berendsen: '
c
c     calculate velocity scaling factor: Berendsen scaling
      ratio_temp = targ_temp/inst_temp
      velscale = dsqrt(1.d0 + (dt/tau)*(ratio_temp - 1.d0))
      do i = 1,nat
        do j = 1,3
          v(j,i) = velscale*v(j,i)
        end do
      end do
c
      return
      end
