c
c     berendsen thermostat
c
      subroutine qmd_berendsen(nat,dt,v,ekin)
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
      integer i,j
      double precision inst_temp
      double precision ratio_temp
      double precision scalefac
c
      character*32 pname
      pname='qmd_berendsen: '
c
c     instantaneous temperature
      inst_temp = 2.0d0*ekin/kb/(1.d0*ndof)
c
c     velocity scaling factor: Berendsen scaling
      ratio_temp = targ_temp/inst_temp
      scalefac = dsqrt(1.d0 + (dt/tau)*(ratio_temp - 1.d0))
      do i = 1,nat
        do j = 1,3
          v(j,i) = scalefac*v(j,i)
        end do
      end do
c
      return
      end
