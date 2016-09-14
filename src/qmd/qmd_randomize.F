c
c     qmd_randomize: randomize velocities
c
      subroutine qmd_randomize(n,m,v)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "util_params.fh"
c
#include "qmd_common.fh"
c
      integer n                  ! number of atoms
      double precision m(n)    ! atomic masses
      double precision v(3,n)  ! atomic velocities
c
      integer i,j
c
      external qmd_gasdev
      double precision qmd_gasdev
c
      character*32 pname
      pname='qmd_randomize: '
c
      do i = 1,n
       do j = 1,3
         v(j,i) = dsqrt(targ_temp*kb_au/(m(i)*camu2au))*qmd_gasdev(idum)
       end do
      end do
c
      return
      end
