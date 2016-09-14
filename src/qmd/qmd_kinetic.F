c
c     calculate kinetic energy 
c
      subroutine qmd_kinetic(n,m,v,ekin)
c
      implicit none
c
#include "util_params.fh"
#include "qmd_common.fh"
c
      integer n
      double precision m(n)
      double precision v(3,n)
      double precision ekin
c
      integer i
c
c     total kinetic energy
      ekin = 0.d0
      do i = 1,n
       ekin = ekin +  m(i)*(v(1,i)*v(1,i)+v(2,i)*v(2,i)+v(3,i)*v(3,i))
      end do
      ekin = 0.5d0*ekin*camu2au
c
      return
      end
c $Id$
