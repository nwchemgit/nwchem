      Subroutine getNxyz(lambda,ic,Nxyz)
c $Id: getNxyz.F,v 1.1 1994-09-28 17:56:58 d3e129 Exp $
      implicit none

      integer lambda, ic, Nxyz(3)

#include "sh_order.fh"

      integer l, l3, nc

c
c Get the angular momentum exponents for the "ic" component of a shell of
c Cartesian GTFs of angular momentum "lambda".
c
c The order of components in shells of increasing angular momentum was defined
c by the routine "defNxyz" as follows:
c
c     Data s/ 0,0,0 /
c
c     Data p/ 1,0,0,
c    &        0,1,0,
c    &        0,0,1 /
c
c     Data d/ 2,0,0,
c    &        1,1,0,
c    &        1,0,1,
c    &        0,2,0,
c    &        0,1,1,
c    &        0,0,2 /
c
c     Data f/ 3,0,0,
c    &        2,1,0,
c    &        2,0,1,
c    &        1,2,0,
c    &        1,1,1,
c    &        1,0,2,
c    &        0,3,0,
c    &        0,2,1,
c    &        0,1,2,
c    &        0,0,3 /
c
c     etc.
c
c******************************************************************************

c Define the angular momentum indices for the given shell component.

      if( lambda.le.linit )then

       l = lambda - 1
       l3 = ((l+1)*(l+2)*(l+3))/6

       nc = l3 + ic

       Nxyz(1) = Ixyz(1,nc)
       Nxyz(2) = Ixyz(2,nc)
       Nxyz(3) = Ixyz(3,nc)

      else

       write(*,*)
       write(*,*) 'GETNXYZ:  Shell component ordering was defined for'
       write(*,*) '          shells up through L = ',linit,'.'
       write(*,*) '          Initialization (defNxyz) should include '
       write(*,*) '          shells up through L = ',lambda,'.'
       stop

      end if

      end
