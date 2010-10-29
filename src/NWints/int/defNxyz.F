      Subroutine defNxyz(lambda)
c $Id$
      implicit none

      integer lambda

#include "sh_order.fh"

      integer l0, l1, l2, lr, nc
      external bdIxyz
c
c Initialize the "Ixyz" data structure up through 4 times the given angular 
c momentum. This is necessary to calculate a 4-center ERI comprising four
c shells of this maximum angular momentum.
c
c The order of components in shells of increasing angular momentum is defined
c as follows:
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

c Multiply by a factor of 4 to intialize for 4-ctr ERIs.

      if (lambda.eq.0) then
        lR = 4    ! minimum specification for s (s functions multi-pole=4)
      else
        lR = 4*lambda
      endif

c Determine whether initialization up through the specified angular momentum
c has already occurred.

      if( lR.le.linit ) return

      if( lR.gt.lRmax )then
       write(*,*)
       write(*,*) 'DEFNXYZ:  ''lRmax'' must be increased to define the'
       write(*,*) '          shell component ordering of higher angular'
       write(*,*) '          momentum basis functions.'
       write(*,*)
       write(*,*) '          For integral evaluation only:'
       write(*,*)
       write(*,*) '          lRmax = 4*(highest ang. mom. in basis) '
       write(*,*)
       stop
      end if

      nc = 0

      do 20 l0 = 0,lR

       do 10 l1 = l0,0,-1
       do 10 l2 = l0-l1,0,-1

        nc = nc + 1

        Ixyz(1,nc) = l1
        Ixyz(2,nc) = l2
        Ixyz(3,nc) = l0 - (l1 + l2)

   10  continue

   20 continue

      linit = lR

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      Block data bdIxyz
      implicit none

#include "sh_order.fh"

      Data Ixyz, linit / lRmax3*0, lRmax3*0, lRmax3*0, -375 /

      end
