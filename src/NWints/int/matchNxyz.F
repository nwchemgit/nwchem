      Subroutine matchNxyz(lambda,ic,Nxyz)
c $Id: matchNxyz.F,v 1.1 1994-10-29 16:41:23 d3g681 Exp $
      implicit none

      integer lambda, ic, Nxyz(3)

#include "sh_order.fh"

      integer l, l3, nc

c
c Given the angular momentum exponents for the "ic" component of a shell of
c Cartesian GTFs of angular momentum "lambda" return ic. Hard error if no
c match is found.
c
c 

      if( lambda.le.linit )then

       l = lambda - 1
       l3 = ((l+1)*(l+2)*(l+3))/6

       do ic = 1, ((l+2)*(l+3))/2
          nc = l3 + ic
          if (
     $         Nxyz(1) .eq. Ixyz(1,nc) .and.
     $         Nxyz(2) .eq. Ixyz(2,nc) .and.
     $         Nxyz(3) .eq. Ixyz(3,nc))
     $         return
       enddo
       call errquit('matchNxyz: invalid shell indices', 
     $      1000*lambda+100*nxyz(1)+10*nxyz(2)+nxyz(3))
      else
       write(*,*)
       write(*,*) 'MATCHNXYZ:  Shell component ordering was defined for'
       write(*,*) '            shells up through L = ',linit,'.'
       write(*,*) '            Initialization (defNxyz) should include '
       write(*,*) '            shells up through L = ',lambda,'.'
       stop

      end if

      end
