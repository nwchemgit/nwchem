      Subroutine dftg_wderiv(ictr, nat, nq_orig, nq, 
     .     iandex, ncenters, func, 
     &                       d_qwght, force, oactive)
*
* $Id$
*
      Implicit none
c
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "tcgmsg.fh"
c
c     compute the product \delta w(j)*func(j)
c
      integer ictr
      integer nat, nq,nq_orig
      logical oactive(nat)
      integer iandex(*), ncenters
      double precision func(*)          ! func value [input]
      double precision d_qwght(3,nq_orig,*)  ! w deriv [input]
c     
c     force vector
c     
      double precision force(3,nat)
c
c     local
c
      integer iat, i
      double precision out(3)
c
      do i = 1, ncenters
         iat = iandex(i)
         if(iat.ne.ictr)then
            call dgemv('n', 3, nq, 1.d0, d_qwght(1,1,i), 3, 
     &           func, 1, 0d0, out, 1)
           if (oactive(iat)) then
              force(1,iat)=force(1,iat)+out(1)
              force(2,iat)=force(2,iat)+out(2)
              force(3,iat)=force(3,iat)+out(3)
           endif
           if (oactive(ictr)) then
              force(1,ictr)=force(1,ictr)-out(1)
              force(2,ictr)=force(2,ictr)-out(2)
              force(3,ictr)=force(3,ictr)-out(3)
           endif
         endif
      enddo
c
      return 
      end 
