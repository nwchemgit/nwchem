      Subroutine dftg_wderiv(ictr, nat, nq, iandex, ncenters, func, 
     &                       d_qwght, force, oactive,  nactive)
*
* $Id: dftg_wderiv.F,v 1.4 1998-03-24 22:59:09 d3h449 Exp $
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
      integer nat, nq
      Integer nactive
      logical oactive(nat)
      integer iandex(*), ncenters
      double precision func(*)          ! func value [input]
      double precision d_qwght(3,nq,*)  ! w deriv [input]
c     
c     force vector
c     
      double precision force(3,nat)
c
c     local
c
      integer iat, i
c
      do i = 1, ncenters
         iat = iandex(i)
         if(iat.ne.ictr)then
           if (oactive(iat)) 
     &          call dgemv('n', 3, nq, 1.d0, d_qwght(1,1,i), 3, 
     &          func, 1, 1.d0, force(1,iat), 1)
           if (oactive(ictr)) 
     &          call dgemv('n', 3, nq, -1.d0, d_qwght(1,1,i), 3, 
     &          func, 1, 1.d0, force(1,ictr), 1)
         endif
      enddo
c
      return 
      end 
