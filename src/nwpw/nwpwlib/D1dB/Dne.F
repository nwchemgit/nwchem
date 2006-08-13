* $Id: Dne.F,v 1.1 2006-08-13 01:03:24 bylaska Exp $
*


*     ***********************************
*     *					*
*     *	        Dne_init		*	
*     *					*
*     ***********************************

      subroutine Dne_init(ispin_in,ne_in,map_in)
      implicit none
      integer ispin_in
      integer ne_in(2)
      integer map_in

#include "Dne.fh"

      ispin   = ispin_in
      ne(1)   = ne_in(1)
      ne(2)   = ne_in(2)

      neq(1) = 0
      neq(2) = 0
      call D1dB_init(1,ne(1),map_in)
      call D1dB_nq(1,neq(1))
      if (ispin.eq.2) then
        call D1dB_init(2,ne(2),map_in)
        call D1dB_nq(2,neq(2))
      end if

      return
      end


*     ***********************************
*     *					*
*     *	          Dne_end   		*	
*     *					*
*     ***********************************

      subroutine Dne_end()
      implicit none

#include "Dne.fh"

      call D1dB_end(1)
      if (ispin.eq.2) call D1dB_end(2)

      return
      end
         

*     ***********************************
*     *					*
*     *	          Dneall_ntoqp 		*	
*     *					*
*     ***********************************

      subroutine Dneall_ntoqp(n,q,p)      
      implicit none
      integer n,q,p

#include "Dne.fh"

      if (n.le.ne(1)) then
        call D1dB_ktoqp(1,n,q,p)
      else
        call D1dB_ktoqp(2,n,q,p)
        q = q + neq(1)
      end if

      return
      end



*     ***********************************
*     *					*
*     *	        Dneall_neq		*	
*     *					*
*     ***********************************

      subroutine Dneall_neq(nqtmp)
      implicit none
      integer nqtmp(2)

#include "Dne.fh"

      nqtmp(1) = neq(1)
      nqtmp(2) = neq(2)
      return 
      end

