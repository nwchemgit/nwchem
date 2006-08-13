* $Id: Dnex.F,v 1.1 2006-08-13 01:03:24 bylaska Exp $
*


*     ***********************************
*     *					*
*     *	        Dnex_init		*	
*     *					*
*     ***********************************

      subroutine Dnex_init(ispin_in,ne_in,map_in)
      implicit none
      integer ispin_in
      integer ne_in(2)
      integer map_in

#include "Dnex.fh"

      ispin   = ispin_in
      ne(1)   = ne_in(1)
      ne(2)   = ne_in(2)

      neq(1) = 0
      neq(2) = 0
      call D1dB_init(3,ne(1),map_in)
      call D1dB_nq(3,neq(1))
      if (ispin.eq.2) then
        call D1dB_init(4,ne(2),map_in)
        call D1dB_nq(4,neq(2))
      end if

      return
      end


*     ***********************************
*     *					*
*     *	          Dnex_end   		*	
*     *					*
*     ***********************************

      subroutine Dnex_end()
      implicit none

#include "Dnex.fh"

      call D1dB_end(3)
      if (ispin.eq.2) call D1dB_end(4)

      return
      end
         

*     ***********************************
*     *					*
*     *	         Dnexall_ntoqp 		*	
*     *					*
*     ***********************************

      subroutine Dnexall_ntoqp(n,q,p)      
      implicit none
      integer n,q,p

#include "Dnex.fh"

      if (n.le.ne(1)) then
        call D1dB_ktoqp(3,n,q,p)
      else
        call D1dB_ktoqp(4,n,q,p)
        q = q + neq(1)
      end if

      return
      end



*     ***********************************
*     *					*
*     *	        Dnexall_neq		*	
*     *					*
*     ***********************************

      subroutine Dnexall_neq(nqtmp)
      implicit none
      integer nqtmp(2)

#include "Dnex.fh"

      nqtmp(1) = neq(1)
      nqtmp(2) = neq(2)
      return 
      end

