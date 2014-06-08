c     This function allows the user to request printing of a specific
c     item or the printing of all quantities relevant to the TDDFT
c     gradients. Its only purpose is convenience in programming.
c
c     The type definition of this function is stored in 
c     tddft_grad_util.fh which also includes util.fh to have everything
c     together in one place.
c
c     Written by Huub van Dam, Oct 2007.
c
      logical function tddft_grad_util_print(name,level)
c
      implicit none
c
#include "util.fh"
c
c     Input:
c
      character*(*) name  ! the item for which the print level is to
                          ! be checked.
      integer       level ! the print level that needs to be exceeded
                          ! for printing to be activated.
c
c     Code:
c
      tddft_grad_util_print = util_print('tddft grad',level).or.
     +                        util_print(name,level)
c
      return
      end
c $Id$
