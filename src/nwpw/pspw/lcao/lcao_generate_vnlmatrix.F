*
* $Id$
*
      subroutine lcao_generate_vnlmatrix(vnlmatrix,phi1,phi2)
      implicit none
      integer vnlmatrix
      complex*16 phi1(*)
      complex*16 phi2(*)

#include "bafdecls.fh"
#include "global.fh"

*     **** local variables ****
      logical move,frac
      integer i,j,nbasis,npack1
      integer ne(2),ispin
      real*8  sum,fion(3),occ(1)

*     ***** external functions ****
      integer  aorbs_nbasis
      external aorbs_nbasis

      
      call Pack_npack(1,npack1)
      nbasis = aorbs_nbasis()
      ispin = 1
      ne(1) = 1
      ne(2) = 0
      move =.false.
      frac = .false.
      do j=1,nbasis

*        **** get phi2 ****
         call aorbs_aorb(j,phi1)
         call Pack_c_Zero(1,phi2)
         call v_nonlocal(ispin,ne,phi1,phi2,move,fion,frac,occ)

         do i=1,j

*           **** get phi1 ****
            call aorbs_aorb(i,phi1)
            call Pack_cc_dot(1,phi1,phi2,sum)
            sum = -sum
            call ga_fill_patch(vnlmatrix,i,i,j,j,sum)
            call ga_fill_patch(vnlmatrix,j,j,i,i,sum)
         end do
      end do

      return
      end


