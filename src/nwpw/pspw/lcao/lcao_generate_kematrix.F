*
* $Id$
*
      subroutine lcao_generate_kematrix(kematrix,phi1,phi2)
      implicit none
      integer kematrix
      complex*16 phi1(*)
      complex*16 phi2(*)

#include "bafdecls.fh"
#include "global.fh"

*     **** local variables ****
      integer i,j,nbasis,npack1
      integer ispin,ne(2)
      real*8  sum

*     ***** external functions ****
      integer  aorbs_nbasis
      external aorbs_nbasis

      
      ispin = 1
      ne(1) = 1
      ne(2) = 0
      call Pack_npack(1,npack1)
      nbasis = aorbs_nbasis()
      do j=1,nbasis

*        **** get phi2 ****
         call aorbs_aorb(j,phi1)
         call ke(ispin,ne,phi1,phi2)

         do i=1,j

*           **** get phi1 ****
            call aorbs_aorb(i,phi1)
            call Pack_cc_dot(1,phi1,phi2,sum)

            sum = -sum
            call ga_fill_patch(kematrix,i,i,j,j,sum)
            call ga_fill_patch(kematrix,j,j,i,i,sum)
         end do
      end do

      return
      end


