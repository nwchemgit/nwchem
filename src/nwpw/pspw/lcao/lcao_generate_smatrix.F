*
* $Id$
*
      subroutine lcao_generate_smatrix(smatrix,phi1,phi2)
      implicit none
      integer smatrix
      complex*16 phi1(*)
      complex*16 phi2(*)

#include "bafdecls.fh"
#include "global.fh"

*     **** local variables ****
      integer i,j,nbasis,npack1
      real*8  sum

*     ***** external functions ****
      integer  aorbs_nbasis
      external aorbs_nbasis

      
      call Pack_npack(1,npack1)
      nbasis = aorbs_nbasis()
      do i=1,nbasis

*        **** get phi1 ****
         call aorbs_aorb(i,phi1)
         do j=1,i

*           **** get phi2 ****
            call aorbs_aorb(j,phi2)
            call Pack_cc_dot(1,phi1,phi2,sum)

            call ga_fill_patch(smatrix,i,i,j,j,sum)
            call ga_fill_patch(smatrix,j,j,i,i,sum)
         end do
      end do

      return
      end


