*
* $Id: lcao_write_matrix.F,v 1.1 2001-11-30 23:52:05 bylaska Exp $
*
      subroutine lcao_write_matrix(unit,matrix)
      implicit none
      integer unit
      integer matrix

#include "mafdecls.fh"
#include "global.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)
      integer i,j,nbasis,ld
      real*8  sum
      integer  aorbs_nbasis
      external aorbs_nbasis

      call Parallel_taskid(taskid)
      
      nbasis = aorbs_nbasis()
      do i=1,nbasis

         do j=1,nbasis
         call ga_get(matrix,i,i,j,j,sum,ld)
         if (taskid.eq.MASTER) then
         write(unit,*) i,j,sum
         end if

         end do
      end do

      return
      end


