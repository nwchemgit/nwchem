*
* $Id$
*
      subroutine lcao_write_matrix_seq(unit,size,matrix)
      implicit none
      integer unit,size
      real*8 matrix(size,size)

      integer i,j
      do i=1,size
         write(unit,1000) (matrix(i,j),j=1,size)
      end do
      return 
 1000 format(100E10.3)
      end

      subroutine lcao_write_matrix(unit,matrix)
      implicit none
#include "errquit.fh"
      integer unit
      integer matrix

#include "bafdecls.fh"
#include "global.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)
      logical value
      integer i,j,nbasis
      integer coef(2)
      integer  aorbs_nbasis
      external aorbs_nbasis

      call Parallel_taskid(taskid)
      
      nbasis = aorbs_nbasis()
      value = BA_push_get(mt_dbl,(nbasis),'coef',coef(2),coef(1))
      if (.not. value)
     >  call errquit('lcao_write_matrix:out of stack memory',0, MA_ERR)


      do i=1,nbasis

         call ga_get(matrix,i,i,1,nbasis,dbl_mb(coef(1)),1)
         if (taskid.eq.MASTER) then
         write(unit,1000) (dbl_mb(coef(1)+j-1),j=1,nbasis)
         end if

      end do


      value = BA_pop_stack(coef(2))
      if (.not. value)
     >  call errquit('lcao_write_matrix:popping stack memory',0, MA_ERR)

      return
 1000 format(100E10.3)
      end


