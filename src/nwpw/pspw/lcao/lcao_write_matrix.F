*
* $Id: lcao_write_matrix.F,v 1.2 2001-12-02 19:38:14 bylaska Exp $
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
      logical value
      integer i,j,nbasis
      integer coef(2)
      integer  aorbs_nbasis
      external aorbs_nbasis

      call Parallel_taskid(taskid)
      
      nbasis = aorbs_nbasis()
      value = MA_push_get(mt_dbl,(nbasis),'coef',coef(2),coef(1))
      if (.not. value)
     >  call errquit('lcao_write_matrix:out of stack memory',0)


      do i=1,nbasis

         call ga_get(matrix,i,i,1,nbasis,dbl_mb(coef(1)),1)
         if (taskid.eq.MASTER) then
         write(unit,1000) (dbl_mb(coef(1)+j-1),j=1,nbasis)
         end if

      end do


      value = MA_pop_stack(coef(2))
      if (.not. value)
     >  call errquit('lcao_write_matrix:popping stack memory',0)

      return
 1000 format(100E10.3)
      end


