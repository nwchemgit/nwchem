
*
* $Id: generate_map_indexes.F,v 1.3 2004-04-24 21:07:11 bylaska Exp $
*

*     ***********************************
*     *                                 *
*     *        generate_map_indexes     *
*     *                                 *
*     ***********************************

      subroutine generate_map_indexes(taskid,np,ny,nz,p_map,q_map,nq)
      implicit none
      integer taskid,np
      integer ny,nz
      integer p_map(*)
      integer q_map(*)
      integer nq

#include "mafdecls.fh"
#include "errquit.fh"


      !**** local variables ****
      logical value
      integer i,j,k,q,p
      integer nq1,nq2,rmdr1
      integer indx_proc(2),indx_q(2)

      !**** allocate temporary space ****
      value =           MA_push_get(mt_int,(ny*nz),'indx_proc',
     >                    indx_proc(2),
     >                    indx_proc(1))
      value = value.and.MA_push_get(mt_int,(ny*nz),'indx_q',
     >                    indx_q(2),
     >                    indx_q(1))
      if (.not. value) call errquit('out of stack memory',0,MA_ERR)


      nq1   = (ny*nz)/np
      rmdr1 = mod(ny*nz,np)
      nq2 = nq1
      if (rmdr1.gt.0) nq2 = nq2+1
      nq = 0
      p=0
      q=1
      do i=1,(ny*nz)
        int_mb(indx_proc(1)+i-1) = p
        int_mb(indx_q(1)+i-1)    = q

        if (taskid.eq.p) nq = nq + 1
        q = q+1
        if (q.gt.nq2) then
           q = 1
           p = p + 1
           p = mod(p,np)
           if (p.gt.rmdr1) nq2 = nq1
        end if
      end do


      do k=1,nz
      do j=1,ny
        i = p_map(j+(k-1)*ny)+1
        p = int_mb(indx_proc(1)+i-1)
        q = int_mb(indx_q(1)+i-1)
        p_map(j+(k-1)*ny) = p
        q_map(j+(k-1)*ny) = q
      end do
      end do

      !**** deallocate temporary space  ****
      value =           MA_pop_stack(indx_q(2))
      value = value.and.MA_pop_stack(indx_proc(2))
      if (.not. value) call errquit('error popping stack',0,MA_ERR)

      return
      end
