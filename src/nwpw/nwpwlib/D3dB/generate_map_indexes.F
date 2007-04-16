
*
* $Id: generate_map_indexes.F,v 1.7 2007-04-16 23:28:25 bylaska Exp $
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
           if (p.ge.rmdr1) nq2 = nq1
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

*     ***********************************
*     *                                 *
*     *        expand_hilbert2d         *
*     *                                 *
*     ***********************************

      subroutine expand_hilbert2d(np,ny,nz, 
     >                            pmap_in, qmap_in,
     >                            pmap_out,qmap_out)
      implicit none
      integer np,ny,nz
      integer pmap_in(ny,nz)
      integer qmap_in(ny,nz)
      integer pmap_out(2*ny,2*nz) 
      integer qmap_out(2*ny,2*nz)

      !**** local variables ****
      integer j,k,p,nqp


      !*** set pmap_out ***
      do k=1,nz
      do j=1,ny
         pmap_out(j,k)       = pmap_in(j,k)
         pmap_out(j+ny,k)    = pmap_in(j,k)
         pmap_out(j,k+nz)    = pmap_in(j,k)
         pmap_out(j+ny,k+nz) = pmap_in(j,k)
      end do
      end do


      !*** set qmap_out ***
      do p=0,np-1

         nqp=0
         do k=1,nz
         do j=1,ny
           if (pmap_in(j,k).eq.p) then
              if (qmap_in(j,k).gt.nqp) nqp = qmap_in(j,k)
           end if
         end do
         end do

         do k=1,nz
         do j=1,ny
           if (pmap_in(j,k).eq.p) then
              qmap_out(j,k)       = qmap_in(j,k)
              qmap_out(j+ny,k)    = qmap_in(j,k) + nqp
              qmap_out(j,k+nz)    = qmap_in(j,k) + 2*nqp
              qmap_out(j+ny,k+nz) = qmap_in(j,k) + 3*nqp
           end if
         end do
         end do

      end do

      return
      end
         



*     ***********************************
*     *                                 *
*     *        generate_map_indexes_2   *
*     *                                 *
*     ***********************************

      subroutine generate_map_indexes_2(taskid,np,ny,nz,p_map,q_map,nq)
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


      nq = 0
      p=0
      q=1
      do i=1,(ny*nz)
        int_mb(indx_proc(1)+i-1) = p
        int_mb(indx_q(1)+i-1)    = q
        if (taskid.eq.p) nq = nq + 1

        p = p + 1
        if (p.ge.np) then
          p = 0
          q = q+1
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


*     ***********************************
*     *                                 *
*     *        generate_map_indexes_3   *
*     *                                 *
*     ***********************************

      subroutine generate_map_indexes_3(taskid,np,ny,nz,p_map,q_map,nq)
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
      integer i,j,k,q,p,jj,kk,p1,q1,p2,q2
      integer nq1,nq2,rmdr1


      nq = 0
      p=0
      q=1
      do k=1,nz
      do j=1,ny
        p_map(j+(k-1)*ny) = p
        q_map(j+(k-1)*ny) = q
        if (taskid.eq.p) nq = nq + 1

        p = p + 1
        if (p.ge.np) then
          p = 0
          q = q+1
        end if
      end do
      end do

      do k=0,(nz/2-1)
        i = k+2
        j = nz-k
        p1 = p_map(1+(i-1)*ny)
        q1 = q_map(1+(i-1)*ny)
        p2 = p_map(1+(j-1)*ny)
        q2 = q_map(1+(j-1)*ny)

        value = .false.
        jj=2
        kk=1
        do while (.not.value)
          p = p_map(jj+(kk-1)*ny)
          q = q_map(jj+(kk-1)*ny)

          value = (p.eq.p1)
          if (value) then
             p_map(jj+(kk-1)*ny) = p2
             q_map(jj+(kk-1)*ny) = q2
             p_map(1+(j-1)*ny)   = p
             q_map(1+(j-1)*ny)   = q
          else
             kk = kk+1
             if (kk.gt.nz) then
               kk=1
               jj=jj+1
             end if
          end if
        end do

      end do

      return
      end


