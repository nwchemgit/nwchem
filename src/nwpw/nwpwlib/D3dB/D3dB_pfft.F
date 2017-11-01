
#define NBLOCKS 3

*
* $Id$
*

*     ***********************************************************
*     *								*
*     *   		   D3dB_pfft library	 		*
*     *			(MPI implemenation)			*
*     *								*
*     *   Author - Eric Bylaska					*
*     *   date   - 3/23/96					*
*     *								*
*     ***********************************************************

*     ***********************************
*     *					*
*     *	        D3dB_pfft_init		*
*     *					*
*     ***********************************

      subroutine D3dB_pfft_init()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "D3dB.fh"
#include "D3dB_pfft.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      double precision eps
      parameter (eps=1d-12)

      logical value,yrow,zrow,yzslab,has_small
      integer nxh,nyh,nzh,q,p
      integer k1,k2,k3,ii,nb
      integer i,j,k,id,nbmax
      integer index,index2
      double precision ggcut,g1,g2,g3,gg,unitg(3,3)
      integer zero_arow2(2),zero_arow3(2)

*     **** external functions ***
      integer  control_pfft3_qsize
      real*8   lattice_ggcut,lattice_wggcut,lattice_unitg
      real*8   lattice_unitg_small
      external control_pfft3_qsize
      external lattice_ggcut,lattice_wggcut,lattice_unitg
      external lattice_unitg_small
      logical  lattice_has_small
      external lattice_has_small



      call Parallel2d_taskid_i(taskid)
      nxh = nx(1)/2
      nyh = ny(1)/2
      nzh = nz(1)/2

      has_small = lattice_has_small()

      if (mapping.eq.1) then
         nxh = nx(1)/2
         nyh = ny(1)/2
         nzh = nz(1)/2
         value =           BA_alloc_get(mt_log,(nxh+1)*nq(1),
     >                     'zero_row3',zero_row3(2,0),zero_row3(1,0))
         value = value.and.BA_alloc_get(mt_log,(nxh+1)*nq(1),
     >                     'zero_row3',zero_row3(2,1),zero_row3(1,1))

         value = value.and.BA_alloc_get(mt_log,(nxh+1)*nq(1),
     >                     'zero_row2',zero_row2(2,0),zero_row2(1,0))
         value = value.and.BA_alloc_get(mt_log,(nxh+1)*nq(1),
     >                     'zero_row2',zero_row2(2,1),zero_row2(1,1))

         value = value.and.BA_alloc_get(mt_log,(nxh+1),
     >               'zero_slab23',zero_slab23(2,0),zero_slab23(1,0))
         value = value.and.BA_alloc_get(mt_log,(nxh+1),
     >               'zero_slab23',zero_slab23(2,1),zero_slab23(1,1))
         if (.not. value) call errquit('out of heap memory',0, MA_ERR)


         value = value.and.BA_push_get(mt_log,(nxh+1)*ny(1),
     >                     'zero_arow3',zero_arow3(2),zero_arow3(1))
         if (.not. value) call errquit('out of stack memory',0,MA_ERR)



         if (has_small) then
         nxh = nx(3)/2
         nyh = ny(3)/2
         nzh = nz(3)/2
         value = value.and.BA_alloc_get(mt_log,(nxh+1)*nq(3),
     >                     'zero_row3',zero_row3(2,2),zero_row3(1,2))
         value = value.and.BA_alloc_get(mt_log,(nxh+1)*nq(3),
     >                     'zero_row3',zero_row3(2,3),zero_row3(1,3))

         value = value.and.BA_alloc_get(mt_log,(nxh+1)*nq(3),
     >                     'zero_row2',zero_row2(2,2),zero_row2(1,2))
         value = value.and.BA_alloc_get(mt_log,(nxh+1)*nq(3),
     >                     'zero_row2',zero_row2(2,3),zero_row2(1,3))

         value = value.and.BA_alloc_get(mt_log,(nxh+1),
     >               'zero_slab23',zero_slab23(2,2),zero_slab23(1,2))
         value = value.and.BA_alloc_get(mt_log,(nxh+1),
     >               'zero_slab23',zero_slab23(2,3),zero_slab23(1,3))
         if (.not. value) call errquit('out of heap memory',0, MA_ERR)
         end if

         if (has_small) then
            nbmax = 3
         else
            nbmax = 1
         end if

         do nb=0,nbmax

         if (nb.lt.2) then
            id = 1
            nxh = nx(1)/2
            nyh = ny(1)/2
            nzh = nz(1)/2
            unitg(1,1) = lattice_unitg(1,1)
            unitg(2,1) = lattice_unitg(2,1)
            unitg(3,1) = lattice_unitg(3,1)
            unitg(1,2) = lattice_unitg(1,2)
            unitg(2,2) = lattice_unitg(2,2)
            unitg(3,2) = lattice_unitg(3,2)
            unitg(1,3) = lattice_unitg(1,3)
            unitg(2,3) = lattice_unitg(2,3)
            unitg(3,3) = lattice_unitg(3,3)
         else
            id = 3
            nxh = nx(3)/2
            nyh = ny(3)/2
            nzh = nz(3)/2
            unitg(1,1) = lattice_unitg_small(1,1)
            unitg(2,1) = lattice_unitg_small(2,1)
            unitg(3,1) = lattice_unitg_small(3,1)
            unitg(1,2) = lattice_unitg_small(1,2)
            unitg(2,2) = lattice_unitg_small(2,2)
            unitg(3,2) = lattice_unitg_small(3,2)
            unitg(1,3) = lattice_unitg_small(1,3)
            unitg(2,3) = lattice_unitg_small(2,3)
            unitg(3,3) = lattice_unitg_small(3,3)
         end if

         if ((nb.eq.0).or.(nb.eq.2)) then
            ggcut = lattice_ggcut()
         else
            ggcut = lattice_wggcut()
         end if

*        **** find zero_row3 - (i,j,*) rows that are zero ****
         do q=1,(nxh+1)*nq(id)
             log_mb(zero_row3(1,nb)+q-1) =.true.
         end do
         do q=1,(nxh+1)*ny(id)
           log_mb(zero_arow3(1)+q-1) =.true.
         end do


         do k2 = -nyh+1, nyh-1
         do k1 = 0,nxh-1
          i=k1
          j=k2
          if (i .lt. 0) i = i + nx(id)
          if (j .lt. 0) j = j + ny(id)
          i=i+1
          j=j+1
          zrow = .true.
          do k3 = -nzh+1, nzh-1
            g1 = k1*unitg(1,1)
     >         + k3*unitg(1,2)
     >         + k2*unitg(1,3)
            g2 = k1*unitg(2,1)
     >         + k3*unitg(2,2)
     >         + k2*unitg(2,3)
            g3 = k1*unitg(3,1)
     >         + k3*unitg(3,2)
     >         + k2*unitg(3,3)
            gg = g1*g1 + g2*g2 + g3*g3
            gg= gg-ggcut
            if (gg.lt.-eps) zrow = .false.
          end do
          if (.not.zrow) then
            log_mb(zero_arow3(1)+i-1+(nxh+1)*(j-1)) =.false.
            call D3dB_ktoqp(id,j,q,p)
            if (p.eq.taskid) then
              index2 = i + (nxh+1)*(q-1)
              log_mb(zero_row3(1,nb)+index2-1) =.false.
            end if
          end if

         end do
         end do

         call D3dB_c_ptranspose_jk_init(nb,log_mb(zero_arow3(1)))


*        **** find zero_slab23 - (i,*,*) slabs that are zero ****
         do k1 = 1,nxh+1
           log_mb(zero_slab23(1,nb)+k1-1) =.true.
         end do

         do k1 = 0,nxh-1
          i=k1
          if (i .lt. 0) i = i + nx(id)
          i=i+1
          yzslab = .true.
          do k3 = -nzh+1, nzh-1
          do k2 = -nyh+1, nyh-1
            g1 = k1*unitg(1,1)
     >         + k2*unitg(1,2)
     >         + k3*unitg(1,3)
            g2 = k1*unitg(2,1)
     >         + k2*unitg(2,2)
     >         + k3*unitg(2,3)
            g3 = k1*unitg(3,1)
     >         + k2*unitg(3,2)
     >         + k3*unitg(3,3)
            gg = g1*g1 + g2*g2 + g3*g3
            gg= gg-ggcut
            if (gg.lt.-eps) yzslab = .false.
          end do
          end do
          if (.not.yzslab) then
            log_mb(zero_slab23(1,nb)+i-1) =.false.
          end if

         end do

*        **** find zero_row2 - (i,*,k) rows that are zero after fft of (i,j,*) ****
         do k3 = 1,nz(id)
         do k1 = 1,nxh+1
           call D3dB_ktoqp(id,k3,q,p)
           if (p.eq.taskid) then
            index2 = k1 + (nxh+1)*(q-1)
            log_mb(zero_row2(1,nb)+index2-1) 
     >       = log_mb(zero_slab23(1,nb)+k1-1)
           end if
         end do
         end do

        
         end do !*nb*
         value = BA_pop_stack(zero_arow3(2))
         if (.not. value) call errquit('error freeing stack',0,MA_ERR)

      !*** mapping == 2 ***
      else
         nxh = nx(1)/2
         nyh = ny(1)/2
         nzh = nz(1)/2
         value =           BA_alloc_get(mt_log,nq3(1),
     >                     'zero_row3',zero_row3(2,0),zero_row3(1,0))
         value = value.and.BA_alloc_get(mt_log,nq3(1),
     >                     'zero_row3',zero_row3(2,1),zero_row3(1,1))

         value = value.and.BA_alloc_get(mt_log,nq2(1),
     >                     'zero_row2',zero_row2(2,0),zero_row2(1,0))
         value = value.and.BA_alloc_get(mt_log,nq2(1),
     >                     'zero_row2',zero_row2(2,1),zero_row2(1,1))

         value = value.and.BA_alloc_get(mt_log,(nxh+1),
     >               'zero_slab23',zero_slab23(2,0),zero_slab23(1,0))
         value = value.and.BA_alloc_get(mt_log,(nxh+1),
     >               'zero_slab23',zero_slab23(2,1),zero_slab23(1,1))
         if (.not. value) call errquit('out of heap memory',0, MA_ERR)

         value = value.and.BA_push_get(mt_log,(nxh+1)*nz(1),
     >                     'zero_arow2',zero_arow2(2),zero_arow2(1))
         value = value.and.BA_push_get(mt_log,(nxh+1)*ny(1),
     >                     'zero_arow3',zero_arow3(2),zero_arow3(1))
         if (.not. value) call errquit('out of stack memory',0,MA_ERR)


         if (has_small) then
         nxh = nx(3)/2
         nyh = ny(3)/2
         nzh = nz(3)/2
         value =           BA_alloc_get(mt_log,nq3(3),
     >                     'zero_row3',zero_row3(2,2),zero_row3(1,2))
         value = value.and.BA_alloc_get(mt_log,nq3(3),
     >                     'zero_row3',zero_row3(2,3),zero_row3(1,3))

         value = value.and.BA_alloc_get(mt_log,nq2(3),
     >                     'zero_row2',zero_row2(2,2),zero_row2(1,2))
         value = value.and.BA_alloc_get(mt_log,nq2(3),
     >                     'zero_row2',zero_row2(2,3),zero_row2(1,3))

         value = value.and.BA_alloc_get(mt_log,(nxh+1),
     >               'zero_slab23',zero_slab23(2,2),zero_slab23(1,2))
         value = value.and.BA_alloc_get(mt_log,(nxh+1),
     >               'zero_slab23',zero_slab23(2,3),zero_slab23(1,3))
         if (.not. value) call errquit('out of heap memory',0, MA_ERR)

         end if


         if (has_small) then
            nbmax = 3
         else
            nbmax = 1
         end if


         do nb=0,nbmax

         if (nb.lt.2) then
            id = 1
            nxh = nx(1)/2
            nyh = ny(1)/2
            nzh = nz(1)/2
            unitg(1,1) = lattice_unitg(1,1)
            unitg(2,1) = lattice_unitg(2,1)
            unitg(3,1) = lattice_unitg(3,1)
            unitg(1,2) = lattice_unitg(1,2)
            unitg(2,2) = lattice_unitg(2,2)
            unitg(3,2) = lattice_unitg(3,2)
            unitg(1,3) = lattice_unitg(1,3)
            unitg(2,3) = lattice_unitg(2,3)
            unitg(3,3) = lattice_unitg(3,3)
         else
            id = 3
            nxh = nx(3)/2
            nyh = ny(3)/2
            nzh = nz(3)/2
            unitg(1,1) = lattice_unitg_small(1,1)
            unitg(2,1) = lattice_unitg_small(2,1)
            unitg(3,1) = lattice_unitg_small(3,1)
            unitg(1,2) = lattice_unitg_small(1,2)
            unitg(2,2) = lattice_unitg_small(2,2)
            unitg(3,2) = lattice_unitg_small(3,2)
            unitg(1,3) = lattice_unitg_small(1,3)
            unitg(2,3) = lattice_unitg_small(2,3)
            unitg(3,3) = lattice_unitg_small(3,3)
         end if


         if ((nb.eq.0).or.(nb.eq.2)) then
            ggcut = lattice_ggcut()
         else
            ggcut = lattice_wggcut()
         end if

*        **** find zero_row3 - (i,j,*) rows that are zero ****
         do q = 1,nq3(id)
           log_mb(zero_row3(1,nb)+q-1) =.true.
         end do
         do q = 1,(nxh+1)*ny(id)
           log_mb(zero_arow3(1)+q-1) =.true.
         end do

         do k2 = -nyh+1, nyh-1
         do k1 = 0,nxh-1
          i=k1
          j=k2
          if (i .lt. 0) i = i + nx(id)
          if (j .lt. 0) j = j + ny(id)
          i=i+1
          j=j+1
          zrow = .true.
          do k3 = -nzh+1, nzh-1
            g1 = k1*unitg(1,1)
     >         + k2*unitg(1,2)
     >         + k3*unitg(1,3)
            g2 = k1*unitg(2,1)
     >         + k2*unitg(2,2)
     >         + k3*unitg(2,3)
            g3 = k1*unitg(3,1)
     >         + k2*unitg(3,2)
     >         + k3*unitg(3,3)
            gg = g1*g1 + g2*g2 + g3*g3
            gg= gg-ggcut
            if (gg.lt.-eps) zrow = .false.
          end do
          if (.not.zrow) then
            log_mb(zero_arow3(1)+i-1+(nxh+1)*(j-1)) =.false.
            q = int_mb(q_map3(1,id)+(i-1)+(j-1)*(nx(id)/2+1))
            p = int_mb(p_map3(1,id)+(i-1)+(j-1)*(nx(id)/2+1))
            if (p.eq.taskid) then
              log_mb(zero_row3(1,nb)+q-1) =.false.
            end if
          end if

         end do
         end do

*        **** find zero_slab23 - (i,*,*) slabs that are zero ****
         do k1 = 1,nxh+1
           log_mb(zero_slab23(1,nb)+k1-1) =.true.
         end do

         do k1 = 0,nxh-1
          i=k1
          if (i .lt. 0) i = i + nx(id)
          i=i+1
          yzslab = .true.
          do k3 = -nzh+1, nzh-1
          do k2 = -nyh+1, nyh-1
            g1 = k1*unitg(1,1)
     >         + k2*unitg(1,2)
     >         + k3*unitg(1,3)
            g2 = k1*unitg(2,1)
     >         + k2*unitg(2,2)
     >         + k3*unitg(2,3)
            g3 = k1*unitg(3,1)
     >         + k2*unitg(3,2)
     >         + k3*unitg(3,3)
            gg = g1*g1 + g2*g2 + g3*g3
            gg= gg-ggcut
            if (gg.lt.-eps) yzslab = .false.
          end do
          end do
          if (.not.yzslab) then
            log_mb(zero_slab23(1,nb)+i-1) =.false.
          end if

         end do
         

*        **** find zero_row2 - (i,*,k) rows that are zero after fft of (i,j,*) ****
         do k = 1,nz(id)
         do i = 1,nxh+1
           q = int_mb(q_map2(1,id)+(k-1)+(i-1)*(nz(id)))
           p = int_mb(p_map2(1,id)+(k-1)+(i-1)*(nz(id)))

           log_mb(zero_arow2(1)+i-1+(nxh+1)*(k-1))
     >       = log_mb(zero_slab23(1,nb)+i-1)
           if (p.eq.taskid) then
            log_mb(zero_row2(1,nb)+q-1)
     >       = log_mb(zero_slab23(1,nb)+i-1)
           end if
         end do
         end do
         
         call D3dB_c_ptranspose_ijk_init(nb,
     >                                   log_mb(zero_arow2(1)),
     >                                   log_mb(zero_arow3(1)))

         
         end do

       value =           BA_pop_stack(zero_arow3(2))
       value = value.and.BA_pop_stack(zero_arow2(2))
       if (.not. value) call errquit('error pop stack',0,MA_ERR)


      end if

      call D3dB_pfft3_queue_init(control_pfft3_qsize())

      return
      end




*     ***********************************
*     *                                 *
*     *         D3dB_pfft_end           *
*     *                                 *
*     ***********************************

      subroutine D3dB_pfft_end()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "D3dB.fh"
#include "D3dB_pfft.fh"



*     **** indexing variables ****
      integer iq_to_i1(2,0:3)
      integer iq_to_i2(2,0:3)
      integer iz_to_i2(2,0:3)
      integer i1_start(2,0:3)
      integer i2_start(2,0:3)
      common / ptrans_blk1 / iq_to_i1,iq_to_i2,iz_to_i2,
     >                       i1_start,i2_start

*     **** indexing variables ****
      integer jq_to_i1(2,0:3)
      integer jq_to_i2(2,0:3)
      integer jz_to_i2(2,0:3)
      integer j1_start(2,0:3)
      integer j2_start(2,0:3)
      common / ptrans_blk2 / jq_to_i1,jq_to_i2,jz_to_i2,
     >                       j1_start,j2_start


*     *** hilbert tranpose data structure ****
      integer h_iq_to_i1(2,6,0:3)
      integer h_iq_to_i2(2,6,0:3)
      integer h_iz_to_i2(2,6,0:3)
      integer h_iz_to_i2_count(6,0:3)
      integer h_i1_start(2,6,0:3)
      integer h_i2_start(2,6,0:3)
      common / ptrans_blk_ijk / h_iq_to_i1,
     >                         h_iq_to_i2,
     >                         h_iz_to_i2,
     >                         h_iz_to_i2_count,
     >                         h_i1_start,
     >                         h_i2_start


      logical value,has_small
      integer i
      logical  lattice_has_small
      external lattice_has_small

      has_small = lattice_has_small()

      value =           BA_free_heap(zero_row2(2,0))
      value = value.and.BA_free_heap(zero_row2(2,1))
      value = value.and.BA_free_heap(zero_row3(2,0))
      value = value.and.BA_free_heap(zero_row3(2,1))
      value = value.and.BA_free_heap(zero_slab23(2,0))
      value = value.and.BA_free_heap(zero_slab23(2,1))

      if (has_small) then
      value = value.and.BA_free_heap(zero_row2(2,2))
      value = value.and.BA_free_heap(zero_row2(2,3))
      value = value.and.BA_free_heap(zero_row3(2,2))
      value = value.and.BA_free_heap(zero_row3(2,3))
      value = value.and.BA_free_heap(zero_slab23(2,2))
      value = value.and.BA_free_heap(zero_slab23(2,3))
      end if

      if (mapping.eq.1) then
      value = value.and.BA_free_heap(iq_to_i1(2,0))
      value = value.and.BA_free_heap(iq_to_i1(2,1))
      value = value.and.BA_free_heap(iq_to_i2(2,0))
      value = value.and.BA_free_heap(iq_to_i2(2,1))
      value = value.and.BA_free_heap(iz_to_i2(2,0))
      value = value.and.BA_free_heap(iz_to_i2(2,1))
      value = value.and.BA_free_heap(i1_start(2,0))
      value = value.and.BA_free_heap(i1_start(2,1))
      value = value.and.BA_free_heap(i2_start(2,0))
      value = value.and.BA_free_heap(i2_start(2,1))

      value = value.and.BA_free_heap(jq_to_i1(2,0))
      value = value.and.BA_free_heap(jq_to_i1(2,1))
      value = value.and.BA_free_heap(jq_to_i2(2,0))
      value = value.and.BA_free_heap(jq_to_i2(2,1))
      value = value.and.BA_free_heap(jz_to_i2(2,0))
      value = value.and.BA_free_heap(jz_to_i2(2,1))
      value = value.and.BA_free_heap(j1_start(2,0))
      value = value.and.BA_free_heap(j1_start(2,1))
      value = value.and.BA_free_heap(j2_start(2,0))
      value = value.and.BA_free_heap(j2_start(2,1))

      if (has_small) then
      value = value.and.BA_free_heap(iq_to_i1(2,2))
      value = value.and.BA_free_heap(iq_to_i1(2,3))
      value = value.and.BA_free_heap(iq_to_i2(2,2))
      value = value.and.BA_free_heap(iq_to_i2(2,3))
      value = value.and.BA_free_heap(iz_to_i2(2,2))
      value = value.and.BA_free_heap(iz_to_i2(2,3))
      value = value.and.BA_free_heap(i1_start(2,2))
      value = value.and.BA_free_heap(i1_start(2,3))
      value = value.and.BA_free_heap(i2_start(2,2))
      value = value.and.BA_free_heap(i2_start(2,3))

      value = value.and.BA_free_heap(jq_to_i1(2,2))
      value = value.and.BA_free_heap(jq_to_i1(2,3))
      value = value.and.BA_free_heap(jq_to_i2(2,2))
      value = value.and.BA_free_heap(jq_to_i2(2,3))
      value = value.and.BA_free_heap(jz_to_i2(2,2))
      value = value.and.BA_free_heap(jz_to_i2(2,3))
      value = value.and.BA_free_heap(j1_start(2,2))
      value = value.and.BA_free_heap(j1_start(2,3))
      value = value.and.BA_free_heap(j2_start(2,2))
      value = value.and.BA_free_heap(j2_start(2,3))
      end if

      end if

      if (mapping.eq.2) then
      do i=1,6
      value = value.and.BA_free_heap(h_i1_start(2,i,0))
      value = value.and.BA_free_heap(h_i2_start(2,i,0))
      value = value.and.BA_free_heap(h_iq_to_i1(2,i,0))
      value = value.and.BA_free_heap(h_iq_to_i2(2,i,0))
      value = value.and.BA_free_heap(h_iz_to_i2(2,i,0))
      value = value.and.BA_free_heap(h_i1_start(2,i,1))
      value = value.and.BA_free_heap(h_i2_start(2,i,1))
      value = value.and.BA_free_heap(h_iq_to_i1(2,i,1))
      value = value.and.BA_free_heap(h_iq_to_i2(2,i,1))
      value = value.and.BA_free_heap(h_iz_to_i2(2,i,1))

      if (has_small) then
      value = value.and.BA_free_heap(h_i1_start(2,i,2))
      value = value.and.BA_free_heap(h_i2_start(2,i,2))
      value = value.and.BA_free_heap(h_iq_to_i1(2,i,2))
      value = value.and.BA_free_heap(h_iq_to_i2(2,i,2))
      value = value.and.BA_free_heap(h_iz_to_i2(2,i,2))
      value = value.and.BA_free_heap(h_i1_start(2,i,3))
      value = value.and.BA_free_heap(h_i2_start(2,i,3))
      value = value.and.BA_free_heap(h_iq_to_i1(2,i,3))
      value = value.and.BA_free_heap(h_iq_to_i2(2,i,3))
      value = value.and.BA_free_heap(h_iz_to_i2(2,i,3))
      end if

      end do
      end if

      if (.not. value) 
     > call errquit('D3dB_pfft_end:error freeing heap',0,MA_ERR)

      call D3dB_pfft3_queue_end()
      return
      end




*     ***********************************
*     *					*
*     *	        D3dB_cr_pfft3b		*
*     *					*
*     ***********************************

      subroutine D3dB_cr_pfft3b(nb,nbb,A)

*****************************************************
*                                                   *
*      This routine performs the operation of       *
*      a three dimensional complex to complex       *
*      inverse fft                                  *
*           A(nx,ny(nb),nz(nb)) <- FFT3^(-1)[A(kx,ky,kz)]   * 
*                                                   *
*      Entry - 					    *
*              A: a column distribuded 3d block     *
*              tmp: tempory work space must be at   *
*                    least the size of (complex)    *
*                    (nfft*nfft + 1) + 10*nfft      * 
*                                                   *
*       Exit - A is transformed and the imaginary   *
*              part of A is set to zero             *
*       uses - D3dB_c_ptranspose_jk, dcopy           *
*                                                   *
*****************************************************

      implicit none
      integer nb,nbb
      complex*16  A(*)

#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"

      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / D3dB_fft / tmpx,tmpy,tmpz


*     *** local variables ***
      integer i,j,k,q,indx,ierr
      integer nxh,nxhy,nxhz,indx0,indx1,indx2

      
      !integer tmp1(2),tmp2(2),tmp3(2)
      integer tmp2(2),tmp3(2)
      logical value
      integer tid,nthr,offset
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads


      call nwpw_timing_start(1)
      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

*     ***** allocate temporary space ****
      !call D3dB_nfft3d(nb,nfft3d)
      value = BA_push_get(mt_dcpl,(nfft3d(nb)),'ffttmp2',
     >                    tmp2(2),tmp2(1))
      value = value.and.
     >      BA_push_get(mt_dbl,(n2ft3d(nb)),'ffttmp3',tmp3(2),tmp3(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

       nxh = (nx(nb)/2+1)
       nxhz = nxh*nz(nb)
       nxhy = nxh*ny(nb)
      
      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     ********************************************
*     ***         Do a transpose of A          ***
*     ***      A(kx,kz,ky) <- A(kx,ky,kz)      ***
*     ********************************************
c     call D3dB_c_transpose_jk(nb,A,dcpl_mb(tmp2(1)),dbl_mb(tmp3(1)))

*     *************************************************
*     ***     do fft along kz dimension             ***
*     ***   A(kx,nz(nb),ky) <- fft1d^(-1)[A(kx,kz,ky)]  ***
*     *************************************************
#ifdef MLIB
      !call z1dfft(dbl_mb(tmp3(1)),nz(nb),dcpl_mb(tmpz(1)),-3,ierr)
      indx0 = 0
      indx2 = 0
      do q=1,nq(nb)
      do i=1,nxh
       indx2 = indx2 + 1
       if (.not.log_mb(zero_row3(1,nbb)+indx2-1)) then

         indx = i+indx0
         indx1 = indx
         do k=1,nz(nb)
            dcpl_mb(tmp2(1)+k-1) = A(indx)
            indx = indx + nxh
         end do
         call z1dfft(dcpl_mb(tmp2(1)),nz(nb),
     >               dcpl_mb(tmpz(1,nb)),-2,ierr)
         do k=1,nz(nb)
            A(indx1) = dcpl_mb(tmp2(1)+k-1)
            indx1 = indx1 + nxh
         end do
       end if
      end do
      indx0 = indx0+nxhz
      end do
      !call dscal((nx(nb)+2)*ny(nb)*nq(nb),dble(nz(nb)),A,1)

#else
      !call dcffti(nz(nb),dcpl_mb(tmp1(1)))
      indx0 = 0
      indx2 = 0
      do q=1,nq(nb)
      do i=1,nxh
       indx2 = indx2 + 1
       if (.not.log_mb(zero_row3(1,nbb)+indx2-1)) then

#ifdef FFTW3
         indx  = i + indx0
         call dfftw_execute_dft(plans(7,nb),A(indx),A(indx))

#else
         indx  = i + indx0
         indx1 = indx
         do k=1,nz(nb)
            dcpl_mb(tmp2(1)+k-1) = A(indx)
            indx = indx + nxh
         end do
         call dcfftb(nz(nb),dcpl_mb(tmp2(1)),dcpl_mb(tmpz(1,nb)))
         do k=1,nz(nb)
            A(indx1) = dcpl_mb(tmp2(1)+k-1)
            indx1 = indx1 + nxh
         end do
#endif

       end if
      end do
      indx0 = indx0 + nxhz
      end do
#endif

*     ********************************************
*     ***         Do a ptranspose of A          ***
*     ***      A(kx,ky,nz(nb)) <- A(kx,nz(nb),ky)      ***
*     ********************************************
      call D3dB_c_ptranspose1_jk(nbb,A,dcpl_mb(tmp2(1)),dbl_mb(tmp3(1)))

*     *************************************************
*     ***     do fft along ky dimension             ***
*     ***   A(kx,ny(nb),nz(nb)) <- fft1d^(-1)[A(kx,ky,nz(nb))]  ***
*     *************************************************
#ifdef MLIB
     
      indx0 = 0
      indx2 = 0
      do q=1,nq(nb)
      do i=1,nxh
        indx2 = indx2 + 1
        if (.not.log_mb(zero_row2(1,nbb)+indx2-1)) then

          indx  = i + indx0
          indx1 = indx
          do j=1,ny(nb)
            dcpl_mb(tmp2(1)+j-1) = A(indx)
            indx = indx + nxh
          end do
          call z1dfft(dcpl_mb(tmp2(1)),ny(nb),
     >               dcpl_mb(tmpy(1,nb)),-2,ierr)
          do j=1,ny(nb)
            A(indx1) = dcpl_mb(tmp2(1)+j-1)
            indx1 = indx1 + nxh
          end do
        end if
      end do
      indx0 = indx0 + nxhy
      end do
      
#else
    
      indx0 = 0
      indx2 = 0
      do q=1,nq(nb)
      do i=1,nxh
        indx2 = indx2 + 1

        if (.not.log_mb(zero_row2(1,nbb)+indx2-1)) then

#ifdef FFTW3
         indx  = i + indx0
         call dfftw_execute_dft(plans(8,nb),A(indx),A(indx))
#else
         indx  = i + indx0
         indx1 = indx
         do j=1,ny(nb)
            dcpl_mb(tmp2(1)+j-1) = A(indx)
            indx = indx + nxh
         end do
         call dcfftb(ny(nb),dcpl_mb(tmp2(1)),dcpl_mb(tmpy(1,nb)))
         do j=1,ny(nb)
            A(indx1) = dcpl_mb(tmp2(1)+j-1)
            indx1 = indx1 + nxh
         end do
#endif

        end if
      end do
      indx0 = indx0 + nxhy
      end do
#endif

*     *************************************************
*     ***     do fft along kx dimension             ***
*     ***   A(nx(nb),ny(nb),nz(nb)) <- fft1d^(-1)[A(kx,ny(nb),nz(nb))]  ***
*     *************************************************
#ifdef MLIB
      !call drc1ft (dbl_mb(tmp3(1)),nx(nb),dcpl_mb(tmp1(1)),-3,ierr)
      indx = 1
      do q=1,nq(nb)
      do j=1,ny(nb)
         call drc1ft(A(indx),nx(nb),dcpl_mb(tmpx(1,nb)),-2,ierr)
         indx = indx + nxh
      end do
      end do

#else
#ifdef FFTW3
      call dfftw_execute_dft_c2r(plans(3,nb),A,A)

#else
      call cshift1_fftb(nx(nb),ny(nb),nq(nb),1,A)
      indx = 1
      do q=1,nq(nb)
      do j=1,ny(nb)
         !indx = 1 + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
         call drfftb(nx(nb),A(indx),dcpl_mb(tmpx(1,nb)))
         indx = indx + nxh
      end do
      end do
      call zeroend_fftb(nx(nb),ny(nb),nq(nb),1,A)
#endif
#endif


      !*************************
      !**** hilbert mapping ****
      !*************************
      else


*     *************************************************
*     ***     do fft along kz dimension             ***
*     ***   A(nz(nb),kx,ky) <- fft1d^(-1)[A(kz,kx,ky)]  ***
*     *************************************************
#ifdef MLIB
      indx = 1
      do q=1,nq3(nb)
         if (.not.log_mb(zero_row3(1,nbb)+q-1)) then
         call z1dfft(A(indx),nz(nb),dcpl_mb(tmpz(1,nb)),-2,ierr)
         end if
         indx = indx + nz(nb)
      end do
#else
#ifdef FFTW3
      indx = 1
      do q=1,nq3(nb)
         if (.not.log_mb(zero_row3(1,nbb)+q-1)) then
            call dfftw_execute_dft(plans(17,nb),A(indx),A(indx))
         end if
         indx = indx + nz(nb)
      end do
#else

      !indx = 1
      offset = tid*(2*nz(nb)+15)
      do q=tid+1,nq3(nb),nthr
         if (.not.log_mb(zero_row3(1,nbb)+q-1)) then
            call dcfftb(nz(nb),A(1+(q-1)*nz(nb)),
     >                  dcpl_mb(tmpz(1,nb)+offset))
         end if
         !indx = indx + nz(nb)
      end do
!$OMP BARRIER

#endif
#endif

      call D3dB_c_ptranspose_ijk(nbb,3,A,
     >                           dcpl_mb(tmp2(1)),
     >                           dbl_mb(tmp3(1)))

*     *************************************************
*     ***     do fft along ky dimension             ***
*     ***   A(ny(nb),nz(nb),kx) <- fft1d^(-1)[A(ky,nz(nb),kx)]  ***
*     *************************************************
#ifdef MLIB
      indx = 1
      do q=1,nq2(nb)
         if (.not.log_mb(zero_row2(1,nbb)+q-1)) then
         call z1dfft(A(indx),ny(nb),dcpl_mb(tmpy(1,nb)),-2,ierr)
         end if
         indx = indx + ny(nb)
      end do
#else
#ifdef FFTW3
      indx = 1
      do q=1,nq2(nb)
         if (.not.log_mb(zero_row2(1,nbb)+q-1)) then
            call dfftw_execute_dft(plans(18,nb),A(indx),A(indx))
         end if
         indx = indx + ny(nb)
      end do
#else

      !indx = 1
      offset = tid*(2*ny(nb)+15)
      do q=tid+1,nq2(nb),nthr
         if (.not.log_mb(zero_row2(1,nbb)+q-1)) then
            call dcfftb(ny(nb),A(1+(q-1)*ny(nb)),
     >                  dcpl_mb(tmpy(1,nb)+offset))
         end if
         !indx = indx + ny(nb)
      end do
!$OMP BARRIER

#endif
#endif

      call D3dB_c_ptranspose_ijk(nbb,4,A,
     >                           dcpl_mb(tmp2(1)),
     >                           dbl_mb(tmp3(1)))

*     *************************************************
*     ***     do fft along kx dimension             ***
*     ***   A(nx(nb),ny(nb),nz(nb)) <- fft1d^(-1)[A(kx,ny(nb),nz(nb))]  ***
*     *************************************************
#ifdef MLIB
      indx = 1
      do q=1,nq1(nb)
         !indx = 1 + (q-1)*(nx(nb)/2+1)
         call drc1ft(A(indx),nx(nb),dcpl_mb(tmpx(1,nb)),-2,ierr)
         indx = indx + nxh
      end do
#else

#ifdef FFTW3
         call dfftw_execute_dft_c2r(plans(13,nb),A,A)
#else

      !indx = 1
      offset = tid*(2*nx(nb)+15)
      call cshift1_fftb(nx(nb),nq1(nb),1,1,A)
      do q=tid+1,nq1(nb),nthr
         call drfftb(nx(nb),A(1+(q-1)*nxh),
     >              dcpl_mb(tmpx(1,nb)+offset))
         !indx = indx + nxh
      end do
!$OMP BARRIER
      call zeroend_fftb(nx(nb),nq1(nb),1,1,A)

#endif
#endif

      end if
    
*     **** deallocate temporary space  ****
      value = BA_pop_stack(tmp3(2))
      value = BA_pop_stack(tmp2(2))
      !value = BA_pop_stack(tmp1(2))

      call nwpw_timing_end(1)
      return
      end




*     ***********************************
*     *					*
*     *	        D3dB_rc_pfft3f		*
*     *					*
*     ***********************************

      subroutine D3dB_rc_pfft3f(nb,nbb,A)

*****************************************************
*                                                   *
*      This routine performs the operation of       *
*      a three dimensional complex to complex fft   *
*           A(kx,ky,kz) <- FFT3[A(nx(nb),ny(nb),nz(nb))]        * 
*                                                   *
*      Entry - 					    *
*              A: a column distribuded 3d block     *
*              tmp: tempory work space must be at   *
*                    least the size of (complex)    *
*                    (nfft*nfft + 1) + 10*nfft      * 
*                                                   *
*       Exit - A is transformed                     *
*                                                   *
*       uses - ptranspose1 subroutine                *
*                                                   *
*****************************************************

      implicit none
      integer nb,nbb
      complex*16  A(*)

#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"

      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / D3dB_fft / tmpx,tmpy,tmpz


*     *** local variables ***
      integer i,j,k,q,indx,indx1,ierr,indx2,indx0
      integer nxh,nxhy,nxhz

      !integer tmp1(2),tmp2(2),tmp3(2)
      integer tmp2(2),tmp3(2)
      logical value
      integer tid,nthr,offset
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads


      call nwpw_timing_start(1)

      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

*     ***** allocate temporary space ****
      !call D3dB_nfft3d(nb,nfft3d)
      value = BA_push_get(mt_dcpl,(nfft3d(nb)),'tmp2',tmp2(2),tmp2(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(n2ft3d(nb)),'tmp3',tmp3(2),tmp3(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      nxh = (nx(nb)/2+1)
      nxhz = nxh*nz(nb)
      nxhy = nxh*ny(nb)

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     ********************************************
*     ***     do fft along nx(nb) dimension        ***
*     ***   A(kx,ny(nb),nz(nb)) <- fft1d[A(nx(nb),ny(nb),nz(nb))]  ***
*     ********************************************
#ifdef MLIB
      call drcfts(A,nx(nb),1,ny(nb)*nq(nb),
     >                  nx(nb)+2,1,ierr)

#else

#ifdef FFTW3
      call dfftw_execute_dft_r2c(plans(4,nb),A,A)

#else
      !call drffti(nx(nb),dcpl_mb(tmp1(1)))
      indx = 1
      do q=1,nq(nb)
      do j=1,ny(nb)
         call drfftf(nx(nb),A(indx),dcpl_mb(tmpx(1,nb)))
         indx = indx + nxh
      end do
      end do
      call cshift_fftf(nx(nb),ny(nb),nq(nb),1,A)
#endif
#endif


*     ********************************************
*     ***     do fft along ny(nb) dimension        ***
*     ***   A(kx,ky,nz(nb)) <- fft1d[A(kx,ny(nb),nz(nb))]  ***
*     ********************************************

#ifdef MLIB
      
      indx0 = 0
      indx2 = 0
      do q=1,nq(nb)
      do i=1,nxh
        indx2 = indx2 + 1
        if (.not.log_mb(zero_row2(1,nbb)+indx2-1)) then

         indx  = i + indx0
         indx1 = indx
         do j=1,ny(nb)
            dcpl_mb(tmp2(1)+j-1) = A(indx)
            indx = indx + nxh
         end do
         call z1dfft(dcpl_mb(tmp2(1)),ny(nb),dcpl_mb(tmpy(1,nb)),1,ierr)
         do j=1,ny(nb)
            A(indx1) = dcpl_mb(tmp2(1)+j-1)
            indx1 = indx1 + nxh
         end do
        
        end if
      end do
      indx0 = indx0 + nxhy
      end do
#else

       indx0 = 0
       indx2 = 0
       do q=1,nq(nb)
       do i=1,nxh
        indx2 = indx2 + 1
        if (.not.log_mb(zero_row2(1,nbb)+indx2-1)) then

#ifdef FFTW3
          indx  = i + indx0
          call dfftw_execute_dft(plans(9,nb),A(indx),A(indx))

#else
          indx  = i + indx0
          indx1 = indx
          do k=1,ny(nb)
             dcpl_mb(tmp2(1)+k-1) = A(indx)
             indx = indx + nxh
          end do
          call dcfftf(ny(nb),dcpl_mb(tmp2(1)),dcpl_mb(tmpy(1,nb)))
          do k=1,ny(nb)
             A(indx1) = dcpl_mb(tmp2(1)+k-1)
             indx1 = indx1 + nxh
          end do
#endif

        end if
       end do
       indx0 = indx0 + nxhy
       end do

#endif


*     ********************************************
*     ***         Do a ptranspose of A          ***
*     ***      A(ky,nz(nb),ky) <- A(kx,ky,nz(nb))      ***
*     ********************************************
      call D3dB_c_ptranspose2_jk(nbb,A,dcpl_mb(tmp2(1)),dbl_mb(tmp3(1)))
c      call D3dB_c_transpose_jk(nbb,A,dcpl_mb(tmp2(1)),dbl_mb(tmp3(1)))


*     ********************************************
*     ***     do fft along nz(nb) dimension        ***
*     ***   A(kx,kz,ky) <- fft1d[A(kx,nz(nb),ky)]  ***
*     ********************************************
#ifdef MLIB
      
      indx0 = 0
      indx2 = 0
      do q=1,nq(nb)
      do i=1,nxh
        indx2 = indx2 + 1
        if (.not.log_mb(zero_row3(1,nbb)+indx2-1)) then

         indx  = i + indx0
         indx1 = indx
         do k=1,nz(nb)
            dcpl_mb(tmp2(1)+k-1) = A(indx)
            indx = indx + nxh
         end do
         call z1dfft(dcpl_mb(tmp2(1)),nz(nb),dcpl_mb(tmpz(1,nb)),1,ierr)
         do k=1,nz(nb)
            A(indx1) = dcpl_mb(tmp2(1)+k-1)
            indx1 = indx1 + nxh
         end do
        end if
      end do
      indx0 = indx0 + nxhz
      end do
#else

       indx0 = 0
       indx2 = 0
       do q=1,nq(nb)
       do i=1,nxh
        indx2 = indx2 + 1
        if (.not.log_mb(zero_row3(1,nbb)+indx2-1)) then
 
#ifdef FFTW3
          indx  = i + indx0
          call dfftw_execute_dft(plans(10,nb),A(indx),A(indx))

#else
          indx  = i + indx0
          indx1 = indx
          do k=1,nz(nb)
             dcpl_mb(tmp2(1)+k-1) = A(indx)
             indx = indx + nxh
          end do
          call dcfftf(nz(nb),dcpl_mb(tmp2(1)),dcpl_mb(tmpz(1,nb)))
          do k=1,nz(nb)
             A(indx1) = dcpl_mb(tmp2(1)+k-1)
             indx1 = indx1 + nxh
          end do
#endif
 
        end if
       end do
       indx0 = indx0 + nxhz
       end do

#endif

*     ********************************************
*     ***         Do a transpose of A          ***
*     ***      A(kx,ky,kz) <- A(kx,kz,ky)      ***
*     ********************************************
c     call D3dB_c_transpose_jk(nb,A,dcpl_mb(tmp2(1)),dbl_mb(tmp3(1)))



      !*************************
      !**** hilbert mapping ****
      !*************************
      else

*     ********************************************
*     ***     do fft along nx(nb) dimension        ***
*     ***   A(kx,ny(nb),nz(nb)) <- fft1d[A(nx(nb),ny(nb),nz(nb))]  ***
*     ********************************************
#ifdef MLIB
      call drcfts(A,nx(nb),1,nq1(nb),
     >                  nx(nb)+2,1,ierr)
#else

#ifdef FFTW3
          call dfftw_execute_dft_r2c(plans(14,nb),A,A)

#else
      !indx = 1
      offset = tid*(2*nx(nb)+15)
      do q=tid+1,nq1(nb),nthr
         call drfftf(nx(nb),A(1+(q-1)*nxh),dcpl_mb(tmpx(1,nb)+offset))
         !indx = indx + nxh
      end do
      call cshift_fftf(nx(nb),nq1(nb),1,1,A)

#endif
#endif

      call D3dB_c_ptranspose_ijk(nbb,1,A,
     >                           dcpl_mb(tmp2(1)),
     >                           dbl_mb(tmp3(1)))

*     ********************************************
*     ***     do fft along ny(nb) dimension        ***
*     ***   A(ky,nz(nb),kx) <- fft1d[A(ny(nb),nz(nb),kx)]  ***
*     ********************************************
#ifdef MLIB
      indx = 1
      do q=1,nq2(nb)
         if (.not.log_mb(zero_row2(1,nbb)+q-1)) then
         call z1dfft(A(indx),ny(nb),dcpl_mb(tmpy(1,nb)),1,ierr)
         end if
         indx = indx + ny(nb)
      end do
#else
#ifdef FFTW3
      indx = 1
      do q=1,nq2(nb)
         if (.not.log_mb(zero_row2(1,nbb)+q-1)) then
            call dfftw_execute_dft(plans(19,nb),A(indx),A(indx))
         end if
         indx = indx + ny(nb)
      end do
#else

      !indx = 1
      offset = tid*(2*ny(nb)+15)
      do q=tid+1,nq2(nb),nthr
         if (.not.log_mb(zero_row2(1,nbb)+q-1)) then
            call dcfftf(ny(nb),A(1+(q-1)*ny(nb)),
     >                  dcpl_mb(tmpy(1,nb)+offset))
         end if
         !indx = indx + ny(nb)
      end do
!$OMP BARRIER

#endif
#endif

      call D3dB_c_ptranspose_ijk(nbb,2,A,
     >                          dcpl_mb(tmp2(1)),
     >                          dbl_mb(tmp3(1)))

*     ********************************************
*     ***     do fft along nz(nb) dimension        ***
*     ***   A(kz,kx,ky) <- fft1d[A(nz(nb),kx,ky)]  ***
*     ********************************************
#ifdef MLIB
      indx = 1
      do q=1,nq3(nb)
         if (.not.log_mb(zero_row3(1,nbb)+q-1)) then
         call z1dfft(A(indx),nz(nb),dcpl_mb(tmpz(1,nb)),1,ierr)
         end if
         indx = indx + nz(nb)
      end do
#else
#ifdef FFTW3
      indx = 1
      do q=1,nq3(nb)
         if (.not.log_mb(zero_row3(1,nbb)+q-1)) then
            call dfftw_execute_dft(plans(20,nb),A(indx),A(indx))
         end if
         indx = indx + nz(nb)
      end do
#else

      !indx = 1
      offset = tid*(2*nz(nb)+15)
      do q=tid+1,nq3(nb),nthr
         if (.not.log_mb(zero_row3(1,nbb)+q-1)) then
            call dcfftf(nz(nb),A(1+(q-1)*nz(nb)),
     >                  dcpl_mb(tmpz(1,nb)+offset))
         end if
         !indx = indx + nz(nb)
      end do
!$OMP BARRIER

#endif
#endif

      end if

*     **** deallocate temporary space  ****
      value = BA_pop_stack(tmp3(2))
      value = BA_pop_stack(tmp2(2))
      !value = BA_pop_stack(tmp1(2))
      call nwpw_timing_end(1)
      return
      end





*     ***********************************
*     *					*
*     *	   D3dB_c_ptranspose_jk_init	*
*     *					*
*     ***********************************

      subroutine D3dB_c_ptranspose_jk_init(nbb,zero_arow3)
      implicit none
      integer nbb
      logical zero_arow3(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "D3dB.fh"


*     **** indexing variables ****
      integer iq_to_i1(2,0:3)
      integer iq_to_i2(2,0:3)
      integer iz_to_i2(2,0:3)
      integer i1_start(2,0:3)
      integer i2_start(2,0:3)
      common / ptrans_blk1 / iq_to_i1,iq_to_i2,iz_to_i2,
     >                       i1_start,i2_start

*     **** indexing variables ****
      integer jq_to_i1(2,0:3)
      integer jq_to_i2(2,0:3)
      integer jz_to_i2(2,0:3)
      integer j1_start(2,0:3)
      integer j2_start(2,0:3)
      common / ptrans_blk2 / jq_to_i1,jq_to_i2,jz_to_i2,
     >                       j1_start,j2_start



*     **** local variables ****
      integer proc_to,proc_from
      integer pto,qto,np,taskid
      integer pfrom,qfrom
      integer phere,qhere
      integer index1,index2,index3
      integer jndex1,jndex2,jndex3
      integer itmp,ii,jj
      integer i,j,k,it,id
      logical value

      if (nbb.lt.2) then
         id = 1
      else
         id = 3
      end if

*     **** allocate ptrans_blk1 and ptrans_blk2 common block ****
      value = BA_alloc_get(mt_int,((nx(id)/2+1)*ny(id)*nq(id)),
     >                     'piq_to_i1',iq_to_i1(2,nbb),iq_to_i1(1,nbb))
      value=value.and.BA_alloc_get(mt_int,((nx(id)/2+1)*ny(id)*nq(id)),
     >                     'piq_to_i2',iq_to_i2(2,nbb),iq_to_i2(1,nbb))
      value=value.and.BA_alloc_get(mt_int,((nx(id)/2+1)*ny(id)*nq(id)),
     >                     'piz_to_i2',iz_to_i2(2,nbb),iz_to_i2(1,nbb))

      value = value.and.BA_alloc_get(mt_int,(nz(id)+1),
     >                     'pi1_start',i1_start(2,nbb),i1_start(1,nbb))
      value = value.and.BA_alloc_get(mt_int,(nz(id)+1),
     >                     'pi2_start',i2_start(2,nbb),i2_start(1,nbb))


      value = BA_alloc_get(mt_int,((nx(id)/2+1)*ny(id)*nq(id)),
     >                     'riq_to_i1',jq_to_i1(2,nbb),jq_to_i1(1,nbb))
      value=value.and.BA_alloc_get(mt_int,((nx(id)/2+1)*ny(id)*nq(id)),
     >                     'riq_to_i2',jq_to_i2(2,nbb),jq_to_i2(1,nbb))
      value=value.and.BA_alloc_get(mt_int,((nx(id)/2+1)*ny(id)*nq(id)),
     >                     'riz_to_i2',jz_to_i2(2,nbb),jz_to_i2(1,nbb))

      value = value.and.BA_alloc_get(mt_int,(nz(id)+1),
     >                     'ri1_start',j1_start(2,nbb),j1_start(1,nbb))
      value = value.and.BA_alloc_get(mt_int,(nz(id)+1),
     >                     'ri2_start',j2_start(2,nbb),j2_start(1,nbb))

      if (.not. value) 
     > call errquit('D3dB_ptranspose_jk_init:out of heap',0,MA_ERR)

      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)


      index1 = 1 
      index2 = 1
      index3 = 1
      jndex1 = 1 
      jndex2 = 1
      jndex3 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
         int_mb(i1_start(1,nbb)+it) = index1
         int_mb(i2_start(1,nbb)+it) = index2
         int_mb(j1_start(1,nbb)+it) = jndex1
         int_mb(j2_start(1,nbb)+it) = jndex2

         do k=1,nz(id)
         do j=1,ny(id)

*           **** packing scheme **** 
            call D3dB_ktoqp(id,k,qhere,phere)
            call D3dB_ktoqp(id,j,qto,pto)
            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
               do i=1,(nx(id)/2+1)
                  ii = i + (nx(id)/2+1)*(k-1)
                  jj = i + (nx(id)/2+1)*(j-1)
                  itmp = i + (j-1)*(nx(id)/2+1) 
     >                     + (qhere-1)*(nx(id)/2+1)*ny(id)

                  if (.not.zero_arow3(ii)) then
                  int_mb(iq_to_i1(1,nbb)+index1-1) = itmp
                  index1 = index1 + 1
                  end if

                  if (.not.zero_arow3(jj)) then
                  int_mb(jq_to_i1(1,nbb)+jndex1-1) = itmp
                  jndex1 = jndex1 + 1
                  end if

               end do
            end if
             
*           **** unpacking scheme ****
            call D3dB_ktoqp(id,j,qhere,phere)
            call D3dB_ktoqp(id,k,qfrom,pfrom)
            if ((phere.eq.taskid).and.(pfrom.eq.proc_from)) then
               do i=1,(nx(id)/2+1)
                  ii = i + (nx(id)/2+1)*(k-1)
                  jj = i + (nx(id)/2+1)*(j-1)
                  itmp = i + (k-1)*(nx(id)/2+1) 
     >                     + (qhere-1)*(nx(id)/2+1)*ny(id)
                  if (zero_arow3(ii)) then
                  int_mb(iz_to_i2(1,nbb)+index3-1) = itmp
                  index3 = index3 + 1
                  else
                  int_mb(iq_to_i2(1,nbb)+index2-1) = itmp
                  index2 = index2 + 1
                  end if

                  if (zero_arow3(jj)) then
                  int_mb(jz_to_i2(1,nbb)+jndex3-1) = itmp
                  jndex3 = jndex3 + 1
                  else
                  int_mb(jq_to_i2(1,nbb)+jndex2-1) = itmp
                  jndex2 = jndex2 + 1
                  end if
               end do
            end if
         end do
         end do
      end do
      int_mb(i1_start(1,nbb)+np) = index1
      int_mb(i2_start(1,nbb)+np) = index2
      int_mb(j1_start(1,nbb)+np) = jndex1
      int_mb(j2_start(1,nbb)+np) = jndex2


      return
      end





*     ***********************************
*     *					*
*     *	   D3dB_c_ptranspose_ijk_init	*
*     *					*
*     ***********************************

      subroutine D3dB_c_ptranspose_ijk_init(nbb,zero_arow2,zero_arow3)
      implicit none
      integer nbb
      logical zero_arow2(*)
      logical zero_arow3(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "D3dB.fh"


*     *** hilbert tranpose data structure ****
      integer h_iq_to_i1(2,6,0:3)
      integer h_iq_to_i2(2,6,0:3)
      integer h_iz_to_i2(2,6,0:3)
      integer h_iz_to_i2_count(6,0:3)
      integer h_i1_start(2,6,0:3)
      integer h_i2_start(2,6,0:3)
      common / ptrans_blk_ijk / h_iq_to_i1,
     >                         h_iq_to_i2,
     >                         h_iz_to_i2,
     >                         h_iz_to_i2_count,
     >                         h_i1_start,
     >                         h_i2_start

*     **** local variables ****
      logical value,iszero
      integer proc_to,proc_from
      integer pto,qto,np,taskid
      integer phere,qhere
      integer index1,index2,index3,itmp
      integer i,j,k,it,id

      if (nbb.lt.2) then
         id = 1
      else
         id = 3
      end if

      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)



*     ********************************************************
*     **** map1to2 mapping - done - tranpose operation #1  ***
*     ****   (ny,nz,nx/2+1)  <-- (nx/2+1,ny,nz)            ***
*     ****   use zero_arow2                                ***
*     ********************************************************

*     **** allocate trans_blk_ijk common block ****
      value = BA_alloc_get(mt_int,((nx(id)/2+1)*nq1(id)),
     >                     'h_iq_to_i1_1',
     >                      h_iq_to_i1(2,1,nbb),
     >                      h_iq_to_i1(1,1,nbb))
      value = value.and.
     >        BA_alloc_get(mt_int,(ny(id)*nq2(id)),
     >                     'h_iq_to_i2_1',
     >                      h_iq_to_i2(2,1,nbb),
     >                      h_iq_to_i2(1,1,nbb))
      value = value.and.
     >        BA_alloc_get(mt_int,(ny(id)*nq2(id)),
     >                     'h_iz_to_i2_1',
     >                      h_iz_to_i2(2,1,nbb),
     >                      h_iz_to_i2(1,1,nbb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i1_start_1',
     >                      h_i1_start(2,1,nbb),
     >                      h_i1_start(1,1,nbb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i2_start_1',
     >                      h_i2_start(2,1,nbb),
     >                      h_i2_start(1,1,nbb))
      if (.not.value)
     > call errquit('D3dB_transpose_ijk_initt:out of heap',0,MA_ERR)

      index1 = 1 
      index2 = 1
      index3 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
         int_mb(h_i1_start(1,1,nbb)+it) = index1
         int_mb(h_i2_start(1,1,nbb)+it) = index2

         do k=1,nz(id)
         do j=1,ny(id)
         do i=1,(nx(id)/2+1)
            iszero = zero_arow2(i+(k-1)*(nx(id)/2+1)) 

*           **** packing scheme **** 
            phere = int_mb(p_map1(1,id)+(j-1)+(k-1)*ny(id))
            qhere = int_mb(q_map1(1,id)+(j-1)+(k-1)*ny(id))

            pto   = int_mb(p_map2(1,id)+(k-1)+(i-1)*nz(id))
            qto   = int_mb(q_map2(1,id)+(k-1)+(i-1)*nz(id))
            
           
            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
               itmp = i + (qhere-1)*(nx(id)/2+1)
               if (.not.iszero) then
               int_mb(h_iq_to_i1(1,1,nbb)+index1-1) = itmp
               index1 = index1 + 1
               end if
            end if
             
*           **** unpacking scheme ****
            if ((pto.eq.taskid).and.(phere.eq.proc_from)) then
               itmp = j + (qto-1)*ny(id)
               if (.not.iszero) then
               int_mb(h_iq_to_i2(1,1,nbb)+index2-1) =  itmp
               index2 = index2 + 1
               else
               int_mb(h_iz_to_i2(1,1,nbb)+index3-1) =  itmp
               index3 = index3 + 1
               end if
            end if

         end do
         end do
         end do

      end do
      int_mb(h_i1_start(1,1,nbb)+np) = index1
      int_mb(h_i2_start(1,1,nbb)+np) = index2
      h_iz_to_i2_count(1,nbb) = index3 - 1



*     *********************************************************
*     **** map2to3 mapping - done - transpose operation #2 ****
*     ****    (nz,nx/2+1,ny)  <-- (ny,nz,nx/2+1)           ****
*     ****    use zero_arow3                               ****
*     *********************************************************

*     **** allocate trans_blk_ijk common block ****
      value = BA_alloc_get(mt_int,(ny(id)*nq2(id)),
     >                     'h_iq_to_i1_2',
     >                      h_iq_to_i1(2,2,nbb),
     >                      h_iq_to_i1(1,2,nbb))
      value = value.and.
     >        BA_alloc_get(mt_int,(nz(id)*nq3(id)),
     >                     'h_iq_to_i2_2',
     >                      h_iq_to_i2(2,2,nbb),
     >                      h_iq_to_i2(1,2,nbb))
      value = value.and.
     >        BA_alloc_get(mt_int,(nz(id)*nq3(id)),
     >                     'h_iz_to_i2_2',
     >                      h_iz_to_i2(2,2,nbb),
     >                      h_iz_to_i2(1,2,nbb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i1_start_2',
     >                      h_i1_start(2,2,nbb),
     >                      h_i1_start(1,2,nbb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i2_start_2',
     >                      h_i2_start(2,2,nbb),
     >                      h_i2_start(1,2,nbb))
      if (.not.value)
     > call errquit('D3dB_transpose_ijk_init:out of heap',0,MA_ERR)

      index1 = 1
      index2 = 1
      index3 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
         int_mb(h_i1_start(1,2,nbb)+it) = index1
         int_mb(h_i2_start(1,2,nbb)+it) = index2

         do k=1,nz(id)
         do j=1,ny(id)
         do i=1,(nx(id)/2+1)

            iszero = zero_arow3(i+(j-1)*(nx(id)/2+1))

*           **** packing scheme ****
            phere = int_mb(p_map2(1,id)+(k-1)+(i-1)*nz(id))
            qhere = int_mb(q_map2(1,id)+(k-1)+(i-1)*nz(id))

            pto   = int_mb(p_map3(1,id)+(i-1)+(j-1)*(nx(id)/2+1))
            qto   = int_mb(q_map3(1,id)+(i-1)+(j-1)*(nx(id)/2+1))


            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
               itmp = j + (qhere-1)*ny(id)
               if (.not.iszero) then
               int_mb(h_iq_to_i1(1,2,nbb)+index1-1) = itmp
               index1 = index1 + 1
               end if
            end if

*           **** unpacking scheme ****
            if ((pto.eq.taskid).and.(phere.eq.proc_from)) then
               itmp = k + (qto-1)*nz(id)
               if (.not.iszero) then
               int_mb(h_iq_to_i2(1,2,nbb)+index2-1) = itmp
               index2 = index2 + 1
               else  
               int_mb(h_iz_to_i2(1,2,nbb)+index3-1) = itmp
               index3 = index3 + 1
               end if
            end if


         end do
         end do
         end do

      end do
      int_mb(h_i1_start(1,2,nbb)+np) = index1
      int_mb(h_i2_start(1,2,nbb)+np) = index2
      h_iz_to_i2_count(2,nbb) = index3 - 1



*     ********************************************************
*     **** map3to2 mapping - done - tranpose operation #3 ****
*     ****  (ny,nz,nx/2+1)  <-- (nz,nx/2+1,ny)            ****
*     ****   use zero_arow3                               ****
*     ********************************************************

*     **** allocate trans_blk_ijk common block ****
      value = BA_alloc_get(mt_int,(nz(id)*nq3(id)),
     >                     'h_iq_to_i1_3',
     >                      h_iq_to_i1(2,3,nbb),
     >                      h_iq_to_i1(1,3,nbb))
      value = value.and.
     >        BA_alloc_get(mt_int,(ny(id)*nq2(id)),
     >                     'h_iq_to_i2_3',
     >                      h_iq_to_i2(2,3,nbb),
     >                      h_iq_to_i2(1,3,nbb))
      value = value.and.
     >        BA_alloc_get(mt_int,(ny(id)*nq2(id)),
     >                     'h_iz_to_i2_3',
     >                      h_iz_to_i2(2,3,nbb),
     >                      h_iz_to_i2(1,3,nbb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i1_start_3',
     >                      h_i1_start(2,3,nbb),
     >                      h_i1_start(1,3,nbb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i2_start_3',
     >                      h_i2_start(2,3,nbb),
     >                      h_i2_start(1,3,nbb))
      if (.not.value)
     > call errquit('D3dB_transpose_ijk_init:out of heap',0,MA_ERR)

      index1 = 1
      index2 = 1
      index3 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
         int_mb(h_i1_start(1,3,nbb)+it) = index1
         int_mb(h_i2_start(1,3,nbb)+it) = index2

         do k=1,nz(id)
         do j=1,ny(id)
         do i=1,(nx(id)/2+1)

            iszero = zero_arow3(i+(j-1)*(nx(id)/2+1))

*           **** packing scheme ****
            phere = int_mb(p_map3(1,id)+(i-1)+(j-1)*(nx(id)/2+1))
            qhere = int_mb(q_map3(1,id)+(i-1)+(j-1)*(nx(id)/2+1))

            pto   = int_mb(p_map2(1,id)+(k-1)+(i-1)*nz(id))
            qto   = int_mb(q_map2(1,id)+(k-1)+(i-1)*nz(id))


            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
               itmp = k + (qhere-1)*nz(id)
               if (.not.iszero) then
               int_mb(h_iq_to_i1(1,3,nbb)+index1-1) = itmp
               index1 = index1 + 1
               end if
            end if

*           **** unpacking scheme ****
            if ((pto.eq.taskid).and.(phere.eq.proc_from)) then
               itmp = j + (qto-1)*ny(id)
               if (.not.iszero) then
               int_mb(h_iq_to_i2(1,3,nbb)+index2-1) = itmp
               index2 = index2 + 1
               else
               int_mb(h_iz_to_i2(1,3,nbb)+index3-1) = itmp
               index3 = index3 + 1
               end if
            end if


         end do
         end do
         end do

      end do
      int_mb(h_i1_start(1,3,nbb)+np) = index1
      int_mb(h_i2_start(1,3,nbb)+np) = index2
      h_iz_to_i2_count(3,nbb) = index3 - 1


*     ********************************************************
*     **** map2to1 mapping - done - tranpose operation #4 ****
*     ****  (nx/2+1,ny,nz)  <-- (ny,nz,nx/2+1)            ****
*     ****  use zero_arow2                                ****
*     ********************************************************

*     **** allocate trans_blk_ijk common block ****
      value = BA_alloc_get(mt_int,(ny(id)*nq2(id)),
     >                     'h_iq_to_i1_4',
     >                      h_iq_to_i1(2,4,nbb),
     >                      h_iq_to_i1(1,4,nbb))
      value = value.and.
     >        BA_alloc_get(mt_int,((nx(id)/2+1)*nq1(id)),
     >                     'h_iq_to_i2_4',
     >                      h_iq_to_i2(2,4,nbb),
     >                      h_iq_to_i2(1,4,nbb))
      value = value.and.
     >        BA_alloc_get(mt_int,((nx(id)/2+1)*nq1(id)),
     >                     'h_iz_to_i2_4',
     >                      h_iz_to_i2(2,4,nbb),
     >                      h_iz_to_i2(1,4,nbb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i1_start_4',
     >                      h_i1_start(2,4,nbb),
     >                      h_i1_start(1,4,nbb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i2_start_4',
     >                      h_i2_start(2,4,nbb),
     >                      h_i2_start(1,4,nbb))
      if (.not.value)
     > call errquit('D3dB_transpose_ijk_init:out of heap',0,MA_ERR)

      index1 = 1
      index2 = 1
      index3 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
         int_mb(h_i1_start(1,4,nbb)+it) = index1
         int_mb(h_i2_start(1,4,nbb)+it) = index2

         do k=1,nz(id)
         do j=1,ny(id)
         do i=1,(nx(id)/2+1)

            iszero = zero_arow2(i+(k-1)*(nx(id)/2+1))

*           **** packing scheme ****
            phere = int_mb(p_map2(1,id)+(k-1)+(i-1)*nz(id))
            qhere = int_mb(q_map2(1,id)+(k-1)+(i-1)*nz(id))

            pto   = int_mb(p_map1(1,id)+(j-1)+(k-1)*ny(id))
            qto   = int_mb(q_map1(1,id)+(j-1)+(k-1)*ny(id))


            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
               itmp = j + (qhere-1)*ny(id)
               if (.not.iszero) then
               int_mb(h_iq_to_i1(1,4,nbb)+index1-1) = itmp
               index1 = index1 + 1
               end if
            end if

*           **** unpacking scheme ****
            if ((pto.eq.taskid).and.(phere.eq.proc_from)) then
               itmp = i + (qto-1)*(nx(id)/2+1)
               if (.not.iszero) then
               int_mb(h_iq_to_i2(1,4,nbb)+index2-1) = itmp
               index2 = index2 + 1
               else
               int_mb(h_iz_to_i2(1,4,nbb)+index3-1) = itmp
               index3 = index3 + 1
               end if
            end if


         end do
         end do
         end do

      end do
      int_mb(h_i1_start(1,4,nbb)+np) = index1
      int_mb(h_i2_start(1,4,nbb)+np) = index2
      h_iz_to_i2_count(4,nbb) = index3 - 1



*     **********************************************************
*     **** map1to3 mapping  - done - tranpose operation # 5 ****
*     **********************************************************

*     **** allocate trans_blk_ijk common block ****
      value = BA_alloc_get(mt_int,((nx(id)/2+1)*nq1(id)),
     >                     'h_iq_to_i1_5',
     >                      h_iq_to_i1(2,5,nbb),
     >                      h_iq_to_i1(1,5,nbb))
      value = value.and.
     >        BA_alloc_get(mt_int,(nz(id)*nq3(id)),
     >                     'h_iq_to_i2_5',
     >                      h_iq_to_i2(2,5,nbb),
     >                      h_iq_to_i2(1,5,nbb))
      value = value.and.
     >        BA_alloc_get(mt_int,(nz(id)*nq3(id)),
     >                     'h_iz_to_i2_5',
     >                      h_iz_to_i2(2,5,nbb),
     >                      h_iz_to_i2(1,5,nbb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i1_start_5',
     >                      h_i1_start(2,5,nbb),
     >                      h_i1_start(1,5,nbb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i2_start_5',
     >                      h_i2_start(2,5,nbb),
     >                      h_i2_start(1,5,nbb))
      if (.not.value)
     > call errquit('D3dB_transpose_ijk_init:out of heap',0,MA_ERR)

      index1 = 1
      index2 = 1
      index3 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
         int_mb(h_i1_start(1,5,nbb)+it) = index1
         int_mb(h_i2_start(1,5,nbb)+it) = index2

         do k=1,nz(id)
         do j=1,ny(id)
         do i=1,(nx(id)/2+1)

*           **** packing scheme ****
            phere = int_mb(p_map1(1,id)+(j-1)+(k-1)*ny(id))
            qhere = int_mb(q_map1(1,id)+(j-1)+(k-1)*ny(id))

            pto   = int_mb(p_map3(1,id)+(i-1)+(j-1)*(nx(id)/2+1))
            qto   = int_mb(q_map3(1,id)+(i-1)+(j-1)*(nx(id)/2+1))


            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
               itmp = i + (qhere-1)*(nx(id)/2+1)
               int_mb(h_iq_to_i1(1,5,nbb)+index1-1) = itmp
               index1 = index1 + 1
            end if

*           **** unpacking scheme ****
            if ((pto.eq.taskid).and.(phere.eq.proc_from)) then
               itmp = k + (qto-1)*nz(id)
               int_mb(h_iq_to_i2(1,5,nbb)+index2-1) = itmp
               index2 = index2 + 1
            end if

         end do
         end do
         end do

      end do
      int_mb(h_i1_start(1,5,nbb)+np) = index1
      int_mb(h_i2_start(1,5,nbb)+np) = index2
      h_iz_to_i2_count(5,nbb) = index3 - 1



*     *************************
*     **** map3to1 mapping ****
*     *************************

*     **** allocate trans_blk_ijk common block ****
      value = BA_alloc_get(mt_int,(nz(id)*nq3(id)),
     >                     'h_iq_to_i1_6',
     >                      h_iq_to_i1(2,6,nbb),
     >                      h_iq_to_i1(1,6,nbb))
      value = value.and.
     >        BA_alloc_get(mt_int,((nx(id)/2+1)*nq1(id)),
     >                     'h_iq_to_i2_6',
     >                      h_iq_to_i2(2,6,nbb),
     >                      h_iq_to_i2(1,6,nbb))
      value = value.and.
     >        BA_alloc_get(mt_int,((nx(id)/2+1)*nq1(id)),
     >                     'h_iz_to_i2_6',
     >                      h_iz_to_i2(2,6,nbb),
     >                      h_iz_to_i2(1,6,nbb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i1_start_6',
     >                      h_i1_start(2,6,nbb),
     >                      h_i1_start(1,6,nbb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i2_start_6',
     >                      h_i2_start(2,6,nbb),
     >                      h_i2_start(1,6,nbb))
      if (.not.value)
     > call errquit('D3dB_transpose_ijk_init:out of heap',0,MA_ERR)

      index1 = 1
      index2 = 1
      index3 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
         int_mb(h_i1_start(1,6,nbb)+it) = index1
         int_mb(h_i2_start(1,6,nbb)+it) = index2

         do k=1,nz(id)
         do j=1,ny(id)
         do i=1,(nx(id)/2+1)

*           **** packing scheme ****
            phere = int_mb(p_map3(1,id)+(i-1)+(j-1)*(nx(id)/2+1))
            qhere = int_mb(q_map3(1,id)+(i-1)+(j-1)*(nx(id)/2+1))

            pto   = int_mb(p_map1(1,id)+(j-1)+(k-1)*ny(id))
            qto   = int_mb(q_map1(1,id)+(j-1)+(k-1)*ny(id))


            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
               itmp = k + (qhere-1)*nz(id)
               int_mb(h_iq_to_i1(1,6,nbb)+index1-1) = itmp
               index1 = index1 + 1
            end if

*           **** unpacking scheme ****
            if ((pto.eq.taskid).and.(phere.eq.proc_from)) then
               itmp = i + (qto-1)*(nx(id)/2+1)
               int_mb(h_iq_to_i2(1,6,nbb)+index2-1) = itmp
               index2 = index2 + 1
            end if

         end do
         end do
         end do

      end do
      int_mb(h_i1_start(1,6,nbb)+np) = index1
      int_mb(h_i2_start(1,6,nbb)+np) = index2
      h_iz_to_i2_count(6,nbb) = index3 - 1

      return
      end






*     ***********************************
*     *					*
*     *	        D3dB_cr_mpfft3b		*
*     *					*
*     ***********************************

      subroutine D3dB_cr_mpfft3b(nb,nbb,m,A)

*****************************************************
*                                                   *
*      This routine performs the operation of       *
*      a three dimensional complex to complex       *
*      inverse fft                                  *
*           A(nx,ny(nb),nz(nb)) <- FFT3^(-1)[A(kx,ky,kz)]   * 
*                                                   *
*      Entry - 					    *
*              A: a column distribuded 3d block     *
*              tmp: tempory work space must be at   *
*                    least the size of (complex)    *
*                    (nfft*nfft + 1) + 10*nfft      * 
*                                                   *
*       Exit - A is transformed and the imaginary   *
*              part of A is set to zero             *
*       uses - D3dB_c_ptranspose_jk, dcopy           *
*                                                   *
*****************************************************

      implicit none
      integer nb,nbb,m
      complex*16  A(*)

#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"

      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / D3dB_fft / tmpx,tmpy,tmpz



*     *** local variables ***
      integer i,j,k,q,indx,ierr,ii,shift,shift0
      integer nxh,nxhy,nxhz,indx0,indx1,indx2

      
      !integer tmp1(2),tmp2(2),tmp3(2)
      integer tmp2(2),tmp3(2)
      logical value


      call nwpw_timing_start(1)

*     ***** allocate temporary space ****
      !call D3dB_nfft3d(nb,nfft3d)
      value = BA_push_get(mt_dcpl,(nfft3d(nb)),'ffttmp2',
     >                    tmp2(2),tmp2(1))
      value = value.and.
     >      BA_push_get(mt_dbl,(n2ft3d(nb)),'ffttmp3',tmp3(2),tmp3(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

       nxh = (nx(nb)/2+1)
       nxhz = nxh*nz(nb)
       nxhy = nxh*ny(nb)
      
      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then

      shift = 1
      do ii=1,m
      shift0=shift-1
*     *************************************************
*     ***     do fft along kz dimension             ***
*     ***   A(kx,nz(nb),ky) <- fft1d^(-1)[A(kx,kz,ky)]  ***
*     *************************************************
#ifdef MLIB
      indx0 = 0
      indx2 = 0
      do q=1,nq(nb)
      do i=1,nxh
       indx2 = indx2 + 1
       if (.not.log_mb(zero_row3(1,nbb)+indx2-1)) then

         indx = i+indx0
         indx1 = indx
         do k=1,nz(nb)
            dcpl_mb(tmp2(1)+k-1) = A(indx+shift0)
            indx = indx + nxh
         end do
         call z1dfft(dcpl_mb(tmp2(1)),nz(nb),
     >               dcpl_mb(tmpz(1,nb)),-2,ierr)
         do k=1,nz(nb)
            A(indx1+shift0) = dcpl_mb(tmp2(1)+k-1)
            indx1 = indx1 + nxh
         end do
       end if
      end do
      indx0 = indx0+nxhz
      end do
      !call dscal((nx(nb)+2)*ny(nb)*nq(nb),dble(nz(nb)),A(shift),1)

#else


      !call dcffti(nz(nb),dcpl_mb(tmp1(1)))
      indx0 = 0
      indx2 = 0
      do q=1,nq(nb)
      do i=1,nxh
       indx2 = indx2 + 1
       if (.not.log_mb(zero_row3(1,nbb)+indx2-1)) then

#ifdef FFTW3
         indx = i+indx0+shift0
         call dfftw_execute_dft(plans(7,nb),A(indx),A(indx))

#else
         indx  = i + indx0
         indx1 = indx
         do k=1,nz(nb)
            dcpl_mb(tmp2(1)+k-1) = A(shift0+indx)
            indx = indx + nxh
         end do
         call dcfftb(nz(nb),dcpl_mb(tmp2(1)),dcpl_mb(tmpz(1,nb)))
         do k=1,nz(nb)
            A(shift0+indx1) = dcpl_mb(tmp2(1)+k-1)
            indx1 = indx1 + nxh
         end do
#endif

       end if
      end do
      indx0 = indx0 + nxhz
      end do
#endif


*     ********************************************
*     ***         Do a ptranspose of A          ***
*     ***      A(kx,ky,nz(nb)) <- A(kx,nz(nb),ky)      ***
*     ********************************************
      call D3dB_c_ptranspose1_jk(nbb,A(shift),
     >                           dcpl_mb(tmp2(1)),
     >                           dbl_mb(tmp3(1)))
      shift = shift + nfft3d(nb)
      end do

      shift = 1
      do ii=1,m
      shift0 = shift-1
*     *************************************************
*     ***     do fft along ky dimension             ***
*     ***   A(kx,ny(nb),nz(nb)) <- fft1d^(-1)[A(kx,ky,nz(nb))]  ***
*     *************************************************
#ifdef MLIB
     
      indx0 = 0
      indx2 = 0
      do q=1,nq(nb)
      do i=1,nxh
        indx2 = indx2 + 1
        if (.not.log_mb(zero_row2(1,nbb)+indx2-1)) then

          indx  = i + indx0
          indx1 = indx
          do j=1,ny(nb)
            dcpl_mb(tmp2(1)+j-1) = A(shift0+indx)
            indx = indx + nxh
          end do
          call z1dfft(dcpl_mb(tmp2(1)),ny(nb),
     >               dcpl_mb(tmpy(1,nb)),-2,ierr)
          do j=1,ny(nb)
            A(shift0+indx1) = dcpl_mb(tmp2(1)+j-1)
            indx1 = indx1 + nxh
          end do
        end if
      end do
      indx0 = indx0 + nxhy
      end do
      
#else
    
      indx0 = 0
      indx2 = 0
      do q=1,nq(nb)
      do i=1,nxh
        indx2 = indx2 + 1

        if (.not.log_mb(zero_row2(1,nbb)+indx2-1)) then

#ifdef FFTW3
         indx = i+indx0+shift0
         call dfftw_execute_dft(plans(8,nb),A(indx),A(indx))
#else
         indx  = i + indx0
         indx1 = indx
         do j=1,ny(nb)
            dcpl_mb(tmp2(1)+j-1) = A(shift0+indx)
            indx = indx + nxh
         end do
         call dcfftb(ny(nb),dcpl_mb(tmp2(1)),dcpl_mb(tmpy(1,nb)))
         do j=1,ny(nb)
            A(shift0+indx1) = dcpl_mb(tmp2(1)+j-1)
            indx1 = indx1 + nxh
         end do
#endif

        end if
      end do
      indx0 = indx0 + nxhy
      end do
#endif

*     *************************************************
*     ***     do fft along kx dimension             ***
*     ***   A(nx(nb),ny(nb),nz(nb)) <- fft1d^(-1)[A(kx,ny(nb),nz(nb))]  ***
*     *************************************************
#ifdef MLIB
      !call drc1ft (dbl_mb(tmp3(1)),nx(nb),dcpl_mb(tmp1(1)),-3,ierr)
      indx = 1
      do q=1,nq(nb)
      do j=1,ny(nb)
         call drc1ft(A(shift0+indx),nx(nb),dcpl_mb(tmpx(1,nb)),-2,ierr)
         indx = indx + nxh
      end do
      end do

#else

#ifdef FFTW3
      call dfftw_execute_dft_c2r(plans(3,nb),A(shift0),A(shift0))

#else
      call cshift1_fftb(nx(nb),ny(nb),nq(nb),1,A(shift))
      indx = 1
      do q=1,nq(nb)
      do j=1,ny(nb)
         !indx = 1 + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
         call drfftb(nx(nb),A(shift0+indx),dcpl_mb(tmpx(1,nb)))
         indx = indx + nxh
      end do
      end do
      call zeroend_fftb(nx(nb),ny(nb),nq(nb),1,A(shift))
#endif

#endif

      shift = shift + nfft3d(nb)
      end do

      !*************************
      !**** hilbert mapping ****
      !*************************
      else

      shift = 1
      do ii=1,m
      shift0=shift-1

*     *************************************************
*     ***     do fft along kz dimension             ***
*     ***   A(nz(nb),kx,ky) <- fft1d^(-1)[A(kz,kx,ky)]  ***
*     *************************************************
#ifdef MLIB
      indx = 1
      do q=1,nq3(nb)
         if (.not.log_mb(zero_row3(1,nbb)+q-1)) then
         call z1dfft(A(shift0+indx),nz(nb),dcpl_mb(tmpz(1,nb)),-2,ierr)
         end if
         indx = indx + nz(nb)
      end do
#else
      indx = 1
      do q=1,nq3(nb)
         if (.not.log_mb(zero_row3(1,nbb)+q-1)) then

#ifdef FFTW3
         call dfftw_execute_dft(plans(17,nb),
     >                          A(indx+shift0),A(indx+shift0))
#else
         call dcfftb(nz(nb),A(shift0+indx),dcpl_mb(tmpz(1,nb)))
#endif

         end if
         indx = indx + nz(nb)
      end do
#endif

      call D3dB_c_ptranspose_ijk(nbb,3,A(shift),
     >                           dcpl_mb(tmp2(1)),
     >                           dbl_mb(tmp3(1)))

*     *************************************************
*     ***     do fft along ky dimension             ***
*     ***   A(ny(nb),nz(nb),kx) <- fft1d^(-1)[A(ky,nz(nb),kx)]  ***
*     *************************************************
#ifdef MLIB
      indx = 1
      do q=1,nq2(nb)
         if (.not.log_mb(zero_row2(1,nbb)+q-1)) then
         call z1dfft(A(shift0+indx),ny(nb),dcpl_mb(tmpy(1,nb)),-2,ierr)
         end if
         indx = indx + ny(nb)
      end do
#else
      indx = 1
      do q=1,nq2(nb)
         if (.not.log_mb(zero_row2(1,nbb)+q-1)) then

#ifdef FFTW3
         call dfftw_execute_dft(plans(18,nb),
     >                          A(indx+shift0),A(indx+shift0))
#else
         call dcfftb(ny(nb),A(shift0+indx),dcpl_mb(tmpy(1,nb)))
#endif

         end if
         indx = indx + ny(nb)
      end do
#endif

      call D3dB_c_ptranspose_ijk(nbb,4,A(shift),
     >                           dcpl_mb(tmp2(1)),
     >                           dbl_mb(tmp3(1)))

*     *************************************************
*     ***     do fft along kx dimension             ***
*     ***   A(nx(nb),ny(nb),nz(nb)) <- fft1d^(-1)[A(kx,ny(nb),nz(nb))]  ***
*     *************************************************
#ifdef MLIB
      indx = 1
      do q=1,nq1(nb)
         !indx = 1 + (q-1)*(nx(nb)/2+1)
         call drc1ft(A(shift0+indx),nx(nb),dcpl_mb(tmpx(1,nb)),-2,ierr)
         indx = indx + nxh
      end do
#else

#ifdef FFTW3
      call dfftw_execute_dft_c2r(plans(13,nb),A(shift0),A(shift0))

#else
      call cshift1_fftb(nx(nb),nq1(nb),1,1,A(shift))
      indx = 1
      do q=1,nq1(nb)
         !indx = 1 + (q-1)*(nx(nb)/2+1)
         call drfftb(nx(nb),A(shift0+indx),dcpl_mb(tmpx(1,nb)))
         indx = indx + nxh
      end do
      call zeroend_fftb(nx(nb),nq1(nb),1,1,A(shift))
#endif

#endif

      shift = shift + nfft3d(nb)
      end do


      end if
    
*     **** deallocate temporary space  ****
      value = BA_pop_stack(tmp3(2))
      value = BA_pop_stack(tmp2(2))
      !value = BA_pop_stack(tmp1(2))

      call nwpw_timing_end(1)
      return
      end






*     ***********************************
*     *                                 *
*     *         D3dB_pfftfx             *
*     *                                 *
*     ***********************************
*
*     do fft along nx(id) dimension       
*
*        A(kx,ny(id),nz(id)) <- fft1d[A(nx(id),ny(id),nz(id))] 
*     

      subroutine D3dB_pfftfx(nbb,A,tmp1,tmp2,request,reqcnt)
      implicit none
      integer nbb
      complex*16 A(*)
      complex*16 tmp1(*)
      complex*16 tmp2(*)
      integer    request(*),reqcnt


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / D3dB_fft / tmpx,tmpy,tmpz

      !*** local variables ***
      integer j,q,indx,nxh,id
      integer tid,nthr,offset
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

      if (nbb.lt.2) then
         id = 1
      else
         id = 3
      end if

      nxh = (nx(id)/2+1)

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
      call nwpw_timing_start(31)
#ifdef MLIB
      call drcfts(A,nx(id),1,ny(id)*nq(id),
     >                  nx(id)+2,1,ierr)
#else

#ifdef FFTW3
      call dfftw_execute_dft_r2c(plans(4,id),A,A)

#else

#ifdef NUSE_OPENMP
      call D3dB_pfftfx_sub(ny(id)*nq(id),nx(id),nxh,
     >                    dcpl_mb(tmpx(1,id)),A)
      call cshift_fftf(nx(id),ny(id),nq(id),1,A)

#else
      indx = 1
      do q=1,nq(id)
      do j=1,ny(id)
         call drfftf(nx(id),A(indx),dcpl_mb(tmpx(1,id)))
         indx = indx + nxh
      end do
      end do
      call cshift_fftf(nx(id),ny(id),nq(id),1,A)
#endif
#endif
#endif
      call dcopy((nx(id)+2)*ny(id)*nq(id),A,1,tmp1,1)
      call nwpw_timing_end(31)


      !*************************
      !**** hilbert mapping ****
      !*************************
      else
      call nwpw_timing_start(31)

#ifdef MLIB
      call drcfts(A,nx(id),1,nq1(id),
     >                  nx(id)+2,1,ierr)
#else

#ifdef FFTW3
      call dfftw_execute_dft_r2c(plans(14,id),A,A)

#else

      !indx = 1
      offset = tid*(2*nx(id)+15)
      do q=tid+1,nq1(id),nthr
         indx = 1+(q-1)*nxh
         call drfftf(nx(id),A(indx),dcpl_mb(tmpx(1,id)+offset))
         !indx = indx + nxh
      end do
      call cshift_fftf(nx(id),nq1(id),1,1,A)

#endif
#endif

      call nwpw_timing_end(31)
      call nwpw_timing_start(32)
      call D3dB_c_ptranspose_ijk_start(nbb,1,A,tmp1,tmp2,
     >                                 request,reqcnt,40)
      call nwpw_timing_end(32)

      end if

      return
      end

      subroutine  D3dB_pfftfx_sub(n,nx,nxh,tmpx,A)
      implicit none
      integer n,nx,nxh
      real*8     tmpx(2*nx+15)
      !complex*16 A(nxh,n)
      real*8 A(2*nxh,n)
      integer i,indx

      do i=1,n
         call drfftf(nx,A(1,i),tmpx)
      end do
      return
      end





*     ***********************************
*     *                                 *
*     *         D3dB_pfftfx0            *
*     *                                 *
*     ***********************************
*
*     do fft along nx(id) dimension       
*
*        A(kx,ny(id),nz(id)) <- fft1d[A(nx(id),ny(id),nz(id))] 
*     

      subroutine D3dB_pfftfx0(nbb,tmp1,tmp2,request,reqcnt)
      implicit none
      integer nbb
      complex*16 tmp1(*)
      complex*16 tmp2(*)
      integer    request(*),reqcnt


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / D3dB_fft / tmpx,tmpy,tmpz

      !*** local variables ***
      integer j,q,indx,nxh,id
      integer tid,nthr,offset
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()
      if (nbb.lt.2) then
         id = 1
      else
         id = 3
      end if

      nxh = (nx(id)/2+1)

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
      call nwpw_timing_start(31)
#ifdef MLIB
      call drcfts(tmp1,nx(id),1,ny(id)*nq(id),
     >                  nx(id)+2,1,ierr)
#else

#ifdef FFTW3
      call dfftw_execute_dft_r2c(plans(4,id),tmp1,tmp1)

#else

#ifdef NUSE_OPENMP
      call D3dB_pfftfx_sub(ny(id)*nq(id),nx(id),nxh,
     >                    dcpl_mb(tmpx(1,id)),tmp1)
      call cshift_fftf(nx(id),ny(id),nq(id),1,tmp1)

#else
      indx = 1
      do q=1,nq(id)
      do j=1,ny(id)
         call drfftf(nx(id),tmp1(indx),dcpl_mb(tmpx(1,id)))
         indx = indx + nxh
      end do
      end do
      call cshift_fftf(nx(id),ny(id),nq(id),1,tmp1)
#endif
#endif
#endif
      call nwpw_timing_end(31)

      !*************************
      !**** hilbert mapping ****
      !*************************
      else
      call nwpw_timing_start(31)

#ifdef MLIB
      call drcfts(tmp1,nx(id),1,nq1(id),
     >                  nx(id)+2,1,ierr)
#else

#ifdef FFTW3
      call dfftw_execute_dft_r2c(plans(14,id),tmp1,tmp1)

#else

      !indx = 1
      offset = tid*(2*nx(id)+15)
      do q=tid+1,nq1(id),nthr
         indx = 1+(q-1)*nxh
         call drfftf(nx(id),tmp1(indx),dcpl_mb(tmpx(1,id)+offset))
         !indx = indx + nxh
      end do
      call cshift_fftf(nx(id),nq1(id),1,1,tmp1)

#endif
#endif

      call nwpw_timing_end(31)
      call nwpw_timing_start(32)
      call D3dB_c_ptranspose_ijk_start(nbb,1,tmp1,tmp2,tmp1,
     >                                 request,reqcnt,40)
      call nwpw_timing_end(32)

      end if

      return
      end








*     ***********************************
*     *                                 *
*     *         D3dB_pfftfy             *
*     *                                 *
*     ***********************************
*
*     do fft along ny(id) dimension       
*
*        A(kx,ny(id),nz(id)) <- fft1d[A(nx(id),ny(id),nz(id))] 
*     

      subroutine D3dB_pfftfy(nbb,tmp1,tmp2,request,reqcnt)
      implicit none
      integer nbb
      complex*16 tmp1(*)
      complex*16 tmp2(*)
      integer    request(*),reqcnt


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / D3dB_fft / tmpx,tmpy,tmpz


      !*** local variables ***
      integer i,j,k,indx,indx0,indx1,indx2,q,nxh,nxhy,id
      integer tid,nthr,offset
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      if (nbb.lt.2) then
         id = 1
      else
         id = 3
      end if
      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

      nxh = (nx(id)/2+1)
      nxhy = nxh*ny(id)


      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
      call nwpw_timing_start(31)

*     ********************************************
*     ***     do fft along ny(id) dimension        ***
*     ***   A(kx,ky,nz(id)) <- fft1d[A(kx,ny(id),nz(id))]  ***
*     ********************************************

#ifdef MLIB
      indx0 = 0
      indx2 = 0
      do q=1,nq(id)
      do i=1,nxh
        indx2 = indx2 + 1
        if (.not.log_mb(zero_row2(1,nbb)+indx2-1)) then

         indx  = i + indx0
         indx1 = indx
         do j=1,ny(id)
            tmp2(j) = tmp1(indx)
            indx = indx + nxh
         end do
         call z1dfft(tmp2,ny(id),dcpl_mb(tmpy(1,id)),1,ierr)
         do j=1,ny(id)
            tmp1(indx1) = tmp2(j)
            indx1 = indx1 + nxh
         end do

        end if
      end do
      indx0 = indx0 + nxhy
      end do
#else

#ifdef FFTW3
       indx0 = 0
       indx2 = 0
       do q=1,nq(id)
       do i=1,nxh
        indx2 = indx2 + 1
        if (.not.log_mb(zero_row2(1,nbb)+indx2-1)) then
          indx  = i + indx0
          call dfftw_execute_dft(plans(9,1),tmp1(indx),tmp1(indx))
        end if
       end do
       indx0 = indx0 + nxhy
       end do

#else
       indx0 = 0
       indx2 = 0
       do q=1,nq(id)
       do i=1,nxh
        indx2 = indx2 + 1
        if (.not.log_mb(zero_row2(1,nbb)+indx2-1)) then
          indx  = i + indx0
          indx1 = indx
          do k=1,ny(id)
             tmp2(k) = tmp1(indx)
             indx = indx + nxh
          end do
          call dcfftf(ny(id),tmp2,dcpl_mb(tmpy(1,id)))
          do k=1,ny(id)
             tmp1(indx1) = tmp2(k)
             indx1 = indx1 + nxh
          end do
        end if
       end do
       indx0 = indx0 + nxhy
       end do
#endif
#endif


*     ********************************************
*     ***         Do a ptranspose of A          ***
*     ***      A(ky,nz(id),ky) <- A(kx,ky,nz(id))      ***
*     ********************************************
      call nwpw_timing_end(31)
      call nwpw_timing_start(32)
      call D3dB_c_ptranspose2_jk_start(nbb,tmp1,tmp2,tmp1,
     >                                 request,reqcnt,41)
      call nwpw_timing_end(32)


      !*************************
      !**** hilbert mapping ****
      !*************************
      else
      call nwpw_timing_start(32)
      call D3dB_c_ptranspose_ijk_end(nbb,1,tmp1,tmp2,request,reqcnt)
      call nwpw_timing_end(32)
      call nwpw_timing_start(31)

*     ********************************************
*     ***     do fft along ny(id) dimension        ***
*     ***   A(ky,nz(id),kx) <- fft1d[A(ny(id),nz(id),kx)]  ***
*     ********************************************
#ifdef MLIB
      indx = 1
      do q=1,nq2(id)
         if (.not.log_mb(zero_row2(1,nbb)+q-1)) then
         call z1dfft(tmp1(indx),ny(id),dcpl_mb(tmpy(1,id)),1,ierr)
         end if
         indx = indx + ny(id)
      end do
#else
#ifdef FFTW3
      indx = 1
      do q=1,nq2(id)
         if (.not.log_mb(zero_row2(1,nbb)+q-1)) then
         call dfftw_execute_dft(plans(19,1),tmp1(indx),tmp1(indx))
         end if
         indx = indx + ny(id)
      end do
#else

      !indx = 1
      offset = tid*(2*ny(id)+15)
      do q=tid+1,nq2(id),nthr
         indx = 1+(q-1)*ny(id)
         if (.not.log_mb(zero_row2(1,nbb)+q-1)) then
            call dcfftf(ny(id),tmp1(indx),dcpl_mb(tmpy(1,id)+offset))
         end if
         !indx = indx + ny(id)
      end do
!$OMP BARRIER

#endif
#endif

      call nwpw_timing_end(31)
      call nwpw_timing_start(32)
      call D3dB_c_ptranspose_ijk_start(nbb,2,tmp1,tmp2,tmp1,
     >                                 request,reqcnt,42)
      call nwpw_timing_end(32)


      end if

      return
      end

      subroutine  D3dB_pfftfy_sub2(n,ny,row2,tmpy,A)
      implicit none
      integer n,ny
      logical    row2(n)
      real*8     tmpy(4*ny+15)
      !complex*16 A(ny,n)
      real*8     A(2*ny,n)
      integer i
      do i=1,n
         if (.not.row2(i)) then
            call dcfftf(ny,A(1,i),tmpy)
         end if
      end do
      return
      end







*     ***********************************
*     *                                 *
*     *         D3dB_pfftfz             *
*     *                                 *
*     ***********************************
*

      subroutine D3dB_pfftfz(nbb,tmp1,tmp2,request,reqcnt)
      implicit none
      integer nbb
      complex*16 tmp1(*)
      complex*16 tmp2(*)
      integer    request(*),reqcnt


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / D3dB_fft / tmpx,tmpy,tmpz


      !*** local variables ***
      integer i,k,q,nxh,nxhz,id
      integer indx,indx0,indx1,indx2
      integer tid,nthr,offset
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      if (nbb.lt.2) then
         id = 1
      else
         id = 3
      end if
      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

      nxh = (nx(id)/2+1)
      nxhz = nxh*nz(id)



      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then

      call nwpw_timing_start(32)
      call D3dB_c_ptranspose2_jk_end(nbb,tmp2,tmp1,request,reqcnt)
      call nwpw_timing_end(32)
      call nwpw_timing_start(31)

*     ********************************************
*     ***     do fft along nz(id) dimension        ***
*     ***   A(kx,kz,ky) <- fft1d[A(kx,nz(id),ky)]  ***
*     ********************************************

#ifdef MLIB

      indx0 = 0
      indx2 = 0
      do q=1,nq(id)
      do i=1,nxh
        indx2 = indx2 + 1
        if (.not.log_mb(zero_row3(1,nbb)+indx2-1)) then

         indx  = i + indx0
         indx1 = indx
         do k=1,nz(id)
            tmp1(k) = tmp2(indx)
            indx = indx + nxh
         end do
         call z1dfft(tmp1,nz(id),dcpl_mb(tmpz(1,id)),1,ierr)
         do k=1,nz(id)
            tmp2(indx1) = tmp1(k)
            indx1 = indx1 + nxh
         end do
        end if
      end do
      indx0 = indx0 + nxhz
      end do
#else

       indx0 = 0
       indx2 = 0
       do q=1,nq(id)
       do i=1,nxh
        indx2 = indx2 + 1
        if (.not.log_mb(zero_row3(1,nbb)+indx2-1)) then

#ifdef FFTW3
          indx  = i + indx0
          call dfftw_execute_dft(plans(10,1),tmp2(indx),tmp2(indx))

#else
          indx  = i + indx0
          indx1 = indx
          do k=1,nz(id)
             tmp1(k) = tmp2(indx)
             indx = indx + nxh
          end do
          call dcfftf(nz(id),tmp1,dcpl_mb(tmpz(1,id)))
          do k=1,nz(id)
             tmp2(indx1) = tmp1(k)
             indx1 = indx1 + nxh
          end do
#endif

        end if
       end do
       indx0 = indx0 + nxhz
       end do
      call nwpw_timing_end(31)

#endif


      !*************************
      !**** hilbert mapping ****
      !*************************
      else
      call nwpw_timing_start(32)
      call D3dB_c_ptranspose_ijk_end(nbb,2,tmp2,tmp1,request,reqcnt)
      call nwpw_timing_end(32)
      call nwpw_timing_start(31)

*     ********************************************
*     ***     do fft along nz(id) dimension        ***
*     ***   A(kz,kx,ky) <- fft1d[A(nz(id),kx,ky)]  ***
*     ********************************************
#ifdef MLIB
      indx = 1
      do q=1,nq3(id)
         if (.not.log_mb(zero_row3(1,nbb)+q-1)) then
         call z1dfft(tmp2(indx),nz(id),dcpl_mb(tmpz(1,id)),1,ierr)
         end if
         indx = indx + nz(id)
      end do
#else
#ifdef FFTW3
      indx = 1
      do q=1,nq3(id)
         if (.not.log_mb(zero_row3(1,nbb)+q-1)) then
            call dfftw_execute_dft(plans(20,1),tmp2(indx),tmp2(indx))
         end if
         indx = indx + nz(id)
      end do
#else

      !indx = 1
      offset = tid*(2*nz(id)+15)
      do q=tid+1,nq3(id),nthr
         indx = 1+(q-1)*nz(id)
         if (.not.log_mb(zero_row3(1,nbb)+q-1)) then
            call dcfftf(nz(id),tmp2(indx),dcpl_mb(tmpz(1,id)+offset))
         end if
         !indx = indx + nz(id)
      end do
!$OMP BARRIER

#endif
#endif
      call nwpw_timing_end(31)

      end if

      call nwpw_timing_start(32)
      call Pack_c_pack_start(nbb,tmp2,tmp1,request,reqcnt,43)
      call nwpw_timing_end(32)

      return
      end

      subroutine  D3dB_pfftfz_sub2(n,nz,row3,tmpz,A)
      implicit none
      integer n,nz
      logical    row3(n)
      real*8     tmpz(4*nz+15)
      !complex*16 A(nz,n)
      real*8     A(2*nz,n)
      integer i
      do i=1,n
         if (.not.row3(i)) then
            call dcfftf(nz,A(1,i),tmpz)
         end if
      end do
      return
      end




*     ***********************************
*     *                                 *
*     *         D3dB_pfftf_final        *
*     *                                 *
*     ***********************************
*

      subroutine D3dB_pfftf_final(nbb,tmp1,tmp2,request,reqcnt)
      implicit none
      integer nbb
      complex*16 tmp1(*)
      complex*16 tmp2(*)
      integer    request(*),reqcnt


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      call nwpw_timing_start(32)
      call Pack_c_pack_end(nbb,tmp2,request,reqcnt)
      call nwpw_timing_end(32)

      return
      end




*     ***********************************
*     *                                 *
*     *         D3dB_pfftf_step         *
*     *                                 *
*     ***********************************
      subroutine D3dB_pfftf_step(step,nbb,A,tmp1,tmp2,request,reqcnt)
      implicit none
      integer step
      integer nbb
      complex*16 A(*)
      complex*16 tmp1(*)
      complex*16 tmp2(*)
      integer    request(*),reqcnt

     
      if (step.eq.1) then
         call D3dB_pfftfx(nbb,A,tmp1,tmp2,request,reqcnt)
      else if (step.eq.2) then
         call D3dB_pfftfy(nbb,tmp1,tmp2,request,reqcnt)
      else if (step.eq.3) then
         call D3dB_pfftfz(nbb,tmp1,tmp2,request,reqcnt)
      else if (step.eq.4) then
         call D3dB_pfftf_final(nbb,tmp1,tmp2,request,reqcnt)
      end if
!$OMP BARRIER

      return
      end


*     ***********************************
*     *                                 *
*     *      D3dB_pfftf_step_center     *
*     *                                 *
*     ***********************************
      subroutine D3dB_pfftf_step_center(center,step,nbb,A,tmp1,tmp2,
     >                                  request,reqcnt)
      implicit none
      integer center(3)
      integer step
      integer nbb
      complex*16 A(*)
      complex*16 tmp1(*)
      complex*16 tmp2(*)
      integer    request(*),reqcnt

      if (step.eq.1) then
         call D3dB_r_1to3_Contract_start3(center,A,tmp2,request,reqcnt)
      else if (step.eq.2) then
         call D3dB_r_1to3_Contract_mid3(tmp2,tmp1,request,reqcnt)
      else if (step.eq.3) then
         call D3dB_r_1to3_Contract_final3(tmp1,tmp2,request,reqcnt)

      else if (step.eq.4) then
         call D3dB_pfftfx0(nbb,tmp2,tmp1,request,reqcnt)
      else if (step.eq.5) then
         call D3dB_pfftfy(nbb,tmp1,tmp2,request,reqcnt)
      else if (step.eq.6) then
         call D3dB_pfftfz(nbb,tmp1,tmp2,request,reqcnt)
      else if (step.eq.7) then
         call D3dB_pfftf_final(nbb,tmp1,tmp2,request,reqcnt)
      end if

!$OMP BARRIER

      return
      end






*     ***********************************
*     *					*
*     *	        D3dB_pfftbz		*
*     *					*
*     ***********************************

*                                                  
*      This routine performs the operation of      
*      a three dimensional complex to complex     
*      inverse fft                                
*           A(nx,ny(id),nz(id)) <- FFT3^(-1)[A(kx,ky,kz)]  
*                                                  
*      Entry -                                     
*              A: a column distribuded 3d block    
*              tmp: tempory work space must be at  
*                    least the size of (complex)   
*                    (nfft*nfft + 1) + 10*nfft     
*                                               
*       Exit - A is transformed and the imaginary
*              part of A is set to zero          
*       uses - D3dB_c_ptranspose_jk, dcopy        
*                                                  


      subroutine D3dB_pfftbz(nbb,tmp1,tmp2,request,reqcnt)
      implicit none
      integer nbb
      complex*16  tmp1(*)
      complex*16  tmp2(*)
      integer     request(*),reqcnt
      

#include "bafdecls.fh"
#include "errquit.fh"
#include "D3dB.fh"
#include "D3dB_pfft.fh"

      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / D3dB_fft / tmpx,tmpy,tmpz



*     *** local variables ***
      integer i,k,q,indx,ierr,id
      integer nxh,nxhz,indx0,indx1,indx2,tid,nthr
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      if (nbb.lt.2) then
         id = 1
      else
         id = 3
      end if
      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

       nxh = (nx(id)/2+1)
       nxhz = nxh*nz(id)
      
      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then

      call nwpw_timing_start(31)
*     *************************************************
*     ***     do fft along kz dimension             ***
*     ***   A(kx,nz(id),ky) <- fft1d^(-1)[A(kx,kz,ky)]  ***
*     *************************************************
#ifdef MLIB
      indx0 = 0
      indx2 = 0
      do q=1,nq(id)
      do i=1,nxh
       indx2 = indx2 + 1
       if (.not.log_mb(zero_row3(1,nbb)+indx2-1)) then

         indx = i+indx0
         indx1 = indx
         do k=1,nz(id)
            tmp2(k) = tmp1(indx)
            indx = indx + nxh
         end do
         call z1dfft(tmp2,nz(id),
     >               dcpl_mb(tmpz(1,id)),-2,ierr)
         do k=1,nz(id)
            tmp1(indx1) = tmp2(k)
            indx1 = indx1 + nxh
         end do
       end if
      end do
      indx0 = indx0+nxhz
      end do

#else
      indx0 = 0
      indx2 = 0
      do q=1,nq(id)
      do i=1,nxh
       indx2 = indx2 + 1
       if (.not.log_mb(zero_row3(1,nbb)+indx2-1)) then

#ifdef FFTW3
         indx  = i + indx0
         call dfftw_execute_dft(plans(7,1),tmp1(indx),tmp1(indx))
#else

         indx  = i + indx0
         indx1 = indx
         do k=1,nz(id)
            tmp2(k) = tmp1(indx)
            indx = indx + nxh
         end do
         call dcfftb(nz(id),tmp2,dcpl_mb(tmpz(1,id)))
         do k=1,nz(id)
            tmp1(indx1) = tmp2(k)
            indx1 = indx1 + nxh
         end do
#endif

       end if
      end do
      indx0 = indx0 + nxhz
      end do
#endif

*     ********************************************
*     ***         Do a ptranspose of A          ***
*     ***      A(kx,ky,nz(id)) <- A(kx,nz(id),ky)      ***
*     ********************************************
      call nwpw_timing_end(31)
      call nwpw_timing_start(32)
      call D3dB_c_ptranspose1_jk_start(nbb,tmp1,tmp2,tmp1,
     >                                 request,reqcnt,44)
      call nwpw_timing_end(32)


      !*************************
      !**** hilbert mapping ****
      !*************************
      else


      call nwpw_timing_start(31)
*     *************************************************
*     ***     do fft along kz dimension             ***
*     ***   A(nz(id),kx,ky) <- fft1d^(-1)[A(kz,kx,ky)]  ***
*     *************************************************
#ifdef MLIB
      indx = 1
      do q=1,nq3(id)
         if (.not.log_mb(zero_row3(1,nbb)+q-1)) then
         call z1dfft(tmp1(indx),nz(id),dcpl_mb(tmpz(1,id)),-2,ierr)
         end if
         indx = indx + nz(id)
      end do
#else
#ifdef FFTW3
      indx = 1
      do q=1,nq3(id)
         if (.not.log_mb(zero_row3(1,nbb)+q-1)) then
            call dfftw_execute_dft(plans(17,1),tmp1(indx),tmp1(indx))
         end if
         indx = indx + nz(id)
      end do
#else


      !indx = 1
      do q=1+tid,nq3(id),nthr
         indx = 1 + (q-1)*nz(id)
         if (.not.log_mb(zero_row3(1,nbb)+q-1)) then
            call dcfftb(nz(id),tmp1(indx),
     >                  dcpl_mb(tmpz(1,id)+tid*(2*nz(id)+15)))
         end if
         !indx = indx + nz(id)
      end do
!$OMP BARRIER

#endif
#endif

      call nwpw_timing_end(31)
      call nwpw_timing_start(32)
      call D3dB_c_ptranspose_ijk_start(nbb,3,tmp1,tmp2,tmp1,
     >                                 request,reqcnt,45)
      call nwpw_timing_end(32)

      end if
      return
      end

      subroutine  D3dB_pfftbz_sub2(n,nz,row3,tmpz,A)
      implicit none
      integer n,nz
      logical    row3(n)
      real*8     tmpz(4*nz+15)
      !complex*16 A(nz,n)
      real*8     A(2*nz,n)
      integer i
      do i=1,n
         if (.not.row3(i)) then
            call dcfftb(nz,A(1,i),tmpz)
         end if
      end do
      return
      end





*     ***********************************
*     *					*
*     *	        D3dB_pfftby		*
*     *					*
*     ***********************************

      subroutine D3dB_pfftby(nbb,tmp1,tmp2,request,reqcnt)
      implicit none
      integer nbb
      complex*16  tmp1(*)
      complex*16  tmp2(*)
      integer     request(*),reqcnt

#include "bafdecls.fh"
#include "errquit.fh"
#include "D3dB.fh"
#include "D3dB_pfft.fh"

      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / D3dB_fft / tmpx,tmpy,tmpz


*     *** local variables ***
      integer i,j,q,indx,ierr,id
      integer nxh,nxhy,indx0,indx1,indx2,tid,nthr
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

       if (nbb.lt.2) then
          id = 1
       else
          id = 3
       end if
      
       tid  = Parallel_threadid()
       nthr = Parallel_nthreads()

       nxh = (nx(id)/2+1)
       nxhy = nxh*ny(id)
      
      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then


*     ********************************************
*     ***         Do a ptranspose of A          ***
*     ***      A(kx,ky,nz(id)) <- A(kx,nz(id),ky)      ***
*     ********************************************
      call nwpw_timing_start(32)
      call D3dB_c_ptranspose1_jk_end(nbb,tmp2,tmp1,request,reqcnt)
      call nwpw_timing_end(32)
      call nwpw_timing_start(31)

*     *************************************************
*     ***     do fft along ky dimension             ***
*     ***   A(kx,ny(id),nz(id)) <- fft1d^(-1)[A(kx,ky,nz(id))]  ***
*     *************************************************
#ifdef MLIB
     
      indx0 = 0
      indx2 = 0
      do q=1,nq(id)
      do i=1,nxh
        indx2 = indx2 + 1
        if (.not.log_mb(zero_row2(1,nbb)+indx2-1)) then

          indx  = i + indx0
          indx1 = indx
          do j=1,ny(id)
            tmp1(j) = tmp2(indx)
            indx = indx + nxh
          end do
          call z1dfft(tmp1,ny(id),
     >               dcpl_mb(tmpy(1,id)),-2,ierr)
          do j=1,ny(id)
            tmp2(indx1) = tmp1(j)
            indx1 = indx1 + nxh
          end do
        end if
      end do
      indx0 = indx0 + nxhy
      end do
      
#else
    
      indx0 = 0
      indx2 = 0
      do q=1,nq(id)
      do i=1,nxh
        indx2 = indx2 + 1

        if (.not.log_mb(zero_row2(1,nbb)+indx2-1)) then

#ifdef FFTW3
         indx  = i + indx0
         call dfftw_execute_dft(plans(8,id),tmp2(indx),tmp2(indx))

#else
         indx  = i + indx0
         indx1 = indx
         do j=1,ny(id)
            tmp1(j) = tmp2(indx)
            indx = indx + nxh
         end do
         call dcfftb(ny(id),tmp1,dcpl_mb(tmpy(1,id)))
         do j=1,ny(id)
            tmp2(indx1) = tmp1(j)
            indx1 = indx1 + nxh
         end do
#endif

        end if
      end do
      indx0 = indx0 + nxhy
      end do
#endif

      call nwpw_timing_end(31)


      !*************************
      !**** hilbert mapping ****
      !*************************
      else

      call nwpw_timing_start(32)
      call D3dB_c_ptranspose_ijk_end(nbb,3,tmp2,tmp1,request,reqcnt)

      call nwpw_timing_end(32)
      call nwpw_timing_start(31)
*     *************************************************
*     ***     do fft along ky dimension             ***
*     ***   A(ny(id),nz(id),kx) <- fft1d^(-1)[A(ky,nz(id),kx)]  ***
*     *************************************************
#ifdef MLIB
      indx = 1
      do q=1,nq2(id)
         if (.not.log_mb(zero_row2(1,nbb)+q-1)) then
         call z1dfft(tmp2(indx),ny(id),dcpl_mb(tmpy(1,id)),-2,ierr)
         end if
         indx = indx + ny(id)
      end do
#else
#ifdef FFTW3
      indx = 1
      do q=1,nq2(id)
         if (.not.log_mb(zero_row2(1,nbb)+q-1)) then
            call dfftw_execute_dft(plans(18,id),tmp2(indx),tmp2(indx))
         end if
         indx = indx + ny(id)
      end do
#else
      !indx = 1
      do q=1+tid,nq2(id),nthr
         indx = 1 + (q-1)*ny(id)
         if (.not.log_mb(zero_row2(1,nbb)+q-1)) then
            call dcfftb(ny(id),tmp2(indx),
     >                  dcpl_mb(tmpy(1,id)+tid*(2*ny(id)+15)))
         end if
         !indx = indx + ny(id)
      end do
!$OMP BARRIER

#endif
#endif

      call nwpw_timing_end(31)
      call nwpw_timing_start(32)
      call D3dB_c_ptranspose_ijk_start(nbb,4,tmp2,tmp1,tmp2,
     >                                 request,reqcnt,46)
      call nwpw_timing_end(32)


      end if

      return
      end

      subroutine  D3dB_pfftby_sub2(n,ny,row2,tmpy,A)
      implicit none
      integer n,ny
      logical    row2(n)
      real*8     tmpy(4*ny+15)
      !complex*16 A(ny,n)
      real*8     A(2*ny,n)
      integer i
      do i=1,n
         if (.not.row2(i)) then
            call dcfftb(ny,A(1,i),tmpy)
         end if
      end do
      return
      end





*     ***********************************
*     *					*
*     *	        D3dB_pfftbx		*
*     *					*
*     ***********************************

      subroutine D3dB_pfftbx(nbb,tmp1,tmp2,request,reqcnt)
      implicit none
      integer nbb
      complex*16  tmp1(*)
      complex*16  tmp2(*)
      integer     request(*),reqcnt

#include "bafdecls.fh"
#include "errquit.fh"
#include "D3dB.fh"
#include "D3dB_pfft.fh"

      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / D3dB_fft / tmpx,tmpy,tmpz



*     *** local variables ***
      integer j,q,indx,ierr,id
      integer nxh,tid,nthr

*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      if (nbb.lt.2) then
         id = 1
      else
         id = 3
      end if
      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

      nxh = (nx(id)/2+1)
      
      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then

      call nwpw_timing_start(31)
*     *************************************************
*     ***     do fft along kx dimension             ***
*     ***   A(nx(id),ny(id),nz(id)) <- fft1d^(-1)[A(kx,ny(id),nz(id))]  ***
*     *************************************************
#ifdef MLIB
      indx = 1
      do q=1,nq(id)
      do j=1,ny(id)
         call drc1ft(tmp2(indx),nx(id),dcpl_mb(tmpx(1,id)),-2,ierr)
         indx = indx + nxh
      end do
      end do

#else

#ifdef FFTW3
      call dfftw_execute_dft_c2r(plans(3,id),tmp2,tmp2)

#else
#ifdef NUSE_OPENMP
      call cshift1_fftb(nx(id),ny(id),nq(id),1,tmp2)
      call D3dB_pfftbx_sub(nq(id)*ny(id),nx(id),nxh,
     >                     dcpl_mb(tmpx(1,id)),tmp2)
      call zeroend_fftb(nx(id),ny(id),nq(id),1,tmp2)
#else
      call cshift1_fftb(nx(id),ny(id),nq(id),1,tmp2)
      indx = 1
      do q=1,nq(id)
      do j=1,ny(id)
         call drfftb(nx(id),tmp2(indx),dcpl_mb(tmpx(1,id)))
         indx = indx + nxh
      end do
      end do
      call zeroend_fftb(nx(id),ny(id),nq(id),1,tmp2)
#endif
#endif
#endif
      call dcopy((nx(id)+2)*ny(id)*nq(id),tmp2,1,tmp1,1)



      call nwpw_timing_end(31)
*     *************************************************
      !*************************
      !**** hilbert mapping ****
      !*************************
      else



      call nwpw_timing_start(32)
      call D3dB_c_ptranspose_ijk_end(nbb,4,tmp1,tmp2,request,reqcnt)
      call nwpw_timing_end(32)
      call nwpw_timing_start(31)

*     *************************************************
*     ***     do fft along kx dimension             ***
*     ***   A(nx(id),ny(id),nz(id)) <- fft1d^(-1)[A(kx,ny(id),nz(id))]  ***
*     *************************************************
#ifdef MLIB
      indx = 1
      do q=1,nq1(id)
         call drc1ft(tmp1(indx),nx(id),dcpl_mb(tmpx(1,id)),-2,ierr)
         indx = indx + nxh
      end do
#else

#ifdef FFTW3
      call dfftw_execute_dft_c2r(plans(13,id),tmp1,tmp1)

#else


      !call cshift1_fftb(nx(id),nq1(id),1,1,tmp1)
      indx = 1
      do q=1+tid,nq1(id),nthr
         indx = 1+(q-1)*nxh
         call cshift1_fftb1(nx(id),tmp1(indx))
         call drfftb(nx(id),tmp1(indx),
     >               dcpl_mb(tmpx(1,id)+tid*(2*nx(id)+15))) 
         call zeroend_fftb1(nx(id),tmp1(indx))
         !indx = indx + nxh
      end do
      !call zeroend_fftb(nx(id),nq1(id),1,1,tmp1)

#endif
#endif

      call nwpw_timing_end(31)
      end if
    
      return
      end


      subroutine  D3dB_pfftbx_sub(n,nx,nxh,tmpx,A)
      implicit none
      integer n,nx,nxh
      real*8     tmpx(2*nx+15)
      !complex*16 A(nxh,n)
      real*8     A(2*nxh,n)
      integer i
      do i=1,n
         call drfftb(nx,A(1,i),tmpx)
      end do

      return
      end



*     ***********************************
*     *                                 *
*     *         D3dB_pfftb_step         *
*     *                                 *
*     ***********************************
*

      subroutine D3dB_pfftb_step(step,nbb,A,tmp1,tmp2,request,reqcnt)
      implicit none
      integer step
      integer nbb
      complex*16 A(*)
      complex*16 tmp1(*)
      complex*16 tmp2(*)
      integer    request(*),reqcnt


#include "bafdecls.fh"

      integer tmp1z(2),tmp2z(2)
      common / pfft_queuez_common / tmp1z,tmp2z

c      integer  Parallel_threadid,tid
c      external Parallel_threadid

c      tid = Parallel_threadid()

      if (step.eq.1) then
         call nwpw_timing_start(32)
         call Pack_c_unpack_start(nbb,A,tmp1,request,reqcnt,47)
         call nwpw_timing_end(32)
      else if (step.eq.2) then
         call nwpw_timing_start(32)
         call Pack_c_unpack_mid(nbb,tmp1,tmp2,
     >                          dcpl_mb(tmp1z(1)),dcpl_mb(tmp2z(1)),
     >                          request,reqcnt,48)
         call nwpw_timing_end(32)
      else if (step.eq.3) then
         call nwpw_timing_start(32)
         call Pack_c_unpack_end(nbb,tmp1,
     >                          dcpl_mb(tmp1z(1)),dcpl_mb(tmp2z(1)),
     >                          request,reqcnt)
         call nwpw_timing_end(32)
      else if (step.eq.4) then
         call D3dB_pfftbz(nbb,tmp1,tmp2,request,reqcnt)
      else if (step.eq.5) then
         call D3dB_pfftby(nbb,tmp1,tmp2,request,reqcnt)
      else if (step.eq.6) then
         call D3dB_pfftbx(nbb,tmp1,tmp2,request,reqcnt)
      end if

      return
      end



*     ***********************************
*     *                                 *
*     *      D3dB_pfftb_step_center     *
*     *                                 *
*     ***********************************

      subroutine D3dB_pfftb_step_center(center,
     >                            step,nbb,A,tmp1,tmp2,request,reqcnt)
      implicit none
      integer center(3)
      integer step
      integer nbb
      complex*16 A(*)
      complex*16 tmp1(*)
      complex*16 tmp2(*)
      integer    request(*),reqcnt


#include "bafdecls.fh"

      integer tmp1z(2),tmp2z(2)
      common / pfft_queuez_common / tmp1z,tmp2z

c      integer  Parallel_threadid,tid
c      external Parallel_threadid

c      tid = Parallel_threadid()

      if (step.eq.1) then
         call nwpw_timing_start(32)
         call Pack_c_unpack_start(nbb,A,tmp1,request,reqcnt,47)
         call nwpw_timing_end(32)
      else if (step.eq.2) then
         call nwpw_timing_start(32)
         call Pack_c_unpack_mid(nbb,tmp1,tmp2,
     >                          dcpl_mb(tmp1z(1)),dcpl_mb(tmp2z(1)),
     >                          request,reqcnt,48)
         call nwpw_timing_end(32)
      else if (step.eq.3) then
         call nwpw_timing_start(32)
         call Pack_c_unpack_end(nbb,tmp1,
     >                          dcpl_mb(tmp1z(1)),dcpl_mb(tmp2z(1)),
     >                          request,reqcnt)
         call nwpw_timing_end(32)
      else if (step.eq.4) then
         call D3dB_pfftbz(nbb,tmp1,tmp2,request,reqcnt)
      else if (step.eq.5) then
         call D3dB_pfftby(nbb,tmp1,tmp2,request,reqcnt)
      else if (step.eq.6) then
         call D3dB_pfftbx(nbb,tmp1,tmp2,request,reqcnt)

      else if (step.eq.7) then
         call D3dB_r_3to1_Expand_start3(center,tmp1,tmp2,request,reqcnt)
      else if (step.eq.8) then
         call D3dB_r_3to1_Expand_mid3(tmp2,tmp1,request,reqcnt)
      else if (step.eq.9) then
         call D3dB_r_3to1_Expand_final3(tmp1,tmp2,request,reqcnt)
      end if

      return
      end


*     ***********************************
*     *                                 *
*     *     D3dB_pfft3_queue_init       *
*     *                                 *
*     ***********************************

      subroutine D3dB_pfft3_queue_init(qmax_in)
      implicit none
      integer qmax_in


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      !**** pfft_queues ****
      integer aqindx(2),aqstatus(2)
      integer atmp(2),arequest(2),areqcnt(2),acenter(2)
      integer aqmax,aqsize,alast_index
      common / pfft_aqueue_common / aqindx,aqstatus,atmp,arequest,
     >                              areqcnt,acenter,
     >                              aqmax,aqsize,alast_index

      integer bqindx(2),bqstatus(2)
      integer btmp(2),brequest(2),breqcnt(2),bcenter(2)
      integer bqmax,bqsize,blast_index
      common / pfft_bqueue_common / bqindx,bqstatus,btmp,brequest,
     >                              breqcnt,bcenter,
     >                              bqmax,bqsize,blast_index
      
      integer tmp1z(2),tmp2z(2)
      common / pfft_queuez_common / tmp1z,tmp2z

      logical value
      integer np,size,zplane_sizetmp

      call Parallel2d_np_i(np)

c     **** allocate aqueue ****
      aqmax       = qmax_in
      aqsize      = 0
      alast_index = aqmax
      size  = nfft3d(1)*aqmax*2
      value = BA_alloc_get(mt_dcpl,size,'atmp',atmp(2),atmp(1))
      size  = np*aqmax*2
      value = value.and.
     >      BA_alloc_get(mt_int,2*size,'arequest',
     >                   arequest(2),arequest(1))
      size  = aqmax
      value = value.and.
     >        BA_alloc_get(mt_int,size,'aqindx',aqindx(2),aqindx(1))
      value = value.and.
     >      BA_alloc_get(mt_int,size,'aqstatus',aqstatus(2),aqstatus(1))
      value = value.and.
     >        BA_alloc_get(mt_int,size,'areqcnt',areqcnt(2),areqcnt(1))
      size  = 3*aqmax
      value = value.and.
     >        BA_alloc_get(mt_int,size,'acenter',acenter(2),acenter(1))
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)

      call dcopy(4*nfft3d(1)*aqmax,0.0d0,0,dcpl_mb(atmp(1)),1)


c     **** allocate bqueue ****
      bqmax       = qmax_in
      bqsize      = 0
      blast_index = bqmax
      size  = nfft3d(1)*bqmax*2
      value = BA_alloc_get(mt_dcpl,size,'btmp',btmp(2),btmp(1))
      size  = np*bqmax*2
      value = value.and.
     >      BA_alloc_get(mt_int,2*size,'brequest',
     >                   brequest(2),brequest(1))
      size  = bqmax
      value = value.and.
     >        BA_alloc_get(mt_int,size,'bqindx',bqindx(2),bqindx(1))
      value = value.and.
     >      BA_alloc_get(mt_int,size,'bqstatus',bqstatus(2),bqstatus(1))
      value = value.and.
     >        BA_alloc_get(mt_int,size,'breqcnt',breqcnt(2),breqcnt(1))
      size  = 3*bqmax
      value = value.and.
     >        BA_alloc_get(mt_int,size,'bcenter',bcenter(2),bcenter(1))
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)
      call dcopy(4*nfft3d(1)*bqmax,0.0d0,0,dcpl_mb(btmp(1)),1)



      call D3dB_zplane_size(1,zplane_sizetmp)
      value = value.and.
     >        BA_alloc_get(mt_dcpl,(zplane_sizetmp),
     >                    'tmp1z',tmp1z(2),tmp1z(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,(zplane_sizetmp),
     >                    'tmp2z',tmp2z(2),tmp2z(1))
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)
      


      return
      end




*     ***********************************
*     *                                 *
*     *     D3dB_pfft3_queue_end        *
*     *                                 *
*     ***********************************

      subroutine D3dB_pfft3_queue_end()
      implicit none


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      !**** pfft_queues ****
      integer aqindx(2),aqstatus(2)
      integer atmp(2),arequest(2),areqcnt(2),acenter(2)
      integer aqmax,aqsize,alast_index
      common / pfft_aqueue_common / aqindx,aqstatus,atmp,arequest,
     >                              areqcnt,acenter,
     >                              aqmax,aqsize,alast_index

      integer bqindx(2),bqstatus(2)
      integer btmp(2),brequest(2),breqcnt(2),bcenter(2)
      integer bqmax,bqsize,blast_index
      common / pfft_bqueue_common / bqindx,bqstatus,btmp,brequest,
     >                              breqcnt,bcenter,
     >                              bqmax,bqsize,blast_index

      integer tmp1z(2),tmp2z(2)
      common / pfft_queuez_common / tmp1z,tmp2z

      logical value

      value =           BA_free_heap(atmp(2))
      value = value.and.BA_free_heap(arequest(2))
      value = value.and.BA_free_heap(aqindx(2))
      value = value.and.BA_free_heap(aqstatus(2))
      value = value.and.BA_free_heap(areqcnt(2))
      value = value.and.BA_free_heap(acenter(2))

      value = value.and.BA_free_heap(btmp(2))
      value = value.and.BA_free_heap(brequest(2))
      value = value.and.BA_free_heap(bqindx(2))
      value = value.and.BA_free_heap(bqstatus(2))
      value = value.and.BA_free_heap(breqcnt(2))
      value = value.and.BA_free_heap(bcenter(2))

      value = value.and.BA_free_heap(tmp1z(2))
      value = value.and.BA_free_heap(tmp2z(2))
      if (.not. value) 
     > call errquit('D3dB_pfft3_queue_end:error freeing heap',0,MA_ERR)

      return
      end







*     ***********************************
*     *                                 *
*     *     D3dB_cr_pfft3_queue_filled  *
*     *                                 *
*     ***********************************

      logical function D3dB_cr_pfft3_queue_filled()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      !**** pfft_queue ****
      integer qindx(2),qstatus(2)
      integer tmp(2),request(2),reqcnt(2),qcenter(2)
      integer qmax,qsize,last_index
      common / pfft_aqueue_common / qindx,qstatus,tmp,request,reqcnt,
     >                              qcenter,qmax,qsize,last_index

      D3dB_cr_pfft3_queue_filled = (qsize.ge.qmax)
      return
      end



*     ***********************************
*     *                                 *
*     *     D3dB_rc_pfft3_queue_filled  *
*     *                                 *
*     ***********************************

      logical function D3dB_rc_pfft3_queue_filled()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      !**** pfft_queue ****
      integer qindx(2),qstatus(2)
      integer tmp(2),request(2),reqcnt(2),qcenter(2)
      integer qmax,qsize,last_index
      common / pfft_bqueue_common / qindx,qstatus,tmp,request,reqcnt,
     >                              qcenter,qmax,qsize,last_index

      D3dB_rc_pfft3_queue_filled = (qsize.ge.qmax)
      return
      end




*     ***********************************
*     *                                 *
*     *     D3dB_rc_pfft3f_queuein      *
*     *                                 *
*     ***********************************

      subroutine D3dB_rc_pfft3f_queuein(nbb,A)
      implicit none
      integer nbb
      complex*16 A(*)


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      !**** pfft_queue ****
      integer qindx(2),qstatus(2)
      integer tmp(2),request(2),reqcnt(2),qcenter(2)
      integer qmax,qsize,last_index
      common / pfft_bqueue_common / qindx,qstatus,tmp,request,reqcnt,
     >                              qcenter,qmax,qsize,last_index

      !*** local variables ***
      integer shift1,shift2,shift3,id
      integer q,indx,status,np

c      call nwpw_timing_start(1)
      call nwpw_timing_start(30)
      if (nbb.lt.2) then
         id = 1
      else
         id = 3
      end if
      call Parallel2d_np_i(np)

      do q=1,qsize
         indx   = int_mb(qindx(1)+q-1)
         !int_mb(qstatus(1)+indx-1) = int_mb(qstatus(1)+indx-1) + 1
         status = int_mb(qstatus(1)+indx-1)+1
         shift1=nfft3d(id)*(2*(indx-1))
         shift2=nfft3d(id)*(2*(indx-1)+1)
         shift3=2*np*(indx-1)
         call D3dB_pfftf_step(status,nbb,A,
     >                        dcpl_mb(tmp(1)+shift1),
     >                        dcpl_mb(tmp(1)+shift2),
     >                        int_mb(request(1)+shift3),
     >                        int_mb(reqcnt(1)+indx-1))
         call Parallel_shared_iinc(.true.,1,int_mb(qstatus(1)+indx-1))
      end do

c      qsize = qsize + 1
c      last_index = last_index+1
c      if (last_index.gt.qmax) last_index = 1
c      int_mb(qindx(1)+qsize-1)        = last_index
c      int_mb(qstatus(1)+last_index-1) = 1
      call Parallel_shared_irot(.true.,qmax,last_index)
      call Parallel_shared_iinc(.false.,1,qsize)
      call Parallel_shared_icopy(.false.,last_index,
     >                           int_mb(qindx(1)+qsize-1))
      call Parallel_shared_icopy(.true.,1,
     >                           int_mb(qstatus(1)+last_index-1))

      status = 1
      shift1=nfft3d(id)*(2*(last_index-1))
      shift2=nfft3d(id)*(2*(last_index-1)+1)
      shift3=2*np*(last_index-1)

      call D3dB_pfftf_step(status,nbb,A,
     >                     dcpl_mb(tmp(1)+shift1),
     >                     dcpl_mb(tmp(1)+shift2),
     >                     int_mb(request(1)+shift3),
     >                     int_mb(reqcnt(1)+last_index-1))

c      call nwpw_timing_end(1)
      call nwpw_timing_end(30)
      return
      end




*     ***********************************
*     *                                 *
*     *     D3dB_rc_pfft3f_queueout     *
*     *                                 *
*     ***********************************

      subroutine D3dB_rc_pfft3f_queueout(nbb,A)
      implicit none
      integer nbb
      complex*16 A(*)


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      !**** pfft_queue ****
      integer qindx(2),qstatus(2)
      integer tmp(2),request(2),reqcnt(2),qcenter(2)
      integer qmax,qsize,last_index
      common / pfft_bqueue_common / qindx,qstatus,tmp,request,reqcnt,
     >                              qcenter,qmax,qsize,last_index

      !*** local variables ***
      integer shift1,shift2,shift3,id
      integer q,indx,indx1,status,np

      if (nbb.lt.2) then
         id = 1
      else
         id = 3
      end if
c      call nwpw_timing_start(1)
      call nwpw_timing_start(30)
      call Parallel2d_np_i(np)

c      do q=1,qsize
c         indx   = int_mb(qindx(1)+q-1)
c         int_mb(qstatus(1)+indx-1) = int_mb(qstatus(1)+indx-1) + 1
c         status = int_mb(qstatus(1)+indx-1)
c         shift1=nfft3d(1)*(2*(indx-1))
c         shift2=nfft3d(1)*(2*(indx-1)+1)
c         shift3=2*np*(indx-1)
c         call D3dB_pfftf_step(status,nbb,A,
c     >                        dcpl_mb(tmp(1)+shift1),
c     >                        dcpl_mb(tmp(1)+shift2),
c     >                        int_mb(request(1)+shift3),
c     >                        int_mb(reqcnt(1)+indx-1))
c      end do

      indx1   = int_mb(qindx(1))
      do while (int_mb(qstatus(1)+indx1-1).lt.4)
         do q=1,qsize
            indx   = int_mb(qindx(1)+q-1)
            !int_mb(qstatus(1)+indx-1) = int_mb(qstatus(1)+indx-1) + 1
            status = int_mb(qstatus(1)+indx-1)+1
            shift1=nfft3d(id)*(2*(indx-1))
            shift2=nfft3d(id)*(2*(indx-1)+1)
            shift3=2*np*(indx-1)
            call D3dB_pfftf_step(status,nbb,A,
     >                        dcpl_mb(tmp(1)+shift1),
     >                        dcpl_mb(tmp(1)+shift2),
     >                        int_mb(request(1)+shift3),
     >                        int_mb(reqcnt(1)+indx-1))
            call Parallel_shared_iinc(.true.,1,
     >                                int_mb(qstatus(1)+indx-1))
         end do
      end do

c      qsize = qsize -1
c      do q=1,qsize
c        int_mb(qindx(1)+q-1) = int_mb(qindx(1)+q)
c      end do

      shift2=nfft3d(id)*(2*(indx1-1)+1)
      call Pack_c_Copy(nbb,dcpl_mb(tmp(1)+shift2),A)

      call Parallel_shared_iinc(.true.,-1,qsize)
      call Parallel_shared_lshift(.true.,qsize,int_mb(qindx(1)))

c      call nwpw_timing_end(1)
      call nwpw_timing_end(30)
      return
      end





*     ***********************************
*     *                                 *
*     *     D3dB_cr_pfft3b_queuein      *
*     *                                 *
*     ***********************************

      subroutine D3dB_cr_pfft3b_queuein(nbb,A)
      implicit none
      integer nbb
      complex*16 A(*)


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      !**** pfft_queue ****
      integer qindx(2),qstatus(2)
      integer tmp(2),request(2),reqcnt(2),qcenter(2)
      integer qmax,qsize,last_index
      common / pfft_aqueue_common / qindx,qstatus,tmp,request,reqcnt,
     >                              qcenter,qmax,qsize,last_index

      !*** local variables ***
      integer shift1,shift2,shift3,id
      integer q,indx,status,np
      integer  Parallel_threadid,Parallel_nthreads,tid,nthr
      external Parallel_threadid,Parallel_nthreads

      if (nbb.lt.2) then
         id = 1
      else
         id = 3
      end if
c      call nwpw_timing_start(1)
      call nwpw_timing_start(30)
      call Parallel2d_np_i(np)
c      tid  = Parallel_threadid()
c      nthr = Parallel_nthreads()


      do q=1,qsize
         indx   = int_mb(qindx(1)+q-1)
         !int_mb(qstatus(1)+indx-1) = int_mb(qstatus(1)+indx-1) + 1
         status = int_mb(qstatus(1)+indx-1)+1
         shift1=nfft3d(id)*(2*(indx-1))
         shift2=nfft3d(id)*(2*(indx-1)+1)
         shift3=2*np*(indx-1)
         call D3dB_pfftb_step(status,nbb,A,
     >                        dcpl_mb(tmp(1)+shift1),
     >                        dcpl_mb(tmp(1)+shift2),
     >                        int_mb(request(1)+shift3),
     >                        int_mb(reqcnt(1)+indx-1))
      
        call Parallel_shared_iinc(.true.,1,int_mb(qstatus(1)+indx-1))
      end do

c      qsize = qsize + 1
c      last_index = last_index+1
c      if (last_index.gt.qmax) last_index = 1
c      int_mb(qindx(1)+qsize-1)        = last_index
c      int_mb(qstatus(1)+last_index-1) = 1

      call Parallel_shared_irot(.true.,qmax,last_index)
      call Parallel_shared_iinc(.false.,1,qsize)
      call Parallel_shared_icopy(.false.,last_index,
     >                           int_mb(qindx(1)+qsize-1))
      call Parallel_shared_icopy(.true.,1,
     >                           int_mb(qstatus(1)+last_index-1))

      status = 1
      shift1=nfft3d(id)*(2*(last_index-1))
      shift2=nfft3d(id)*(2*(last_index-1)+1)
      shift3=2*np*(last_index-1)

      call D3dB_pfftb_step(status,nbb,A,
     >                     dcpl_mb(tmp(1)+shift1),
     >                     dcpl_mb(tmp(1)+shift2),
     >                     int_mb(request(1)+shift3),
     >                     int_mb(reqcnt(1)+last_index-1))

!$OMP BARRIER
c      call nwpw_timing_end(1)
      call nwpw_timing_end(30)

      return
      end




*     ***********************************
*     *                                 *
*     *     D3dB_cr_pfft3b_queueout     *
*     *                                 *
*     ***********************************

      subroutine D3dB_cr_pfft3b_queueout(nbb,A)
      implicit none
      integer nbb
      complex*16 A(*)


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      !**** pfft_queue ****
      integer qindx(2),qstatus(2)
      integer tmp(2),request(2),reqcnt(2),qcenter(2)
      integer qmax,qsize,last_index
      common / pfft_aqueue_common / qindx,qstatus,tmp,request,reqcnt,
     >                              qcenter, qmax,qsize,last_index

      !*** local variables ***
      integer shift1,shift2,shift3
      integer q,indx,indx1,status,np,id
c      integer  Parallel_threadid,tid
c      external Parallel_threadid


      if (nbb.lt.2) then
         id = 1
      else
         id = 3
      end if
c      tid = Parallel_threadid()

      call Parallel2d_np_i(np)

c      call nwpw_timing_start(1)
      call nwpw_timing_start(30)

c      do q=1,qsize
c         indx   = int_mb(qindx(1)+q-1)
c         int_mb(qstatus(1)+indx-1) = int_mb(qstatus(1)+indx-1) + 1
c         status = int_mb(qstatus(1)+indx-1)
c         shift1=nfft3d(1)*(2*(indx-1))
c         shift2=nfft3d(1)*(2*(indx-1)+1)
c         shift3=2*np*(indx-1)
c         call D3dB_pfftb_step(status,nbb,A,
c     >                        dcpl_mb(tmp(1)+shift1),
c     >                        dcpl_mb(tmp(1)+shift2),
c     >                        int_mb(request(1)+shift3),
c     >                        int_mb(reqcnt(1)+indx-1))
c      end do


      indx1   = int_mb(qindx(1))
      do while (int_mb(qstatus(1)+indx1-1).lt.6)
         do q=1,qsize
            indx   = int_mb(qindx(1)+q-1)
            !int_mb(qstatus(1)+indx-1) = int_mb(qstatus(1)+indx-1) + 1
            status = int_mb(qstatus(1)+indx-1) + 1
            shift1=nfft3d(id)*(2*(indx-1))
            shift2=nfft3d(id)*(2*(indx-1)+1)
            shift3=2*np*(indx-1)
            call D3dB_pfftb_step(status,nbb,A,
     >                        dcpl_mb(tmp(1)+shift1),
     >                        dcpl_mb(tmp(1)+shift2),
     >                        int_mb(request(1)+shift3),
     >                        int_mb(reqcnt(1)+indx-1))

            call Parallel_shared_iinc(.true.,1,
     >                                int_mb(qstatus(1)+indx-1))
         end do
      end do
!$OMP BARRIER

c      call Parallel_shared_iinc(.true.,-1,qsize)
c      call Parallel_shared_lshift(.true.,qsize,int_mb(qindx(1)))

      shift1=nfft3d(id)*(2*(indx1-1))
      call Parallel_shared_vector_copy(.true.,2*nfft3d(id),
     >                                 dcpl_mb(tmp(1)+shift1),
     >                                 A)

      call Parallel_shared_iinc(.true.,-1,qsize)
      call Parallel_shared_lshift(.true.,qsize,int_mb(qindx(1)))

c      call nwpw_timing_end(1)
      call nwpw_timing_end(30)
      return
      end





cccccccccccccccccccccc  queue center routines   ccccccccccccccccccccccccccccccccccccc



*     ******************************************
*     *                                        *
*     *     D3dB_rc_pfft3f_queuein_center      *
*     *                                        *
*     ******************************************

      subroutine D3dB_rc_pfft3f_queuein_center(nbb,center,A)
      implicit none
      integer nbb
      integer center(3)
      complex*16 A(*)


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      !**** pfft_queue ****
      integer qindx(2),qstatus(2)
      integer tmp(2),request(2),reqcnt(2),qcenter(2)
      integer qmax,qsize,last_index
      common / pfft_bqueue_common / qindx,qstatus,tmp,request,reqcnt,
     >                              qcenter,qmax,qsize,last_index

      !*** local variables ***
      integer shift1,shift2,shift3
      integer q,indx,status,np

c      call nwpw_timing_start(1)
      call nwpw_timing_start(30)
      !if (nbb.lt.2) then
      !   id = 1
      !else
      !   id = 3
      !end if
      call Parallel2d_np_i(np)

      do q=1,qsize
         indx   = int_mb(qindx(1)+q-1)
         !int_mb(qstatus(1)+indx-1) = int_mb(qstatus(1)+indx-1) + 1
         status = int_mb(qstatus(1)+indx-1)+1
         shift1=nfft3d(1)*(2*(indx-1))
         shift2=nfft3d(1)*(2*(indx-1)+1)
         shift3=4*np*(indx-1)
         call D3dB_pfftf_step_center(int_mb(qcenter(1)+3*(indx-1)),
     >                        status,nbb,A,
     >                        dcpl_mb(tmp(1)+shift1),
     >                        dcpl_mb(tmp(1)+shift2),
     >                        int_mb(request(1)+shift3),
     >                        int_mb(reqcnt(1)+indx-1))
         call Parallel_shared_iinc(.true.,1,int_mb(qstatus(1)+indx-1))
      end do

c      qsize = qsize + 1
c      last_index = last_index+1
c      if (last_index.gt.qmax) last_index = 1
c      int_mb(qindx(1)+qsize-1)        = last_index
c      int_mb(qstatus(1)+last_index-1) = 1
      call Parallel_shared_irot(.true.,qmax,last_index)
      call Parallel_shared_iinc(.false.,1,qsize)
      call Parallel_shared_icopy(.false.,last_index,
     >                           int_mb(qindx(1)+qsize-1))
      call Parallel_shared_vector_icopy(.false.,3,center,
     >                           int_mb(qcenter(1)+3*(last_index-1)))
      call Parallel_shared_icopy(.true.,1,
     >                           int_mb(qstatus(1)+last_index-1))

      status = 1
      shift1=nfft3d(1)*(2*(last_index-1))
      shift2=nfft3d(1)*(2*(last_index-1)+1)
      shift3=4*np*(last_index-1)

      call D3dB_pfftf_step_center(int_mb(qcenter(1)+3*(last_index-1)),
     >                     status,nbb,A,
     >                     dcpl_mb(tmp(1)+shift1),
     >                     dcpl_mb(tmp(1)+shift2),
     >                     int_mb(request(1)+shift3),
     >                     int_mb(reqcnt(1)+last_index-1))

c      call nwpw_timing_end(1)
      call nwpw_timing_end(30)
      return
      end




*     ******************************************
*     *                                        *
*     *     D3dB_rc_pfft3f_queueout_center     *
*     *                                        *
*     ******************************************

      subroutine D3dB_rc_pfft3f_queueout_center(nbb,A)
      implicit none
      integer nbb
      complex*16 A(*)


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      !**** pfft_queue ****
      integer qindx(2),qstatus(2)
      integer tmp(2),request(2),reqcnt(2),qcenter(2)
      integer qmax,qsize,last_index
      common / pfft_bqueue_common / qindx,qstatus,tmp,request,reqcnt,
     >                              qcenter,qmax,qsize,last_index

      !*** local variables ***
      integer shift1,shift2,shift3
      integer q,indx,indx1,status,np

      !if (nbb.lt.2) then
      !   id = 1
      !else
      !   id = 3
      !end if
c      call nwpw_timing_start(1)
      call nwpw_timing_start(30)
      call Parallel2d_np_i(np)

c      do q=1,qsize
c         indx   = int_mb(qindx(1)+q-1)
c         int_mb(qstatus(1)+indx-1) = int_mb(qstatus(1)+indx-1) + 1
c         status = int_mb(qstatus(1)+indx-1)
c         shift1=nfft3d(1)*(2*(indx-1))
c         shift2=nfft3d(1)*(2*(indx-1)+1)
c         shift3=2*np*(indx-1)
c         call D3dB_pfftf_step(status,nbb,A,
c     >                        dcpl_mb(tmp(1)+shift1),
c     >                        dcpl_mb(tmp(1)+shift2),
c     >                        int_mb(request(1)+shift3),
c     >                        int_mb(reqcnt(1)+indx-1))
c      end do

      indx1   = int_mb(qindx(1))
      do while (int_mb(qstatus(1)+indx1-1).lt.7)
         do q=1,qsize
            indx   = int_mb(qindx(1)+q-1)
            !int_mb(qstatus(1)+indx-1) = int_mb(qstatus(1)+indx-1) + 1
            status = int_mb(qstatus(1)+indx-1)+1
            shift1=nfft3d(1)*(2*(indx-1))
            shift2=nfft3d(1)*(2*(indx-1)+1)
            shift3=4*np*(indx-1)
            call D3dB_pfftf_step_center(int_mb(qcenter(1)+3*(indx-1)),
     >                        status,nbb,A,
     >                        dcpl_mb(tmp(1)+shift1),
     >                        dcpl_mb(tmp(1)+shift2),
     >                        int_mb(request(1)+shift3),
     >                        int_mb(reqcnt(1)+indx-1))
            call Parallel_shared_iinc(.true.,1,
     >                                int_mb(qstatus(1)+indx-1))
         end do
      end do

c      qsize = qsize -1
c      do q=1,qsize
c        int_mb(qindx(1)+q-1) = int_mb(qindx(1)+q)
c      end do

      shift2=nfft3d(1)*(2*(indx1-1)+1)
      call Pack_c_Copy(nbb,dcpl_mb(tmp(1)+shift2),A)

      call Parallel_shared_iinc(.true.,-1,qsize)
      call Parallel_shared_lshift(.true.,qsize,int_mb(qindx(1)))

c      call nwpw_timing_end(1)
      call nwpw_timing_end(30)
      return
      end





*     ******************************************
*     *                                        *
*     *     D3dB_cr_pfft3b_queuein_center      *
*     *                                        *
*     ******************************************

      subroutine D3dB_cr_pfft3b_queuein_center(nbb,center,A)
      implicit none
      integer nbb
      integer center(3)
      complex*16 A(*)


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      !**** pfft_queue ****
      integer qindx(2),qstatus(2)
      integer tmp(2),request(2),reqcnt(2),qcenter(2)
      integer qmax,qsize,last_index
      common / pfft_aqueue_common / qindx,qstatus,tmp,request,reqcnt,
     >                              qcenter,qmax,qsize,last_index

      !*** local variables ***
      integer shift1,shift2,shift3
      integer q,indx,status,np
      integer  Parallel_threadid,Parallel_nthreads,tid,nthr
      external Parallel_threadid,Parallel_nthreads

      !if (nbb.lt.2) then
      !   id = 1
      !else
      !   id = 3
      !end if
      !id = 1
c      call nwpw_timing_start(1)
      call nwpw_timing_start(30)
      call Parallel2d_np_i(np)
c      tid  = Parallel_threadid()
c      nthr = Parallel_nthreads()


      do q=1,qsize
         indx   = int_mb(qindx(1)+q-1)
         !int_mb(qstatus(1)+indx-1) = int_mb(qstatus(1)+indx-1) + 1
         status = int_mb(qstatus(1)+indx-1)+1
         shift1=nfft3d(1)*(2*(indx-1))
         shift2=nfft3d(1)*(2*(indx-1)+1)
         shift3=4*np*(indx-1)
         call D3dB_pfftb_step_center(int_mb(qcenter(1)+3*(indx-1)),
     >                        status,nbb,A,
     >                        dcpl_mb(tmp(1)+shift1),
     >                        dcpl_mb(tmp(1)+shift2),
     >                        int_mb(request(1)+shift3),
     >                        int_mb(reqcnt(1)+indx-1))
      
        call Parallel_shared_iinc(.true.,1,int_mb(qstatus(1)+indx-1))
      end do

c      qsize = qsize + 1
c      last_index = last_index+1
c      if (last_index.gt.qmax) last_index = 1
c      int_mb(qindx(1)+qsize-1)        = last_index
c      int_mb(qstatus(1)+last_index-1) = 1

      call Parallel_shared_irot(.true.,qmax,last_index)
      call Parallel_shared_iinc(.false.,1,qsize)
      call Parallel_shared_icopy(.false.,last_index,
     >                           int_mb(qindx(1)+qsize-1))
      call Parallel_shared_vector_icopy(.false.,3,center,
     >                           int_mb(qcenter(1)+3*(last_index-1)))
      call Parallel_shared_icopy(.true.,1,
     >                           int_mb(qstatus(1)+last_index-1))

      status = 1
      shift1=nfft3d(1)*(2*(last_index-1))
      shift2=nfft3d(1)*(2*(last_index-1)+1)
      shift3=4*np*(last_index-1)

      call D3dB_pfftb_step_center(int_mb(qcenter(1)+3*(last_index-1)),
     >                     status,nbb,A,
     >                     dcpl_mb(tmp(1)+shift1),
     >                     dcpl_mb(tmp(1)+shift2),
     >                     int_mb(request(1)+shift3),
     >                     int_mb(reqcnt(1)+last_index-1))

!$OMP BARRIER
c      call nwpw_timing_end(1)
      call nwpw_timing_end(30)

      return
      end




*     ******************************************
*     *                                        *
*     *     D3dB_cr_pfft3b_queueout_center     *
*     *                                        *
*     ******************************************

      subroutine D3dB_cr_pfft3b_queueout_center(nbb,A)
      implicit none
      integer nbb
      complex*16 A(*)


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      !**** pfft_queue ****
      integer qindx(2),qstatus(2)
      integer tmp(2),request(2),reqcnt(2),qcenter(2)
      integer qmax,qsize,last_index
      common / pfft_aqueue_common / qindx,qstatus,tmp,request,reqcnt,
     >                              qcenter,qmax,qsize,last_index

      !*** local variables ***
      integer shift1,shift2,shift3
      integer q,indx,indx1,status,np
c      integer  Parallel_threadid,tid
c      external Parallel_threadid


      !if (nbb.lt.2) then
      !   id = 1
      !else
      !   id = 3
      !end if
      !id = 1
c      tid = Parallel_threadid()

      call Parallel2d_np_i(np)

c      call nwpw_timing_start(1)
      call nwpw_timing_start(30)

c      do q=1,qsize
c         indx   = int_mb(qindx(1)+q-1)
c         int_mb(qstatus(1)+indx-1) = int_mb(qstatus(1)+indx-1) + 1
c         status = int_mb(qstatus(1)+indx-1)
c         shift1=nfft3d(1)*(2*(indx-1))
c         shift2=nfft3d(1)*(2*(indx-1)+1)
c         shift3=2*np*(indx-1)
c         call D3dB_pfftb_step(status,nbb,A,
c     >                        dcpl_mb(tmp(1)+shift1),
c     >                        dcpl_mb(tmp(1)+shift2),
c     >                        int_mb(request(1)+shift3),
c     >                        int_mb(reqcnt(1)+indx-1))
c      end do


      indx1   = int_mb(qindx(1))
      do while (int_mb(qstatus(1)+indx1-1).lt.9)
         do q=1,qsize
            indx   = int_mb(qindx(1)+q-1)
            !int_mb(qstatus(1)+indx-1) = int_mb(qstatus(1)+indx-1) + 1
            status = int_mb(qstatus(1)+indx-1) + 1
            shift1=nfft3d(1)*(2*(indx-1))
            shift2=nfft3d(1)*(2*(indx-1)+1)
            shift3=4*np*(indx-1)
            call D3dB_pfftb_step_center(int_mb(qcenter(1)+3*(indx-1)),
     >                        status,nbb,A,
     >                        dcpl_mb(tmp(1)+shift1),
     >                        dcpl_mb(tmp(1)+shift2),
     >                        int_mb(request(1)+shift3),
     >                        int_mb(reqcnt(1)+indx-1))

            call Parallel_shared_iinc(.true.,1,
     >                                int_mb(qstatus(1)+indx-1))
         end do
      end do
!$OMP BARRIER

c      call Parallel_shared_iinc(.true.,-1,qsize)
c      call Parallel_shared_lshift(.true.,qsize,int_mb(qindx(1)))

      shift2=nfft3d(1)*(2*(indx1-1)+1)
      call Parallel_shared_vector_copy(.true.,2*nfft3d(1),
     >                                 dcpl_mb(tmp(1)+shift2),
     >                                 A)

      call Parallel_shared_iinc(.true.,-1,qsize)
      call Parallel_shared_lshift(.true.,qsize,int_mb(qindx(1)))

c      call nwpw_timing_end(1)
      call nwpw_timing_end(30)
      return
      end


