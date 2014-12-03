
#define NBLOCKS 2

*
* $Id$
*

*     ***********************************************************
*     *								*
*     *   		   C3dB_pfft library	 		*
*     *			(MPI implemenation)			*
*     *								*
*     *   Author - Eric Bylaska					*
*     *   date   - 3/23/96					*
*     *								*
*     ***********************************************************

*     ***********************************
*     *					*
*     *	        C3dB_pfft_init		*
*     *					*
*     ***********************************

      subroutine C3dB_pfft_init()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "C3dB.fh"
#include "C3dB_pfft.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      double precision eps
      parameter (eps=1d-12)

      logical value,yrow,zrow,yzslab
      integer nxh,nyh,nzh,q,p
      integer k1,k2,k3,fb,nbb
      integer i,j,k
      integer index2
      double precision ggcut,g1,g2,g3,gg
      integer zero_arow2(2),zero_arow3(2)

*     **** external functions ***
      integer  control_pfft3_qsize,brillioun_nbrillq
      real*8   lattice_ggcut,lattice_wggcut,lattice_unitg
      real*8   brillioun_k
      external control_pfft3_qsize,brillioun_nbrillq
      external lattice_ggcut,lattice_wggcut,lattice_unitg
      external brillioun_k

      call Parallel3d_taskid_i(taskid)
      nxh = nx(1)/2
      nyh = ny(1)/2
      nzh = nz(1)/2

      if (mapping.eq.1) then
         
         value =           BA_alloc_get(mt_log,nx(1)*nq(1),
     >                     'zero_row3',zero_row3(2,0),zero_row3(1,0))
         value = value.and.BA_alloc_get(mt_log,nx(1)*nq(1),
     >                     'zero_row3',zero_row3(2,1),zero_row3(1,1))

         value = value.and.BA_alloc_get(mt_log,nx(1)*nq(1),
     >                     'zero_row2',zero_row2(2,0),zero_row2(1,0))
         value = value.and.BA_alloc_get(mt_log,nx(1)*nq(1),
     >                     'zero_row2',zero_row2(2,1),zero_row2(1,1))

         value = value.and.BA_alloc_get(mt_log,nx(1),
     >               'zero_slab23',zero_slab23(2,0),zero_slab23(1,0))
         value = value.and.BA_alloc_get(mt_log,nx(1),
     >               'zero_slab23',zero_slab23(2,1),zero_slab23(1,1))
         if (.not. value) call errquit('out of heap memory',0, MA_ERR)

         value = value.and.BA_push_get(mt_log,(nx(1)*ny(1)),
     >                     'zero_arow3',zero_arow3(2),zero_arow3(1))
         if (.not. value) call errquit('out of stack memory',0,MA_ERR)


         !*** fb==0 - density, fb=1 - all wavefunctions/brillioun zones ***
         do fb=0,1

         if (fb.eq.0) then
            ggcut = lattice_ggcut()
         else
            ggcut = lattice_wggcut()
         end if

*        **** find zero_row3 - (i,*,k) rows that are zero ****
*        **** note the crazy notation used below ****
         do q=1,nx(1)*nq(1)
             log_mb(zero_row3(1,fb)+q-1) =.true.
         end do
         do q=1,nx(1)*ny(1)
           log_mb(zero_arow3(1)+q-1) =.true.
         end do


         do k2 = -nyh+1, nyh-1
         do k1 = -nxh+1, nxh-1
          i=k1
          j=k2
          if (i .lt. 0) i = i + nx(1)
          if (j .lt. 0) j = j + ny(1)
          i=i+1
          j=j+1
          zrow = .true.
          do k3 = -nzh+1, nzh-1
            if (fb.eq.0) then
               g1 = k1*lattice_unitg(1,1)
     >            + k3*lattice_unitg(1,2)
     >            + k2*lattice_unitg(1,3)
               g2 = k1*lattice_unitg(2,1)
     >            + k3*lattice_unitg(2,2)
     >            + k2*lattice_unitg(2,3)
               g3 = k1*lattice_unitg(3,1)
     >            + k3*lattice_unitg(3,2)
     >            + k2*lattice_unitg(3,3)
               gg = g1*g1 + g2*g2 + g3*g3
               gg= gg-ggcut
               if (gg.lt.-eps) zrow = .false.
            else
               do nbb=1,brillioun_nbrillq()
                  g1 = k1*lattice_unitg(1,1)
     >               + k3*lattice_unitg(1,2)
     >               + k2*lattice_unitg(1,3)
     >               + brillioun_k(1,nbb)
                  g2 = k1*lattice_unitg(2,1)
     >               + k3*lattice_unitg(2,2)
     >               + k2*lattice_unitg(2,3)
     >               + brillioun_k(2,nbb)
                  g3 = k1*lattice_unitg(3,1)
     >               + k3*lattice_unitg(3,2)
     >               + k2*lattice_unitg(3,3)
     >               + brillioun_k(3,nbb)
                  gg = g1*g1 + g2*g2 + g3*g3
                  gg= gg-ggcut
                  if (gg.lt.-eps) zrow = .false.
               end do
            endif

          end do
          if (.not.zrow) then
            log_mb(zero_arow3(1)+(i-1)+nx(1)*(j-1)) =.false.
            call C3dB_ktoqp(1,j,q,p)
            if (p.eq.taskid) then
              index2 = i + nx(1)*(q-1)
              log_mb(zero_row3(1,fb)+index2-1) =.false.
            end if
          end if

         end do
         end do

         call C3dB_c_ptranspose_jk_init(fb,log_mb(zero_arow3(1)))


*        **** find zero_slab23 - (i,*,*) slabs that are zero ****
         do k1 = 1,nx(1)
           log_mb(zero_slab23(1,fb)+k1-1) =.true.
         end do

         do k1 = -nxh+1, nxh-1
          i=k1
          if (i .lt. 0) i = i + nx(1)
          i=i+1
          yzslab = .true.
          do k3 = -nzh+1, nzh-1
          do k2 = -nyh+1, nyh-1
            if (fb.eq.0) then
               g1 = k1*lattice_unitg(1,1)
     >            + k2*lattice_unitg(1,2)
     >            + k3*lattice_unitg(1,3)
               g2 = k1*lattice_unitg(2,1)
     >            + k2*lattice_unitg(2,2)
     >            + k3*lattice_unitg(2,3)
               g3 = k1*lattice_unitg(3,1)
     >            + k2*lattice_unitg(3,2)
     >            + k3*lattice_unitg(3,3)
               gg = g1*g1 + g2*g2 + g3*g3
               gg= gg-ggcut
               if (gg.lt.-eps) yzslab = .false.
            else
               do nbb=1,brillioun_nbrillq()
                  g1 = k1*lattice_unitg(1,1)
     >               + k2*lattice_unitg(1,2)
     >               + k3*lattice_unitg(1,3)
     >               + brillioun_k(1,nbb)
                  g2 = k1*lattice_unitg(2,1)
     >               + k2*lattice_unitg(2,2)
     >               + k3*lattice_unitg(2,3)
     >               + brillioun_k(2,nbb)
                  g3 = k1*lattice_unitg(3,1)
     >               + k2*lattice_unitg(3,2)
     >               + k3*lattice_unitg(3,3)
     >               + brillioun_k(3,nbb)
                  gg = g1*g1 + g2*g2 + g3*g3
                  gg= gg-ggcut
                  if (gg.lt.-eps) yzslab = .false.
               end do
            end if
          end do
          end do
          if (.not.yzslab) then
            log_mb(zero_slab23(1,fb)+i-1) =.false.
          end if

         end do

*        **** find zero_row2 - (i,*,k) rows that are zero after fft of (i,j,*) ****
         do k3 = 1,nz(1)
         do k1 = 1,nx(1)
           call C3dB_ktoqp(1,k3,q,p)
           if (p.eq.taskid) then
            index2 = k1 + nx(1)*(q-1)
            log_mb(zero_row2(1,fb)+index2-1) 
     >       = log_mb(zero_slab23(1,fb)+k1-1)
           end if
         end do
         end do

        
         end do !*fb*
         value = BA_pop_stack(zero_arow3(2))
         if (.not. value) call errquit('error freeing stack',0,MA_ERR)

      !*** mapping == 2 ***
      else
         value =           BA_alloc_get(mt_log,nq3(1),
     >                     'zero_row3',zero_row3(2,0),zero_row3(1,0))
         value = value.and.BA_alloc_get(mt_log,nq3(1),
     >                     'zero_row3',zero_row3(2,1),zero_row3(1,1))

         value = value.and.BA_alloc_get(mt_log,nq2(1),
     >                     'zero_row2',zero_row2(2,0),zero_row2(1,0))
         value = value.and.BA_alloc_get(mt_log,nq2(1),
     >                     'zero_row2',zero_row2(2,1),zero_row2(1,1))

         value = value.and.BA_alloc_get(mt_log,nx(1),
     >               'zero_slab23',zero_slab23(2,0),zero_slab23(1,0))
         value = value.and.BA_alloc_get(mt_log,nx(1),
     >               'zero_slab23',zero_slab23(2,1),zero_slab23(1,1))
         if (.not. value) call errquit('out of heap memory',0, MA_ERR)

         value = value.and.BA_push_get(mt_log,nx(1)*nz(1),
     >                     'zero_arow2',zero_arow2(2),zero_arow2(1))
         value = value.and.BA_push_get(mt_log,nx(1)*ny(1),
     >                     'zero_arow3',zero_arow3(2),zero_arow3(1))
         if (.not. value) call errquit('out of stack memory',0,MA_ERR)

         !*** fb==0 - density, fb=1 - all wavefunctions/brillioun zones ***
         do fb=0,1

         if (fb.eq.0) then
            ggcut = lattice_ggcut()
         else
            ggcut = lattice_wggcut()
         end if

*        **** find zero_row3 - (i,j,*) rows that are zero ****
         do q = 1,nq3(1)
           log_mb(zero_row3(1,fb)+q-1) =.true.
         end do
         do q = 1,nx(1)*ny(1)
           log_mb(zero_arow3(1)+q-1) =.true.
         end do

         do k2 = -nyh+1, nyh-1
         do k1 = -nxh+1, nxh-1
          i=k1
          j=k2
          if (i .lt. 0) i = i + nx(1)
          if (j .lt. 0) j = j + ny(1)
          i=i+1
          j=j+1
          zrow = .true.
          do k3 = -nzh+1, nzh-1
            if (fb.eq.0) then
               g1 = k1*lattice_unitg(1,1)
     >            + k2*lattice_unitg(1,2)
     >            + k3*lattice_unitg(1,3)
               g2 = k1*lattice_unitg(2,1)
     >            + k2*lattice_unitg(2,2)
     >            + k3*lattice_unitg(2,3)
               g3 = k1*lattice_unitg(3,1)
     >            + k2*lattice_unitg(3,2)
     >            + k3*lattice_unitg(3,3)
               gg = g1*g1 + g2*g2 + g3*g3
               gg= gg-ggcut
               if (gg.lt.-eps) zrow = .false.
            else
               do nbb=1,brillioun_nbrillq()
                  g1 = k1*lattice_unitg(1,1)
     >               + k2*lattice_unitg(1,2)
     >               + k3*lattice_unitg(1,3)
     >               + brillioun_k(1,nbb)
                  g2 = k1*lattice_unitg(2,1)
     >               + k2*lattice_unitg(2,2)
     >               + k3*lattice_unitg(2,3)
     >               + brillioun_k(2,nbb)
                  g3 = k1*lattice_unitg(3,1)
     >               + k2*lattice_unitg(3,2)
     >               + k3*lattice_unitg(3,3)
     >               + brillioun_k(3,nbb)
                  gg = g1*g1 + g2*g2 + g3*g3
                  gg= gg-ggcut
                  if (gg.lt.-eps) zrow = .false.
               end do
            end if
          end do
          if (.not.zrow) then
            log_mb(zero_arow3(1)+(i-1)+nx(1)*(j-1)) =.false.
            q = int_mb(q_map3(1,1)+(i-1)+(j-1)*(nx(1)))
            p = int_mb(p_map3(1,1)+(i-1)+(j-1)*(nx(1)))
            if (p.eq.taskid) then
              log_mb(zero_row3(1,fb)+q-1) =.false.
            end if
          end if

         end do
         end do

*        **** find zero_slab23 - (i,*,*) slabs that are zero ****
         do k1 = 1,nx(1)
           log_mb(zero_slab23(1,fb)+k1-1) =.true.
         end do

         do k1 = -nxh+1, nxh-1
          i=k1
          if (i .lt. 0) i = i + nx(1)
          i=i+1
          yzslab = .true.
          do k3 = -nzh+1, nzh-1
          do k2 = -nyh+1, nyh-1
            if (fb.eq.0) then
               g1 = k1*lattice_unitg(1,1)
     >            + k2*lattice_unitg(1,2)
     >            + k3*lattice_unitg(1,3)
               g2 = k1*lattice_unitg(2,1)
     >            + k2*lattice_unitg(2,2)
     >            + k3*lattice_unitg(2,3)
               g3 = k1*lattice_unitg(3,1)
     >            + k2*lattice_unitg(3,2)
     >            + k3*lattice_unitg(3,3)
               gg = g1*g1 + g2*g2 + g3*g3
               gg= gg-ggcut
               if (gg.lt.-eps) yzslab = .false.
            else
               do nbb=1,brillioun_nbrillq()
                  g1 = k1*lattice_unitg(1,1)
     >               + k2*lattice_unitg(1,2)
     >               + k3*lattice_unitg(1,3)
     >               + brillioun_k(1,nbb)
                  g2 = k1*lattice_unitg(2,1)
     >               + k2*lattice_unitg(2,2)
     >               + k3*lattice_unitg(2,3)
     >               + brillioun_k(2,nbb)
                  g3 = k1*lattice_unitg(3,1)
     >               + k2*lattice_unitg(3,2)
     >               + k3*lattice_unitg(3,3)
     >               + brillioun_k(3,nbb)
                  gg = g1*g1 + g2*g2 + g3*g3
                  gg= gg-ggcut
                  if (gg.lt.-eps) yzslab = .false.
               end do
            end if
          end do
          end do
          if (.not.yzslab) then
            log_mb(zero_slab23(1,fb)+i-1) =.false.
          end if

         end do


*        **** find zero_row2 - (i,*,k) rows that are zero after fft of (i,j,*) ****
         do k = 1,nz(1)
         do i = 1,nx(1)
           q = int_mb(q_map2(1,1)+(k-1)+(i-1)*(nz(1)))
           p = int_mb(p_map2(1,1)+(k-1)+(i-1)*(nz(1)))

           log_mb(zero_arow2(1)+i-1+nx(1)*(k-1))
     >       = log_mb(zero_slab23(1,fb)+i-1)
           if (p.eq.taskid) then
            log_mb(zero_row2(1,fb)+q-1)
     >       = log_mb(zero_slab23(1,fb)+i-1)
           end if
         end do
         end do

         call C3dB_c_ptranspose_ijk_init(fb,
     >                                   log_mb(zero_arow2(1)),
     >                                   log_mb(zero_arow3(1)))


         end do

       value =           BA_pop_stack(zero_arow3(2))
       value = value.and.BA_pop_stack(zero_arow2(2))
       if (.not. value) call errquit('error pop stack',0,MA_ERR)


      end if

      call C3dB_pfft3_queue_init(control_pfft3_qsize())

      return
      end




*     ***********************************
*     *                                 *
*     *         C3dB_pfft_end           *
*     *                                 *
*     ***********************************

      subroutine C3dB_pfft_end()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "C3dB.fh"
#include "C3dB_pfft.fh"



*     **** indexing variables ****
      integer iq_to_i1(2,0:1)
      integer iq_to_i2(2,0:1)
      integer iz_to_i2(2,0:1)
      integer i1_start(2,0:1)
      integer i2_start(2,0:1)
      common / c_ptrans_blk1 / iq_to_i1,iq_to_i2,iz_to_i2,
     >                       i1_start,i2_start

*     **** indexing variables ****
      integer jq_to_i1(2,0:1)
      integer jq_to_i2(2,0:1)
      integer jz_to_i2(2,0:1)
      integer j1_start(2,0:1)
      integer j2_start(2,0:1)
      common / c_ptrans_blk2 / jq_to_i1,jq_to_i2,jz_to_i2,
     >                       j1_start,j2_start


*     *** hilbert tranpose data structure ****
      integer h_iq_to_i1(2,6,0:1)
      integer h_iq_to_i2(2,6,0:1)
      integer h_iz_to_i2(2,6,0:1)
      integer h_iz_to_i2_count(6,0:1)
      integer h_i1_start(2,6,0:1)
      integer h_i2_start(2,6,0:1)
      common / c_ptrans_blk_ijk / h_iq_to_i1,
     >                         h_iq_to_i2,
     >                         h_iz_to_i2,
     >                         h_iz_to_i2_count,
     >                         h_i1_start,
     >                         h_i2_start


      logical value
      integer i

      value =           BA_free_heap(zero_row2(2,0))
      value = value.and.BA_free_heap(zero_row2(2,1))
      value = value.and.BA_free_heap(zero_row3(2,0))
      value = value.and.BA_free_heap(zero_row3(2,1))
      value = value.and.BA_free_heap(zero_slab23(2,0))
      value = value.and.BA_free_heap(zero_slab23(2,1))

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
      end do
      end if

      if (.not. value) call errquit('error freeing heap',0,MA_ERR)

      call C3dB_pfft3_queue_end()
      return
      end




*     ***********************************
*     *					*
*     *	        C3dB_cr_pfft3b		*
*     *					*
*     ***********************************

      subroutine C3dB_cr_pfft3b(nb,nbb,A)

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
*       uses - C3dB_c_ptranspose_jk, dcopy           *
*                                                   *
*****************************************************

      implicit none
      integer nb,nbb
      complex*16  A(*)

#include "bafdecls.fh"
#include "errquit.fh"

#include "C3dB.fh"
#include "C3dB_pfft.fh"

      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / C3dB_fft / tmpx,tmpy,tmpz


*     *** local variables ***
      integer i,j,k,q,indx,fb
      integer indx2
#ifdef MLIB
      integer ierr
#endif

      
      !integer tmp1(2),tmp2(2),tmp3(2)
      integer tmp2(2),tmp3(2)
      logical value


      call nwpw_timing_start(1)

      !*** set fb ***
      if (nbb.eq.0) then
         fb = 0 
      else
         fb = 1 
      end if
 

*     ***** allocate temporary space ****
      !call C3dB_nfft3d(nb,nfft3d)
      value = BA_push_get(mt_dcpl,(nfft3d(nb)),'ffttmp2',
     >                    tmp2(2),tmp2(1))
      value = value.and.
     >      BA_push_get(mt_dcpl,(nfft3d(nb)),'ffttmp3',tmp3(2),tmp3(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      
      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then

*     ****************************************************
*     ***     do fft along ky dimension                ***
*     ***  A(kx,ny(nb),kz) <- fft1d^(-1)[A(kx,ky,kz)]  ***
*     ****************************************************
#ifdef MLIB
      !call z1dfft(dbl_mb(tmp3(1)),nz(nb),dcpl_mb(tmpz(1)),-3,ierr)
      indx2 = 0
      do q=1,nq(nb)
      do i=1,nx(nb)
       indx2 = indx2 + 1
       if (.not.log_mb(zero_row3(1,fb)+indx2-1)) then

         indx = i + (q-1)*nx(nb)*nz(nb)
         call zcopy(nz(nb),A(indx),nx(nb),dcpl_mb(tmp2(1)),1)
         call z1dfft(dcpl_mb(tmp2(1)),nz(nb),
     >               dcpl_mb(tmpz(1,nb)),-2,ierr)
         call zcopy(nz(nb),dcpl_mb(tmp2(1)),1,A(indx),nx(nb))

       end if
      end do
      end do

#else
      !call dcffti(nz(nb),dcpl_mb(tmp1(1)))
      indx2 = 0
      do q=1,nq(nb)
      do i=1,nx(nb)
       indx2 = indx2 + 1
       if (.not.log_mb(zero_row3(1,fb)+indx2-1)) then
         indx = i + (q-1)*nx(nb)*nz(nb)
         call zcopy(nz(nb),A(indx),nx(nb),dcpl_mb(tmp2(1)),1)
         call dcfftb(nz(nb),dcpl_mb(tmp2(1)),dcpl_mb(tmpz(1,nb)))
         call zcopy(nz(nb),dcpl_mb(tmp2(1)),1,A(indx),nx(nb))
       end if
      end do
      end do
#endif

*     ********************************************
*     ***         Do a ptranspose of A          ***
*     ***   A(kx,kz,ny(nb)) <- A(kx,ny(nb),kz)  ***
*     ********************************************
      call C3dB_c_ptranspose1_jk(fb,A,dcpl_mb(tmp2(1)),dcpl_mb(tmp3(1)))

*     *************************************************
*     ***     do fft along kz dimension             ***
*     ***   A(kx,nz(nb),ny(nb)) <- fft1d^(-1)[A(kx,kz,ny(nb))]  ***
*     *************************************************
#ifdef MLIB
     
      indx2 = 0
      do q=1,nq(nb)
      do i=1,nx(nb)
        indx2 = indx2 + 1
        if (.not.log_mb(zero_row2(1,fb)+indx2-1)) then
          indx = i + (q-1)*nx(nb)*ny(nb)
          call zcopy(ny(nb),A(indx),nx(nb),dcpl_mb(tmp2(1)),1)
          call z1dfft(dcpl_mb(tmp2(1)),ny(nb),
     >               dcpl_mb(tmpy(1,nb)),-2,ierr)
          call zcopy(ny(nb),dcpl_mb(tmp2(1)),1,A(indx),nx(nb))
        end if
      end do
      end do
      
#else
      indx2 = 0
      do q=1,nq(nb)
      do i=1,nx(nb)
        indx2 = indx2 + 1
        if (.not.log_mb(zero_row2(1,fb)+indx2-1)) then
         indx = i + (q-1)*nx(nb)*ny(nb)
         call zcopy(ny(nb),A(indx),nx(nb),dcpl_mb(tmp2(1)),1)
         call dcfftb(ny(nb),dcpl_mb(tmp2(1)),dcpl_mb(tmpy(1,nb)))
         call zcopy(ny(nb),dcpl_mb(tmp2(1)),1,A(indx),nx(nb))
        end if
      end do
      end do
#endif

*     *********************************************************************
*     ***     do fft along kx dimension                                 ***
*     ***   A(nx(nb),nz(nb),ny(nb)) <- fft1d^(-1)[A(kx,nz(nb),ny(nb))]  ***
*     *********************************************************************
#ifdef MLIB
      !call drc1ft (dbl_mb(tmp3(1)),nx(nb),dcpl_mb(tmp1(1)),-3,ierr)
      indx = 1
      do q=1,nq(nb)
      do j=1,ny(nb)
         call z1dfft(A(indx),nx(nb),
     >               dcpl_mb(tmpx(1,nb)),-2,ierr)
         indx = indx + nx(nb)
      end do
      end do
#else
      indx = 1
      do q=1,nq(nb)
      do j=1,ny(nb)
         call dcfftb(nx(nb),A(indx),dcpl_mb(tmpx(1,nb)))
         indx = indx + nx(nb)
      end do
      end do
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
         if (.not.log_mb(zero_row3(1,fb)+q-1)) then
         call z1dfft(A(indx),nz(nb),dcpl_mb(tmpz(1,nb)),-2,ierr)
         end if
         indx = indx + nz(nb)
      end do
#else
      indx = 1
      do q=1,nq3(nb)
         if (.not.log_mb(zero_row3(1,fb)+q-1)) then
         call dcfftb(nz(nb),A(indx),dcpl_mb(tmpz(1,nb)))
         end if
         indx = indx + nz(nb)
      end do
#endif

      call C3dB_c_ptranspose_ijk(fb,3,A,
     >                           dcpl_mb(tmp2(1)),
     >                           dcpl_mb(tmp3(1)))

*     *************************************************
*     ***     do fft along ky dimension             ***
*     ***   A(ny(nb),nz(nb),kx) <- fft1d^(-1)[A(ky,nz(nb),kx)]  ***
*     *************************************************
#ifdef MLIB
      indx = 1
      do q=1,nq2(nb)
         if (.not.log_mb(zero_row2(1,fb)+q-1)) then
         call z1dfft(A(indx),ny(nb),dcpl_mb(tmpy(1,nb)),-2,ierr)
         end if
         indx = indx + ny(nb)
      end do
#else
      indx = 1
      do q=1,nq2(nb)
         if (.not.log_mb(zero_row2(1,fb)+q-1)) then
         call dcfftb(ny(nb),A(indx),dcpl_mb(tmpy(1,nb)))
         end if
         indx = indx + ny(nb)
      end do
#endif

      call C3dB_c_ptranspose_ijk(fb,4,A,
     >                           dcpl_mb(tmp2(1)),
     >                           dcpl_mb(tmp3(1)))

*     *************************************************
*     ***     do fft along kx dimension             ***
*     ***   A(nx(nb),ny(nb),nz(nb)) <- fft1d^(-1)[A(kx,ny(nb),nz(nb))]  ***
*     *************************************************
#ifdef MLIB
      indx = 1
      do q=1,nq1(nb)
         call z1dfft(A(indx),nx(nb),dcpl_mb(tmpx(1,nb)),-2,ierr)
         indx = indx + nx(nb)
      end do
#else
      indx = 1
      do q=1,nq1(nb)
         call dcfftb(nx(nb),A(indx),dcpl_mb(tmpx(1,nb)))
         indx = indx + nx(nb)
      end do
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
*     *	        C3dB_rc_pfft3f		*
*     *					*
*     ***********************************

      subroutine C3dB_rc_pfft3f(nb,nbb,A)

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

#include "C3dB.fh"
#include "C3dB_pfft.fh"

      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / C3dB_fft / tmpx,tmpy,tmpz

*     *** local variables ***
      integer i,j,k,q,indx,indx1,indx2,fb
#ifdef MLIB
      integer ierr
#endif

      !integer tmp1(2),tmp2(2),tmp3(2)
      integer tmp2(2),tmp3(2)
      logical value


      call nwpw_timing_start(1)

      !*** set fb ***
      if (nbb.eq.0) then
         fb = 0
      else
         fb = 1
      end if

*     ***** allocate temporary space ****
      !call C3dB_nfft3d(nb,nfft3d)
      value = BA_push_get(mt_dcpl,(nfft3d(nb)),'tmp2',tmp2(2),tmp2(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(nfft3d(nb)),'tmp3',tmp3(2),tmp3(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     ********************************************
*     ***     do fft along nx(nb) dimension        ***
*     ***   A(kx,nz(nb),ny(nb)) <- fft1d[A(nx(nb),nz(nb),ny(nb))]  ***
*     ********************************************
#ifdef MLIB
      indx = 1
      do q=1,nq(nb)
      do j=1,ny(nb)
         call z1dfft(A(indx),nx(nb),dcpl_mb(tmpx(1,nb)),1,ierr)
         indx = indx + nx(nb)
      end do
      end do
#else
      indx = 1
      do q=1,nq(nb)
      do j=1,ny(nb)
         call dcfftf(nx(nb),A(indx),dcpl_mb(tmpx(1,nb)))
         indx = indx + nx(nb)
      end do
      end do
#endif


*     ********************************************
*     ***     do fft along nz(nb) dimension        ***
*     ***   A(kx,kz,nz(nb)) <- fft1d[A(kx,nz(nb),ny(nb))]  ***
*     ********************************************

#ifdef MLIB
      
      indx2 = 0
      do q=1,nq(nb)
      do i=1,nx(nb)
        indx2 = indx2 + 1
        if (.not.log_mb(zero_row2(1,fb)+indx2-1)) then

         indx = i + (q-1)*nx(nb)*ny(nb)
         call zcopy(ny(nb),A(indx),nx(nb),dcpl_mb(tmp2(1)),1)
         call z1dfft(dcpl_mb(tmp2(1)),ny(nb),dcpl_mb(tmpy(1,nb)),1,ierr)
         call zcopy(ny(nb),dcpl_mb(tmp2(1)),1,A(indx),nx(nb))

        end if
      end do
      end do
#else

       indx2 = 0
       do q=1,nq(nb)
       do i=1,nx(nb)
        indx2 = indx2 + 1
        if (.not.log_mb(zero_row2(1,fb)+indx2-1)) then

          indx = i + (q-1)*nx(nb)*ny(nb)
          call zcopy(ny(nb),A(indx),nx(nb),dcpl_mb(tmp2(1)),1)
          call dcfftf(ny(nb),dcpl_mb(tmp2(1)),dcpl_mb(tmpy(1,nb)))
          call zcopy(ny(nb),dcpl_mb(tmp2(1)),1,A(indx),nx(nb))
        end if
       end do
       end do

#endif

*     ********************************************
*     ***         Do a ptranspose of A          ***
*     ***      A(ky,ny(nb),kz) <- A(kx,kz,ny(nb))      ***
*     ********************************************
      call C3dB_c_ptranspose2_jk(fb,A,dcpl_mb(tmp2(1)),dcpl_mb(tmp3(1)))


*     ********************************************
*     ***     do fft along ny(nb) dimension        ***
*     ***   A(kx,ky,kz) <- fft1d[A(kx,ny(nb),kz)]  ***
*     ********************************************
#ifdef MLIB
      
      indx2 = 0
      do q=1,nq(nb)
      do i=1,nx(nb)
        indx2 = indx2 + 1
        if (.not.log_mb(zero_row3(1,fb)+indx2-1)) then
         indx = i + (q-1)*nx(nb)*ny(nb)
         call zcopy(nz(nb),A(indx),nx(nb),dcpl_mb(tmp2(1)),1)
         call z1dfft(dcpl_mb(tmp2(1)),nz(nb),dcpl_mb(tmpz(1,nb)),1,ierr)
         call zcopy(nz(nb),dcpl_mb(tmp2(1)),1,A(indx),nx(nb))

        end if
      end do
      end do
#else

       indx2 = 0
       do q=1,nq(nb)
       do i=1,nx(nb)
        indx2 = indx2 + 1
        if (.not.log_mb(zero_row3(1,fb)+indx2-1)) then
          indx = i + (q-1)*nx(nb)*ny(nb)
          call zcopy(nz(nb),A(indx),nx(nb),dcpl_mb(tmp2(1)),1)
          call dcfftf(nz(nb),dcpl_mb(tmp2(1)),dcpl_mb(tmpz(1,nb)))
          call zcopy(nz(nb),dcpl_mb(tmp2(1)),1,A(indx),nx(nb))
        end if
       end do
       end do

#endif


      !*************************
      !**** hilbert mapping ****
      !*************************
      else

*     ********************************************
*     ***     do fft along nx(nb) dimension        ***
*     ***   A(kx,ny(nb),nz(nb)) <- fft1d[A(nx(nb),ny(nb),nz(nb))]  ***
*     ********************************************
#ifdef MLIB
      indx = 1
      do q=1,nq1(nb)
         call z1dfft(A(indx),nx(nb),dcpl_mb(tmpx(1,nb)),1,ierr)
         indx = indx + nx(nb)
      end do
#else
      indx = 1
      do q=1,nq1(nb)
         call dcfftf(nx(nb),A(indx),dcpl_mb(tmpx(1,nb)))
         indx = indx + nx(nb)
      end do
#endif

      call C3dB_c_ptranspose_ijk(fb,1,A,
     >                           dcpl_mb(tmp2(1)),
     >                           dcpl_mb(tmp3(1)))

*     ********************************************
*     ***     do fft along ny(nb) dimension        ***
*     ***   A(ky,nz(nb),kx) <- fft1d[A(ny(nb),nz(nb),kx)]  ***
*     ********************************************
#ifdef MLIB
      indx = 1
      do q=1,nq2(nb)
         if (.not.log_mb(zero_row2(1,fb)+q-1)) then
         call z1dfft(A(indx),ny(nb),dcpl_mb(tmpy(1,nb)),1,ierr)
         end if
         indx = indx + ny(nb)
      end do
#else
      indx = 1
      do q=1,nq2(nb)
         if (.not.log_mb(zero_row2(1,fb)+q-1)) then
         call dcfftf(ny(nb),A(indx),dcpl_mb(tmpy(1,nb)))
         end if
         indx = indx + ny(nb)
      end do
#endif

      call C3dB_c_ptranspose_ijk(fb,2,A,
     >                          dcpl_mb(tmp2(1)),
     >                          dcpl_mb(tmp3(1)))

*     ********************************************
*     ***     do fft along nz(nb) dimension        ***
*     ***   A(kz,kx,ky) <- fft1d[A(nz(nb),kx,ky)]  ***
*     ********************************************
#ifdef MLIB
      indx = 1
      do q=1,nq3(nb)
         if (.not.log_mb(zero_row3(1,fb)+q-1)) then
         call z1dfft(A(indx),nz(nb),dcpl_mb(tmpz(1,nb)),1,ierr)
         end if
         indx = indx + nz(nb)
      end do
#else
      indx = 1
      do q=1,nq3(nb)
         if (.not.log_mb(zero_row3(1,fb)+q-1)) then
         call dcfftf(nz(nb),A(indx),dcpl_mb(tmpz(1,nb)))
         end if
         indx = indx + nz(nb)
      end do
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
*     *	   C3dB_c_ptranspose_jk_init	*
*     *					*
*     ***********************************

      subroutine C3dB_c_ptranspose_jk_init(fb,zero_arow3)
      implicit none
      integer fb
      logical zero_arow3(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "C3dB.fh"


*     **** indexing variables ****
      integer iq_to_i1(2,0:1)
      integer iq_to_i2(2,0:1)
      integer iz_to_i2(2,0:1)
      integer i1_start(2,0:1)
      integer i2_start(2,0:1)
      common / c_ptrans_blk1 / iq_to_i1,iq_to_i2,iz_to_i2,
     >                       i1_start,i2_start

*     **** indexing variables ****
      integer jq_to_i1(2,0:1)
      integer jq_to_i2(2,0:1)
      integer jz_to_i2(2,0:1)
      integer j1_start(2,0:1)
      integer j2_start(2,0:1)
      common / c_ptrans_blk2 / jq_to_i1,jq_to_i2,jz_to_i2,
     >                       j1_start,j2_start



*     **** local variables ****
      integer proc_to,proc_from
      integer pto,qto,np,taskid
      integer pfrom,qfrom
      integer phere,qhere
      integer index1,index2,index3
      integer jndex1,jndex2,jndex3
      integer itmp,ii,jj
      integer i,j,k,it
      logical value



*     **** allocate ptrans_blk1 and ptrans_blk2 common block ****
      value = BA_alloc_get(mt_int,(nx(1)*ny(1)*nq(1)),
     >                     'piq_to_i1',iq_to_i1(2,fb),iq_to_i1(1,fb))
      value=value.and.BA_alloc_get(mt_int,(nx(1)*ny(1)*nq(1)),
     >                     'piq_to_i2',iq_to_i2(2,fb),iq_to_i2(1,fb))
      value=value.and.BA_alloc_get(mt_int,(nx(1)*ny(1)*nq(1)),
     >                     'piz_to_i2',iz_to_i2(2,fb),iz_to_i2(1,fb))

      value = value.and.BA_alloc_get(mt_int,(nz(1)+1),
     >                     'pi1_start',i1_start(2,fb),i1_start(1,fb))
      value = value.and.BA_alloc_get(mt_int,(nz(1)+1),
     >                     'pi2_start',i2_start(2,fb),i2_start(1,fb))


      value = BA_alloc_get(mt_int,(nx(1)*ny(1)*nq(1)),
     >                     'riq_to_i1',jq_to_i1(2,fb),jq_to_i1(1,fb))
      value=value.and.BA_alloc_get(mt_int,(nx(1)*ny(1)*nq(1)),
     >                     'riq_to_i2',jq_to_i2(2,fb),jq_to_i2(1,fb))
      value=value.and.BA_alloc_get(mt_int,(nx(1)*ny(1)*nq(1)),
     >                     'riz_to_i2',jz_to_i2(2,fb),jz_to_i2(1,fb))

      value = value.and.BA_alloc_get(mt_int,(nz(1)+1),
     >                     'ri1_start',j1_start(2,fb),j1_start(1,fb))
      value = value.and.BA_alloc_get(mt_int,(nz(1)+1),
     >                     'ri2_start',j2_start(2,fb),j2_start(1,fb))

      if (.not. value) 
     > call errquit('C3dB_ptranspose_jk_init:out of heap',0,MA_ERR)

      call Parallel3d_taskid_i(taskid)
      call Parallel3d_np_i(np)


      index1 = 1 
      index2 = 1
      index3 = 1
      jndex1 = 1 
      jndex2 = 1
      jndex3 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
         int_mb(i1_start(1,fb)+it) = index1
         int_mb(i2_start(1,fb)+it) = index2
         int_mb(j1_start(1,fb)+it) = jndex1
         int_mb(j2_start(1,fb)+it) = jndex2

         do k=1,nz(1)
         do j=1,ny(1)

*           **** packing scheme **** 
            call C3dB_ktoqp(1,k,qhere,phere)
            call C3dB_ktoqp(1,j,qto,pto)
            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
               do i=1,nx(1)
                  ii = i + nx(1)*(k-1)
                  jj = i + nx(1)*(j-1)
                  itmp = i + (j-1)*nx(1) 
     >                     + (qhere-1)*nx(1)*ny(1)

                  if (.not.zero_arow3(ii)) then
                  int_mb(iq_to_i1(1,fb)+index1-1) = itmp
                  index1 = index1 + 1
                  end if

                  if (.not.zero_arow3(jj)) then
                  int_mb(jq_to_i1(1,fb)+jndex1-1) = itmp
                  jndex1 = jndex1 + 1
                  end if

               end do
            end if
             
*           **** unpacking scheme ****
            call C3dB_ktoqp(1,j,qhere,phere)
            call C3dB_ktoqp(1,k,qfrom,pfrom)
            if ((phere.eq.taskid).and.(pfrom.eq.proc_from)) then
               do i=1,nx(1)
                  ii = i + nx(1)*(k-1)
                  jj = i + nx(1)*(j-1)
                  itmp = i + (k-1)*nx(1)
     >                     + (qhere-1)*nx(1)*ny(1)
                  if (zero_arow3(ii)) then
                  int_mb(iz_to_i2(1,fb)+index3-1) = itmp
                  index3 = index3 + 1
                  else
                  int_mb(iq_to_i2(1,fb)+index2-1) = itmp
                  index2 = index2 + 1
                  end if

                  if (zero_arow3(jj)) then
                  int_mb(jz_to_i2(1,fb)+jndex3-1) = itmp
                  jndex3 = jndex3 + 1
                  else
                  int_mb(jq_to_i2(1,fb)+jndex2-1) = itmp
                  jndex2 = jndex2 + 1
                  end if
               end do
            end if
         end do
         end do
      end do
      int_mb(i1_start(1,fb)+np) = index1
      int_mb(i2_start(1,fb)+np) = index2
      int_mb(j1_start(1,fb)+np) = jndex1
      int_mb(j2_start(1,fb)+np) = jndex2


      return
      end





*     ***********************************
*     *					*
*     *	   C3dB_c_ptranspose_ijk_init	*
*     *					*
*     ***********************************

      subroutine C3dB_c_ptranspose_ijk_init(fb,zero_arow2,zero_arow3)
      implicit none
      integer fb
      logical zero_arow2(*)
      logical zero_arow3(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "C3dB.fh"


*     *** hilbert tranpose data structure ****
      integer h_iq_to_i1(2,6,0:1)
      integer h_iq_to_i2(2,6,0:1)
      integer h_iz_to_i2(2,6,0:1)
      integer h_iz_to_i2_count(6,0:1)
      integer h_i1_start(2,6,0:1)
      integer h_i2_start(2,6,0:1)
      common / c_ptrans_blk_ijk / h_iq_to_i1,
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
      integer i,j,k,it


      call Parallel3d_taskid_i(taskid)
      call Parallel3d_np_i(np)


*     ********************************************************
*     **** map1to2 mapping - done - tranpose operation #1  ***
*     ****   (ny,nz,nx)  <-- (nx,ny,nz)                    ***
*     ****   use zero_arow2                                ***
*     ********************************************************

*     **** allocate trans_blk_ijk common block ****
      value = BA_alloc_get(mt_int,(nx(1)*nq1(1)),
     >                     'h_iq_to_i1_1',
     >                      h_iq_to_i1(2,1,fb),
     >                      h_iq_to_i1(1,1,fb))
      value = value.and.
     >        BA_alloc_get(mt_int,(ny(1)*nq2(1)),
     >                     'h_iq_to_i2_1',
     >                      h_iq_to_i2(2,1,fb),
     >                      h_iq_to_i2(1,1,fb))
      value = value.and.
     >        BA_alloc_get(mt_int,(ny(1)*nq2(1)),
     >                     'h_iz_to_i2_1',
     >                      h_iz_to_i2(2,1,fb),
     >                      h_iz_to_i2(1,1,fb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i1_start_1',
     >                      h_i1_start(2,1,fb),
     >                      h_i1_start(1,1,fb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i2_start_1',
     >                      h_i2_start(2,1,fb),
     >                      h_i2_start(1,1,fb))
      if (.not.value)
     > call errquit('C3dB_transpose_ijk_initt:out of heap',0,MA_ERR)

      index1 = 1 
      index2 = 1
      index3 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
         int_mb(h_i1_start(1,1,fb)+it) = index1
         int_mb(h_i2_start(1,1,fb)+it) = index2

         do k=1,nz(1)
         do j=1,ny(1)
         do i=1,nx(1)
            iszero = zero_arow2(i+(k-1)*nx(1)) 

*           **** packing scheme **** 
            phere = int_mb(p_map1(1,1)+(j-1)+(k-1)*ny(1))
            qhere = int_mb(q_map1(1,1)+(j-1)+(k-1)*ny(1))

            pto   = int_mb(p_map2(1,1)+(k-1)+(i-1)*nz(1))
            qto   = int_mb(q_map2(1,1)+(k-1)+(i-1)*nz(1))
            
           
            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
               itmp = i + (qhere-1)*nx(1)
               if (.not.iszero) then
               int_mb(h_iq_to_i1(1,1,fb)+index1-1) = itmp
               index1 = index1 + 1
               end if
            end if
             
*           **** unpacking scheme ****
            if ((pto.eq.taskid).and.(phere.eq.proc_from)) then
               itmp = j + (qto-1)*ny(1)
               if (.not.iszero) then
               int_mb(h_iq_to_i2(1,1,fb)+index2-1) =  itmp
               index2 = index2 + 1
               else
               int_mb(h_iz_to_i2(1,1,fb)+index3-1) =  itmp
               index3 = index3 + 1
               end if
            end if

         end do
         end do
         end do

      end do
      int_mb(h_i1_start(1,1,fb)+np) = index1
      int_mb(h_i2_start(1,1,fb)+np) = index2
      h_iz_to_i2_count(1,fb) = index3 - 1





*     *********************************************************
*     **** map2to3 mapping - done - transpose operation #2 ****
*     ****    (nz,nx,ny)  <-- (ny,nz,nx)                   ****
*     ****    use zero_arow3                               ****
*     *********************************************************

*     **** allocate trans_blk_ijk common block ****
      value = BA_alloc_get(mt_int,(ny(1)*nq2(1)),
     >                     'h_iq_to_i1_2',
     >                      h_iq_to_i1(2,2,fb),
     >                      h_iq_to_i1(1,2,fb))
      value = value.and.
     >        BA_alloc_get(mt_int,(nz(1)*nq3(1)),
     >                     'h_iq_to_i2_2',
     >                      h_iq_to_i2(2,2,fb),
     >                      h_iq_to_i2(1,2,fb))
      value = value.and.
     >        BA_alloc_get(mt_int,(nz(1)*nq3(1)),
     >                     'h_iz_to_i2_2',
     >                      h_iz_to_i2(2,2,fb),
     >                      h_iz_to_i2(1,2,fb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i1_start_2',
     >                      h_i1_start(2,2,fb),
     >                      h_i1_start(1,2,fb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i2_start_2',
     >                      h_i2_start(2,2,fb),
     >                      h_i2_start(1,2,fb))
      if (.not.value)
     > call errquit('C3dB_transpose_ijk_init:out of heap',0,MA_ERR)

      index1 = 1
      index2 = 1
      index3 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
         int_mb(h_i1_start(1,2,fb)+it) = index1
         int_mb(h_i2_start(1,2,fb)+it) = index2

         do k=1,nz(1)
         do j=1,ny(1)
         do i=1,nx(1)

            iszero = zero_arow3(i+(j-1)*nx(1))

*           **** packing scheme ****
            phere = int_mb(p_map2(1,1)+(k-1)+(i-1)*nz(1))
            qhere = int_mb(q_map2(1,1)+(k-1)+(i-1)*nz(1))

            pto   = int_mb(p_map3(1,1)+(i-1)+(j-1)*nx(1))
            qto   = int_mb(q_map3(1,1)+(i-1)+(j-1)*nx(1))


            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
               itmp = j + (qhere-1)*ny(1)
               if (.not.iszero) then
               int_mb(h_iq_to_i1(1,2,fb)+index1-1) = itmp
               index1 = index1 + 1
               end if
            end if

*           **** unpacking scheme ****
            if ((pto.eq.taskid).and.(phere.eq.proc_from)) then
               itmp = k + (qto-1)*nz(1)
               if (.not.iszero) then
               int_mb(h_iq_to_i2(1,2,fb)+index2-1) = itmp
               index2 = index2 + 1
               else  
               int_mb(h_iz_to_i2(1,2,fb)+index3-1) = itmp
               index3 = index3 + 1
               end if
            end if


         end do
         end do
         end do

      end do
      int_mb(h_i1_start(1,2,fb)+np) = index1
      int_mb(h_i2_start(1,2,fb)+np) = index2
      h_iz_to_i2_count(2,fb) = index3 - 1






*     ********************************************************
*     **** map3to2 mapping - done - tranpose operation #3 ****
*     ****  (ny,nz,nx)  <-- (nz,nx,ny)                    ****
*     ****   use zero_arow3                               ****
*     ********************************************************

*     **** allocate trans_blk_ijk common block ****
      value = BA_alloc_get(mt_int,(nz(1)*nq3(1)),
     >                     'h_iq_to_i1_3',
     >                      h_iq_to_i1(2,3,fb),
     >                      h_iq_to_i1(1,3,fb))
      value = value.and.
     >        BA_alloc_get(mt_int,(ny(1)*nq2(1)),
     >                     'h_iq_to_i2_3',
     >                      h_iq_to_i2(2,3,fb),
     >                      h_iq_to_i2(1,3,fb))
      value = value.and.
     >        BA_alloc_get(mt_int,(ny(1)*nq2(1)),
     >                     'h_iz_to_i2_3',
     >                      h_iz_to_i2(2,3,fb),
     >                      h_iz_to_i2(1,3,fb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i1_start_3',
     >                      h_i1_start(2,3,fb),
     >                      h_i1_start(1,3,fb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i2_start_3',
     >                      h_i2_start(2,3,fb),
     >                      h_i2_start(1,3,fb))
      if (.not.value)
     > call errquit('C3dB_transpose_ijk_init:out of heap',0,MA_ERR)

      index1 = 1
      index2 = 1
      index3 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
         int_mb(h_i1_start(1,3,fb)+it) = index1
         int_mb(h_i2_start(1,3,fb)+it) = index2

         do k=1,nz(1)
         do j=1,ny(1)
         do i=1,nx(1)

            iszero = zero_arow3(i+(j-1)*nx(1))

*           **** packing scheme ****
            phere = int_mb(p_map3(1,1)+(i-1)+(j-1)*nx(1))
            qhere = int_mb(q_map3(1,1)+(i-1)+(j-1)*nx(1))

            pto   = int_mb(p_map2(1,1)+(k-1)+(i-1)*nz(1))
            qto   = int_mb(q_map2(1,1)+(k-1)+(i-1)*nz(1))


            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
               itmp = k + (qhere-1)*nz(1)
               if (.not.iszero) then
               int_mb(h_iq_to_i1(1,3,fb)+index1-1) = itmp
               index1 = index1 + 1
               end if
            end if

*           **** unpacking scheme ****
            if ((pto.eq.taskid).and.(phere.eq.proc_from)) then
               itmp = j + (qto-1)*ny(1)
               if (.not.iszero) then
               int_mb(h_iq_to_i2(1,3,fb)+index2-1) = itmp
               index2 = index2 + 1
               else
               int_mb(h_iz_to_i2(1,3,fb)+index3-1) = itmp
               index3 = index3 + 1
               end if
            end if


         end do
         end do
         end do

      end do
      int_mb(h_i1_start(1,3,fb)+np) = index1
      int_mb(h_i2_start(1,3,fb)+np) = index2
      h_iz_to_i2_count(3,fb) = index3 - 1




*     ********************************************************
*     **** map2to1 mapping - done - tranpose operation #4 ****
*     ****  (nx,ny,nz)  <-- (ny,nz,nx)                    ****
*     ****  use zero_arow2                                ****
*     ********************************************************

*     **** allocate trans_blk_ijk common block ****
      value = BA_alloc_get(mt_int,(ny(1)*nq2(1)),
     >                     'h_iq_to_i1_4',
     >                      h_iq_to_i1(2,4,fb),
     >                      h_iq_to_i1(1,4,fb))
      value = value.and.
     >        BA_alloc_get(mt_int,(nx(1)*nq1(1)),
     >                     'h_iq_to_i2_4',
     >                      h_iq_to_i2(2,4,fb),
     >                      h_iq_to_i2(1,4,fb))
      value = value.and.
     >        BA_alloc_get(mt_int,(nx(1)*nq1(1)),
     >                     'h_iz_to_i2_4',
     >                      h_iz_to_i2(2,4,fb),
     >                      h_iz_to_i2(1,4,fb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i1_start_4',
     >                      h_i1_start(2,4,fb),
     >                      h_i1_start(1,4,fb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i2_start_4',
     >                      h_i2_start(2,4,fb),
     >                      h_i2_start(1,4,fb))
      if (.not.value)
     > call errquit('C3dB_transpose_ijk_init:out of heap',0,MA_ERR)

      index1 = 1
      index2 = 1
      index3 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
         int_mb(h_i1_start(1,4,fb)+it) = index1
         int_mb(h_i2_start(1,4,fb)+it) = index2

         do k=1,nz(1)
         do j=1,ny(1)
         do i=1,nx(1)

            iszero = zero_arow2(i+(k-1)*nx(1))

*           **** packing scheme ****
            phere = int_mb(p_map2(1,1)+(k-1)+(i-1)*nz(1))
            qhere = int_mb(q_map2(1,1)+(k-1)+(i-1)*nz(1))

            pto   = int_mb(p_map1(1,1)+(j-1)+(k-1)*ny(1))
            qto   = int_mb(q_map1(1,1)+(j-1)+(k-1)*ny(1))


            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
               itmp = j + (qhere-1)*ny(1)
               if (.not.iszero) then
               int_mb(h_iq_to_i1(1,4,fb)+index1-1) = itmp
               index1 = index1 + 1
               end if
            end if

*           **** unpacking scheme ****
            if ((pto.eq.taskid).and.(phere.eq.proc_from)) then
               itmp = i + (qto-1)*nx(1)
               if (.not.iszero) then
               int_mb(h_iq_to_i2(1,4,fb)+index2-1) = itmp
               index2 = index2 + 1
               else
               int_mb(h_iz_to_i2(1,4,fb)+index3-1) = itmp
               index3 = index3 + 1
               end if
            end if


         end do
         end do
         end do

      end do
      int_mb(h_i1_start(1,4,fb)+np) = index1
      int_mb(h_i2_start(1,4,fb)+np) = index2
      h_iz_to_i2_count(4,fb) = index3 - 1





*     **********************************************************
*     **** map1to3 mapping  - done - tranpose operation # 5 ****
*     **********************************************************

*     **** allocate trans_blk_ijk common block ****
      value = BA_alloc_get(mt_int,(nx(1)*nq1(1)),
     >                     'h_iq_to_i1_5',
     >                      h_iq_to_i1(2,5,fb),
     >                      h_iq_to_i1(1,5,fb))
      value = value.and.
     >        BA_alloc_get(mt_int,(nz(1)*nq3(1)),
     >                     'h_iq_to_i2_5',
     >                      h_iq_to_i2(2,5,fb),
     >                      h_iq_to_i2(1,5,fb))
      value = value.and.
     >        BA_alloc_get(mt_int,(nz(1)*nq3(1)),
     >                     'h_iz_to_i2_5',
     >                      h_iz_to_i2(2,5,fb),
     >                      h_iz_to_i2(1,5,fb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i1_start_5',
     >                      h_i1_start(2,5,fb),
     >                      h_i1_start(1,5,fb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i2_start_5',
     >                      h_i2_start(2,5,fb),
     >                      h_i2_start(1,5,fb))
      if (.not.value)
     > call errquit('C3dB_transpose_ijk_init:out of heap',0,MA_ERR)

      index1 = 1
      index2 = 1
      index3 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
         int_mb(h_i1_start(1,5,fb)+it) = index1
         int_mb(h_i2_start(1,5,fb)+it) = index2

         do k=1,nz(1)
         do j=1,ny(1)
         do i=1,nx(1)

*           **** packing scheme ****
            phere = int_mb(p_map1(1,1)+(j-1)+(k-1)*ny(1))
            qhere = int_mb(q_map1(1,1)+(j-1)+(k-1)*ny(1))

            pto   = int_mb(p_map3(1,1)+(i-1)+(j-1)*nx(1))
            qto   = int_mb(q_map3(1,1)+(i-1)+(j-1)*nx(1))


            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
               itmp = i + (qhere-1)*nx(1)
               int_mb(h_iq_to_i1(1,5,fb)+index1-1) = itmp
               index1 = index1 + 1
            end if

*           **** unpacking scheme ****
            if ((pto.eq.taskid).and.(phere.eq.proc_from)) then
               itmp = k + (qto-1)*nz(1)
               int_mb(h_iq_to_i2(1,5,fb)+index2-1) = itmp
               index2 = index2 + 1
            end if

         end do
         end do
         end do

      end do
      int_mb(h_i1_start(1,5,fb)+np) = index1
      int_mb(h_i2_start(1,5,fb)+np) = index2
      h_iz_to_i2_count(5,fb) = index3 - 1






*     *************************
*     **** map3to1 mapping ****
*     *************************

*     **** allocate trans_blk_ijk common block ****
      value = BA_alloc_get(mt_int,(nz(1)*nq3(1)),
     >                     'h_iq_to_i1_6',
     >                      h_iq_to_i1(2,6,fb),
     >                      h_iq_to_i1(1,6,fb))
      value = value.and.
     >        BA_alloc_get(mt_int,(nx(1)*nq1(1)),
     >                     'h_iq_to_i2_6',
     >                      h_iq_to_i2(2,6,fb),
     >                      h_iq_to_i2(1,6,fb))
      value = value.and.
     >        BA_alloc_get(mt_int,(nx(1)*nq1(1)),
     >                     'h_iz_to_i2_6',
     >                      h_iz_to_i2(2,6,fb),
     >                      h_iz_to_i2(1,6,fb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i1_start_6',
     >                      h_i1_start(2,6,fb),
     >                      h_i1_start(1,6,fb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i2_start_6',
     >                      h_i2_start(2,6,fb),
     >                      h_i2_start(1,6,fb))
      if (.not.value)
     > call errquit('C3dB_transpose_ijk_init:out of heap',0,MA_ERR)

      index1 = 1
      index2 = 1
      index3 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
         int_mb(h_i1_start(1,6,fb)+it) = index1
         int_mb(h_i2_start(1,6,fb)+it) = index2

         do k=1,nz(1)
         do j=1,ny(1)
         do i=1,nx(1)

*           **** packing scheme ****
            phere = int_mb(p_map3(1,1)+(i-1)+(j-1)*nx(1))
            qhere = int_mb(q_map3(1,1)+(i-1)+(j-1)*nx(1))

            pto   = int_mb(p_map1(1,1)+(j-1)+(k-1)*ny(1))
            qto   = int_mb(q_map1(1,1)+(j-1)+(k-1)*ny(1))


            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
               itmp = k + (qhere-1)*nz(1)
               int_mb(h_iq_to_i1(1,6,fb)+index1-1) = itmp
               index1 = index1 + 1
            end if

*           **** unpacking scheme ****
            if ((pto.eq.taskid).and.(phere.eq.proc_from)) then
               itmp = i + (qto-1)*nx(1)
               int_mb(h_iq_to_i2(1,6,fb)+index2-1) = itmp
               index2 = index2 + 1
            end if

         end do
         end do
         end do

      end do
      int_mb(h_i1_start(1,6,fb)+np) = index1
      int_mb(h_i2_start(1,6,fb)+np) = index2
      h_iz_to_i2_count(6,fb) = index3 - 1


      return
      end



*     ***********************************
*     *                                 *
*     *         C3dB_pfftfx             *
*     *                                 *
*     ***********************************
*
*     do fft along nx(1) dimension       
*
*        A(kx,ny(1),nz(1)) <- fft1d[A(nx(1),ny(1),nz(1))] 
*     

      subroutine C3dB_pfftfx(nbb,A,tmp1,tmp2,request,reqcnt)
      implicit none
      integer nbb
      complex*16 A(*)
      complex*16 tmp1(*)
      complex*16 tmp2(*)
      integer    request(*),reqcnt


#include "bafdecls.fh"
#include "errquit.fh"

#include "C3dB.fh"
#include "C3dB_pfft.fh"


      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / C3dB_fft / tmpx,tmpy,tmpz

      !*** local variables ***
      integer j,q,indx,fb

      if (nbb.eq.0) then
         fb = 0
      else
         fb = 1
      end if

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
      call nwpw_timing_start(31)
#ifdef MLIB
      indx = 1
      do q=1,nq(1)
      do j=1,ny(1)
         call z1dfft(A(indx),nx(1),dcpl_mb(tmpx(1,1)),1,ierr)
         indx = indx + nx(1)
      end do
      end do

#else
      indx = 1
      do q=1,nq(1)
      do j=1,ny(1)
         call dcfftf(nx(1),A(indx),dcpl_mb(tmpx(1,1)))
         indx = indx + nx(1)
      end do
      end do
#endif
      call dcopy(2*nx(1)*ny(1)*nq(1),A,1,tmp1,1)
      call nwpw_timing_end(31)


      !*************************
      !**** hilbert mapping ****
      !*************************
      else
      call nwpw_timing_start(31)

#ifdef MLIB
      indx = 1
      do q=1,nq1(1)
         call z1dfft(A(indx),nx(1),dcpl_mb(tmpx(1,1)),1,ierr)
         indx = indx + nx(1)
      end do
#else
      indx = 1
      do q=1,nq1(1)
         call dcfftf(nx(1),A(indx),dcpl_mb(tmpx(1,1)))
         indx = indx + nx(1)
      end do
#endif

      call nwpw_timing_end(31)
      call nwpw_timing_start(32)
      call C3dB_c_ptranspose_ijk_start(fb,1,A,tmp1,tmp2,
     >                                 request,reqcnt,40)
      call nwpw_timing_end(32)

      end if

      return
      end


*     ***********************************
*     *                                 *
*     *         C3dB_pfftfy             *
*     *                                 *
*     ***********************************
*
*     do fft along ny(1) dimension       
*
*        A(kx,ny(1),nz(1)) <- fft1d[A(nx(1),ny(1),nz(1))] 
*     

      subroutine C3dB_pfftfy(nbb,tmp1,tmp2,request,reqcnt)
      implicit none
      integer nbb
      complex*16 tmp1(*)
      complex*16 tmp2(*)
      integer    request(*),reqcnt


#include "bafdecls.fh"
#include "errquit.fh"

#include "C3dB.fh"
#include "C3dB_pfft.fh"


      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / C3dB_fft / tmpx,tmpy,tmpz


      !*** local variables ***
      integer i,j,k,indx,indx2,q,fb

      if (nbb.eq.0) then
         fb = 0
      else
         fb = 1
      end if

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
      call nwpw_timing_start(31)

*     ********************************************
*     ***     do fft along ny(1) dimension        ***
*     ***   A(kx,ky,nz(1)) <- fft1d[A(kx,ny(1),nz(1))]  ***
*     ********************************************

#ifdef MLIB

      indx2 = 0
      do q=1,nq(1)
      do i=1,nx(1)
        indx2 = indx2 + 1
        if (.not.log_mb(zero_row2(1,fb)+indx2-1)) then
         indx = i + (q-1)*nx(1)*ny(1)
         call zcopy(ny(1),tmp1(indx),nx(1),tmp2,1)
         call z1dfft(tmp2,ny(1),dcpl_mb(tmpy(1,1)),1,ierr)
         call zcopy(ny(1),tmp2,1,tmp1(indx),nx(1))
        end if
      end do
      end do
#else
       indx2 = 0
       do q=1,nq(1)
       do i=1,nx(1)
        indx2 = indx2 + 1
        if (.not.log_mb(zero_row2(1,fb)+indx2-1)) then
          indx = i + (q-1)*nx(1)*ny(1)
          call zcopy(ny(1),tmp1(indx),nx(1),tmp2,1)
          call dcfftf(ny(1),tmp2,dcpl_mb(tmpy(1,1)))
          call zcopy(ny(1),tmp2,1,tmp1(indx),nx(1))
        end if
       end do
       end do
#endif


*     ********************************************
*     ***         Do a ptranspose of A          ***
*     ***      A(ky,nz(1),ky) <- A(kx,ky,nz(1))      ***
*     ********************************************
      call nwpw_timing_end(31)
      call nwpw_timing_start(32)
      call C3dB_c_ptranspose2_jk_start(fb,tmp1,tmp2,tmp1,
     >                                 request,reqcnt,41)
      call nwpw_timing_end(32)




      !*************************
      !**** hilbert mapping ****
      !*************************
      else
      call nwpw_timing_start(32)
      call C3dB_c_ptranspose_ijk_end(fb,1,tmp1,tmp2,request,reqcnt)
      call nwpw_timing_end(32)
      call nwpw_timing_start(31)

*     ********************************************
*     ***     do fft along ny(1) dimension        ***
*     ***   A(ky,nz(1),kx) <- fft1d[A(ny(1),nz(1),kx)]  ***
*     ********************************************
#ifdef MLIB
      indx = 1
      do q=1,nq2(1)
         if (.not.log_mb(zero_row2(1,fb)+q-1)) then
         call z1dfft(tmp1(indx),ny(1),dcpl_mb(tmpy(1,1)),1,ierr)
         end if
         indx = indx + ny(1)
      end do
#else
      indx = 1
      do q=1,nq2(1)
         if (.not.log_mb(zero_row2(1,fb)+q-1)) then
         call dcfftf(ny(1),tmp1(indx),dcpl_mb(tmpy(1,1)))
         end if
         indx = indx + ny(1)
      end do
#endif

      call nwpw_timing_end(31)
      call nwpw_timing_start(32)
      call C3dB_c_ptranspose_ijk_start(fb,2,tmp1,tmp2,tmp1,
     >                                 request,reqcnt,42)
      call nwpw_timing_end(32)


      end if

      return
      end





*     ***********************************
*     *                                 *
*     *         C3dB_pfftfz             *
*     *                                 *
*     ***********************************
*

      subroutine C3dB_pfftfz(nbb,tmp1,tmp2,request,reqcnt)
      implicit none
      integer nbb
      complex*16 tmp1(*)
      complex*16 tmp2(*)
      integer    request(*),reqcnt


#include "bafdecls.fh"
#include "errquit.fh"

#include "C3dB.fh"
#include "C3dB_pfft.fh"


      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / C3dB_fft / tmpx,tmpy,tmpz


      !*** local variables ***
      integer i,k,q,fb
      integer indx,indx2

#ifdef MLIB
      integer ierr
#endif

      if (nbb.eq.0) then
         fb = 0
      else
         fb = 1
      end if

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then

      call nwpw_timing_start(32)
      call C3dB_c_ptranspose2_jk_end(fb,tmp2,tmp1,request,reqcnt)
      call nwpw_timing_end(32)
      call nwpw_timing_start(31)

*     ********************************************
*     ***     do fft along nz(1) dimension        ***
*     ***   A(kx,kz,ky) <- fft1d[A(kx,nz(1),ky)]  ***
*     ********************************************

#ifdef MLIB
      indx2 = 0
      do q=1,nq(1)
      do i=1,nx(1)
        indx2 = indx2 + 1
        if (.not.log_mb(zero_row3(1,fb)+indx2-1)) then
         indx = i + (q-1)*nx(1)*nz(1)
         call zcopy(nz(1),tmp2(indx),nx(1),tmp1,1)
         call z1dfft(tmp1,nz(1),dcpl_mb(tmpz(1,1)),1,ierr)
         call zcopy(nz(1),tmp1,1,tmp2(indx),nx(1))
        end if
      end do
      end do
#else
      indx2 = 0
      do q=1,nq(1)
      do i=1,nx(1)
        indx2 = indx2 + 1
        if (.not.log_mb(zero_row3(1,fb)+indx2-1)) then
          indx = i + (q-1)*nx(1)*nz(1)
          call zcopy(nz(1),tmp2(indx),nx(1),tmp1,1)
          call dcfftf(nz(1),tmp1,dcpl_mb(tmpz(1,1)))
          call zcopy(nz(1),tmp1,1,tmp2(indx),nx(1))
        end if
      end do
      end do
      call nwpw_timing_end(31)

#endif


      !*************************
      !**** hilbert mapping ****
      !*************************
      else
      call nwpw_timing_start(32)
      call C3dB_c_ptranspose_ijk_end(fb,2,tmp2,tmp1,request,reqcnt)
      call nwpw_timing_end(32)
      call nwpw_timing_start(31)

*     ********************************************
*     ***     do fft along nz(1) dimension        ***
*     ***   A(kz,kx,ky) <- fft1d[A(nz(1),kx,ky)]  ***
*     ********************************************
#ifdef MLIB
      indx = 1
      do q=1,nq3(1)
         if (.not.log_mb(zero_row3(1,fb)+q-1)) then
         call z1dfft(tmp2(indx),nz(1),dcpl_mb(tmpz(1,1)),1,ierr)
         end if
         indx = indx + nz(1)
      end do
#else
      indx = 1
      do q=1,nq3(1)
         if (.not.log_mb(zero_row3(1,fb)+q-1)) then
         call dcfftf(nz(1),tmp2(indx),dcpl_mb(tmpz(1,1)))
         end if
         indx = indx + nz(1)
      end do
#endif
      call nwpw_timing_end(31)

      end if



      call nwpw_timing_start(32)
      call Cram_c_pack_start(nbb,tmp2,tmp1,request,reqcnt,43)
      call nwpw_timing_end(32)

      return
      end



*     ***********************************
*     *                                 *
*     *         C3dB_pfftf_final        *
*     *                                 *
*     ***********************************
*

      subroutine C3dB_pfftf_final(nbb,tmp1,tmp2,request,reqcnt)
      implicit none
      integer nbb
      complex*16 tmp1(*)
      complex*16 tmp2(*)
      integer    request(*),reqcnt

#include "bafdecls.fh"
#include "errquit.fh"

#include "C3dB.fh"
#include "C3dB_pfft.fh"

      call nwpw_timing_start(32)
      call Cram_c_pack_end(nbb,tmp2,request,reqcnt)
      call nwpw_timing_end(32)

      return
      end




*     ***********************************
*     *                                 *
*     *         C3dB_pfftf_step         *
*     *                                 *
*     ***********************************
*

      subroutine C3dB_pfftf_step(step,nbb,A,tmp1,tmp2,request,reqcnt)
      implicit none
      integer step
      integer nbb
      complex*16 A(*)
      complex*16 tmp1(*)
      complex*16 tmp2(*)
      integer    request(*),reqcnt


      if (step.eq.1) then
         call C3dB_pfftfx(nbb,A,tmp1,tmp2,request,reqcnt)
      else if (step.eq.2) then
         call C3dB_pfftfy(nbb,tmp1,tmp2,request,reqcnt)
      else if (step.eq.3) then
         call C3dB_pfftfz(nbb,tmp1,tmp2,request,reqcnt)
      else if (step.eq.4) then
         call C3dB_pfftf_final(nbb,tmp1,tmp2,request,reqcnt)
      end if

      return
      end




*     ***********************************
*     *					*
*     *	        C3dB_pfftbz		*
*     *					*
*     ***********************************

*                                                  
*      This routine performs the operation of      
*      a three dimensional complex to complex     
*      inverse fft                                
*           A(nx,ny(1),nz(1)) <- FFT3^(-1)[A(kx,ky,kz)]  
*                                                  
*      Entry -                                     
*              A: a column distribuded 3d block    
*              tmp: tempory work space must be at  
*                    least the size of (complex)   
*                    (nfft*nfft + 1) + 10*nfft     
*                                               
*       Exit - A is transformed and the imaginary
*              part of A is set to zero          
*       uses - C3dB_c_ptranspose_jk, dcopy        
*                                                  


      subroutine C3dB_pfftbz(nbb,tmp1,tmp2,request,reqcnt)
      implicit none
      integer nbb
      complex*16  tmp1(*)
      complex*16  tmp2(*)
      integer     request(*),reqcnt
      

#include "bafdecls.fh"
#include "errquit.fh"
#include "C3dB.fh"
#include "C3dB_pfft.fh"

      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / C3dB_fft / tmpx,tmpy,tmpz



*     *** local variables ***
      integer i,k,q,indx,ierr
      integer indx2,fb

      if (nbb.eq.0) then
         fb = 0
      else
         fb = 1
      end if

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then

      call nwpw_timing_start(31)
*     *************************************************
*     ***     do fft along kz dimension             ***
*     ***   A(kx,nz(1),ky) <- fft1d^(-1)[A(kx,kz,ky)]  ***
*     *************************************************
#ifdef MLIB
      indx2 = 0
      do q=1,nq(1)
      do i=1,nx(1)
       indx2 = indx2 + 1
       if (.not.log_mb(zero_row3(1,fb)+indx2-1)) then
         indx = i + (q-1)*nx(1)*nz(1)
         call zcopy(nz(1),tmp1(indx),nx(1),tmp2,1)
         call z1dfft(tmp2,nz(1),dcpl_mb(tmpz(1,1)),-2,ierr)
         call zcopy(nz(1),tmp2,1,tmp1(indx),nx(1))
       end if
      end do
      end do
#else
      indx2 = 0
      do q=1,nq(1)
      do i=1,nx(1)
       indx2 = indx2 + 1
       if (.not.log_mb(zero_row3(1,fb)+indx2-1)) then
         indx = i + (q-1)*nx(1)*nz(1)
         call zcopy(nz(1),tmp1(indx),nx(1),tmp2,1)
         call dcfftb(nz(1),tmp2,dcpl_mb(tmpz(1,1)))
         call zcopy(nz(1),tmp2,1,tmp1(indx),nx(1))
       end if
      end do
      end do
#endif

*     ********************************************
*     ***         Do a ptranspose of A          ***
*     ***      A(kx,ky,nz(1)) <- A(kx,nz(1),ky)      ***
*     ********************************************
      call nwpw_timing_end(31)
      call nwpw_timing_start(32)
      call C3dB_c_ptranspose1_jk_start(fb,tmp1,tmp2,tmp1,
     >                                 request,reqcnt,44)
      call nwpw_timing_end(32)


      !*************************
      !**** hilbert mapping ****
      !*************************
      else


      call nwpw_timing_start(31)
*     *************************************************
*     ***     do fft along kz dimension             ***
*     ***   A(nz(1),kx,ky) <- fft1d^(-1)[A(kz,kx,ky)]  ***
*     *************************************************
#ifdef MLIB
      indx = 1
      do q=1,nq3(1)
         if (.not.log_mb(zero_row3(1,fb)+q-1)) then
         call z1dfft(tmp1(indx),nz(1),dcpl_mb(tmpz(1,1)),-2,ierr)
         end if
         indx = indx + nz(1)
      end do
#else
      indx = 1
      do q=1,nq3(1)
         if (.not.log_mb(zero_row3(1,fb)+q-1)) then
         call dcfftb(nz(1),tmp1(indx),dcpl_mb(tmpz(1,1)))
         end if
         indx = indx + nz(1)
      end do
#endif

      call nwpw_timing_end(31)
      call nwpw_timing_start(32)
      call C3dB_c_ptranspose_ijk_start(fb,3,tmp1,tmp2,tmp1,
     >                                 request,reqcnt,45)
      call nwpw_timing_end(32)


      end if
    
      return
      end





*     ***********************************
*     *					*
*     *	        C3dB_pfftby		*
*     *					*
*     ***********************************

      subroutine C3dB_pfftby(nbb,tmp1,tmp2,request,reqcnt)
      implicit none
      integer nbb
      complex*16  tmp1(*)
      complex*16  tmp2(*)
      integer     request(*),reqcnt

#include "bafdecls.fh"
#include "errquit.fh"
#include "C3dB.fh"
#include "C3dB_pfft.fh"

      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / C3dB_fft / tmpx,tmpy,tmpz


*     *** local variables ***
      integer i,j,q,indx,ierr
      integer indx2,fb

      if (nbb.eq.0) then
         fb = 0
      else
         fb = 1
      end if
      
      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then


*     ********************************************
*     ***         Do a ptranspose of A          ***
*     ***      A(kx,ky,nz(1)) <- A(kx,nz(1),ky)      ***
*     ********************************************
      call nwpw_timing_start(32)
      call C3dB_c_ptranspose1_jk_end(fb,tmp2,tmp1,request,reqcnt)
      call nwpw_timing_end(32)
      call nwpw_timing_start(31)

*     *************************************************
*     ***     do fft along ky dimension             ***
*     ***   A(kx,ny(1),nz(1)) <- fft1d^(-1)[A(kx,ky,nz(1))]  ***
*     *************************************************
#ifdef MLIB
     
      indx2 = 0
      do q=1,nq(1)
      do i=1,nx(1)
        indx2 = indx2 + 1
        if (.not.log_mb(zero_row2(1,fb)+indx2-1)) then
          indx = i + (q-1)*nx(1)*nz(1)
          call zcopy(ny(1),tmp2(indx),nx(1),tmp1,1)
          call z1dfft(tmp1,ny(1),dcpl_mb(tmpy(1,1)),-2,ierr)
          call zcopy(ny(1),tmp1,1,tmp2(indx),nx(1))
        end if
      end do
      end do
#else
      indx2 = 0
      do q=1,nq(1)
      do i=1,nx(1)
        indx2 = indx2 + 1
        if (.not.log_mb(zero_row2(1,fb)+indx2-1)) then
         indx = i + (q-1)*nx(1)*nz(1)
         call zcopy(ny(1),tmp2(indx),nx(1),tmp1,1)
         call dcfftb(ny(1),tmp1,dcpl_mb(tmpy(1,1)))
         call zcopy(ny(1),tmp1,1,tmp2(indx),nx(1))
        end if
      end do
      end do
#endif

      call nwpw_timing_end(31)


      !*************************
      !**** hilbert mapping ****
      !*************************
      else

      call nwpw_timing_start(32)
      call C3dB_c_ptranspose_ijk_end(fb,3,tmp2,tmp1,request,reqcnt)

      call nwpw_timing_end(32)
      call nwpw_timing_start(31)
*     *************************************************
*     ***     do fft along ky dimension             ***
*     ***   A(ny(1),nz(1),kx) <- fft1d^(-1)[A(ky,nz(1),kx)]  ***
*     *************************************************
#ifdef MLIB
      indx = 1
      do q=1,nq2(1)
         if (.not.log_mb(zero_row2(1,fb)+q-1)) then
         call z1dfft(tmp2(indx),ny(1),dcpl_mb(tmpy(1,1)),-2,ierr)
         end if
         indx = indx + ny(1)
      end do
#else
      indx = 1
      do q=1,nq2(1)
         if (.not.log_mb(zero_row2(1,fb)+q-1)) then
         call dcfftb(ny(1),tmp2(indx),dcpl_mb(tmpy(1,1)))
         end if
         indx = indx + ny(1)
      end do
#endif

      call nwpw_timing_end(31)
      call nwpw_timing_start(32)
      call C3dB_c_ptranspose_ijk_start(fb,4,tmp2,tmp1,tmp2,
     >                                 request,reqcnt,46)
      call nwpw_timing_end(32)


      end if

      return
      end




*     ***********************************
*     *					*
*     *	        C3dB_pfftbx		*
*     *					*
*     ***********************************

      subroutine C3dB_pfftbx(nbb,tmp1,tmp2,request,reqcnt)
      implicit none
      integer nbb
      complex*16  tmp1(*)
      complex*16  tmp2(*)
      integer     request(*),reqcnt

#include "bafdecls.fh"
#include "errquit.fh"
#include "C3dB.fh"
#include "C3dB_pfft.fh"

      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / C3dB_fft / tmpx,tmpy,tmpz

*     *** local variables ***
      integer j,q,indx,ierr,fb


      if (nbb.eq.0) then
         fb = 0
      else
         fb = 1
      endif
      
      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then

      call nwpw_timing_start(31)
*     *************************************************
*     ***     do fft along kx dimension             ***
*     ***   A(nx(1),ny(1),nz(1)) <- fft1d^(-1)[A(kx,ny(1),nz(1))]  ***
*     *************************************************
#ifdef MLIB
      indx = 1
      do q=1,nq(1)
      do j=1,ny(1)
         call z1dfft(tmp2(indx),nx(1),dcpl_mb(tmpx(1,1)),-2,ierr)
         indx = indx + nx(1)
      end do
      end do

#else
      indx = 1
      do q=1,nq(1)
      do j=1,ny(1)
         call dcfftb(nx(1),tmp2(indx),dcpl_mb(tmpx(1,1)))
         indx = indx + nx(1)
      end do
      end do
#endif
      call dcopy(2*nx(1)*ny(1)*nq(1),tmp2,1,tmp1,1)



      call nwpw_timing_end(31)
*     *************************************************
      !*************************
      !**** hilbert mapping ****
      !*************************
      else

      call nwpw_timing_start(32)
      call C3dB_c_ptranspose_ijk_end(fb,4,tmp1,tmp2,request,reqcnt)
      call nwpw_timing_end(32)
      call nwpw_timing_start(31)

*     *************************************************
*     ***     do fft along kx dimension             ***
*     ***   A(nx(1),ny(1),nz(1)) <- fft1d^(-1)[A(kx,ny(1),nz(1))]  ***
*     *************************************************
#ifdef MLIB
      indx = 1
      do q=1,nq1(1)
         call z1dfft(tmp1(indx),nx(1),dcpl_mb(tmpx(1,1)),-2,ierr)
         indx = indx + nx(1)
      end do
#else
      indx = 1
      do q=1,nq1(1)
         call dcfftb(nx(1),tmp1(indx),dcpl_mb(tmpx(1,1)))
         indx = indx + nx(1)
      end do
#endif
      call nwpw_timing_end(31)

      end if
    
      return
      end


*     ***********************************
*     *                                 *
*     *         C3dB_pfftb_step         *
*     *                                 *
*     ***********************************
*

      subroutine C3dB_pfftb_step(step,nbb,A,tmp1,tmp2,request,reqcnt)
      implicit none
      integer step
      integer nbb
      complex*16 A(*)
      complex*16 tmp1(*)
      complex*16 tmp2(*)
      integer    request(*),reqcnt


#include "bafdecls.fh"


      if (step.eq.1) then
         call nwpw_timing_start(32)
         call Cram_c_unpack_start(nbb,A,tmp1,request,reqcnt,47)
         call nwpw_timing_end(32)
      else if (step.eq.2) then
         call nwpw_timing_start(32)
         call Cram_c_unpack_end(nbb,tmp1,tmp2,request,reqcnt)
         call nwpw_timing_end(32)
      else if (step.eq.3) then
         call C3dB_pfftbz(nbb,tmp1,tmp2,request,reqcnt)
      else if (step.eq.4) then
         call C3dB_pfftby(nbb,tmp1,tmp2,request,reqcnt)
      else if (step.eq.5) then
         call C3dB_pfftbx(nbb,tmp1,tmp2,request,reqcnt)
      end if

      return
      end


*     ***********************************
*     *                                 *
*     *     C3dB_pfft3_queue_init       *
*     *                                 *
*     ***********************************

      subroutine C3dB_pfft3_queue_init(qmax_in)
      implicit none
      integer qmax_in


#include "bafdecls.fh"
#include "errquit.fh"

#include "C3dB.fh"
#include "C3dB_pfft.fh"


      !**** pfft_queues ****
      integer aqindx(2),aqstatus(2)
      integer atmp(2),arequest(2),areqcnt(2),aqnbb(2)
      integer aqmax,aqsize,alast_index
      common / c_pfft_aqueue_common / aqindx,aqstatus,atmp,arequest,
     >                              areqcnt,aqnbb,
     >                              aqmax,aqsize,alast_index
      integer bqindx(2),bqstatus(2)
      integer btmp(2),brequest(2),breqcnt(2),bqnbb(2)
      integer bqmax,bqsize,blast_index
      common / c_pfft_bqueue_common / bqindx,bqstatus,btmp,brequest,
     >                              breqcnt,bqnbb,
     >                              bqmax,bqsize,blast_index
      

      logical value
      integer np,size

      call Parallel3d_np_i(np)

c     **** allocate aqueue ****
      aqmax       = qmax_in
      aqsize      = 0
      alast_index = aqmax
      size  = nfft3d(1)*aqmax*2
      value = BA_alloc_get(mt_dcpl,size,'atmp',atmp(2),atmp(1))
      size  = np*aqmax*2
      value = value.and.
     >      BA_alloc_get(mt_int,size,'arequest',arequest(2),arequest(1))
      size  = aqmax
      value = value.and.
     >        BA_alloc_get(mt_int,size,'aqindx',aqindx(2),aqindx(1))
      value = value.and.
     >      BA_alloc_get(mt_int,size,'aqstatus',aqstatus(2),aqstatus(1))
      value = value.and.
     >        BA_alloc_get(mt_int,size,'areqcnt',areqcnt(2),areqcnt(1))
      value = value.and.
     >        BA_alloc_get(mt_int,size,'aqnbb',aqnbb(2),aqnbb(1))


c     **** allocate bqueue ****
      bqmax       = qmax_in
      bqsize      = 0
      blast_index = bqmax
      size  = nfft3d(1)*bqmax*2
      value = BA_alloc_get(mt_dcpl,size,'btmp',btmp(2),btmp(1))
      size  = np*bqmax*2
      value = value.and.
     >      BA_alloc_get(mt_int,size,'brequest',brequest(2),brequest(1))
      size  = bqmax
      value = value.and.
     >        BA_alloc_get(mt_int,size,'bqindx',bqindx(2),bqindx(1))
      value = value.and.
     >      BA_alloc_get(mt_int,size,'bqstatus',bqstatus(2),bqstatus(1))
      value = value.and.
     >        BA_alloc_get(mt_int,size,'breqcnt',breqcnt(2),breqcnt(1))
      value = value.and.
     >        BA_alloc_get(mt_int,size,'bqnbb',bqnbb(2),bqnbb(1))

      return
      end




*     ***********************************
*     *                                 *
*     *     C3dB_pfft3_queue_end        *
*     *                                 *
*     ***********************************

      subroutine C3dB_pfft3_queue_end()
      implicit none


#include "bafdecls.fh"
#include "errquit.fh"

#include "C3dB.fh"
#include "C3dB_pfft.fh"



      !**** pfft_queues ****
      integer aqindx(2),aqstatus(2)
      integer atmp(2),arequest(2),areqcnt(2),aqnbb(2)
      integer aqmax,aqsize,alast_index
      common / c_pfft_aqueue_common / aqindx,aqstatus,atmp,arequest,
     >                              areqcnt,aqnbb,
     >                              aqmax,aqsize,alast_index
      integer bqindx(2),bqstatus(2)
      integer btmp(2),brequest(2),breqcnt(2),bqnbb(2)
      integer bqmax,bqsize,blast_index
      common / c_pfft_bqueue_common / bqindx,bqstatus,btmp,brequest,
     >                              breqcnt,bqnbb,
     >                              bqmax,bqsize,blast_index



      logical value

      value =           BA_free_heap(atmp(2))
      value = value.and.BA_free_heap(arequest(2))
      value = value.and.BA_free_heap(aqindx(2))
      value = value.and.BA_free_heap(aqstatus(2))
      value = value.and.BA_free_heap(areqcnt(2))
      value = value.and.BA_free_heap(aqnbb(2))

      value = value.and.BA_free_heap(btmp(2))
      value = value.and.BA_free_heap(brequest(2))
      value = value.and.BA_free_heap(bqindx(2))
      value = value.and.BA_free_heap(bqstatus(2))
      value = value.and.BA_free_heap(breqcnt(2))
      value = value.and.BA_free_heap(bqnbb(2))

      if (.not. value) call errquit('error freeing heap',0,MA_ERR)

      return
      end







*     ***********************************
*     *                                 *
*     *     C3dB_cr_pfft3_queue_filled  *
*     *                                 *
*     ***********************************

      logical function C3dB_cr_pfft3_queue_filled()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

#include "C3dB.fh"
#include "C3dB_pfft.fh"




      !**** pfft_queue ****
      integer qindx(2),qstatus(2)
      integer tmp(2),request(2),reqcnt(2),qnbb(2)
      integer qmax,qsize,last_index
      common / c_pfft_aqueue_common / qindx,qstatus,tmp,request,reqcnt,
     >                             qnbb,qmax,qsize,last_index

      C3dB_cr_pfft3_queue_filled = (qsize.ge.qmax)
      return
      end



*     ***********************************
*     *                                 *
*     *     C3dB_rc_pfft3_queue_filled  *
*     *                                 *
*     ***********************************

      logical function C3dB_rc_pfft3_queue_filled()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

#include "C3dB.fh"
#include "C3dB_pfft.fh"


      !**** pfft_queue ****
      integer qindx(2),qstatus(2)
      integer tmp(2),request(2),reqcnt(2),qnbb(2)
      integer qmax,qsize,last_index
      common / c_pfft_bqueue_common / qindx,qstatus,tmp,request,reqcnt,
     >                             qnbb,qmax,qsize,last_index

      C3dB_rc_pfft3_queue_filled = (qsize.ge.qmax)
      return
      end




*     ***********************************
*     *                                 *
*     *     C3dB_rc_pfft3f_queuein      *
*     *                                 *
*     ***********************************

      subroutine C3dB_rc_pfft3f_queuein(nbb,A)
      implicit none
      integer nbb
      complex*16 A(*)


#include "bafdecls.fh"
#include "errquit.fh"

#include "C3dB.fh"
#include "C3dB_pfft.fh"


      !**** pfft_queue ****
      integer qindx(2),qstatus(2)
      integer tmp(2),request(2),reqcnt(2),qnbb(2)
      integer qmax,qsize,last_index
      common / c_pfft_bqueue_common / qindx,qstatus,tmp,request,reqcnt,
     >                             qnbb,qmax,qsize,last_index

      !*** local variables ***
      integer shift1,shift2,shift3
      integer q,indx,status,np

      call nwpw_timing_start(30)
      call Parallel3d_np_i(np)

      do q=1,qsize
         indx   = int_mb(qindx(1)+q-1)
         int_mb(qstatus(1)+indx-1) = int_mb(qstatus(1)+indx-1) + 1
         status = int_mb(qstatus(1)+indx-1)
         shift1=nfft3d(1)*(2*(indx-1))
         shift2=nfft3d(1)*(2*(indx-1)+1)
         shift3=2*np*(indx-1)
         call C3dB_pfftf_step(status,
     >                        int_mb(qnbb(1)+indx-1),
     >                        A,
     >                        dcpl_mb(tmp(1)+shift1),
     >                        dcpl_mb(tmp(1)+shift2),
     >                        int_mb(request(1)+shift3),
     >                        int_mb(reqcnt(1)+indx-1))
      end do

      qsize = qsize + 1
      last_index = last_index+1
      if (last_index.gt.qmax) last_index = 1
      int_mb(qindx(1)+qsize-1)        = last_index
      int_mb(qstatus(1)+last_index-1) = 1
      status = 1
      int_mb(qnbb(1)+last_index-1) = nbb
      shift1=nfft3d(1)*(2*(last_index-1))
      shift2=nfft3d(1)*(2*(last_index-1)+1)
      shift3=2*np*(last_index-1)

      call C3dB_pfftf_step(status,nbb,A,
     >                     dcpl_mb(tmp(1)+shift1),
     >                     dcpl_mb(tmp(1)+shift2),
     >                     int_mb(request(1)+shift3),
     >                     int_mb(reqcnt(1)+last_index-1))


      call nwpw_timing_end(30)
      return
      end




*     ***********************************
*     *                                 *
*     *     C3dB_rc_pfft3f_queueout     *
*     *                                 *
*     ***********************************

      subroutine C3dB_rc_pfft3f_queueout(nbb,A)
      implicit none
      integer nbb
      complex*16 A(*)


#include "bafdecls.fh"
#include "errquit.fh"

#include "C3dB.fh"
#include "C3dB_pfft.fh"


      !**** pfft_queue ****
      integer qindx(2),qstatus(2)
      integer tmp(2),request(2),reqcnt(2),qnbb(2)
      integer qmax,qsize,last_index
      common / c_pfft_bqueue_common / qindx,qstatus,tmp,request,reqcnt,
     >                             qnbb,qmax,qsize,last_index

      !*** local variables ***
      integer shift1,shift2,shift3
      integer q,indx,indx1,status,np

      call nwpw_timing_start(30)
      call Parallel3d_np_i(np)

      indx1   = int_mb(qindx(1))
      do while (int_mb(qstatus(1)+indx1-1).lt.4)
         do q=1,qsize
            indx   = int_mb(qindx(1)+q-1)
            int_mb(qstatus(1)+indx-1) = int_mb(qstatus(1)+indx-1) + 1
            status = int_mb(qstatus(1)+indx-1)
            shift1=nfft3d(1)*(2*(indx-1))
            shift2=nfft3d(1)*(2*(indx-1)+1)
            shift3=2*np*(indx-1)
            call C3dB_pfftf_step(status,
     >                        int_mb(qnbb(1)+indx-1),
     >                        A,
     >                        dcpl_mb(tmp(1)+shift1),
     >                        dcpl_mb(tmp(1)+shift2),
     >                        int_mb(request(1)+shift3),
     >                        int_mb(reqcnt(1)+indx-1))
         end do
      end do

      qsize = qsize -1
      do q=1,qsize
        int_mb(qindx(1)+q-1) = int_mb(qindx(1)+q)
      end do

      shift2=nfft3d(1)*(2*(indx1-1)+1)
      call Cram_c_Copy(nbb,dcpl_mb(tmp(1)+shift2),A)

      call nwpw_timing_end(30)
      return
      end





*     ***********************************
*     *                                 *
*     *     C3dB_cr_pfft3b_queuein      *
*     *                                 *
*     ***********************************

      subroutine C3dB_cr_pfft3b_queuein(nbb,A)
      implicit none
      integer nbb
      complex*16 A(*)


#include "bafdecls.fh"
#include "errquit.fh"

#include "C3dB.fh"
#include "C3dB_pfft.fh"


      !**** pfft_queue ****
      integer qindx(2),qstatus(2)
      integer tmp(2),request(2),reqcnt(2),qnbb(2)
      integer qmax,qsize,last_index
      common / c_pfft_aqueue_common / qindx,qstatus,tmp,request,reqcnt,
     >                             qnbb,qmax,qsize,last_index

      !*** local variables ***
      integer shift1,shift2,shift3
      integer q,indx,status,np

      call nwpw_timing_start(30)
      call Parallel3d_np_i(np)

      do q=1,qsize
         indx   = int_mb(qindx(1)+q-1)
         int_mb(qstatus(1)+indx-1) = int_mb(qstatus(1)+indx-1) + 1
         status = int_mb(qstatus(1)+indx-1)
         shift1=nfft3d(1)*(2*(indx-1))
         shift2=nfft3d(1)*(2*(indx-1)+1)
         shift3=2*np*(indx-1)
         call C3dB_pfftb_step(status,
     >                        int_mb(qnbb(1)+indx-1),
     >                        A,
     >                        dcpl_mb(tmp(1)+shift1),
     >                        dcpl_mb(tmp(1)+shift2),
     >                        int_mb(request(1)+shift3),
     >                        int_mb(reqcnt(1)+indx-1))
      end do

      qsize = qsize + 1
      last_index = last_index+1
      if (last_index.gt.qmax) last_index = 1
      int_mb(qindx(1)+qsize-1)        = last_index
      int_mb(qstatus(1)+last_index-1) = 1
      status = 1
      int_mb(qnbb(1)+last_index-1) = nbb
      shift1=nfft3d(1)*(2*(last_index-1))
      shift2=nfft3d(1)*(2*(last_index-1)+1)
      shift3=2*np*(last_index-1)

      call C3dB_pfftb_step(status,nbb,A,
     >                     dcpl_mb(tmp(1)+shift1),
     >                     dcpl_mb(tmp(1)+shift2),
     >                     int_mb(request(1)+shift3),
     >                     int_mb(reqcnt(1)+last_index-1))

      call nwpw_timing_end(30)
      return
      end




*     ***********************************
*     *                                 *
*     *     C3dB_cr_pfft3b_queueout     *
*     *                                 *
*     ***********************************

      subroutine C3dB_cr_pfft3b_queueout(nbb,A)
      implicit none
      integer nbb
      complex*16 A(*)


#include "bafdecls.fh"
#include "errquit.fh"

#include "C3dB.fh"
#include "C3dB_pfft.fh"


      !**** pfft_queue ****
      integer qindx(2),qstatus(2)
      integer tmp(2),request(2),reqcnt(2),qnbb(2)
      integer qmax,qsize,last_index
      common / c_pfft_aqueue_common / qindx,qstatus,tmp,request,reqcnt,
     >                             qnbb,qmax,qsize,last_index

      !*** local variables ***
      integer shift1,shift2,shift3
      integer q,indx,indx1,status,np

      call Parallel3d_np_i(np)

      call nwpw_timing_start(30)

      indx1   = int_mb(qindx(1))
      do while (int_mb(qstatus(1)+indx1-1).lt.5)
         do q=1,qsize
            indx   = int_mb(qindx(1)+q-1)
            int_mb(qstatus(1)+indx-1) = int_mb(qstatus(1)+indx-1) + 1
            status = int_mb(qstatus(1)+indx-1)
            shift1=nfft3d(1)*(2*(indx-1))
            shift2=nfft3d(1)*(2*(indx-1)+1)
            shift3=2*np*(indx-1)
            call C3dB_pfftb_step(status,
     >                        int_mb(qnbb(1)+indx-1),
     >                        A,
     >                        dcpl_mb(tmp(1)+shift1),
     >                        dcpl_mb(tmp(1)+shift2),
     >                        int_mb(request(1)+shift3),
     >                        int_mb(reqcnt(1)+indx-1))
         end do
      end do

      qsize = qsize -1
      do q=1,qsize
        int_mb(qindx(1)+q-1) = int_mb(qindx(1)+q)
      end do

      shift1=nfft3d(1)*(2*(indx1-1))
      call dcopy(2*nfft3d(1),dcpl_mb(tmp(1)+shift1),1,A,1)

      call nwpw_timing_end(30)
      return
      end

