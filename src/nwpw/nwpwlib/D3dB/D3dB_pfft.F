
#define NBLOCKS 2

*
* $Id: D3dB_pfft.F,v 1.1 2005-06-24 19:52:38 bylaska Exp $
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

#include "mafdecls.fh"
#include "errquit.fh"
#include "D3dB.fh"
#include "D3dB_pfft.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      double precision eps
      parameter (eps=1d-12)

      logical value,yrow,zrow,yzslab
      integer nxh,nyh,nzh,q,p
      integer k1,k2,k3,ii,nb
      integer i,j,k
      integer index,index2
      double precision ggcut,g1,g2,g3,gg
      integer zero_arow3(2)

*     **** external functions ***
      real*8   lattice_ggcut,lattice_wggcut,lattice_unitg
      external lattice_ggcut,lattice_wggcut,lattice_unitg



      call Parallel_taskid(taskid)
      nxh = nx(1)/2
      nyh = ny(1)/2
      nzh = nz(1)/2

      if (mapping.eq.1) then
         
         value =           MA_alloc_get(mt_log,(nxh+1)*nq(1),
     >                     'zero_row3',zero_row3(2,0),zero_row3(1,0))
         value = value.and.MA_alloc_get(mt_log,(nxh+1)*nq(1),
     >                     'zero_row3',zero_row3(2,1),zero_row3(1,1))

         value = value.and.MA_alloc_get(mt_log,(nxh+1)*nq(1),
     >                     'zero_row2',zero_row2(2,0),zero_row2(1,0))
         value = value.and.MA_alloc_get(mt_log,(nxh+1)*nq(1),
     >                     'zero_row2',zero_row2(2,1),zero_row2(1,1))

         value = value.and.MA_alloc_get(mt_log,(nxh+1),
     >               'zero_slab23',zero_slab23(2,0),zero_slab23(1,0))
         value = value.and.MA_alloc_get(mt_log,(nxh+1),
     >               'zero_slab23',zero_slab23(2,1),zero_slab23(1,1))

         value = value.and.MA_alloc_get(mt_log,(nxh+1)*ny(1),
     >                     'zero_arow3',zero_arow3(2),zero_arow3(1))

         if (.not. value) call errquit('out of heap memory',0, MA_ERR)


         do nb=0,1



         if (nb.eq.0) then
            ggcut = lattice_ggcut()
         else
            ggcut = lattice_wggcut()
         end if

*        **** find zero_row3 - (i,j,*) rows that are zero ****
         do k2 = 1,ny(1)
         do k1 = 1,nxh+1
           log_mb(zero_arow3(1)+k1-1+(nxh+1)*(k2-1)) =.true.
           !log_mb(zero_arow3(1)+k1-1+(nxh+1)*(k2-1)) =.false.
           call D3dB_ktoqp(1,k2,q,p)
           if (p.eq.taskid) then
             index2 = k1 + (nxh+1)*(q-1)
             log_mb(zero_row3(1,nb)+index2-1) =.true.
           end if
         end do
         end do

         do k2 = -nyh+1, nyh-1
         do k1 = 0,nxh-1
          i=k1
          j=k2
          if (i .lt. 0) i = i + nx(1)
          if (j .lt. 0) j = j + ny(1)
          i=i+1
          j=j+1
          zrow = .true.
          do k3 = -nzh+1, nzh-1
            g1 = k1*lattice_unitg(1,1)
     >         + k2*lattice_unitg(1,2)
     >         + k3*lattice_unitg(1,3)
            g2 = k1*lattice_unitg(2,1)
     >         + k2*lattice_unitg(2,2)
     >         + k3*lattice_unitg(2,3)
            g3 = k1*lattice_unitg(3,1)
     >         + k2*lattice_unitg(3,2)
     >         + k3*lattice_unitg(3,3)
            gg = g1*g1 + g2*g2 + g3*g3
            gg= gg-ggcut
            if (gg.lt.-eps) zrow = .false.
          end do
          if (.not.zrow) then
            log_mb(zero_arow3(1)+i-1+(nxh+1)*(j-1)) =.false.
            call D3dB_ktoqp(1,j,q,p)
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
          if (i .lt. 0) i = i + nx(1)
          i=i+1
          yzslab = .true.
          do k3 = -nzh+1, nzh-1
          do k2 = -nyh+1, nyh-1
            g1 = k1*lattice_unitg(1,1)
     >         + k2*lattice_unitg(1,2)
     >         + k3*lattice_unitg(1,3)
            g2 = k1*lattice_unitg(2,1)
     >         + k2*lattice_unitg(2,2)
     >         + k3*lattice_unitg(2,3)
            g3 = k1*lattice_unitg(3,1)
     >         + k2*lattice_unitg(3,2)
     >         + k3*lattice_unitg(3,3)
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
         do k3 = 1,nz(1)
         do k1 = 1,nxh+1
           call D3dB_ktoqp(1,k3,q,p)
           if (p.eq.taskid) then
            index2 = k1 + (nxh+1)*(q-1)
            log_mb(zero_row2(1,nb)+index2-1) 
     >       = log_mb(zero_slab23(1,nb)+k1-1)
           end if
         end do
         end do



        
         end do
         value = MA_free_heap(zero_arow3(2))
         if (.not. value) call errquit('error freeing heap',0,MA_ERR)

      !*** mapping == 2 ***
      else
         value =           MA_alloc_get(mt_log,nq3(1),
     >                     'zero_row3',zero_row3(2,0),zero_row3(1,0))
         value = value.and.MA_alloc_get(mt_log,nq3(1),
     >                     'zero_row3',zero_row3(2,1),zero_row3(1,1))

         value = value.and.MA_alloc_get(mt_log,nq2(1),
     >                     'zero_row2',zero_row2(2,0),zero_row2(1,0))
         value = value.and.MA_alloc_get(mt_log,nq2(1),
     >                     'zero_row2',zero_row2(2,1),zero_row2(1,1))

         value = value.and.MA_alloc_get(mt_log,(nxh+1),
     >               'zero_slab23',zero_slab23(2,0),zero_slab23(1,0))
         value = value.and.MA_alloc_get(mt_log,(nxh+1),
     >               'zero_slab23',zero_slab23(2,1),zero_slab23(1,1))
         if (.not. value) call errquit('out of heap memory',0, MA_ERR)


         do nb=0,1
         if (nb.eq.0) then
            ggcut = lattice_ggcut()
         else
            ggcut = lattice_wggcut()
         end if

*        **** find zero_row3 - (i,j,*) rows that are zero ****
         do q = 1,nq3(1)
           log_mb(zero_row3(1,nb)+q-1) =.true.
         end do

         do k2 = -nyh+1, nyh-1
         do k1 = 0,nxh-1
          i=k1
          j=k2
          if (i .lt. 0) i = i + nx(1)
          if (j .lt. 0) j = j + ny(1)
          i=i+1
          j=j+1
          zrow = .true.
          do k3 = -nzh+1, nzh-1
            g1 = k1*lattice_unitg(1,1)
     >         + k2*lattice_unitg(1,2)
     >         + k3*lattice_unitg(1,3)
            g2 = k1*lattice_unitg(2,1)
     >         + k2*lattice_unitg(2,2)
     >         + k3*lattice_unitg(2,3)
            g3 = k1*lattice_unitg(3,1)
     >         + k2*lattice_unitg(3,2)
     >         + k3*lattice_unitg(3,3)
            gg = g1*g1 + g2*g2 + g3*g3
            gg= gg-ggcut
            if (gg.lt.-eps) zrow = .false.
          end do
          if (.not.zrow) then
            q = int_mb(q_map3(1,1)+(i-1)+(j-1)*(nx(1)/2+1))
            p = int_mb(p_map3(1,1)+(i-1)+(j-1)*(nx(1)/2+1))
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
          if (i .lt. 0) i = i + nx(1)
          i=i+1
          yzslab = .true.
          do k3 = -nzh+1, nzh-1
          do k2 = -nyh+1, nyh-1
            g1 = k1*lattice_unitg(1,1)
     >         + k2*lattice_unitg(1,2)
     >         + k3*lattice_unitg(1,3)
            g2 = k1*lattice_unitg(2,1)
     >         + k2*lattice_unitg(2,2)
     >         + k3*lattice_unitg(2,3)
            g3 = k1*lattice_unitg(3,1)
     >         + k2*lattice_unitg(3,2)
     >         + k3*lattice_unitg(3,3)
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
         do k = 1,nz(1)
         do i = 1,nxh+1
           q = int_mb(q_map2(1,1)+(k-1)+(i-1)*(nz(1)))
           p = int_mb(p_map2(1,1)+(k-1)+(i-1)*(nz(1)))

           if (p.eq.taskid) then
            log_mb(zero_row2(1,nb)+q-1)
     >       = log_mb(zero_slab23(1,nb)+i-1)
           end if
         end do
         end do



         end do


      end if

      return
      end




*     ***********************************
*     *                                 *
*     *         D3dB_pfft_end           *
*     *                                 *
*     ***********************************

      subroutine D3dB_pfft_end()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "D3dB_pfft.fh"



*     **** indexing variables ****
      integer iq_to_i1(2,0:1)
      integer iq_to_i2(2,0:1)
      integer iz_to_i2(2,0:1)
      integer i1_start(2,0:1)
      integer i2_start(2,0:1)
      common / ptrans_blk1 / iq_to_i1,iq_to_i2,iz_to_i2,
     >                       i1_start,i2_start

*     **** indexing variables ****
      integer jq_to_i1(2,0:1)
      integer jq_to_i2(2,0:1)
      integer jz_to_i2(2,0:1)
      integer j1_start(2,0:1)
      integer j2_start(2,0:1)
      common / ptrans_blk2 / jq_to_i1,jq_to_i2,jz_to_i2,
     >                       j1_start,j2_start


      logical value

      value =           MA_free_heap(zero_row2(2,0))
      value = value.and.MA_free_heap(zero_row2(2,1))
      value = value.and.MA_free_heap(zero_row3(2,0))
      value = value.and.MA_free_heap(zero_row3(2,1))
      value = value.and.MA_free_heap(zero_slab23(2,0))
      value = value.and.MA_free_heap(zero_slab23(2,1))

      value = value.and.MA_free_heap(iq_to_i1(2,0))
      value = value.and.MA_free_heap(iq_to_i1(2,1))
      value = value.and.MA_free_heap(iq_to_i2(2,0))
      value = value.and.MA_free_heap(iq_to_i2(2,1))
      value = value.and.MA_free_heap(iz_to_i2(2,0))
      value = value.and.MA_free_heap(iz_to_i2(2,1))
      value = value.and.MA_free_heap(i1_start(2,0))
      value = value.and.MA_free_heap(i1_start(2,1))
      value = value.and.MA_free_heap(i2_start(2,0))
      value = value.and.MA_free_heap(i2_start(2,1))

      value = value.and.MA_free_heap(jq_to_i1(2,0))
      value = value.and.MA_free_heap(jq_to_i1(2,1))
      value = value.and.MA_free_heap(jq_to_i2(2,0))
      value = value.and.MA_free_heap(jq_to_i2(2,1))
      value = value.and.MA_free_heap(jz_to_i2(2,0))
      value = value.and.MA_free_heap(jz_to_i2(2,1))
      value = value.and.MA_free_heap(j1_start(2,0))
      value = value.and.MA_free_heap(j1_start(2,1))
      value = value.and.MA_free_heap(j2_start(2,0))
      value = value.and.MA_free_heap(j2_start(2,1))

      if (.not. value) call errquit('error freeing heap',0,MA_ERR)

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
*       uses - D3dB_c_transpose_jk, dcopy           *
*                                                   *
*****************************************************

      implicit none
      integer nb,nbb
      complex*16  A(*)

#include "mafdecls.fh"
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


      call nwpw_timing_start(1)

*     ***** allocate temporary space ****
      !call D3dB_nfft3d(nb,nfft3d)
      value = MA_push_get(mt_dcpl,(nfft3d(nb)),'ffttmp2',
     >                    tmp2(2),tmp2(1))
      value = value.and.
     >      MA_push_get(mt_dbl,(n2ft3d(nb)),'ffttmp3',tmp3(2),tmp3(1))
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

       end if
      end do
      indx0 = indx0 + nxhz
      end do
#endif

*     ********************************************
*     ***         Do a transpose of A          ***
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
      indx = 1
      do q=1,nq3(nb)
         if (.not.log_mb(zero_row3(1,nbb)+q-1)) then
         call dcfftb(nz(nb),A(indx),dcpl_mb(tmpz(1,nb)))
         end if
         indx = indx + nz(nb)
      end do
#endif

      call D3dB_c_transpose_ijk(nb,3,A,dcpl_mb(tmp2(1)),dbl_mb(tmp3(1)))

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
      indx = 1
      do q=1,nq2(nb)
         if (.not.log_mb(zero_row2(1,nbb)+q-1)) then
         call dcfftb(ny(nb),A(indx),dcpl_mb(tmpy(1,nb)))
         end if
         indx = indx + ny(nb)
      end do
#endif

      call D3dB_c_transpose_ijk(nb,4,A,dcpl_mb(tmp2(1)),dbl_mb(tmp3(1)))

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
      call cshift1_fftb(nx(nb),nq1(nb),1,1,A)
      indx = 1
      do q=1,nq1(nb)
         !indx = 1 + (q-1)*(nx(nb)/2+1)
         call drfftb(nx(nb),A(indx),dcpl_mb(tmpx(1,nb)))
         indx = indx + nxh
      end do
      call zeroend_fftb(nx(nb),nq1(nb),1,1,A)
#endif



      end if
    
*     **** deallocate temporary space  ****
      value = MA_pop_stack(tmp3(2))
      value = MA_pop_stack(tmp2(2))
      !value = MA_pop_stack(tmp1(2))

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
*       uses - transpose1 subroutine                *
*                                                   *
*****************************************************

      implicit none
      integer nb,nbb
      complex*16  A(*)

#include "mafdecls.fh"
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


      call nwpw_timing_start(1)


*     ***** allocate temporary space ****
      !call D3dB_nfft3d(nb,nfft3d)
      value = MA_push_get(mt_dcpl,(nfft3d(nb)),'tmp2',tmp2(2),tmp2(1))
      value = value.and.
     >        MA_push_get(mt_dbl,(n2ft3d(nb)),'tmp3',tmp3(2),tmp3(1))
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

        end if
       end do
       indx0 = indx0 + nxhy
       end do

#endif


*     ********************************************
*     ***         Do a transpose of A          ***
*     ***      A(ky,nz(nb),ky) <- A(kx,ky,nz(nb))      ***
*     ********************************************
      call D3dB_c_ptranspose2_jk(nbb,A,dcpl_mb(tmp2(1)),dbl_mb(tmp3(1)))


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
      indx = 1
      do q=1,nq1(nb)
         !indx = 1 + (q-1)*(nx(nb)/2+1)
         call drfftf(nx(nb),A(indx),dcpl_mb(tmpx(1,nb)))
         indx = indx + nxh
      end do
      call cshift_fftf(nx(nb),nq1(nb),1,1,A)
#endif

      call D3dB_c_transpose_ijk(nb,1,A,dcpl_mb(tmp2(1)),dbl_mb(tmp3(1)))

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
      indx = 1
      do q=1,nq2(nb)
         if (.not.log_mb(zero_row2(1,nbb)+q-1)) then
         call dcfftf(ny(nb),A(indx),dcpl_mb(tmpy(1,nb)))
         end if
         indx = indx + ny(nb)
      end do
#endif

      call D3dB_c_transpose_ijk(nb,2,A,dcpl_mb(tmp2(1)),dbl_mb(tmp3(1)))

*     ********************************************
*     ***     do fft along ny(nb) dimension        ***
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
      indx = 1
      do q=1,nq3(nb)
         if (.not.log_mb(zero_row3(1,nbb)+q-1)) then
         call dcfftf(nz(nb),A(indx),dcpl_mb(tmpz(1,nb)))
         end if
         indx = indx + nz(nb)
      end do
#endif

      end if


*     **** deallocate temporary space  ****
      value = MA_pop_stack(tmp3(2))
      value = MA_pop_stack(tmp2(2))
      !value = MA_pop_stack(tmp1(2))

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

#include "mafdecls.fh"
#include "errquit.fh"
#include "D3dB.fh"


*     **** indexing variables ****
      integer iq_to_i1(2,0:1)
      integer iq_to_i2(2,0:1)
      integer iz_to_i2(2,0:1)
      integer i1_start(2,0:1)
      integer i2_start(2,0:1)
      common / ptrans_blk1 / iq_to_i1,iq_to_i2,iz_to_i2,
     >                       i1_start,i2_start

*     **** indexing variables ****
      integer jq_to_i1(2,0:1)
      integer jq_to_i2(2,0:1)
      integer jz_to_i2(2,0:1)
      integer j1_start(2,0:1)
      integer j2_start(2,0:1)
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
      integer i,j,k,it
      logical value



*     **** allocate ptrans_blk1 and ptrans_blk2 common block ****
      value = MA_alloc_get(mt_int,((nx(1)/2+1)*ny(1)*nq(1)),
     >                     'piq_to_i1',iq_to_i1(2,nbb),iq_to_i1(1,nbb))
      value=value.and.MA_alloc_get(mt_int,((nx(1)/2+1)*ny(1)*nq(1)),
     >                     'piq_to_i2',iq_to_i2(2,nbb),iq_to_i2(1,nbb))
      value=value.and.MA_alloc_get(mt_int,((nx(1)/2+1)*ny(1)*nq(1)),
     >                     'piz_to_i2',iz_to_i2(2,nbb),iz_to_i2(1,nbb))

      value = value.and.MA_alloc_get(mt_int,(nz(1)+1),
     >                     'pi1_start',i1_start(2,nbb),i1_start(1,nbb))
      value = value.and.MA_alloc_get(mt_int,(nz(1)+1),
     >                     'pi2_start',i2_start(2,nbb),i2_start(1,nbb))


      value = MA_alloc_get(mt_int,((nx(1)/2+1)*ny(1)*nq(1)),
     >                     'riq_to_i1',jq_to_i1(2,nbb),jq_to_i1(1,nbb))
      value=value.and.MA_alloc_get(mt_int,((nx(1)/2+1)*ny(1)*nq(1)),
     >                     'riq_to_i2',jq_to_i2(2,nbb),jq_to_i2(1,nbb))
      value=value.and.MA_alloc_get(mt_int,((nx(1)/2+1)*ny(1)*nq(1)),
     >                     'riz_to_i2',jz_to_i2(2,nbb),jz_to_i2(1,nbb))

      value = value.and.MA_alloc_get(mt_int,(nz(1)+1),
     >                     'ri1_start',j1_start(2,nbb),j1_start(1,nbb))
      value = value.and.MA_alloc_get(mt_int,(nz(1)+1),
     >                     'ri2_start',j2_start(2,nbb),j2_start(1,nbb))

      if (.not. value) 
     > call errquit('D3dB_ptranspose_jk_init:out of heap',0,MA_ERR)

      call Parallel_taskid(taskid)
      call Parallel_np(np)


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

         do k=1,nz(1)
         do j=1,ny(1)

*           **** packing scheme **** 
            call D3dB_ktoqp(1,k,qhere,phere)
            call D3dB_ktoqp(1,j,qto,pto)
            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
               do i=1,(nx(1)/2+1)
                  ii = i + (nx(1)/2+1)*(k-1)
                  jj = i + (nx(1)/2+1)*(j-1)
                  itmp = i + (j-1)*(nx(1)/2+1) 
     >                     + (qhere-1)*(nx(1)/2+1)*ny(1)

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
            call D3dB_ktoqp(1,j,qhere,phere)
            call D3dB_ktoqp(1,k,qfrom,pfrom)
            if ((phere.eq.taskid).and.(pfrom.eq.proc_from)) then
               do i=1,(nx(1)/2+1)
                  ii = i + (nx(1)/2+1)*(k-1)
                  jj = i + (nx(1)/2+1)*(j-1)
                  itmp = i + (k-1)*(nx(1)/2+1) 
     >                     + (qhere-1)*(nx(1)/2+1)*ny(1)
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


