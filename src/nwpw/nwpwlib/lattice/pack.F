

*     ***********************************
*     *					*
*     *	          Pack_Init		*	
*     *					*
*     ***********************************

      subroutine Pack_Init()
      implicit none


#include "mafdecls.fh"

*     **** common block for masker ****
      integer masker(2,0:1),nwave(0:1)
      common / mask_block / masker,nwave


*     **** common block for pack ****
      integer nida(0:1),nidb2(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb2
      integer nidb(0:1)
      common / pack2_blk / nidb

      integer nwave_all(0:1)
      common / pack_blk3x / nwave_all


*     **** local variables ****
      logical value
      integer taskid,nfft3d
      integer i,j,k
      integer k1,k2,k3
      integer q,p,indx
      integer nb
      integer nx,ny,nz,nq
      integer nxh,nyh,nzh      



*     **** allocate pack memory ****
      call D3dB_nfft3d(1,nfft3d)

      value = MA_alloc_get(mt_int,nfft3d,
     >                     'pack',pack(2,0),pack(1,0))
      value = value.and.
     >        MA_alloc_get(mt_int,nfft3d,
     >                     'pack',pack(2,1),pack(1,1))
      if (.not. value) call errquit('Pack_init: out of heap memory',0)


      call Parallel_taskid(taskid)
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      call D3dB_nq(1,nq)
      nxh = nx/2
      nyh = ny/2
      nzh = nz/2
      
      do nb=0,1
         nida(nb) = 0
         nidb(nb) = 0

*        **** k=(0,0,0)  ****
         k1=0
         k2=0
         k3=0
         call D3dB_ktoqp(1,k3+1,q,p)
         if (p.eq.taskid) then
            indx = (q-1)*(nxh+1)*ny + k2*(nxh+1) + k1 + 1
            if (.not.log_mb(masker(1,nb)+indx-1)) then
               nida(nb) = nida(nb) + 1
               int_mb(pack(1,nb)+nida(nb)-1) = indx
            end if
         end if
         
*        **** k=(0,0,k3) **** - neglect (0,0,-k3) points
         do k=1,(nzh-1)
            k1=0
            k2=0
            k3=k
            call D3dB_ktoqp(1,k3+1,q,p)
            if (p.eq.taskid) then
               indx = (q-1)*(nxh+1)*ny + k2*(nxh+1) + k1 + 1
               if (.not.log_mb(masker(1,nb)+indx-1)) then
                  nidb(nb) = nidb(nb) + 1
                  int_mb(pack(1,nb)+nida(nb)+nidb(nb)-1) = indx
               end if
            end if
         end do
   
*        **** k=(0,k2,k3) **** - neglect (0,-k2, -k3) points
         do k=(-nzh+1),(nzh-1)
         do j=1,(nyh-1)
            k1=0
            k2=j
            k3=k
            if (k3.lt.0) k3 = k3 + nz

            call D3dB_ktoqp(1,k3+1,q,p)
            if (p.eq.taskid) then
               indx = (q-1)*(nxh+1)*ny + k2*(nxh+1) + k1 + 1
               if (.not.log_mb(masker(1,nb)+indx-1)) then
                  nidb(nb) = nidb(nb) + 1
                  int_mb(pack(1,nb)+nida(nb)+nidb(nb)-1) = indx
               end if
            end if
         end do
         end do


*        **** k=(k1,k2,k3) **** 
         do k=(-nzh+1),(nzh-1)
         do j=(-nyh+1),(nyh-1)
         do i=1,(nxh-1)
            k1=i
            k2=j
            k3=k
            if (k2.lt.0) k2 = k2 + ny
            if (k3.lt.0) k3 = k3 + nz

            call D3dB_ktoqp(1,k3+1,q,p)
            if (p.eq.taskid) then
               indx = (q-1)*(nxh+1)*ny + k2*(nxh+1) + k1 + 1
               if (.not.log_mb(masker(1,nb)+indx-1)) then
                  nidb(nb) = nidb(nb) + 1
                  int_mb(pack(1,nb)+nida(nb)+nidb(nb)-1) = indx
               end if
            end if
         end do
         end do
         end do

         

      end do

      call Balance_Init(nida,nidb,nidb2)

      nwave_all(0) = nida(0) + nidb(0)
      nwave_all(1) = nida(1) + nidb(1)
      call D3dB_ISumAll(nwave_all(0))
      call D3dB_ISumAll(nwave_all(1)) 

       
      return 
      end

*     ***********************************
*     *					*
*     *	          Pack_end		*	
*     *					*
*     ***********************************

      subroutine Pack_end()
      implicit none

#include "mafdecls.fh"

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

      logical value

      value = MA_free_heap(pack(2,0))
      value = value.and.
     >        MA_free_heap(pack(2,1))
      if (.not. value) call errquit('error freeing heap memory',0)
      call Balance_End()
      return
      end



*     ***********************************
*     *					*
*     *	          Pack_c_pack		*	
*     *					*
*     ***********************************

      subroutine Pack_c_pack(nb,A)
      implicit none
      integer    nb
      complex*16 A(*)

#include "mafdecls.fh"

*     **** common block for pack ****
      integer nida(0:1),nidb2(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb2
      integer nidb(0:1)
      common / pack2_blk / nidb

      
*     **** local variables ****
      logical value
      integer i,nfft3d
      integer tmp1(2)

*     ***** timing variables *****
*     1 - fft_time
*     2 - dot_time
*     3 - ABC_time
*     4 - vxc_time
*     5 - vlocal_time
*     6 - vnl_time
*     7 - vh_time
*     8 - strfac_time
*     9 - maskC_time
      real*8 tim1,tim2
      real*8 times(50)
      common / timing / times

      call current_second(tim1)



      call D3dB_nfft3d(1,nfft3d)
      value = MA_push_get(mt_dcpl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not.value) call errquit('out of stack memory',0)

      call dcopy(2*nfft3d,A,1,dcpl_mb(tmp1(1)),1)
      call dcopy(2*nfft3d,0.0d0,0,A,1)

      do i=1,(nida(nb)+nidb(nb))
        A(i) = dcpl_mb(tmp1(1) + int_mb(pack(1,nb)+i-1)-1)
      end do

      value = MA_pop_stack(tmp1(2))
      if (.not.value) call errquit('error popping stack',0)


      call Pack_c_balance(nb,A)

      call current_second(tim2)
      times(9) = times(9) + (tim2-tim1)

      return
      end

*     ***********************************
*     *					*
*     *	          Pack_t_pack		*	
*     *					*
*     ***********************************

      subroutine Pack_t_pack(nb,A)
      implicit none
      integer nb
      real*8  A(*)

#include "mafdecls.fh"

*     **** common block for pack ****
      integer nida(0:1),nidb2(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb2
      integer nidb(0:1)
      common / pack2_blk / nidb

      
*     **** local variables ****
      logical value
      integer i,nfft3d
      integer tmp1(2)

*     ***** timing variables *****
*     1 - fft_time
*     2 - dot_time
*     3 - ABC_time
*     4 - vxc_time
*     5 - vlocal_time
*     6 - vnl_time
*     7 - vh_time
*     8 - strfac_time
*     9 - maskC_time
      real*8 tim1,tim2
      real*8 times(50)
      common / timing / times

      call current_second(tim1)


      call D3dB_nfft3d(1,nfft3d)
      value = MA_push_get(mt_dbl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not.value) call errquit('out of stack memory',0)

      call dcopy(nfft3d,A,1,dbl_mb(tmp1(1)),1)
      call dcopy(nfft3d,0.0d0,0,A,1)
      do i=1,(nida(nb)+nidb(nb))
        A(i) = dbl_mb(tmp1(1) + int_mb(pack(1,nb)+i-1) - 1)
      end do

      value = MA_pop_stack(tmp1(2))
      if (.not.value) call errquit('error popping stack',0)

      call Pack_t_balance(nb,A)

      call current_second(tim2)
      times(9) = times(9) + (tim2-tim1)

      return
      end


*     ***********************************
*     *					*
*     *	          Pack_c_unpack		*	
*     *					*
*     ***********************************

      subroutine Pack_c_unpack(nb,A)
      implicit none
      integer    nb
      complex*16 A(*)

#include "mafdecls.fh"

*     **** common block for pack ****
      integer nida(0:1),nidb2(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb2
      integer nidb(0:1)
      common / pack2_blk / nidb


*     **** local variables ****
      logical value
      integer i,nfft3d,ny,nq
      integer tmp1(2)
      integer tmp2(2)
c      integer indx1,indx2

*     ***** timing variables *****
*     1 - fft_time
*     2 - dot_time
*     3 - ABC_time
*     4 - vxc_time
*     5 - vlocal_time
*     6 - vnl_time
*     7 - vh_time
*     8 - strfac_time
*     9 - maskC_time
      real*8 tim1,tim2
      real*8 times(50)
      common / timing / times

      call current_second(tim1)


      call Pack_c_unbalance(nb,A)

      call D3dB_nfft3d(1,nfft3d)
      
      value = MA_push_get(mt_dcpl,(nida(nb)+nidb(nb)),
     >                    'tmp1',tmp1(2),tmp1(1))
      if (.not.value) call errquit('out of stack memory',0)

      call dcopy(2*(nida(nb)+nidb(nb)),A,1,dcpl_mb(tmp1(1)),1)
      call dcopy(2*nfft3d,0.0d0,0,A,1)
      do i=1,(nida(nb)+nidb(nb))
        A(int_mb(pack(1,nb)+i-1)) = dcpl_mb(tmp1(1)+i-1)
      end do
      value = MA_pop_stack(tmp1(2))
      if (.not.value) call errquit('error popping stack',0)


*     **** make the kx=0 plane complete **** 
      call D3dB_ny(1,ny)
      call D3dB_nq(1,nq)
      value = MA_push_get(mt_dcpl,(ny*nq),'tmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,(ny*nq),'tmp2',tmp2(2),tmp2(1))

      call D3dB_c_timereverse(1,A,dcpl_mb(tmp1(1)),dcpl_mb(tmp2(1)))

      value = MA_pop_stack(tmp2(2))
      value = value.and.
     >        MA_pop_stack(tmp1(2))
      if (.not.value) call errquit('error popping stack',0)

c     call D3dB_nx(1,nx)
c     call D3dB_nz(1,nz)
c     call Check_Real(nx,ny,nz,nfft3d,A)

      call current_second(tim2)
      times(9) = times(9) + (tim2-tim1)
  
      return
      end



*     ***********************************
*     *					*
*     *	         Pack_npack		*	
*     *					*
*     ***********************************

      subroutine Pack_npack(nb,npack)
      implicit none
      integer nb
      integer npack

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

      npack = nida(nb)+nidb(nb)
      return
      end

*     ***********************************
*     *					*
*     *	         Pack_nwave		*	
*     *					*
*     ***********************************

      integer function Pack_nwave(nb)
      implicit none
      integer nb

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

      Pack_nwave = nida(nb)+nidb(nb)
      return
      end

*     ***********************************
*     *					*
*     *	         Pack_nwave_all		*	
*     *					*
*     ***********************************

      integer function Pack_nwave_all(nb)
      implicit none
      integer nb

*     **** common block for pack ****
      integer nwave_all(0:1)
      common / pack_blk3x / nwave_all

      Pack_nwave_all = nwave_all(nb)
      return
      end


*     ***********************************
*     *					*
*     *	         Pack_zero		*	
*     *					*
*     ***********************************

      subroutine Pack_zero(nb,zero,pzero) 
      implicit none
      integer  nb
      integer zero,pzero

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

      integer qzero

*     *********************************************************
*     **** warning this routine assumes a specific packing ****
*     *********************************************************
*     index = (qzero-1)*(nx/2+1)*ny + (j-1)*(nx/2+1) + i 
      zero = 1
      call D3dB_ktoqp(1,1,qzero,pzero)

      return
      end



*     ***********************************
*     *					*
*     *	         Pack_cc_ndot		*	
*     *					*
*     ***********************************

      subroutine Pack_cc_ndot(nb,ne,A,B,sum)
      implicit none
      integer    nb
      integer    ne
      complex*16 A(*)
      complex*16 B(*)
      real*8     sum(*)


*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer n,np
      integer npack, shift

*     **** external functions ****
      real*8   ddot
      external ddot

*     ***** timing variables *****
*     1 - fft_time
*     2 - dot_time
*     3 - ABC_time
*     4 - vxc_time
*     5 - vlocal_time
*     6 - vnl_time
*     7 - vh_time
*     8 - strfac_time
*     9 - maskC_time
      real*8 tim1,tim2
      real*8 times(50)
      common / timing / times

      call current_second(tim1)


      call Parallel_np(np)
      npack= nida(nb)+nidb(nb)

      do n=1,ne
        shift = (n-1)*npack

        sum(n) = ddot(2*nida(nb),A(1+shift),1,B,1)
        sum(n) = sum(n) + 2*ddot(2*nidb(nb),
     >                           A(nida(nb)+1+shift),1,
     >                           B(nida(nb)+1),      1)
      end do

      if (np.gt.1) call D3dB_Vector_SumAll(ne,sum)


      call current_second(tim2)
      times(2) = times(2) + (tim2-tim1)

      return
      end 


*     ***********************************
*     *					*
*     *	         Pack_cc_idot		*	
*     *					*
*     ***********************************

      subroutine Pack_cc_idot(nb,A,B,sum)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      real*8     sum

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer np

*     **** external functions ****
      real*8   ddot
      external ddot

*     ***** timing variables *****
*     1 - fft_time
*     2 - dot_time
*     3 - ABC_time
*     4 - vxc_time
*     5 - vlocal_time
*     6 - vnl_time
*     7 - vh_time
*     8 - strfac_time
*     9 - maskC_time
      real*8 tim1,tim2
      real*8 times(50)
      common / timing / times

      call current_second(tim1)


      call Parallel_np(np)

      sum = ddot(2*nida(nb),A,1,B,1)
      sum = sum + 2*ddot(2*nidb(nb),A(nida(nb)+1),1,B(nida(nb)+1),1)

c     if (np.gt.1) call D3dB_SumAll(sum)

      call current_second(tim2)
      times(2) = times(2) + (tim2-tim1)

      return
      end 




*     ***********************************
*     *					*
*     *	         Pack_cc_dot		*	
*     *					*
*     ***********************************

      subroutine Pack_cc_dot(nb,A,B,sum)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      real*8     sum

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer np

*     **** external functions ****
      real*8   ddot
      external ddot

*     ***** timing variables *****
*     1 - fft_time
*     2 - dot_time
*     3 - ABC_time
*     4 - vxc_time
*     5 - vlocal_time
*     6 - vnl_time
*     7 - vh_time
*     8 - strfac_time
*     9 - maskC_time
      real*8 tim1,tim2
      real*8 times(50)
      common / timing / times

      call current_second(tim1)


      call Parallel_np(np)

      sum = ddot(2*nida(nb),A,1,B,1)
      sum = sum + 2*ddot(2*nidb(nb),A(nida(nb)+1),1,B(nida(nb)+1),1)

      if (np.gt.1) call D3dB_SumAll(sum)

      call current_second(tim2)
      times(2) = times(2) + (tim2-tim1)
      return
      end 


*     ***********************************
*     *					*
*     *	         Pack_tt_idot		*	
*     *					*
*     ***********************************

      subroutine Pack_tt_idot(nb,A,B,sum)
      implicit none
      integer    nb
      real*8  A(*)
      real*8  B(*)
      real*8  sum

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb


*     **** external functions ****
      real*8   ddot
      external ddot

*     ***** timing variables *****
*     1 - fft_time
*     2 - dot_time
*     3 - ABC_time
*     4 - vxc_time
*     5 - vlocal_time
*     6 - vnl_time
*     7 - vh_time
*     8 - strfac_time
*     9 - maskC_time
      real*8 tim1,tim2
      real*8 times(50)
      common / timing / times

      call current_second(tim1)



      sum =         ddot(nida(nb),A,1,B,1)
      sum = sum + 2*ddot(nidb(nb),A(nida(nb)+1),1,B(nida(nb)+1),1)


      call current_second(tim2)
      times(2) = times(2) + (tim2-tim1)
      return
      end 





*     ***********************************
*     *					*
*     *	         Pack_tt_dot		*	
*     *					*
*     ***********************************

      subroutine Pack_tt_dot(nb,A,B,sum)
      implicit none
      integer nb
      real*8  A(*)
      real*8  B(*)
      real*8  sum

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer np

*     **** external functions ****
      real*8   ddot
      external ddot

*     ***** timing variables *****
*     1 - fft_time
*     2 - dot_time
*     3 - ABC_time
*     4 - vxc_time
*     5 - vlocal_time
*     6 - vnl_time
*     7 - vh_time
*     8 - strfac_time
*     9 - maskC_time
      real*8 tim1,tim2
      real*8 times(50)
      common / timing / times

      call current_second(tim1)


      call Parallel_np(np)

      sum =         ddot(nida(nb),A,1,B,1)
      sum = sum + 2*ddot(nidb(nb),A(nida(nb)+1),1,B(nida(nb)+1),1)

      if (np.gt.1) call D3dB_SumAll(sum)

      call current_second(tim2)
      times(2) = times(2) + (tim2-tim1)
      return
      end 


*     ***********************************
*     *					*
*     *	         Pack_t_dsum		*	
*     *					*
*     ***********************************

      subroutine Pack_t_dsum(nb,A,sum)
      implicit none
      integer nb
      real*8  A(*)
      real*8  sum

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer np

*     **** external functions ****
      real*8   dsum
      external dsum

*     ***** timing variables *****
*     1 - fft_time
*     2 - dot_time
*     3 - ABC_time
*     4 - vxc_time
*     5 - vlocal_time
*     6 - vnl_time
*     7 - vh_time
*     8 - strfac_time
*     9 - maskC_time
      real*8 tim1,tim2
      real*8 times(50)
      common / timing / times

      call current_second(tim1)


      call Parallel_np(np)

      sum =         dsum(nida(nb),A,1)
      sum = sum + 2*dsum(nidb(nb),A(nida(nb)+1),1)

      if (np.gt.1) call D3dB_SumAll(sum)

      call current_second(tim2)
      times(2) = times(2) + (tim2-tim1)
      return
      end



*     ***********************************
*     *					*
*     *	         Pack_c_Copy		*	
*     *					*
*     ***********************************

      subroutine Pack_c_Copy(nb,A,B)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb


      call dcopy(2*(nida(nb)+nidb(nb)),A,1,B,1)

      return
      end

*     ***********************************
*     *					*
*     *	         Pack_t_Copy		*	
*     *					*
*     ***********************************

      subroutine Pack_t_Copy(nb,A,B)
      implicit none
      integer nb
      real*8 A(*)
      real*8 B(*)


*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb


      call dcopy((nida(nb)+nidb(nb)),A,1,B,1)

      return
      end

*     ***********************************
*     *					*
*     *	         Pack_c_Zero		*	
*     *					*
*     ***********************************

      subroutine Pack_c_Zero(nb,A)
      implicit none
      integer    nb
      complex*16 A(*)


*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb


      call dcopy(2*(nida(nb)+nidb(nb)),0.0d0,0,A,1)

      return
      end



*     ***********************************
*     *					*
*     *	         Pack_cc_Sum		*	
*     *					*
*     ***********************************

      subroutine Pack_cc_Sum(nb,A,B,C)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      complex*16 C(*)


*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

      do i=1,(nida(nb)+nidb(nb))
        C(i) = A(i) + B(i)
      end do

      return
      end

*     ***********************************
*     *					*
*     *	         Pack_tt_Sum		*	
*     *					*
*     ***********************************

      subroutine Pack_tt_Sum(nb,A,B,C)
      implicit none
      integer  nb
      real*8 A(*)
      real*8 B(*)
      real*8 C(*)

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

      do i=1,(nida(nb)+nidb(nb))
        C(i) = A(i) + B(i)
      end do

      return
      end



*     ***********************************
*     *					*
*     *	         Pack_cc_Sub		*	
*     *					*
*     ***********************************

      subroutine Pack_cc_Sub(nb,A,B,C)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      complex*16 C(*)


*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

      do i=1,(nida(nb)+nidb(nb))
        C(i) = A(i) - B(i)
      end do

      return
      end

*     ***********************************
*     *					*
*     *	         Pack_tt_Sub		*	
*     *					*
*     ***********************************

      subroutine Pack_tt_Sub(nb,A,B,C)
      implicit none
      integer    nb
      real*8 A(*)
      real*8 B(*)
      real*8 C(*)

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

      do i=1,(nida(nb)+nidb(nb))
        C(i) = A(i) - B(i)
      end do

      return
      end

*     ***********************************
*     *									*
*     *	         Pack_tt_Sqrt			*	
*     *									*
*     ***********************************

      subroutine Pack_tt_Sqrt(nb,A,C)
      implicit none
      integer  nb
      real*8   A(*)
      real*8   C(*)

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

      do i=1,(nida(nb)+nidb(nb))
        C(i) = dsqrt(A(i))
      end do

      return
      end



*     ***********************************
*     *					*
*     *	         Pack_ct_Sqr		*	
*     *					*
*     ***********************************

      subroutine Pack_ct_Sqr(nb,A,C)
      implicit none
      integer    nb
      complex*16 A(*)
      real*8     C(*)

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

      do i=1,(nida(nb)+nidb(nb))
        C(i) = dble(A(i))**2 + dimag(A(i))**2
      end do

      return
      end



*     ***********************************
*     *					*
*     *	         Pack_c_SMul		*	
*     *					*
*     ***********************************

      subroutine Pack_c_SMul(nb,alpha,A,C)
      implicit none
      integer    nb
      real*8     alpha
      complex*16 A(*)
      complex*16 C(*)


*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

      do i=1,(nida(nb)+nidb(nb))
        C(i) = alpha*A(i)
      end do

      return
      end

*     ***********************************
*     *					*
*     *	         Pack_c_ZMul		*	
*     *					*
*     ***********************************

      subroutine Pack_c_ZMul(nb,alpha,A,C)
      implicit none
      integer    nb
      complex*16 alpha
      complex*16 A(*)
      complex*16 C(*)


*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

      do i=1,(nida(nb)+nidb(nb))
        C(i) = alpha*A(i)
      end do

      return
      end



*     ***********************************
*     *					*
*     *	         Pack_t_SMul		*	
*     *					*
*     ***********************************

      subroutine Pack_t_SMul(nb,alpha,A,C)
      implicit none
      integer    nb
      real*8 alpha
      real*8 A(*)
      real*8 C(*)


*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

      do i=1,(nida(nb)+nidb(nb))
        C(i) = alpha*A(i)
      end do

      return
      end



*     ***********************************
*     *					*
*     *	         Pack_cc_daxpy		*	
*     *					*
*     ***********************************

      subroutine Pack_cc_daxpy(nb,alpha,A,B)
      implicit none
      integer    nb
      real*8     alpha
      complex*16 A(*)
      complex*16 B(*)


*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb


      call daxpy(2*(nida(nb)+nidb(nb)),alpha,A,1,B,1)

      return
      end



*     ***********************************
*     *					*
*     *	         Pack_tt_daxpy		*	
*     *					*
*     ***********************************

      subroutine Pack_tt_daxpy(nb,alpha,A,C)
      implicit none
      integer    nb
      real*8 alpha
      real*8 A(*)
      real*8 C(*)


*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

      do i=1,(nida(nb)+nidb(nb))
        C(i) = C(i) + alpha*A(i)
      end do

      return
      end



*     ***********************************
*     *					*
*     *	         Pack_cc_zaxpy		*	
*     *					*
*     ***********************************

      subroutine Pack_cc_zaxpy(nb,alpha,A,C)
      implicit none
      integer    nb
      complex*16 alpha
      complex*16 A(*)
      complex*16 C(*)


*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

      do i=1,(nida(nb)+nidb(nb))
        C(i) = C(i) + alpha*A(i)
      end do

      return
      end


*     ***********************************
*     *									*
*     *	         Pack_tc_Mul            *
*     *									*
*     ***********************************

      subroutine Pack_tc_Mul(nb,A,B,C)
      implicit none
      integer    nb
      real*8     A(*)
      complex*16 B(*)
      complex*16 C(*)


*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

      do i=1,(nida(nb)+nidb(nb))
        C(i) = A(i)*B(i)
      end do

      return
      end

*     ***********************************
*     *									*
*     *	         Pack_itc_Mul           *
*     *									*
*     ***********************************

      subroutine Pack_itc_Mul(nb,A,B,C)
      implicit none
      integer    nb
      real*8     A(*)
      complex*16 B(*)
      complex*16 C(*)


*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

      do i=1,(nida(nb)+nidb(nb))
        C(i) = A(i) * dcmplx(-dimag(B(i)),dble(B(i)))
      end do

      return
      end





*     ***********************************
*     *                 *
*     *          Pack_cc_Mul            *
*     *                 *
*     ***********************************

      subroutine Pack_cc_Mul(nb,A,B,C)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      complex*16 C(*)


*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

      do i=1,(nida(nb)+nidb(nb))
        C(i) = A(i)*B(i)
      end do

      return
      end

*     ***********************************
*     *									*
*     *	         Pack_cct_conjgMul		*	
*     *									*
*     ***********************************
      subroutine Pack_cct_conjgMul(nb,A,B,C)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      real*8     C(*)

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

      do i=1,(nida(nb)+nidb(nb))
        C(i) = dble(dconjg(A(i))*B(i))
      end do

      return
      end


*     ***********************************
*     *                				    *
*     *          Pack_cc_conjgMul       *
*     *                                 *
*     ***********************************

      subroutine Pack_cc_conjgMul(nb,A,B,C)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      complex*16 C(*)


*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

      do i=1,(nida(nb)+nidb(nb))
        C(i) = dconjg(A(i))*B(i)
      end do

      return
      end




*     ***********************************
*     *					*
*     *	         Pack_tc_iMul            *
*     *					*
*     ***********************************

      subroutine Pack_tc_iMul(nb,A,B,C)
      implicit none
      integer    nb
      real*8     A(*)
      complex*16 B(*)
      complex*16 C(*)


*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

      do i=1,(nida(nb)+nidb(nb))
        C(i) = dcmplx(0.0d0,A(i)) * B(i)
      end do

      return
      end



*     ***********************************
*     *					*
*     *	         Pack_tt_Mul            *
*     *					*
*     ***********************************

      subroutine Pack_tt_Mul(nb,A,B,C)
      implicit none
      integer    nb
      real*8     A(*)
      real*8     B(*)
      real*8     C(*)


*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

      do i=1,(nida(nb)+nidb(nb))
        C(i) = A(i)*B(i)
      end do

      return
      end

