*
* $Id: pack.F,v 1.33 2009-02-07 03:50:18 bylaska Exp $
*

*     ***********************************
*     *					*
*     *	          Pack_Init		*	
*     *					*
*     ***********************************

      subroutine Pack_Init()
      implicit none

#include "mafdecls.fh"
#include "mask_common.fh"
#include "errquit.fh"

*     **** common block for pack ****
      integer nida(0:1),nidb2(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb2
      integer nidb(0:1)
      common / pack2_blk / nidb

      integer nwave_all(0:1)
      common / pack_blk3x / nwave_all

*     **** common block for G_pack ****
      integer G_pack(3,0:1),G_pack_hndl
      common / G_pack_blk / G_pack,G_pack_hndl

*     **** local variables ****
      logical value
      integer taskid,nfft3d
      integer i,j,k
      integer k1,k2,k3
      integer q,p,indx
      integer nb,nsize,tmp(2)
      integer nx,ny,nz
      integer nxh,nyh,nzh      
      integer nwave_in(0:1),nwave_out(0:1)

*     **** external functions ****
      logical  control_balance
      integer  G_indx
      external control_balance
      external G_indx


*     **** allocate pack memory ****
      call D3dB_nfft3d(1,nfft3d)

      value = MA_alloc_get(mt_int,nfft3d,
     >                     'pack',pack(2,0),pack(1,0))
      value = value.and.
     >        MA_alloc_get(mt_int,nfft3d,
     >                     'pack',pack(2,1),pack(1,1))
      if (.not. value) 
     > call errquit('Pack_init: out of heap memory',0,MA_ERR)


      call Parallel2d_taskid_i(taskid)
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
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
         !call D3dB_ktoqp(1,k3+1,q,p)
         call D3dB_ijktoindexp(1,k1+1,k2+1,k3+1,indx,p)
         if (p.eq.taskid) then
            !indx = (q-1)*(nxh+1)*ny + k2*(nxh+1) + k1 + 1
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
            !call D3dB_ktoqp(1,k3+1,q,p)
            call D3dB_ijktoindexp(1,k1+1,k2+1,k3+1,indx,p)
            if (p.eq.taskid) then
               !indx = (q-1)*(nxh+1)*ny + k2*(nxh+1) + k1 + 1
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

            !call D3dB_ktoqp(1,k3+1,q,p)
            call D3dB_ijktoindexp(1,k1+1,k2+1,k3+1,indx,p)
            if (p.eq.taskid) then
               !indx = (q-1)*(nxh+1)*ny + k2*(nxh+1) + k1 + 1
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

            !call D3dB_ktoqp(1,k3+1,q,p)
            call D3dB_ijktoindexp(1,k1+1,k2+1,k3+1,indx,p)
            if (p.eq.taskid) then
               !indx = (q-1)*(nxh+1)*ny + k2*(nxh+1) + k1 + 1
               if (.not.log_mb(masker(1,nb)+indx-1)) then
                  nidb(nb) = nidb(nb) + 1
                  int_mb(pack(1,nb)+nida(nb)+nidb(nb)-1) = indx
               end if
            end if
         end do
         end do
         end do

         

      end do


      nwave_in(0) = nida(0) + nidb(0)
      nwave_in(1) = nida(1) + nidb(1)
      if (control_balance()) then
        call Balance_Init(2,nwave_in,nwave_out)
      else
        nwave_out(0) = nwave_in(0)
        nwave_out(1) = nwave_in(1)
      end if
      nidb2(0) = nidb(0) + (nwave_out(0)-nwave_in(0))
      nidb2(1) = nidb(1) + (nwave_out(1)-nwave_in(1))

      nwave_all(0) = nida(0) + nidb(0)
      nwave_all(1) = nida(1) + nidb(1)

      call D3dB_ISumAll(nwave_all(0))
      call D3dB_ISumAll(nwave_all(1)) 

c      call Parallel2d_taskid_i(taskid)
c      write(*,*) taskid," nida=",nida
c      write(*,*) taskid," nidb=",nidb,nidb2
c      write(*,*) " nida(0),nidb(0),nidb2(0):",taskid,nida(0),
c     >                                      nidb(0),nidb2(0)
c      write(*,*) " nida(1),nidb(1),nidb2(1):",taskid,nida(1),
c     >                                      nidb(1),nidb2(1)


*     **** allocate and define G_pack ****
      nsize = 3*((nida(0)+nidb2(0)) + (nida(1)+nidb2(1)))
      value = MA_push_get(mt_dbl,nfft3d,'tmp',tmp(2),tmp(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,nsize,
     >                     'G_pack',G_pack_hndl,G_pack(1,0))
      if (.not. value) 
     > call errquit('Pack_init:error pushing stack',0,MA_ERR)
      G_pack(2,0) = G_pack(1,0) + (nida(0)+nidb2(0))
      G_pack(3,0) = G_pack(2,0) + (nida(0)+nidb2(0))
      G_pack(1,1) = G_pack(3,0) + (nida(0)+nidb2(0))
      G_pack(2,1) = G_pack(1,1) + (nida(1)+nidb2(1))
      G_pack(3,1) = G_pack(2,1) + (nida(1)+nidb2(1))
      do nb=0,1
         do i=1,3
            call D3dB_t_Copy(1,dbl_mb(G_indx(i)),dbl_mb(tmp(1))) 
            call Pack_t_pack(nb,dbl_mb(tmp(1)))
            call Pack_t_Copy(nb,dbl_mb(tmp(1)),dbl_mb(G_pack(i,nb)))
         end do
      end do
      value = MA_pop_stack(tmp(2))
      if (.not. value) 
     > call errquit('Pack_init:error popping stack',0,MA_ERR)
      return 
      end

*     ***********************************
*     *                                 *
*     *           Pack_G_indx           *
*     *                                 *
*     ***********************************
      integer function Pack_G_indx(nb,i)
      implicit none
      integer nb,i

c     **** common block for G_pack ****
      integer G_pack(3,0:1),G_pack_hndl
      common / G_pack_blk / G_pack,G_pack_hndl

      Pack_G_indx = G_pack(i,nb)
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
#include "errquit.fh"

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** common block for G_pack ****
      integer G_pack(3,0:1),G_pack_hndl
      common / G_pack_blk / G_pack,G_pack_hndl

      logical value

*     **** external functions ****
      logical  control_balance
      external control_balance

      value =           MA_free_heap(pack(2,0))
      value = value.and.MA_free_heap(pack(2,1))
      value = value.and.MA_free_heap(G_pack_hndl)
      if (.not. value) 
     > call errquit('Pack_end:error freeing heap',0,MA_ERR)
      if (control_balance()) call Balance_End()
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
#include "errquit.fh"

*     **** common block for pack ****
      integer nida(0:1),nidb2(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb2
      integer nidb(0:1)
      common / pack2_blk / nidb

      
*     **** local variables ****
      logical value
      integer nfft3d,n
      integer tmp1(2)

*     **** external functions ****
      logical  control_balance
      external control_balance

      call nwpw_timing_start(9)

      call D3dB_nfft3d(1,nfft3d)
      value = MA_push_get(mt_dcpl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not.value) call errquit('out of stack memory',0, MA_ERR)

      call dcopy(2*nfft3d,A,1,dcpl_mb(tmp1(1)),1)
      call dcopy(2*nfft3d,0.0d0,0,A,1)

c      do i=1,(nida(nb)+nidb(nb))
c        A(i) = dcpl_mb(tmp1(1) + int_mb(pack(1,nb)+i-1)-1)
c      end do
      n = nida(nb)+nidb(nb)
      call Pack_c_indexcopy(n,int_mb(pack(1,nb)),dcpl_mb(tmp1(1)),A)

      value = MA_pop_stack(tmp1(2))
      if (.not.value) call errquit('error popping stack',0, MA_ERR)


      if (control_balance()) call Balance_c_balance(nb,A)

      call nwpw_timing_end(9)

      return
      end



*     ***********************************
*     *                                 *
*     *           Pack_c_mpack          *
*     *                                 *
*     ***********************************

      subroutine Pack_c_mpack(nb,m,A)
      implicit none
      integer    nb,m
      complex*16 A(*)

#include "mafdecls.fh"
#include "errquit.fh"

*     **** common block for pack ****
      integer nida(0:1),nidb2(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb2
      integer nidb(0:1)
      common / pack2_blk / nidb


*     **** local variables ****
      logical value
      integer i,j,nfft3d,n,shift
      integer tmp1(2)

*     **** external functions ****
      logical  control_balance
      external control_balance

      call nwpw_timing_start(9)


      call D3dB_nfft3d(1,nfft3d)
      value = MA_push_get(mt_dcpl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not.value) call errquit('out of stack memory',0, MA_ERR)

      shift=1
      n = nida(nb)+nidb(nb)
      do j=1,m
        call dcopy(2*nfft3d,A(shift),1,dcpl_mb(tmp1(1)),1)
        call dcopy(2*nfft3d,0.0d0,0,A(shift),1)
        call Pack_c_indexcopy(n,int_mb(pack(1,nb)),
     >                        dcpl_mb(tmp1(1)),
     >                        A(shift))
        if (control_balance()) call Balance_c_balance(nb,A(shift))
        shift = shift+nfft3d
      end do

      value = MA_pop_stack(tmp1(2))
      if (.not.value) call errquit('error popping stack',0, MA_ERR)

      call nwpw_timing_end(9)
      return
      end

      subroutine Pack_c_indexcopy(n,indx,A,B)
      implicit none
      integer n
      integer indx(*)
      complex*16 A(*), B(*)

      integer i
      do i=1,n
        B(i) = A(indx(i))
      end do
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
#include "errquit.fh"

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

*     **** external functions ****
      logical  control_balance
      external control_balance

      call nwpw_timing_start(9)

      call D3dB_nfft3d(1,nfft3d)
      value = MA_push_get(mt_dbl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not.value) call errquit('out of stack memory',0, MA_ERR)

      call dcopy(nfft3d,A,1,dbl_mb(tmp1(1)),1)
      call dcopy(nfft3d,0.0d0,0,A,1)
      do i=1,(nida(nb)+nidb(nb))
        A(i) = dbl_mb(tmp1(1) + int_mb(pack(1,nb)+i-1) - 1)
      end do

      value = MA_pop_stack(tmp1(2))
      if (.not.value) call errquit('error popping stack',0, MA_ERR)

      if (control_balance()) call Balance_t_balance(nb,A)

      call nwpw_timing_end(9)

      return
      end


*     ***********************************
*     *					*
*     *	          Pack_i_pack		*	
*     *					*
*     ***********************************

      subroutine Pack_i_pack(nb,A)
      implicit none
      integer    nb
      integer    A(*)

#include "mafdecls.fh"
#include "errquit.fh"

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

*     **** external functions ****
      logical  control_balance
      external control_balance

      call nwpw_timing_start(9)

      call D3dB_nfft3d(1,nfft3d)
      value = MA_push_get(mt_int,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not.value) call errquit('out of stack memory',0, MA_ERR)

      call icopy(nfft3d,A,1,int_mb(tmp1(1)),1)
      call icopy(nfft3d,0,0,A,1)

      do i=1,(nida(nb)+nidb(nb))
        A(i) = int_mb(tmp1(1) + int_mb(pack(1,nb)+i-1)-1)
      end do

      value = MA_pop_stack(tmp1(2))
      if (.not.value) call errquit('error popping stack',0, MA_ERR)


      if (control_balance()) call Balance_i_balance(nb,A)

      call nwpw_timing_end(9)

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
#include "errquit.fh"

*     **** common block for pack ****
      integer nida(0:1),nidb2(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb2
      integer nidb(0:1)
      common / pack2_blk / nidb


*     **** local variables ****
      logical value
      integer i,nfft3d,zplane_size
      integer tmp1(2)
      integer tmp2(2)
c      integer indx1,indx2

*     **** external functions ****
      logical  control_balance
      external control_balance


      call nwpw_timing_start(9)
      if (control_balance()) call Balance_c_unbalance(nb,A)

      call D3dB_nfft3d(1,nfft3d)
      
      value = MA_push_get(mt_dcpl,(nida(nb)+nidb(nb)),
     >                    'tmp1',tmp1(2),tmp1(1))
      if (.not.value) call errquit('out of stack memory',0, MA_ERR)

      call dcopy(2*(nida(nb)+nidb(nb)),A,1,dcpl_mb(tmp1(1)),1)
      call dcopy(2*nfft3d,0.0d0,0,A,1)
      do i=1,(nida(nb)+nidb(nb))
        A(int_mb(pack(1,nb)+i-1)) = dcpl_mb(tmp1(1)+i-1)
      end do
      value = MA_pop_stack(tmp1(2))
      if (.not.value) call errquit('error popping stack',0, MA_ERR)


*     **** make the kx=0 plane complete **** 
      !call D3dB_ny(1,ny)
      !call D3dB_nq(1,nq)
      call D3dB_zplane_size(1,zplane_size)
      value = MA_push_get(mt_dcpl,(zplane_size),'tmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,(zplane_size),'tmp2',tmp2(2),tmp2(1))

      call D3dB_c_timereverse(1,A,dcpl_mb(tmp1(1)),dcpl_mb(tmp2(1)))

      value = MA_pop_stack(tmp2(2))
      value = value.and.
     >        MA_pop_stack(tmp1(2))
      if (.not.value) call errquit('error popping stack',0, MA_ERR)

c     call D3dB_nx(1,nx)
c     call D3dB_nz(1,nz)
c     call Check_Real(nx,ny,nz,nfft3d,A)

      call nwpw_timing_end(9)
  
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
*     *                                 *
*     *          Pack_nida              *
*     *                                 *
*     ***********************************

      subroutine Pack_nida(nb,npack)
      implicit none
      integer nb
      integer npack

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

      npack = nida(nb)
      return
      end


*     ***********************************
*     *                                 *
*     *          Pack_nidb              *
*     *                                 *
*     ***********************************

      subroutine Pack_nidb(nb,npack)
      implicit none
      integer nb
      integer npack

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

      npack = nidb(nb)
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
      !call D3dB_ktoqp(1,1,qzero,pzero)
      call D3dB_ijktoindexp(1,1,1,1,zero,pzero)

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

      call nwpw_timing_start(2)
      call Parallel2d_np_i(np)
      npack= nida(nb)+nidb(nb)

      do n=1,ne
        shift = (n-1)*npack

        sum(n) = ddot(2*nida(nb),A(1+shift),1,B,1)
        sum(n) = sum(n) + 2*ddot(2*nidb(nb),
     >                           A(nida(nb)+1+shift),1,
     >                           B(nida(nb)+1),      1)
      end do

      if (np.gt.1) call D3dB_Vector_SumAll(ne,sum)

      call nwpw_timing_end(2)

      return
      end 

*     ***********************************
*     *					*
*     *	         Pack_ccm_sym_dot	*	
*     *					*
*     ***********************************

      subroutine Pack_ccm_sym_dot(nb,n,A,B,matrix)
      implicit none
      integer    nb,n
      complex*16 A(*)
      complex*16 B(*)
      real*8     matrix(n,n)

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer j,k
      integer np,npack,npack2


      call nwpw_timing_start(2)
      call Parallel2d_np_i(np)

      npack  = (nida(nb)+nidb(nb))
      npack2 = 2*npack

      do k=1,n
        call DGEMM('T','N',k,1,npack2,
     >             2.0d0,
     >             A,npack2, 
     >             B(1+(k-1)*npack),npack2, 
     >             0.0d0,
     >             matrix(1,k),k)
        call DGEMM('T','N',k,1,2*nida(nb),
     >             -1.0d0,
     >             A,npack2, 
     >             B(1+(k-1)*npack),npack2, 
     >             1.0d0,
     >             matrix(1,k),k)
      end do

      do k=1,n
      do j=k+1,n
        matrix(j,k) = matrix(k,j)
      end do
      end do

      if (np.gt.1) call D3dB_Vector_SumAll(n*n,matrix)

      call nwpw_timing_end(2)
      return
      end 


*     ***********************************
*     *                                 *
*     *      Pack_ccm_combo_sym_dot     *
*     *                                 *
*     ***********************************

      subroutine Pack_ccm_combo_sym_dot(nb,n,A,B,matrix)
      implicit none
      integer    nb,n
      complex*16 A(*)
      complex*16 B(*)
      real*8     matrix(n,n,3)

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i,j,k
      integer np,npack,npack2



      call nwpw_timing_start(2)
      call Parallel2d_np_i(np)

      npack  = (nida(nb)+nidb(nb))
      npack2 = 2*npack

      do k=1,n
        call DGEMM('T','N',k,1,npack2,
     >             2.0d0,
     >             A,npack2,
     >             A(1+(k-1)*npack),npack2,
     >             0.0d0,
     >             matrix(1,k,1),k)
        call DGEMM('T','N',k,1,npack2,
     >             2.0d0,
     >             A,npack2,
     >             B(1+(k-1)*npack),npack2,
     >             0.0d0,
     >             matrix(1,k,2),k)
        call DGEMM('T','N',k,1,npack2,
     >             2.0d0,
     >             B,npack2,
     >             B(1+(k-1)*npack),npack2,
     >             0.0d0,
     >             matrix(1,k,3),k)

        call DGEMM('T','N',k,1,2*nida(nb),
     >             -1.0d0,
     >             A,npack2,
     >             A(1+(k-1)*npack),npack2,
     >             1.0d0,
     >             matrix(1,k,1),k)
        call DGEMM('T','N',k,1,2*nida(nb),
     >             -1.0d0,
     >             A,npack2,
     >             B(1+(k-1)*npack),npack2,
     >             1.0d0,
     >             matrix(1,k,2),k)
        call DGEMM('T','N',k,1,2*nida(nb),
     >             -1.0d0,
     >             B,npack2,
     >             B(1+(k-1)*npack),npack2,
     >             1.0d0,
     >             matrix(1,k,3),k)
      end do

      do i=1,3
      do k=1,n
      do j=k+1,n
        matrix(j,k,i) = matrix(k,j,i)
      end do
      end do
      end do

      if (np.gt.1) call D3dB_Vector_SumAll(3*n*n,matrix)

      call nwpw_timing_end(2)
      return
      end




*     ***********************************
*     *					*
*     *	         Pack_ccm_sym_dot2	*	
*     *					*
*     ***********************************

      subroutine Pack_ccm_sym_dot2(nb,m,n,A,B,matrix)
      implicit none
      integer    nb,m,n
      complex*16 A(*)
      complex*16 B(*)
      real*8     matrix(m,m)

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer j,k
      integer np,npack,npack2


      call nwpw_timing_start(2)
      call Parallel2d_np_i(np)

      npack  = (nida(nb)+nidb(nb))
      npack2 = 2*npack

      do k=1,n
        call DGEMM('T','N',k,1,npack2,
     >             2.0d0,
     >             A,npack2, 
     >             B(1+(k-1)*npack),npack2, 
     >             0.0d0,
     >             matrix(1,k),k)
        call DGEMM('T','N',k,1,2*nida(nb),
     >             -1.0d0,
     >             A,npack2, 
     >             B(1+(k-1)*npack),npack2, 
     >             1.0d0,
     >             matrix(1,k),k)
      end do

      do k=1,n
      do j=k+1,n
        matrix(j,k) = matrix(k,j)
      end do
      end do

      if (np.gt.1) call D3dB_Vector_SumAll(m*m,matrix)

      call nwpw_timing_end(2)
      return
      end 



*     ***********************************
*     *					*
*     *	         Pack_ccm_dot		*	
*     *					*
*     ***********************************

      subroutine Pack_ccm_dot(nb,n,A,B,matrix)
      implicit none
      integer    nb,n
      complex*16 A(*)
      complex*16 B(*)
      real*8     matrix(n,n)

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer np,npack


      call nwpw_timing_start(2)
      call Parallel2d_np_i(np)

      npack = 2*(nida(nb)+nidb(nb))

        call DGEMM('T','N',n,n,npack,
     >             2.0d0,
     >             A,npack, 
     >             B,npack, 
     >             0.0d0,
     >             matrix,n)
        call DGEMM('T','N',n,n,2*nida(nb),
     >             -1.0d0,
     >             A,npack, 
     >             B,npack, 
     >             1.0d0,
     >             matrix,n)

      if (np.gt.1) call D3dB_Vector_SumAll(n*n,matrix)

      call nwpw_timing_end(2)
      return
      end 

*     ***********************************
*     *					*
*     *	         Pack_ccmn_dot		*	
*     *					*
*     ***********************************

      subroutine Pack_ccmn_dot(nb,m,n,A,B,matrix)
      implicit none
      integer    nb,m,n
      complex*16 A(*)
      complex*16 B(*)
      real*8     matrix(m,n)

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer np,npack


      call nwpw_timing_start(2)
      call Parallel2d_np_i(np)

      npack = 2*(nida(nb)+nidb(nb))

        call DGEMM('T','N',m,n,npack,
     >             2.0d0,
     >             A,npack, 
     >             B,npack, 
     >             0.0d0,
     >             matrix,m)
        call DGEMM('T','N',m,n,2*nida(nb),
     >             -1.0d0,
     >             A,npack, 
     >             B,npack, 
     >             1.0d0,
     >             matrix,m)

      if (np.gt.1) call D3dB_Vector_SumAll(m*n,matrix)

      call nwpw_timing_end(2)
      return
      end 




*     ***********************************
*     *					*
*     *	         Pack_cc_indot		*	
*     *					*
*     ***********************************

      subroutine Pack_cc_indot(nb,ne,A,B,sum)
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
c     integer npack, shift
      integer npack2

*     **** external functions ****
c     real*8   ddot
c     external ddot

      call nwpw_timing_start(2)
      call Parallel2d_np_i(np)
      npack2= 2*(nida(nb)+nidb(nb))

c     do n=1,ne
c       shift = (n-1)*npack
c
c       sum(n) = ddot(2*nida(nb),A(1+shift),1,B,1)
c       sum(n) = sum(n) + 2*ddot(2*nidb(nb),
c    >                           A(nida(nb)+1+shift),1,
c    >                           B(nida(nb)+1),      1)
c     end do

      call DGEMM('T','N',ne,1,(npack2),
     >           (2.0d0),
     >           A,npack2,
     >           B,npack2,
     >           (0.0d0),
     >           sum,ne)
      call DGEMM('T','N',ne,1,(2*nida(nb)),
     >           (-1.0d0),
     >           A,npack2,
     >           B,npack2,
     >           (1.0d0),
     >           sum,ne)

   

      call nwpw_timing_end(2)
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

      call nwpw_timing_start(2)
c      call Parallel2d_np_i(np)

      sum = ddot(2*nida(nb),A,1,B,1)
      sum = sum + 2*ddot(2*nidb(nb),A(nida(nb)+1),1,B(nida(nb)+1),1)

c     if (np.gt.1) call D3dB_SumAll(sum)

      call nwpw_timing_end(2)

      return
      end 




*     ***********************************
*     *                                 *
*     *          Pack_cc_izdot          *
*     *                                 *
*     ***********************************

      subroutine Pack_cc_izdot(nb,A,B,sum)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      complex*16 sum


*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb


*     **** external functions ****
      complex*16 tzdotc
      external   tzdotc

      call nwpw_timing_start(2)

      sum =         tzdotc(nida(nb),A,1,B,1)
      sum = sum + 2*tzdotc(nidb(nb),A(nida(nb)+1),1,B(nida(nb)+1),1)

      call nwpw_timing_end(2)

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

      call nwpw_timing_start(2)
      call Parallel2d_np_i(np)

      sum = ddot(2*nida(nb),A,1,B,1)
      sum = sum + 2*ddot(2*nidb(nb),A(nida(nb)+1),1,B(nida(nb)+1),1)

      if (np.gt.1) call D3dB_SumAll(sum)

      call nwpw_timing_end(2)
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

      call nwpw_timing_start(2)

      sum =         ddot(nida(nb),A,1,B,1)
      sum = sum + 2*ddot(nidb(nb),A(nida(nb)+1),1,B(nida(nb)+1),1)

      call nwpw_timing_end(2)
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

      call nwpw_timing_start(2)
      call Parallel2d_np_i(np)

      sum =         ddot(nida(nb),A,1,B,1)
      sum = sum + 2*ddot(nidb(nb),A(nida(nb)+1),1,B(nida(nb)+1),1)

      if (np.gt.1) call D3dB_SumAll(sum)

      call nwpw_timing_end(2)
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


      call nwpw_timing_start(2)
      call Parallel2d_np_i(np)

      sum =         dsum(nida(nb),A,1)
      sum = sum + 2*dsum(nidb(nb),A(nida(nb)+1),1)

      if (np.gt.1) call D3dB_SumAll(sum)

      call nwpw_timing_end(2)
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
*     *                                 *
*     *          Pack_cc_Sum2           *
*     *                                 *
*     ***********************************

      subroutine Pack_cc_Sum2(nb,A,B)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)


*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

      do i=1,(nida(nb)+nidb(nb))
        B(i) = B(i) + A(i)
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
*     *                                 *
*     *          Pack_tt_Sum2           *
*     *                                 *
*     ***********************************

      subroutine Pack_tt_Sum2(nb,A,B)
      implicit none
      integer  nb
      real*8 A(*)
      real*8 B(*)

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

      do i=1,(nida(nb)+nidb(nb))
        B(i) = B(i) + A(i)
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
*     *					*
*     *	         Pack_tt_Sqrt		*	
*     *					*
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
*     *                                 *
*     *          Pack_tt_Sqrt1          *
*     *                                 *
*     ***********************************

      subroutine Pack_tt_Sqrt1(nb,A)
      implicit none
      integer  nb
      real*8   A(*)

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

      do i=1,(nida(nb)+nidb(nb))
        A(i) = dsqrt(A(i))
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
*     *                                 *
*     *          Pack_c_SMul1           *
*     *                                 *
*     ***********************************

      subroutine Pack_c_SMul1(nb,alpha,A)
      implicit none
      integer    nb
      real*8     alpha
      complex*16 A(*)

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

      do i=1,(nida(nb)+nidb(nb))
        A(i) = alpha*A(i)
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
*     *                                 *
*     *          Pack_t_SMul1           *
*     *                                 *
*     ***********************************

      subroutine Pack_t_SMul1(nb,alpha,A)
      implicit none
      integer    nb
      real*8 alpha
      real*8 A(*)

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

      do i=1,(nida(nb)+nidb(nb))
        A(i) = alpha*A(i)
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

c      do i=1,(nida(nb)+nidb(nb))
c        C(i) = C(i) + alpha*A(i)
c      end do
      call daxpy((nida(nb)+nidb(nb)),alpha,A,1,C,1)

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
*     *					*
*     *	         Pack_tc_Mul            *
*     *					*
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
        C(i) = B(i)*A(i)
      end do

      return
      end


*     ***********************************
*     *                                 *
*     *          Pack_tc_Mul2           *
*     *                                 *
*     ***********************************

      subroutine Pack_tc_Mul2(nb,A,B)
      implicit none
      integer    nb
      real*8     A(*)
      complex*16 B(*)

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

      do i=1,(nida(nb)+nidb(nb))
        B(i) = B(i)*A(i)
      end do

      return
      end



*     ***********************************
*     *                                 *
*     *          Pack_tc_aMul           *
*     *                                 *
*     ***********************************

      subroutine Pack_tc_aMul(nb,alpha,A,B,C)
      implicit none
      integer    nb
      real*8     alpha
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
        C(i) = alpha*B(i)*A(i)
      end do

      return
      end

*     ***********************************
*     *                                 *
*     *          Pack_tc_aMulAdd        *
*     *                                 *
*     ***********************************

      subroutine Pack_tc_aMulAdd(nb,alpha,A,B,C)
      implicit none
      integer    nb
      real*8     alpha
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
        C(i) = C(i) + alpha*B(i)*A(i)
      end do

      return
      end


*     ***********************************
*     *					*
*     *	         Pack_tc_MulAdd         *
*     *					*
*     ***********************************

      subroutine Pack_tc_MulAdd(nb,A,B,C)
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
        C(i) = C(i) + A(i)*B(i)
      end do

      return
      end



*     ***********************************
*     *					*
*     *	         Pack_itc_Mul           *
*     *					*
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
*     *                                 *
*     *          Pack_cc_Mul            *
*     *                                 *
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
*     *                                 *
*     *          Pack_cc_Mul2           *
*     *                                 *
*     ***********************************

      subroutine Pack_cc_Mul2(nb,A,B)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

      do i=1,(nida(nb)+nidb(nb))
        B(i) = B(i)*A(i)
      end do

      return
      end


*     ***********************************
*     *					*
*     *	         Pack_cct_conjgMul	*	
*     *					*
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
*     *                		        *
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
*     *	         Pack_tc_iMul           *
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

*     ***********************************
*     *                                 *
*     *          Pack_tt_Mul2           *
*     *                                 *
*     ***********************************

      subroutine Pack_tt_Mul2(nb,A,B)
      implicit none
      integer    nb
      real*8     A(*)
      real*8     B(*)

*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

      do i=1,(nida(nb)+nidb(nb))
        B(i) = B(i)*A(i)
      end do

      return
      end



*     ***********************************
*     *                                 *
*     *          Pack_c_setzero         *
*     *                                 *
*     ***********************************

      subroutine Pack_c_setzero(nb,vzero,A)
      implicit none
      integer    nb
      real*8     vzero
      complex*16 A(*)


*     **** common block for pack ****
      integer nida(0:1),nidb(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer indx,p,q,taskid

      call Parallel2d_taskid_i(taskid)
      !call D3dB_ktoqp(1,1,q,p)
      call D3dB_ijktoindexp(1,1,1,1,indx,p)
      if (p.eq.taskid) A(indx) = dcmplx(vzero,0.0d0)

      return
      end






*     ***********************************
*     *					*
*     *	          Pack_c_pack_start	*	
*     *					*
*     ***********************************

      subroutine Pack_c_pack_start(nb,A,tmp1,request,reqcnt,msgtype)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 tmp1(*)
      integer    request(*),reqcnt
      integer    msgtype

#include "mafdecls.fh"
#include "errquit.fh"

*     **** common block for pack ****
      integer nida(0:1),nidb2(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb2
      integer nidb(0:1)
      common / pack2_blk / nidb

      
*     **** local variables ****
      logical value
      integer nfft3d,n

*     **** external functions ****
      logical  control_balance
      external control_balance

      call nwpw_timing_start(9)

      call D3dB_nfft3d(1,nfft3d)

      call dcopy(2*nfft3d,A,1,tmp1,1)
      call dcopy(2*nfft3d,0.0d0,0,A,1)

      n = nida(nb)+nidb(nb)
      call Pack_c_indexcopy(n,int_mb(pack(1,nb)),tmp1,A)


      if (control_balance()) 
     >  call Balance_c_balance_start(nb,A,request,reqcnt,msgtype)

      call nwpw_timing_end(9)
      return
      end






*     ***********************************
*     *                                 *
*     *           Pack_c_pack_end       *
*     *                                 *
*     ***********************************

      subroutine Pack_c_pack_end(nb,tmp1,request,reqcnt)
      implicit none
      integer    nb
      complex*16 tmp1(*)
      integer    request(*),reqcnt

#include "mafdecls.fh"
#include "errquit.fh"

*     **** common block for pack ****
      integer nida(0:1),nidb2(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb2
      integer nidb(0:1)
      common / pack2_blk / nidb


*     **** local variables ****
      logical value
      integer nfft3d,n

*     **** external functions ****
      logical  control_balance
      external control_balance

      call nwpw_timing_start(9)


      if (control_balance()) 
     > call Balance_c_balance_end(nb,tmp1,request,reqcnt)

      call nwpw_timing_end(9)
      return
      end






*     ***********************************
*     *					*
*     *	     Pack_c_unpack_start	*	
*     *					*
*     ***********************************

      subroutine Pack_c_unpack_start(nb,A,tmp1,request,reqcnt,msgtype)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 tmp1(*)
      integer    request(*),reqcnt
      integer    msgtype

#include "mafdecls.fh"
#include "errquit.fh"

*     **** common block for pack ****
      integer nida(0:1),nidb2(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb2
      integer nidb(0:1)
      common / pack2_blk / nidb

*     **** external functions ****
      logical  control_balance
      external control_balance

      call nwpw_timing_start(9)
      call Pack_c_Copy(nb,A,tmp1)
      if (control_balance())
     >   call Balance_c_unbalance_start(nb,tmp1,request,reqcnt,msgtype)

      call nwpw_timing_end(9)
  
      return
      end




*     ***********************************
*     *					*
*     *	          Pack_c_unpack_mid	*	
*     *					*
*     ***********************************

      subroutine Pack_c_unpack_mid(nb,tmp1,tmp2,tmp1z,tmp2z,
     >                         request,reqcnt,msgtype)
      implicit none
      integer    nb
      complex*16 tmp1(*)
      complex*16 tmp2(*)
      complex*16 tmp1z(*)
      complex*16 tmp2z(*)
      integer    request(*),reqcnt
      integer    msgtype

#include "mafdecls.fh"
#include "errquit.fh"

*     **** common block for pack ****
      integer nida(0:1),nidb2(0:1)
      integer pack(2,0:1)
      common / pack_blk / pack,nida,nidb2
      integer nidb(0:1)
      common / pack2_blk / nidb


*     **** local variables ****
      logical value
      integer i,nfft3d

*     **** external functions ****
      logical  control_balance
      external control_balance


      call nwpw_timing_start(9)
      if (control_balance()) 
     >  call Balance_c_unbalance_end(nb,tmp1,request,reqcnt)

      
      call D3dB_nfft3d(1,nfft3d)

      call dcopy(2*(nida(nb)+nidb(nb)),tmp1,1,tmp2,1)
      call dcopy(2*nfft3d,0.0d0,0,tmp1,1)
      do i=1,(nida(nb)+nidb(nb))
        tmp1(int_mb(pack(1,nb)+i-1)) = tmp2(i)
      end do


*     **** make the kx=0 plane complete **** 
      call D3dB_c_timereverse_start(1,tmp1,tmp1z,tmp2z,
     >                              request,reqcnt,msgtype)

      call nwpw_timing_end(9)
  
      return
      end




*     ***********************************
*     *					*
*     *	       Pack_c_unpack_end	*	
*     *					*
*     ***********************************
      subroutine Pack_c_unpack_end(nb,tmp1,tmp1z,tmp2z,request,reqcnt)
      implicit none
      integer    nb
      complex*16 tmp1(*)
      complex*16 tmp1z(*)
      complex*16 tmp2z(*)
      integer    request(*),reqcnt
      

      call nwpw_timing_start(9)

      call D3dB_c_timereverse_end(1,tmp1,tmp1z,tmp2z,request,reqcnt)

      call nwpw_timing_end(9)
  
      return
      end



