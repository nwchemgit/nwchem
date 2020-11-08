*
* $Id$
*

*     ***********************************
*     *					*
*     *	          Pack_Init		*	
*     *					*
*     ***********************************

      subroutine Pack_Init()
      implicit none

#include "bafdecls.fh"
#include "mask_common.fh"
#include "errquit.fh"

*     **** common block for pack ****
      integer nida(0:3),nidb2(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb2
      integer nidb(0:3)
      common / pack2_blk / nidb

      integer nwave_all(0:3)
      common / pack_blk3x / nwave_all

*     **** common block for G_pack ****
      integer G_pack(3,0:3),G_pack_hndl
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
      integer nwave_in(0:3),nwave_out(0:3)

*     **** external functions ****
      logical  control_balance
      integer  G_indx,Gsmall_indx
      external control_balance
      external G_indx,Gsmall_indx


*     **** allocate pack memory ****
      call D3dB_nfft3d(1,nfft3d)

      value = BA_alloc_get(mt_int,nfft3d,
     >                     'pack',pack(2,0),pack(1,0))
      value = value.and.
     >        BA_alloc_get(mt_int,nfft3d,
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
      
!$OMP MASTER
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
!$OMP END MASTER
!$OMP BARRIER



*     ********************************************
*     ***** set up packing for small lattice *****
*     ********************************************
      if (has_small) then
*     **** allocate pack memory ****
      call D3dB_nfft3d(3,nfft3d)

      value = BA_alloc_get(mt_int,nfft3d,
     >                     'pack',pack(2,2),pack(1,2))
      value = value.and.
     >        BA_alloc_get(mt_int,nfft3d,
     >                     'pack',pack(2,3),pack(1,3))
      if (.not. value)
     > call errquit('Pack_init: out of heap memory',1,MA_ERR)


      call Parallel2d_taskid_i(taskid)
      call D3dB_nx(3,nx)
      call D3dB_ny(3,ny)
      call D3dB_nz(3,nz)
      nxh = nx/2
      nyh = ny/2
      nzh = nz/2
!$OMP MASTER
      do nb=2,3
         nida(nb) = 0
         nidb(nb) = 0

*        **** k=(0,0,0)  ****
         k1=0
         k2=0
         k3=0
         !call D3dB_ktoqp(3,k3+1,q,p)
         call D3dB_ijktoindexp(3,k1+1,k2+1,k3+1,indx,p)
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
            !call D3dB_ktoqp(3,k3+1,q,p)
            call D3dB_ijktoindexp(3,k1+1,k2+1,k3+1,indx,p)
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

            !call D3dB_ktoqp(3,k3+1,q,p)
            call D3dB_ijktoindexp(3,k1+1,k2+1,k3+1,indx,p)
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

            !call D3dB_ktoqp(3,k3+1,q,p)
            call D3dB_ijktoindexp(3,k1+1,k2+1,k3+1,indx,p)
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
!$OMP END MASTER
!$OMP BARRIER

      end if !*has_small*


      nwave_in(0) = nida(0) + nidb(0)
      nwave_in(1) = nida(1) + nidb(1)
      if (has_small) then
         nwave_in(2) = nida(2) + nidb(2)
         nwave_in(3) = nida(3) + nidb(3)
      end if
      if (control_balance()) then
        if (has_small) then
           call Balance_Init(4,nwave_in,nwave_out)
        else
           call Balance_Init(2,nwave_in,nwave_out)
        end if
      else
        nwave_out(0) = nwave_in(0)
        nwave_out(1) = nwave_in(1)
        if (has_small) then
           nwave_out(2) = nwave_in(2)
           nwave_out(3) = nwave_in(3)
        end if
      end if
      nidb2(0) = nidb(0) + (nwave_out(0)-nwave_in(0))
      nidb2(1) = nidb(1) + (nwave_out(1)-nwave_in(1))
      if (has_small) then
         nidb2(2) = nidb(2) + (nwave_out(2)-nwave_in(2))
         nidb2(3) = nidb(3) + (nwave_out(3)-nwave_in(3))
      end if

      nwave_all(0) = nida(0) + nidb(0)
      nwave_all(1) = nida(1) + nidb(1)
      if (has_small) then
         nwave_all(2) = nida(2) + nidb(2)
         nwave_all(3) = nida(3) + nidb(3)
      end if
      ! write(*,*) "nwave_all=",nwave_all

      call D3dB_ISumAll(nwave_all(0))
      call D3dB_ISumAll(nwave_all(1)) 
      if (has_small) then
         call D3dB_ISumAll(nwave_all(2)) 
         call D3dB_ISumAll(nwave_all(3)) 
      end if

c      call Parallel2d_taskid_i(taskid)
c      write(*,*) taskid," nida=",nida
c      write(*,*) taskid," nidb=",nidb,nidb2
c      write(*,*) " nida(0),nidb(0),nidb2(0):",taskid,nida(0),
c     >                                      nidb(0),nidb2(0)
c      write(*,*) " nida(1),nidb(1),nidb2(1):",taskid,nida(1),
c     >                                      nidb(1),nidb2(1)


*     **** allocate and define G_pack ****
      nsize = 3*((nida(0)+nidb2(0)) + (nida(1)+nidb2(1)))
      if (has_small) 
     >   nsize = nsize + 3*((nida(2)+nidb2(2)) + (nida(3)+nidb2(3)))

      call D3dB_nfft3d(1,nfft3d)
      value = BA_push_get(mt_dbl,nfft3d,'tmp',tmp(2),tmp(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,nsize,
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

      if (has_small) then
         G_pack(1,2) = G_pack(3,1) + (nida(1)+nidb2(1))
         G_pack(2,2) = G_pack(1,2) + (nida(2)+nidb2(2))
         G_pack(3,2) = G_pack(2,2) + (nida(2)+nidb2(2))
         G_pack(1,3) = G_pack(3,2) + (nida(2)+nidb2(2))
         G_pack(2,3) = G_pack(1,3) + (nida(3)+nidb2(3))
         G_pack(3,3) = G_pack(2,3) + (nida(3)+nidb2(3))
         do nb=2,3
            do i=1,3
               call D3dB_t_Copy(3,dbl_mb(Gsmall_indx(i)),dbl_mb(tmp(1)))
               call Pack_t_pack(nb,dbl_mb(tmp(1)))
               call Pack_t_Copy(nb,dbl_mb(tmp(1)),dbl_mb(G_pack(i,nb)))
            end do
         end do
      end if

      value = BA_pop_stack(tmp(2))
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
      integer G_pack(3,0:3),G_pack_hndl
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

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block for pack ****
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** common block for G_pack ****
      integer G_pack(3,0:3),G_pack_hndl
      common / G_pack_blk / G_pack,G_pack_hndl

      logical value

*     **** external functions ****
      logical  control_balance,lattice_has_small
      external control_balance,lattice_has_small

      value =           BA_free_heap(pack(2,0))
      value = value.and.BA_free_heap(pack(2,1))
      value = value.and.BA_free_heap(G_pack_hndl)
      if (lattice_has_small()) then
         value = value.and.BA_free_heap(pack(2,2))
         value = value.and.BA_free_heap(pack(2,3))
      end if
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

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block for pack ****
      integer nida(0:3),nidb2(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb2
      integer nidb(0:3)
      common / pack2_blk / nidb

      
*     **** local variables ****
      logical value
      integer nfft3d,n,id
      integer tmp1(2)

*     **** external functions ****
      logical  control_balance
      external control_balance

      call nwpw_timing_start(9)

      if (nb.lt.2) then
         id = 1
      else
         id = 3
      end if

      call D3dB_nfft3d(id,nfft3d)
      value = BA_push_get(mt_dcpl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not.value) call errquit('out of stack memory',0, MA_ERR)

c      call dcopy(2*nfft3d,A,1,dcpl_mb(tmp1(1)),1)
c      call dcopy(2*nfft3d,0.0d0,0,A,1)
      call Parallel_shared_vector_copy(.true.,2*nfft3d,A,
     >                                 dcpl_mb(tmp1(1)))
      call Parallel_shared_vector_zero(.true.,2*nfft3d,A)

c      do i=1,(nida(nb)+nidb(nb))
c        A(i) = dcpl_mb(tmp1(1) + int_mb(pack(1,nb)+i-1)-1)
c      end do
      n = nida(nb)+nidb(nb)
      call Pack_c_indexcopy(n,int_mb(pack(1,nb)),dcpl_mb(tmp1(1)),A)

      value = BA_pop_stack(tmp1(2))
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

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block for pack ****
      integer nida(0:3),nidb2(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb2
      integer nidb(0:3)
      common / pack2_blk / nidb


*     **** local variables ****
      logical value
      integer i,j,nfft3d,n,shift,id
      integer tmp1(2)

*     **** external functions ****
      logical  control_balance
      external control_balance

      call nwpw_timing_start(9)

      if (nb.lt.2) then
         id = 1
      else
         id = 3
      end if

      call D3dB_nfft3d(id,nfft3d)
      value = BA_push_get(mt_dcpl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not.value) call errquit('out of stack memory',0, MA_ERR)

      shift=1
      n = nida(nb)+nidb(nb)
      do j=1,m
        !call dcopy(2*nfft3d,A(shift),1,dcpl_mb(tmp1(1)),1)
        !call dcopy(2*nfft3d,0.0d0,0,A(shift),1)
        call Parallel_shared_vector_copy(.true.,2*nfft3d,A(shift),
     >                                   dcpl_mb(tmp1(1)))
        call Parallel_shared_vector_zero(.true.,2*nfft3d,A(shift))

        call Pack_c_indexcopy(n,int_mb(pack(1,nb)),
     >                        dcpl_mb(tmp1(1)),
     >                        A(shift))
        if (control_balance()) call Balance_c_balance(nb,A(shift))
        shift = shift+nfft3d
      end do

      value = BA_pop_stack(tmp1(2))
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
!$OMP DO
      do i=1,n
        B(i) = A(indx(i))
      end do
!$OMP END DO
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

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block for pack ****
      integer nida(0:3),nidb2(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb2
      integer nidb(0:3)
      common / pack2_blk / nidb

      
*     **** local variables ****
      logical value
      integer i,nfft3d,id
      integer tmp1(2)

*     **** external functions ****
      logical  control_balance
      external control_balance

      call nwpw_timing_start(9)

      if (nb.lt.2) then
         id = 1
      else
         id = 3
      end if

      call D3dB_nfft3d(id,nfft3d)
      value = BA_push_get(mt_dbl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not.value) call errquit('out of stack memory',0, MA_ERR)

      !call dcopy(nfft3d,A,1,dbl_mb(tmp1(1)),1)
      !call dcopy(nfft3d,0.0d0,0,A,1)
      call Parallel_shared_vector_copy(.true.,nfft3d,A,
     >                                 dbl_mb(tmp1(1)))
      call Parallel_shared_vector_zero(.true.,nfft3d,A)

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        A(i) = dbl_mb(tmp1(1) + int_mb(pack(1,nb)+i-1) - 1)
      end do
!$OMP END DO

      value = BA_pop_stack(tmp1(2))
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

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block for pack ****
      integer nida(0:3),nidb2(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb2
      integer nidb(0:3)
      common / pack2_blk / nidb

      
*     **** local variables ****
      logical value
      integer i,nfft3d,id
      integer tmp1(2)

*     **** external functions ****
      logical  control_balance
      external control_balance

      call nwpw_timing_start(9)
      if (nb.lt.2) then
         id = 1
      else
         id = 3
      end if

      call D3dB_nfft3d(id,nfft3d)
      value = BA_push_get(mt_int,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not.value) call errquit('out of stack memory',0, MA_ERR)

      call icopy(nfft3d,A,1,int_mb(tmp1(1)),1)
      call icopy(nfft3d,0,0,A,1)

      do i=1,(nida(nb)+nidb(nb))
        A(i) = int_mb(tmp1(1) + int_mb(pack(1,nb)+i-1)-1)
      end do

      value = BA_pop_stack(tmp1(2))
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

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block for pack ****
      integer nida(0:3),nidb2(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb2
      integer nidb(0:3)
      common / pack2_blk / nidb


*     **** local variables ****
      logical value
      integer i,nfft3d,zplane_size,id
      integer tmp1(2)
      integer tmp2(2)
c      integer indx1,indx2

*     **** external functions ****
      logical  control_balance
      external control_balance


      call nwpw_timing_start(9)
      if (nb.lt.2) then
         id = 1
      else
         id = 3
      end if
      if (control_balance()) call Balance_c_unbalance(nb,A)

      call D3dB_nfft3d(id,nfft3d)
      
      value = BA_push_get(mt_dcpl,(nida(nb)+nidb(nb)),
     >                    'tmp1',tmp1(2),tmp1(1))
      if (.not.value) call errquit('out of stack memory',0, MA_ERR)

c      call dcopy(2*(nida(nb)+nidb(nb)),A,1,dcpl_mb(tmp1(1)),1)
c      call dcopy(2*nfft3d,0.0d0,0,A,1)
      call Parallel_shared_vector_copy(.true.,2*(nida(nb)+nidb(nb)),
     >                                 A,dcpl_mb(tmp1(1)))
      call Parallel_shared_vector_zero(.true.,2*nfft3d,A)
!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        A(int_mb(pack(1,nb)+i-1)) = dcpl_mb(tmp1(1)+i-1)
      end do
!$OMP END DO
      value = BA_pop_stack(tmp1(2))
      if (.not.value) call errquit('error popping stack',0, MA_ERR)


*     **** make the kx=0 plane complete **** 
      !call D3dB_ny(1,ny)
      !call D3dB_nq(1,nq)
      call D3dB_zplane_size(id,zplane_size)
      value = BA_push_get(mt_dcpl,(zplane_size),'tmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(zplane_size),'tmp2',tmp2(2),tmp2(1))

      call D3dB_c_timereverse(id,A,dcpl_mb(tmp1(1)),dcpl_mb(tmp2(1)))

      value = BA_pop_stack(tmp2(2))
      value = value.and.
     >        BA_pop_stack(tmp1(2))
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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
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
      integer nwave_all(0:3)
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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

      integer qzero,id

      if (nb.lt.2) then
         id = 1
      else
         id = 3
      end if
*     *********************************************************
*     **** warning this routine assumes a specific packing ****
*     *********************************************************
*     index = (qzero-1)*(nx/2+1)*ny + (j-1)*(nx/2+1) + i 
      zero = 1
      !call D3dB_ktoqp(1,1,qzero,pzero)
      call D3dB_ijktoindexp(id,1,1,1,zero,pzero)

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer j,k
      integer np,npack,npack2

      call nwpw_timing_start(2)
      call Parallel2d_np_i(np)

      npack  = (nida(nb)+nidb(nb))
      npack2 = 2*npack

      do k=1,n
        call DGEMM_OMP('T','N',k,1,npack2,
     >             2.0d0,
     >             A,npack2, 
     >             B(1+(k-1)*npack),npack2, 
     >             0.0d0,
     >             matrix(1,k),k)
        call DGEMM_OMP('T','N',k,1,2*nida(nb),
     >             -1.0d0,
     >             A,npack2, 
     >             B(1+(k-1)*npack),npack2, 
     >             1.0d0,
     >             matrix(1,k),k)
      end do

!$OMP DO
      do k=1,n
      do j=k+1,n
        matrix(j,k) = matrix(k,j)
      end do
      end do
!$OMP END DO

      if (np.gt.1) call D3dB_Vector_SumAll(n*n,matrix)

      call nwpw_timing_end(2)
      return
      end 


*     ***********************************
*     *                                 *
*     *          Pack_ccm_sym_dot_omp   *       
*     *                                 *
*     ***********************************

      subroutine Pack_ccm_sym_dot_omp(nb,n,A,B,matrix,privatetmp)
      implicit none
      integer    nb,n
      complex*16 A(*)
      complex*16 B(*)
      real*8     matrix(n,n)
      real*8     privatetmp(*)

*     **** common block for pack ****
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer j,k
      integer np,npack,npack2
      integer tid,nthr

*     **** external functions
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      call nwpw_timing_start(2)
      call Parallel2d_np_i(np)
      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

      npack  = (nida(nb)+nidb(nb))
      npack2 = 2*npack
      call Parallel_shared_vector_zero(.true.,n*n,matrix)

      call dgemm2c_omp_group(1,tid,nthr,n,n,npack2,2*nida(nb),
     >                       A,B,matrix,privatetmp)
!$OMP BARRIER

!$OMP DO
      do k=1,n
      do j=k+1,n
        matrix(j,k) = matrix(k,j)
      end do
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i,j,k
      integer np,npack,npack2



      call nwpw_timing_start(2)
      call Parallel2d_np_i(np)

      npack  = (nida(nb)+nidb(nb))
      npack2 = 2*npack

      do k=1,n
        call DGEMM_OMP('T','N',k,1,npack2,
     >             2.0d0,
     >             A,npack2,
     >             A(1+(k-1)*npack),npack2,
     >             0.0d0,
     >             matrix(1,k,1),k)
        call DGEMM_OMP('T','N',k,1,npack2,
     >             2.0d0,
     >             A,npack2,
     >             B(1+(k-1)*npack),npack2,
     >             0.0d0,
     >             matrix(1,k,2),k)
        call DGEMM_OMP('T','N',k,1,npack2,
     >             2.0d0,
     >             B,npack2,
     >             B(1+(k-1)*npack),npack2,
     >             0.0d0,
     >             matrix(1,k,3),k)

        call DGEMM_OMP('T','N',k,1,2*nida(nb),
     >             -1.0d0,
     >             A,npack2,
     >             A(1+(k-1)*npack),npack2,
     >             1.0d0,
     >             matrix(1,k,1),k)
        call DGEMM_OMP('T','N',k,1,2*nida(nb),
     >             -1.0d0,
     >             A,npack2,
     >             B(1+(k-1)*npack),npack2,
     >             1.0d0,
     >             matrix(1,k,2),k)
        call DGEMM_OMP('T','N',k,1,2*nida(nb),
     >             -1.0d0,
     >             B,npack2,
     >             B(1+(k-1)*npack),npack2,
     >             1.0d0,
     >             matrix(1,k,3),k)
      end do

!$OMP MASTER
      do i=1,3
      do k=1,n
      do j=k+1,n
        matrix(j,k,i) = matrix(k,j,i)
      end do
      end do
      end do
!$OMP END MASTER
!$OMP BARRIER

      if (np.gt.1) call D3dB_Vector_SumAll(3*n*n,matrix)

      call nwpw_timing_end(2)
      return
      end



*     ***********************************
*     *                                 *
*     *    Pack_ccm_combo_sym_dot_omp   *
*     *                                 *
*     ***********************************
      subroutine Pack_ccm_combo_sym_dot_omp(nb,n,A,B,matrix,privatetmp)
      implicit none
      integer    nb,n
      complex*16 A(*)
      complex*16 B(*)
      real*8     matrix(n,n,3)
      real*8     privatetmp(*)

*     **** common block for pack ****
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i,j,k,ik
      integer np,npack2
      integer tid,nthr,nthr3,nthr3l,ftid1,ftid2,ftid3

*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      call nwpw_timing_start(2)
      call Parallel2d_np_i(np)
      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

      if (nthr.eq.2) then
         nthr3  = 2
         nthr3l = 2
         ftid1 = 0
         ftid2 = 0
         ftid3 = 0
      else
         nthr3  = max(1,(nthr - mod(nthr,3))/3)
         nthr3l = max(1,nthr-2*nthr3)
         ftid1 = 0
         ftid2 = min(nthr-1,nthr3)
         ftid3 = min(nthr-1,2*nthr3)
      end if

      npack2 = 2*(nida(nb)+nidb(nb))
      call Parallel_shared_vector_zero(.true.,3*n*n,matrix)

      if ((tid.ge.ftid1).and.(tid.lt.(ftid1+nthr3))) then
         call dgemm2c_omp_group(1,tid-ftid1,nthr3,
     >                          n,n,npack2,2*nida(nb),
     >                          A,A,matrix(1,1,1),privatetmp)
      end if

      if ((tid.ge.ftid2).and.(tid.lt.(ftid2+nthr3))) then
         call dgemm2c_omp_group(2,tid-ftid2,nthr3,
     >                          n,n,npack2,2*nida(nb),
     >                          B,A,matrix(1,1,2),privatetmp)
      end if

      if ((tid.ge.ftid3).and.(tid.lt.(ftid3+nthr3l))) then
         call dgemm2c_omp_group(3,tid-ftid3,nthr3l,
     >                          n,n,npack2,2*nida(nb),
     >                          B,B,matrix(1,1,3),privatetmp)
      end if
!$OMP BARRIER

!$OMP DO
      do ik=1,3*n
         i = (ik-1)/n + 1
         k = mod(ik-1,n)+1
         do j=k+1,n
            matrix(j,k,i) = matrix(k,j,i)
         end do
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer j,k
      integer np,npack,npack2


      call nwpw_timing_start(2)
      call Parallel2d_np_i(np)

      npack  = (nida(nb)+nidb(nb))
      npack2 = 2*npack

      do k=1,n
        call DGEMM_OMP('T','N',k,1,npack2,
     >             2.0d0,
     >             A,npack2, 
     >             B(1+(k-1)*npack),npack2, 
     >             0.0d0,
     >             matrix(1,k),k)
        call DGEMM_OMP('T','N',k,1,2*nida(nb),
     >             -1.0d0,
     >             A,npack2, 
     >             B(1+(k-1)*npack),npack2, 
     >             1.0d0,
     >             matrix(1,k),k)
      end do

!$OMP MASTER
      do k=1,n
      do j=k+1,n
        matrix(j,k) = matrix(k,j)
      end do
      end do
!$OMP END MASTER
!$OMP BARRIER

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer np,npack


      call nwpw_timing_start(2)
      call Parallel2d_np_i(np)

      npack = 2*(nida(nb)+nidb(nb))

        call DGEMM_OMP('T','N',n,n,npack,
     >             2.0d0,
     >             A,npack, 
     >             B,npack, 
     >             0.0d0,
     >             matrix,n)
        call DGEMM_OMP('T','N',n,n,2*nida(nb),
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
*     *                                 *
*     *          Pack_ccm_dot_omp       *       
*     *                                 *
*     ***********************************

      subroutine Pack_ccm_dot_omp(nb,n,A,B,matrix,privatetmp)
      implicit none
      integer    nb,n
      complex*16 A(*)
      complex*16 B(*)
      real*8     matrix(n,n)
      real*8     privatetmp(*)

*     **** common block for pack ****
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer np,npack
      integer tid,nthr

*     **** external functions
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads


      call nwpw_timing_start(2)
      call Parallel2d_np_i(np)
      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

      npack = 2*(nida(nb)+nidb(nb))
      call Parallel_shared_vector_zero(.true.,n*n,matrix)

      call dgemm2c_omp_group(1,tid,nthr,n,n,npack,2*nida(nb),
     >                       A,B,matrix,privatetmp)
!$OMP BARRIER

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer np,npack


      call nwpw_timing_start(2)
      call Parallel2d_np_i(np)

      npack = 2*(nida(nb)+nidb(nb))

        call DGEMM_OMP('T','N',m,n,npack,
     >             2.0d0,
     >             A,npack, 
     >             B,npack, 
     >             0.0d0,
     >             matrix,m)
        call DGEMM_OMP('T','N',m,n,2*nida(nb),
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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer n,np
c     integer npack, shift
      integer npack2,npack,shift

*     **** external functions ****
      real*8   ddot
      external ddot

      call nwpw_timing_start(2)
      call Parallel2d_np_i(np)
      npack = (nida(nb)+nidb(nb))
      npack2= 2*npack

c     do n=1,ne
c       shift = (n-1)*npack
c
c       sum(n) = ddot(2*nida(nb),A(1+shift),1,B,1)
c       sum(n) = sum(n) + 2*ddot(2*nidb(nb),
c    >                           A(nida(nb)+1+shift),1,
c    >                           B(nida(nb)+1),      1)
c     end do

      call DGEMM_OMP('T','N',ne,1,(npack2),
     >           (2.0d0),
     >           A,npack2,
     >           B,npack2,
     >           (0.0d0),
     >           sum,ne)
c      call DGEMM_OMP('T','N',ne,1,(2*nida(nb)),
c     >           (-1.0d0),
c     >           A,npack2,
c     >           B,npack2,
c     >           (1.0d0),
c     >           sum,ne)
!$OMP DO
      do n=1,ne
         shift = (n-1)*npack
         sum(n) = sum(n) - ddot(2*nida(nb),A(1+shift),1,B,1)
      end do
!$OMP END DO


      call nwpw_timing_end(2)
      return
      end 



*     ***********************************
*     *                                 *
*     *       Pack_cc_indot_thread      *       
*     *                                 *
*     ***********************************

      subroutine Pack_cc_indot_thread(nb,ne,nprj,A,B,sumt)
      implicit none
      integer    nb
      integer    ne,nprj
      complex*16 A(*)
      complex*16 B(*)
      real*8     sumt(*)


*     **** common block for pack ****
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i,j,m,n,np,tid,nthr
c     integer npack, shift
      integer npack1,shift,shift2,r,Itid,Mtid

*     **** external functions ****
      real*8   ddot
      external ddot
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      call nwpw_timing_start(2)
      call Parallel2d_np_i(np)
      npack1 = (nida(nb)+nidb(nb))

      nthr = Parallel_nthreads()
      tid  = Parallel_threadid()
      r    = mod(npack1,nthr)
      if (tid.lt.r) then
         Itid = tid*(npack1/nthr+1) + 1
         Mtid = npack1/nthr + 1
      else
          Itid = r + tid*(npack1/nthr) + 1
          Mtid = npack1/nthr
      end if

      call DGEMM('T','N',ne,nprj,(2*Mtid),
     >           (2.0d0),
     >           A(Itid),2*npack1,
     >           B(Itid),2*npack1,
     >           (0.0d0),
     >           sumt,ne)

!$OMP DO
      do n=1,ne*nprj
         i = mod(n-1,ne)
         j = (n-1-i)/ne
         shift  = i*npack1
         shift2 = j*npack1
         sumt(n) = sumt(n) - ddot(2*nida(nb),A(1+shift),1,B(1+shift2),1)
      end do
!$OMP END DO


      call nwpw_timing_end(2)
      return
      end






*     ***********************************
*     *                                 *
*     *       Pack_conjg_tcc_indot      *       
*     *                                 *
*     ***********************************

      subroutine Pack_conjg_tcc_indot(nb,ne,A,B,C,sm)
      implicit none
      integer    nb
      integer    ne
      real*8     A(*)
      complex*16 B(*)
      complex*16 C(*)
      real*8     sm(*)


*     **** common block for pack ****
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer n,np,i,k,npack
      real*8 sa,sb
      common /pack_conjg_tcc_indot_sasb/ sa,sb

      npack = nida(nb)+nidb(nb)
      do i=1,ne
!$OMP MASTER
         sa = 0.0d0
         sb = 0.0d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO REDUCTION(+:sa)
         do k=1,nida(nb)
            sa = sa + dimag(A(k)*B(k+(i-1)*npack)*dconjg(C(k)))
         end do
!$OMP END DO

!$OMP DO REDUCTION(+:sb)
         do k=nida(nb)+1,nida(nb)+nidb(nb)
            sb = sb + 2.0d0*dimag(A(k)*B(k+(i-1)*npack)*dconjg(C(k)))
         end do
!$OMP END DO

!$OMP MASTER
         sm(i) = sa + sb
!$OMP END MASTER
      end do
!$OMP BARRIER
      return
      end



*     ***********************************
*     *                                 *
*     *         Pack_ccp_idot           *       
*     *                                 *
*     ***********************************

* return variable sum is assumed to be private
      subroutine Pack_ccp_idot(nb,A,B,sum)
      implicit none
      integer    nb
      real*8 A(*)
      real*8 B(*)
      real*8     sum

*     **** common block for pack ****
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i
      real*8 sa,sb
      common /pack_conjg_tcc_indot_sasb/ sa,sb

*     **** external functions ****
c      real*8   ddot
c      external ddot

      call nwpw_timing_start(2)

!$OMP MASTER
      sa = 0.0d0
      sb = 0.0d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO REDUCTION(+:sa)
      do i=1,2*nida(nb)
         sa = sa + A(i)*B(i)
      end do
!$OMP END DO

!$OMP DO REDUCTION(+:sb)
      do i=1,2*nidb(nb)
         sb = sb + 2.0d0*A(2*nida(nb)+i)*B(2*nida(nb)+i)
      end do
!$OMP END DO


      sum = sa + sb

c     if (np.gt.1) call D3dB_SumAll(sum)

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
      real*8 A(*)
      real*8 B(*)
      real*8     sum

*     **** common block for pack ****
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i
      real*8 sa,sb
      common /pack_conjg_tcc_indot_sasb/ sa,sb

*     **** external functions ****
c      real*8   ddot
c      external ddot

      call nwpw_timing_start(2)
c      call Parallel2d_np_i(np)

c      sum = ddot(2*nida(nb),A,1,B,1)
c      sum = sum + 2*ddot(2*nidb(nb),A(nida(nb)+1),1,B(nida(nb)+1),1)

!$OMP MASTER
      sa = 0.0d0
      sb = 0.0d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO REDUCTION(+:sa)
      do i=1,2*nida(nb)
         sa = sa + A(i)*B(i)
      end do
!$OMP END DO

!$OMP DO REDUCTION(+:sb)
      do i=1,2*nidb(nb)
         sb = sb + 2.0d0*A(2*nida(nb)+i)*B(2*nida(nb)+i)
      end do
!$OMP END DO

!$OMP MASTER
      sum = sa + sb
!$OMP END MASTER
!$OMP BARRIER

c     if (np.gt.1) call D3dB_SumAll(sum)

      call nwpw_timing_end(2)

      return
      end 

      subroutine Pack_conjg_tcc_idot(nb,A,B,C,summer)
      implicit none
      integer    nb
      real*8     A(*)
      complex*16 B(*)
      complex*16 C(*)
      real*8     summer

*     **** common block for pack ****
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer k

c*     **** external functions ****
c      real*8   ddot
c      external ddot

      call nwpw_timing_start(2)
c      call Parallel2d_np_i(np)

      summer = 0.0d0
      do k=1,nida(nb)
         summer = summer - dimag(A(k)*dconjg(B(k))*C(k))
      end do
      do k=nida(nb)+1,nida(nb)+nidb(nb)
         summer = summer - 2.0d0*dimag(A(k)*dconjg(B(k))*C(k))
      end do

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
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

*     shared output variable:
*     - sum
*
      subroutine Pack_cc_dot(nb,A,B,sum)
      implicit none
      integer    nb
      real*8     A(*)
      real*8     B(*)
      real*8     sum

*     **** common block for pack ****
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer np,i
      real*8 sa,sb
      common /pack_conjg_tcc_indot_sasb/ sa,sb

*     **** external functions ****
      real*8   ddot
      external ddot

      call nwpw_timing_start(2)
      call Parallel2d_np_i(np)

c      sum = ddot(2*nida(nb),A,1,B,1)
c      sum = sum + 2*ddot(2*nidb(nb),A(nida(nb)+1),1,B(nida(nb)+1),1)
!$OMP MASTER
      sa = 0.0d0
      sb = 0.0d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO REDUCTION(+:sa)
      do i=1,2*nida(nb)
         sa = sa + A(i)*B(i)
      end do
!$OMP END DO

!$OMP DO REDUCTION(+:sb)
      do i=1,2*nidb(nb)
         sb = sb + 2.0d0*A(2*nida(nb)+i)*B(2*nida(nb)+i)
      end do
!$OMP END DO

!$OMP MASTER
      sum = sa + sb
!$OMP END MASTER
!$OMP BARRIER

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb


*     **** local variables ****
      integer i
      real*8 sa,sb
      common /pack_conjg_tcc_indot_sasb/ sa,sb

*     **** external functions ****
c      real*8   ddot
c      external ddot

      call nwpw_timing_start(2)

c      sum =         ddot(nida(nb),A,1,B,1)
c      sum = sum + 2*ddot(nidb(nb),A(nida(nb)+1),1,B(nida(nb)+1),1)
!$OMP MASTER
      sa = 0.0d0
      sb = 0.0d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO REDUCTION(+:sa)
      do i=1,nida(nb)
         sa = sa + A(i)*B(i)
      end do
!$OMP END DO

!$OMP DO REDUCTION(+:sb)
      do i=1,nidb(nb)
         sb = sb + 2.0d0*A(nida(nb)+i)*B(nida(nb)+i)
      end do
!$OMP END DO

!$OMP MASTER
      sum = sa + sb
!$OMP END MASTER
!$OMP BARRIER

      call nwpw_timing_end(2)
      return
      end 


*     ***********************************
*     *                                 *
*     *          Pack_ttp_dot           *       
*     *                                 *
*     ***********************************
      
      subroutine Pack_ttp_dot(nb,A,B,sump)
      implicit none
      integer nb
      real*8  A(*)
      real*8  B(*)
      real*8  sump

*     **** common block for pack ****
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb
      
*     **** local variables ****
      integer np,i
      real*8 sa,sb
      common /pack_conjg_tcc_indot_sasb/ sa,sb

*     **** external functions ****
c      real*8   ddot
c      external ddot

      call nwpw_timing_start(2)
      call Parallel2d_np_i(np)

c      sum =         ddot(nida(nb),A,1,B,1)
c      sum = sum + 2*ddot(nidb(nb),A(nida(nb)+1),1,B(nida(nb)+1),1)
!$OMP MASTER
      sa = 0.0d0
      sb = 0.0d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO REDUCTION(+:sa)
      do i=1,nida(nb)
         sa = sa + A(i)*B(i)
      end do
!$OMP END DO

!$OMP DO REDUCTION(+:sb)
      do i=1,nidb(nb)
         sb = sb + 2.0d0*A(nida(nb)+i)*B(nida(nb)+i)
      end do
!$OMP END DO

!$OMP MASTER
      sa = sa + sb
!$OMP END MASTER
!$OMP BARRIER

      if (np.gt.1) call D3dB_SumAll(sa)
      sump = sa

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer np,i
      real*8 sa,sb
      common /pack_conjg_tcc_indot_sasb/ sa,sb

*     **** external functions ****
c      real*8   ddot
c      external ddot

      call nwpw_timing_start(2)
      call Parallel2d_np_i(np)

c      sum =         ddot(nida(nb),A,1,B,1)
c      sum = sum + 2*ddot(nidb(nb),A(nida(nb)+1),1,B(nida(nb)+1),1)
!$OMP MASTER
      sa = 0.0d0
      sb = 0.0d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO REDUCTION(+:sa)
      do i=1,nida(nb)
         sa = sa + A(i)*B(i)
      end do
!$OMP END DO

!$OMP DO REDUCTION(+:sb)
      do i=1,nidb(nb)
         sb = sb + 2.0d0*A(nida(nb)+i)*B(nida(nb)+i)
      end do
!$OMP END DO

!$OMP MASTER
      sum = sa + sb
!$OMP END MASTER
!$OMP BARRIER

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
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
*     *                                 *
*     *          Pack_tc_Copy           *       
*     *                                 *
*     ***********************************

      subroutine Pack_tc_Copy(nb,A,B)
      implicit none
      integer    nb
      real*8     A(*)
      complex*16 B(*)

*     **** common block for pack ****
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

      integer i

c      call dcopy(2*(nida(nb)+nidb(nb)),A,1,B,1)
!$OMP DO
      do i=1,2*(nida(nb)+nidb(nb))
         B(i) = dcmplx(A(i),0.0d0)
      end do
!$OMP END DO

      return
      end

*     ***********************************
*     *                                 *
*     *          Pack_ct_Copy           *       
*     *                                 *
*     ***********************************

      subroutine Pack_ct_Copy(nb,A,B)
      implicit none
      integer    nb
      complex*16 A(*)
      real*8     B(*)

*     **** common block for pack ****
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

      integer i

!$OMP DO
      do i=1,2*(nida(nb)+nidb(nb))
         B(i) = dble(A(i))
      end do
!$OMP END DO

      return
      end




*     ***********************************
*     *					*
*     *	         Pack_c_Copy1		*	
*     *					*
*     ***********************************

      subroutine Pack_c_Copy1(nb,A,B)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)

*     **** common block for pack ****
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb


c      call dcopy(2*(nida(nb)+nidb(nb)),A,1,B,1)
      call Parallel_shared_vector_copy(.true.,2*(nida(nb)+nidb(nb)),A,B)
      return
      end

*     ***********************************
*     *					*
*     *	         Pack_c_Copy0		*	
*     *					*
*     ***********************************

      subroutine Pack_c_Copy0(nb,A,B)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)

*     **** common block for pack ****
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb


      call dcopy(2*(nida(nb)+nidb(nb)),A,1,B,1)

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb


c      call dcopy(2*(nida(nb)+nidb(nb)),A,1,B,1)
      call Parallel_shared_vector_copy(.true.,2*(nida(nb)+nidb(nb)),A,B)

      return
      end


      subroutine Pack_c_nCopy(nb,n,A,B)
      implicit none
      integer    nb,n
      complex*16 A(*)
      complex*16 B(*)

*     **** common block for pack ****
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb
        
c      call dcopy(n*2*(nida(nb)+nidb(nb)),A,1,B,1)
      call Parallel_shared_vector_copy(.true.,
     >                                 n*2*(nida(nb)+nidb(nb)),
     >                                 A,B)
        
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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb


      !call dcopy((nida(nb)+nidb(nb)),A,1,B,1)
      call Parallel_shared_vector_copy(.true.,(nida(nb)+nidb(nb)),A,B)

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb


      !call dcopy(2*(nida(nb)+nidb(nb)),0.0d0,0,A,1)
      call Parallel_shared_vector_zero(.true.,2*(nida(nb)+nidb(nb)),A)

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        C(i) = A(i) + B(i)
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        B(i) = B(i) + A(i)
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        C(i) = A(i) + B(i)
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        B(i) = B(i) + A(i)
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        C(i) = A(i) - B(i)
      end do
!$OMP END DO

      return
      end


*     ***********************************
*     *                                 *
*     *          Pack_cc_Sub2           *
*     *                                 *
*     ***********************************

      subroutine Pack_cc_Sub2(nb,A,B)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)


*     **** common block for pack ****
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        B(i) = B(i) - A(i)
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        C(i) = A(i) - B(i)
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        C(i) = dsqrt(A(i))
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        A(i) = dsqrt(A(i))
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        C(i) = dble(A(i))**2 + dimag(A(i))**2
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        C(i) = alpha*A(i)
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        A(i) = alpha*A(i)
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        C(i) = alpha*A(i)
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        C(i) = alpha*A(i)
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        A(i) = alpha*A(i)
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

      integer i

      !call daxpy(2*(nida(nb)+nidb(nb)),alpha,A,1,B,1)
!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
         B(i) = B(i) + alpha*A(i)
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

c      do i=1,(nida(nb)+nidb(nb))
c        C(i) = C(i) + alpha*A(i)
c      end do
      call DAXPY_OMP((nida(nb)+nidb(nb)),alpha,A,1,C,1)

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i
!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        C(i) = C(i) + alpha*A(i)
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        C(i) = B(i)*A(i)
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        B(i) = B(i)*A(i)
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        C(i) = alpha*B(i)*A(i)
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        C(i) = C(i) + alpha*B(i)*A(i)
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        C(i) = C(i) + A(i)*B(i)
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        C(i) = A(i) * dcmplx(-dimag(B(i)),dble(B(i)))
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        C(i) = A(i)*B(i)
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        B(i) = B(i)*A(i)
      end do
!$OMP END DO

      return
      end


*     ***********************************
*     *                                 *
*     *          Pack_cct_Mul           *       
*     *                                 *
*     ***********************************
      subroutine Pack_cct_Mul(nb,A,B,C)
      implicit none
      integer    nb            
      complex*16 A(*)
      complex*16 B(*)
      real*8     C(*)
 
*     **** common block for pack ****
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)         
      common / pack_blk / pack,nida,nidb
              
*     **** local variables ****
      integer i                  

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        C(i) = dble(A(i)*B(i))
      end do 
!$OMP END DO
            
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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        C(i) = dble(dconjg(A(i))*B(i))
      end do
!$OMP END DO

      return
      end

*     ***********************************
*     *                                 *
*     *       Pack_cct_iconjgMul        *
*     *                                 *
*     ***********************************
      subroutine Pack_cct_iconjgMul(nb,A,B,C)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      real*8     C(*)

*     **** common block for pack ****
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        C(i) = dimag(dconjg(A(i))*B(i))
      end do
!$OMP END DO

      return
      end

*     ***********************************
*     *                                 *
*     *       Pack_cct_iconjgMulAdd     *
*     *                                 *
*     ***********************************
      subroutine Pack_cct_iconjgMulAdd(nb,A,B,C)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      real*8     C(*)

*     **** common block for pack ****
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        C(i) = C(i) + dimag(dconjg(A(i))*B(i))
      end do
!$OMP END DO

      return
      end

*     ***********************************
*     *                                 *
*     *       Pack_cct_iaconjgMulAdd     *
*     *                                 *
*     ***********************************
      subroutine Pack_cct_iaconjgMulAdd(nb,alpha,A,B,C)
      implicit none
      integer    nb
      real*8     alpha
      complex*16 A(*)
      complex*16 B(*)
      real*8     C(*)

*     **** common block for pack ****
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        C(i) = C(i) + alpha*dimag(A(i)*dconjg(B(i)))
      end do
!$OMP END DO

      return
      end

! (a+ib)*(c-id) = ac+bd + i*(bc-ad)
! (a(1)+i*a(2))*(b(1)-i*b(2)) = ... + i*(a(2)*b(1) - a(1)*b(2))



      subroutine Pack_cct_iconjgMulb(nb,A,B,C)
      implicit none
      integer    nb
      complex*16 A(*)
      complex*16 B(*)
      real*8     C(*)

*     **** common block for pack ****
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        C(i) = dimag(A(i))*dble(B(i)) - dble(A(i))*dimag(B(i))
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        C(i) = dconjg(A(i))*B(i)
      end do
!$OMP END DO

      return
      end

*     ***********************************
*     *					*
*     *	         Pack_ttcc_AddMul       *
*     *					*
*     ***********************************

      subroutine Pack_ttcc_AddMul(nb,A,B,C,D)
      implicit none
      integer nb
      real*8 A(*),B(*)
      complex*16 C(*),D(*)

*     **** common block for pack ****
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
         D(i) = (A(i)+B(i))*C(i)
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        C(i) = dcmplx(0.0d0,A(i)) * B(i)
      end do
!$OMP END DO

      return
      end



*     ***********************************
*     *                                 *
*     *          Pack_tc_iMul2          *
*     *                                 *
*     ***********************************

      subroutine Pack_tc_iMul2(nb,A,B)
      implicit none
      integer    nb
      real*8     A(*)
      complex*16 B(*)


*     **** common block for pack ****
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        B(i) = dcmplx(0.0d0,A(i)) * B(i)
      end do
!$OMP END DO

      return
      end


*     ***********************************
*     *                                 *
*     *          Pack_cc_xyzgradient    *
*     *                                 *
*     ***********************************
      subroutine Pack_cc_xyzgradient(nb,xyz,A,B)
      implicit none
      integer    nb,xyz
      complex*16 A(*)
      complex*16 B(*)

#include "bafdecls.fh"

*     **** common block for pack ****
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** common block for G_pack ****
      integer G_pack(3,0:3),G_pack_hndl
      common / G_pack_blk / G_pack,G_pack_hndl

      call Pack_tc_iMul(nb,dbl_mb(G_pack(xyz,nb)),A,B)
      return
      end


*     ***********************************
*     *                                 *
*     *       Pack_cc_xyzgradient1      *
*     *                                 *
*     ***********************************
      subroutine Pack_cc_xyzgradient1(nb,xyz,B)
      implicit none
      integer    nb,xyz
      complex*16 B(*)

#include "bafdecls.fh"

*     **** common block for pack ****
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** common block for G_pack ****
      integer G_pack(3,0:3),G_pack_hndl
      common / G_pack_blk / G_pack,G_pack_hndl

      call Pack_tc_iMul2(nb,dbl_mb(G_pack(xyz,nb)),B)
      return
      end



*     ***********************************
*     *                                 *
*     *    Pack_cc_multiplegradients    *
*     *                                 *
*     ***********************************
      subroutine Pack_cc_multiplegradients(nb,ii,A,B)
      implicit none
      integer    nb,ii
      complex*16 A(*)
      complex*16 B(*)

c     **** one derivatives ****
      if (ii.eq.1) call Pack_cc_xyzgradient(nb,1,A,B)
      if (ii.eq.2) call Pack_cc_xyzgradient(nb,2,A,B)
      if (ii.eq.3) call Pack_cc_xyzgradient(nb,3,A,B)

c     **** two derivatives ****
      if (ii.eq.4) then
         call Pack_cc_xyzgradient(nb,1,A,B)
         call Pack_cc_xyzgradient1(nb,2,B)
      end if
      if (ii.eq.5) then
         call Pack_cc_xyzgradient(nb,1,A,B)
         call Pack_cc_xyzgradient1(nb,3,B)
      end if
      if (ii.eq.6) then
         call Pack_cc_xyzgradient(nb,2,A,B)
         call Pack_cc_xyzgradient1(nb,3,B)
      end if
      if (ii.eq.7) then
         call Pack_cc_xyzgradient(nb,1,A,B)
         call Pack_cc_xyzgradient1(nb,1,B)
      end if
      if (ii.eq.8) then
         call Pack_cc_xyzgradient(nb,2,A,B)
         call Pack_cc_xyzgradient1(nb,2,B)
      end if
      if (ii.eq.9) then
         call Pack_cc_xyzgradient(nb,3,A,B)
         call Pack_cc_xyzgradient1(nb,3,B)
      end if

c     **** three derivatives ****
      if (ii.eq.10) then
         call Pack_cc_xyzgradient(nb,1,A,B)
         call Pack_cc_xyzgradient1(nb,2,B)
         call Pack_cc_xyzgradient1(nb,3,B)
      end if
      if (ii.eq.11) then
         call Pack_cc_xyzgradient(nb,1,A,B)
         call Pack_cc_xyzgradient1(nb,1,B)
         call Pack_cc_xyzgradient1(nb,2,B)
      end if
      if (ii.eq.12) then
         call Pack_cc_xyzgradient(nb,1,A,B)
         call Pack_cc_xyzgradient1(nb,1,B)
         call Pack_cc_xyzgradient1(nb,3,B)
      end if
      if (ii.eq.13) then
         call Pack_cc_xyzgradient(nb,1,A,B)
         call Pack_cc_xyzgradient1(nb,2,B)
         call Pack_cc_xyzgradient1(nb,2,B)
      end if
      if (ii.eq.14) then
         call Pack_cc_xyzgradient(nb,3,A,B)
         call Pack_cc_xyzgradient1(nb,2,B)
         call Pack_cc_xyzgradient1(nb,2,B)
      end if
      if (ii.eq.15) then
         call Pack_cc_xyzgradient(nb,1,A,B)
         call Pack_cc_xyzgradient1(nb,3,B)
         call Pack_cc_xyzgradient1(nb,3,B)
      end if
      if (ii.eq.16) then
         call Pack_cc_xyzgradient(nb,2,A,B)
         call Pack_cc_xyzgradient1(nb,3,B)
         call Pack_cc_xyzgradient1(nb,3,B)
      end if
      if (ii.eq.17) then
         call Pack_cc_xyzgradient(nb,1,A,B)
         call Pack_cc_xyzgradient1(nb,1,B)
         call Pack_cc_xyzgradient1(nb,1,B)
      end if
      if (ii.eq.18) then
         call Pack_cc_xyzgradient(nb,2,A,B)
         call Pack_cc_xyzgradient1(nb,2,B)
         call Pack_cc_xyzgradient1(nb,2,B)
      end if
      if (ii.eq.19) then
         call Pack_cc_xyzgradient(nb,3,A,B)
         call Pack_cc_xyzgradient1(nb,3,B)
         call Pack_cc_xyzgradient1(nb,3,B)
      end if

      return 
      end





*     ***********************************
*     *                                 *
*     *          Pack_t_dscal1          *
*     *                                 *
*     ***********************************

      subroutine Pack_t_dsca1l(nb,alpha,A)
      implicit none
      integer    nb
      real*8     alpha,A(*)


*     **** common block for pack ****
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
       A(i) = alpha*A(i)
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        C(i) = A(i)*B(i)
      end do
!$OMP END DO

      return
      end


*     ***********************************
*     *                                 *
*     *          Pack_ttt_dzaxpy        *
*     *                                 *
*     ***********************************
      subroutine Pack_ttt_dzaxpy(nb,A,B,C,D)
      implicit none
      integer  nb
      real*8   A(*)
      real*8   B(*)
      real*8   C(*)
      real*8   D(*)

*     **** common block for pack ****
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        D(i) = C(i) + A(i)*B(i)
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer i

!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        B(i) = B(i)*A(i)
      end do
!$OMP END DO

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
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer indx,p,q,taskid,id

      if (nb.lt.2) then
         id = 1
      else
         id = 3
      end if
      call Parallel2d_taskid_i(taskid)
      !call D3dB_ktoqp(1,1,q,p)
!$OMP MASTER
      call D3dB_ijktoindexp(id,1,1,1,indx,p)
      if (p.eq.taskid) A(indx) = dcmplx(vzero,0.0d0)
!$OMP END MASTER

      return
      end

*     ***********************************
*     *                                 *
*     *          Pack_c_addzero         *
*     *                                 *
*     ***********************************

      subroutine Pack_c_addzero(nb,vzero,A)
      implicit none
      integer    nb
      real*8     vzero
      complex*16 A(*)


*     **** common block for pack ****
      integer nida(0:3),nidb(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb

*     **** local variables ****
      integer indx,p,q,taskid,id

      if (nb.lt.2) then
         id = 1
      else
         id = 3
      end if
      call Parallel2d_taskid_i(taskid)
      !call D3dB_ktoqp(1,1,q,p)
!$OMP MASTER
      call D3dB_ijktoindexp(id,1,1,1,indx,p)
      if (p.eq.taskid) A(indx) = A(indx) + dcmplx(vzero,0.0d0)
!$OMP END MASTER

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

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block for pack ****
      integer nida(0:3),nidb2(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb2
      integer nidb(0:3)
      common / pack2_blk / nidb

      
*     **** local variables ****
      logical value
      integer nfft3d,n,id

*     **** external functions ****
      logical  control_balance
      external control_balance

      call nwpw_timing_start(9)
      if (nb.lt.2) then
         id = 1
      else
         id = 3
      end if

      call D3dB_nfft3d(id,nfft3d)

      !call dcopy(2*nfft3d,A,1,tmp1,1)
      !call dcopy(2*nfft3d,0.0d0,0,A,1)
      call Parallel_shared_vector_copy(.true.,2*nfft3d,A,tmp1)
      call Parallel_shared_vector_zero(.true.,2*nfft3d,A)

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

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block for pack ****
      integer nida(0:3),nidb2(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb2
      integer nidb(0:3)
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

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block for pack ****
      integer nida(0:3),nidb2(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb2
      integer nidb(0:3)
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

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block for pack ****
      integer nida(0:3),nidb2(0:3)
      integer pack(2,0:3)
      common / pack_blk / pack,nida,nidb2
      integer nidb(0:3)
      common / pack2_blk / nidb


*     **** local variables ****
      logical value
      integer i,nfft3d,id

*     **** external functions ****
      logical  control_balance
      external control_balance

      call nwpw_timing_start(9)

      if (nb.lt.2) then
         id = 1
      else
         id = 3
      end if

      if (control_balance()) 
     >  call Balance_c_unbalance_end(nb,tmp1,request,reqcnt)

      
      call D3dB_nfft3d(id,nfft3d)

c      call dcopy(2*(nida(nb)+nidb(nb)),tmp1,1,tmp2,1)
c      call dcopy(2*nfft3d,0.0d0,0,tmp1,1)
c      do i=1,(nida(nb)+nidb(nb))
c        tmp1(int_mb(pack(1,nb)+i-1)) = tmp2(i)
c      end do
      call Parallel_shared_vector_copy(.true.,2*(nida(nb)+nidb(nb)),
     >                                 tmp1,tmp2)
      call Parallel_shared_vector_zero(.true.,2*nfft3d,tmp1)
!$OMP DO
      do i=1,(nida(nb)+nidb(nb))
        tmp1(int_mb(pack(1,nb)+i-1)) = tmp2(i)
      end do
!$OMP END DO


*     **** make the kx=0 plane complete **** 
      call D3dB_c_timereverse_start(id,tmp1,tmp1z,tmp2z,
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

      if (nb.lt.2) then
         call D3dB_c_timereverse_end(1,tmp1,tmp1z,tmp2z,request,reqcnt)
      else
         call D3dB_c_timereverse_end(3,tmp1,tmp1z,tmp2z,request,reqcnt)
      end if

      call nwpw_timing_end(9)
  
      return
      end



