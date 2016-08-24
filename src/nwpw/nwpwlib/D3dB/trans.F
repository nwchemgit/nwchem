#define NBLOCKS 3

*
* $Id$
*

*     ***********************************
*     *					*
*     *	   D3dB_c_transpose_ijk_init	*
*     *					*
*     ***********************************

      subroutine D3dB_c_transpose_ijk_init(nb)
      implicit none
      integer nb

#include "bafdecls.fh"
#include "errquit.fh"
#include "D3dB.fh"


*     *** hilbert tranpose data structure ****
      integer h_iq_to_i1(2,6,NBLOCKS)
      integer h_iq_to_i2(2,6,NBLOCKS)
      integer h_i1_start(2,6,NBLOCKS)
      integer h_i2_start(2,6,NBLOCKS)
      common / trans_blk_ijk / h_iq_to_i1,
     >                         h_iq_to_i2,
     >                         h_i1_start,
     >                         h_i2_start

*     **** local variables ****
      logical value
      integer proc_to,proc_from
      integer pto,qto,np,taskid
      integer phere,qhere
      integer index1,index2,itmp
      integer i,j,k,it

      integer tid,Parallel_threadid
      external Parallel_threadid

      call Parallel2d_taskid_i(taskid)
      call Parallel2d_np_i(np)


*     ********************************************************
*     **** map1to2 mapping - done - tranpose operation #1 ****
*     ********************************************************

*     **** allocate trans_blk_ijk common block ****
      value = BA_alloc_get(mt_int,((nx(nb)/2+1)*nq1(nb)),
     >                     'h_iq_to_i1_1',
     >                      h_iq_to_i1(2,1,nb),
     >                      h_iq_to_i1(1,1,nb))
      value = value.and.
     >        BA_alloc_get(mt_int,(ny(nb)*nq2(nb)),
     >                     'h_iq_to_i2_1',
     >                      h_iq_to_i2(2,1,nb),
     >                      h_iq_to_i2(1,1,nb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i1_start_1',
     >                      h_i1_start(2,1,nb),
     >                      h_i1_start(1,1,nb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i2_start_1',
     >                      h_i2_start(2,1,nb),
     >                      h_i2_start(1,1,nb))
      if (.not.value)
     > call errquit('D3dB_transpose_ijk_initt:out of heap',0,MA_ERR)


!$OMP SINGLE
cb      write(*,*) "Thread doing the job ",Parallel_threadid()
      index1 = 1 
      index2 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
         int_mb(h_i1_start(1,1,nb)+it) = index1
         int_mb(h_i2_start(1,1,nb)+it) = index2

         do k=1,nz(nb)
         do j=1,ny(nb)
         do i=1,(nx(nb)/2+1)

*           **** packing scheme **** 
            phere = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            qhere = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))

            pto   = int_mb(p_map2(1,nb)+(k-1)+(i-1)*nz(nb))
            qto   = int_mb(q_map2(1,nb)+(k-1)+(i-1)*nz(nb))
           
            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
               itmp = i + (qhere-1)*(nx(nb)/2+1)
c               if (itmp.ge.((nx(nb)/2+1)*nq1(nb))) then 
c                 call abort()
c               end if
               int_mb(h_iq_to_i1(1,1,nb)+itmp-1) = index1
               index1 = index1 + 1
            end if
*           **** unpacking scheme ****
            if ((pto.eq.taskid).and.(phere.eq.proc_from)) then
               itmp = j + (qto-1)*ny(nb)

c               if (itmp.ge.(ny(nb)*nq2(nb))) then 
c                 call abort()
c               end if
               int_mb(h_iq_to_i2(1,1,nb)+itmp-1) = index2
               index2 = index2 + 1
            end if

         end do
         end do
         end do

      end do
      int_mb(h_i1_start(1,1,nb)+np) = index1
      int_mb(h_i2_start(1,1,nb)+np) = index2
!$OMP END SINGLE


*     *********************************************************
*     **** map2to3 mapping - done - transpose operation #2 ****
*     *********************************************************

*     **** allocate trans_blk_ijk common block ****
      value = BA_alloc_get(mt_int,(ny(nb)*nq2(nb)),
     >                     'h_iq_to_i1_2',
     >                      h_iq_to_i1(2,2,nb),
     >                      h_iq_to_i1(1,2,nb))
      value = value.and.
     >        BA_alloc_get(mt_int,(nz(nb)*nq3(nb)),
     >                     'h_iq_to_i2_2',
     >                      h_iq_to_i2(2,2,nb),
     >                      h_iq_to_i2(1,2,nb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i1_start_2',
     >                      h_i1_start(2,2,nb),
     >                      h_i1_start(1,2,nb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i2_start_2',
     >                      h_i2_start(2,2,nb),
     >                      h_i2_start(1,2,nb))
      if (.not.value)
     > call errquit('D3dB_transpose_ijk_init:out of heap',0,MA_ERR)

!$OMP SINGLE 
      index1 = 1
      index2 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
         int_mb(h_i1_start(1,2,nb)+it) = index1
         int_mb(h_i2_start(1,2,nb)+it) = index2

         do k=1,nz(nb)
         do j=1,ny(nb)
         do i=1,(nx(nb)/2+1)

*           **** packing scheme ****
            phere = int_mb(p_map2(1,nb)+(k-1)+(i-1)*nz(nb))
            qhere = int_mb(q_map2(1,nb)+(k-1)+(i-1)*nz(nb))

            pto   = int_mb(p_map3(1,nb)+(i-1)+(j-1)*(nx(nb)/2+1))
            qto   = int_mb(q_map3(1,nb)+(i-1)+(j-1)*(nx(nb)/2+1))

            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
               itmp = j + (qhere-1)*ny(nb)
               int_mb(h_iq_to_i1(1,2,nb)+itmp-1) = index1
               index1 = index1 + 1
            end if

*           **** unpacking scheme ****
            if ((pto.eq.taskid).and.(phere.eq.proc_from)) then
               itmp = k + (qto-1)*nz(nb)
               int_mb(h_iq_to_i2(1,2,nb)+itmp-1) = index2
               index2 = index2 + 1
            end if

         end do
         end do
         end do

      end do
      int_mb(h_i1_start(1,2,nb)+np) = index1
      int_mb(h_i2_start(1,2,nb)+np) = index2


!$OMP END SINGLE 


*     ********************************************************
*     **** map3to2 mapping - done - tranpose operation #3 ****
*     ********************************************************

*     **** allocate trans_blk_ijk common block ****
      value = BA_alloc_get(mt_int,(nz(nb)*nq3(nb)),
     >                     'h_iq_to_i1_3',
     >                      h_iq_to_i1(2,3,nb),
     >                      h_iq_to_i1(1,3,nb))
      value = value.and.
     >        BA_alloc_get(mt_int,(ny(nb)*nq2(nb)),
     >                     'h_iq_to_i2_3',
     >                      h_iq_to_i2(2,3,nb),
     >                      h_iq_to_i2(1,3,nb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i1_start_3',
     >                      h_i1_start(2,3,nb),
     >                      h_i1_start(1,3,nb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i2_start_3',
     >                      h_i2_start(2,3,nb),
     >                      h_i2_start(1,3,nb))
      if (.not.value)
     > call errquit('D3dB_transpose_ijk_init:out of heap',0,MA_ERR)

!$OMP SINGLE 
      index1 = 1
      index2 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
         int_mb(h_i1_start(1,3,nb)+it) = index1
         int_mb(h_i2_start(1,3,nb)+it) = index2

         do k=1,nz(nb)
         do j=1,ny(nb)
         do i=1,(nx(nb)/2+1)

*           **** packing scheme ****
            phere = int_mb(p_map3(1,nb)+(i-1)+(j-1)*(nx(nb)/2+1))
            qhere = int_mb(q_map3(1,nb)+(i-1)+(j-1)*(nx(nb)/2+1))

            pto   = int_mb(p_map2(1,nb)+(k-1)+(i-1)*nz(nb))
            qto   = int_mb(q_map2(1,nb)+(k-1)+(i-1)*nz(nb))


            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
               itmp = k + (qhere-1)*nz(nb)
               int_mb(h_iq_to_i1(1,3,nb)+itmp-1) = index1
               index1 = index1 + 1
            end if

*           **** unpacking scheme ****
            if ((pto.eq.taskid).and.(phere.eq.proc_from)) then
               itmp = j + (qto-1)*ny(nb)
               int_mb(h_iq_to_i2(1,3,nb)+itmp-1) = index2
               index2 = index2 + 1
            end if

         end do
         end do
         end do

      end do
      int_mb(h_i1_start(1,3,nb)+np) = index1
      int_mb(h_i2_start(1,3,nb)+np) = index2

 
!$OMP END SINGLE 

*     ********************************************************
*     **** map2to1 mapping - done - tranpose operation #4 ****
*     ********************************************************

*     **** allocate trans_blk_ijk common block ****
      value = BA_alloc_get(mt_int,(ny(nb)*nq2(nb)),
     >                     'h_iq_to_i1_4',
     >                      h_iq_to_i1(2,4,nb),
     >                      h_iq_to_i1(1,4,nb))
      value = value.and.
     >        BA_alloc_get(mt_int,((nx(nb)/2+1)*nq1(nb)),
     >                     'h_iq_to_i2_4',
     >                      h_iq_to_i2(2,4,nb),
     >                      h_iq_to_i2(1,4,nb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i1_start_4',
     >                      h_i1_start(2,4,nb),
     >                      h_i1_start(1,4,nb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i2_start_4',
     >                      h_i2_start(2,4,nb),
     >                      h_i2_start(1,4,nb))
      if (.not.value)
     > call errquit('D3dB_transpose_ijk_init:out of heap',0,MA_ERR)

!$OMP SINGLE 
      index1 = 1
      index2 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
         int_mb(h_i1_start(1,4,nb)+it) = index1
         int_mb(h_i2_start(1,4,nb)+it) = index2

         do k=1,nz(nb)
         do j=1,ny(nb)
         do i=1,(nx(nb)/2+1)

*           **** packing scheme ****
            phere = int_mb(p_map2(1,nb)+(k-1)+(i-1)*nz(nb))
            qhere = int_mb(q_map2(1,nb)+(k-1)+(i-1)*nz(nb))

            pto   = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            qto   = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))


            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
               itmp = j + (qhere-1)*ny(nb)
               int_mb(h_iq_to_i1(1,4,nb)+itmp-1) = index1
               index1 = index1 + 1
            end if

*           **** unpacking scheme ****
            if ((pto.eq.taskid).and.(phere.eq.proc_from)) then
               itmp = i + (qto-1)*(nx(nb)/2+1)
               int_mb(h_iq_to_i2(1,4,nb)+itmp-1) = index2
               index2 = index2 + 1
            end if

         end do
         end do
         end do

      end do
      int_mb(h_i1_start(1,4,nb)+np) = index1
      int_mb(h_i2_start(1,4,nb)+np) = index2


!$OMP END SINGLE 


*     **********************************************************
*     **** map1to3 mapping  - done - tranpose operation # 5 ****
*     **********************************************************

*     **** allocate trans_blk_ijk common block ****
      value = BA_alloc_get(mt_int,((nx(nb)/2+1)*nq1(nb)),
     >                     'h_iq_to_i1_5',
     >                      h_iq_to_i1(2,5,nb),
     >                      h_iq_to_i1(1,5,nb))
      value = value.and.
     >        BA_alloc_get(mt_int,(nz(nb)*nq3(nb)),
     >                     'h_iq_to_i2_5',
     >                      h_iq_to_i2(2,5,nb),
     >                      h_iq_to_i2(1,5,nb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i1_start_5',
     >                      h_i1_start(2,5,nb),
     >                      h_i1_start(1,5,nb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i2_start_5',
     >                      h_i2_start(2,5,nb),
     >                      h_i2_start(1,5,nb))
      if (.not.value)
     > call errquit('D3dB_transpose_ijk_init:out of heap',0,MA_ERR)
 
!$OMP SINGLE 
      index1 = 1
      index2 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
         int_mb(h_i1_start(1,5,nb)+it) = index1
         int_mb(h_i2_start(1,5,nb)+it) = index2

         do k=1,nz(nb)
         do j=1,ny(nb)
         do i=1,(nx(nb)/2+1)

*           **** packing scheme ****
            phere = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            qhere = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))

            pto   = int_mb(p_map3(1,nb)+(i-1)+(j-1)*(nx(nb)/2+1))
            qto   = int_mb(q_map3(1,nb)+(i-1)+(j-1)*(nx(nb)/2+1))


            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
               itmp = i + (qhere-1)*(nx(nb)/2+1)
               int_mb(h_iq_to_i1(1,5,nb)+itmp-1) = index1
               index1 = index1 + 1
            end if

*           **** unpacking scheme ****
            if ((pto.eq.taskid).and.(phere.eq.proc_from)) then
               itmp = k + (qto-1)*nz(nb)
               int_mb(h_iq_to_i2(1,5,nb)+itmp-1) = index2
               index2 = index2 + 1
            end if

         end do
         end do
         end do

      end do
      int_mb(h_i1_start(1,5,nb)+np) = index1
      int_mb(h_i2_start(1,5,nb)+np) = index2


!$OMP END SINGLE 


*     *************************
*     **** map3to1 mapping ****
*     *************************

*     **** allocate trans_blk_ijk common block ****
      value = BA_alloc_get(mt_int,(nz(nb)*nq3(nb)),
     >                     'h_iq_to_i1_6',
     >                      h_iq_to_i1(2,6,nb),
     >                      h_iq_to_i1(1,6,nb))
      value = value.and.
     >        BA_alloc_get(mt_int,((nx(nb)/2+1)*nq1(nb)),
     >                     'h_iq_to_i2_6',
     >                      h_iq_to_i2(2,6,nb),
     >                      h_iq_to_i2(1,6,nb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i1_start_6',
     >                      h_i1_start(2,6,nb),
     >                      h_i1_start(1,6,nb))
      value = value.and.
     >        BA_alloc_get(mt_int,(np+1),
     >                     'h_i2_start_6',
     >                      h_i2_start(2,6,nb),
     >                      h_i2_start(1,6,nb))
      if (.not.value)
     > call errquit('D3dB_transpose_ijk_init:out of heap',0,MA_ERR)

!$OMP SINGLE 
      index1 = 1
      index2 = 1
      do it=0,np-1
         proc_to   = mod(taskid+it,np)
         proc_from = mod(taskid-it+np,np)
         int_mb(h_i1_start(1,6,nb)+it) = index1
         int_mb(h_i2_start(1,6,nb)+it) = index2

         do k=1,nz(nb)
         do j=1,ny(nb)
         do i=1,(nx(nb)/2+1)

*           **** packing scheme ****
            phere = int_mb(p_map3(1,nb)+(i-1)+(j-1)*(nx(nb)/2+1))
            qhere = int_mb(q_map3(1,nb)+(i-1)+(j-1)*(nx(nb)/2+1))

            pto   = int_mb(p_map1(1,nb)+(j-1)+(k-1)*ny(nb))
            qto   = int_mb(q_map1(1,nb)+(j-1)+(k-1)*ny(nb))


            if ((phere.eq.taskid).and.(pto.eq.proc_to)) then
               itmp = k + (qhere-1)*nz(nb)
               int_mb(h_iq_to_i1(1,6,nb)+itmp-1) = index1
               index1 = index1 + 1
            end if

*           **** unpacking scheme ****
            if ((pto.eq.taskid).and.(phere.eq.proc_from)) then
               itmp = i + (qto-1)*(nx(nb)/2+1)
               int_mb(h_iq_to_i2(1,6,nb)+itmp-1) = index2
               index2 = index2 + 1
            end if

         end do
         end do
         end do

      end do
      int_mb(h_i1_start(1,6,nb)+np) = index1
      int_mb(h_i2_start(1,6,nb)+np) = index2

!$OMP END SINGLE 

      return
      end
