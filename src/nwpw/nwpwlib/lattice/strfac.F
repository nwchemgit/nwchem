*
* $Id$
*

*     ***********************************
*     *                                 *
*     *         strfac_init             *
*     *                                 *
*     ***********************************

      subroutine strfac_init()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "strfac_common.fh"

      integer nion,nx,ny,nz,nxh,nfft3d
      integer i,j,k,p,q,index,nb,taskid
      logical value

*     ***** external functions ****
      integer  ion_nion,ion_nion_q
      external ion_nion,ion_nion_q

      call Parallel2d_taskid_i(taskid)

*     **** allocate phase factors ****
      nion = ion_nion()+ion_nion_q()
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      nxh = nx/2
    
      value = BA_alloc_get(mt_dcpl,(nx*nion),'ex1',ex1(2),ex1(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,(ny*nion),'ex2',ex2(2),ex2(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,(nz*nion),'ex3',ex3(2),ex3(1))
      if (.not. value) 
     >  call errquit('strfac_init:out of heap memory',0, MA_ERR)

*     **** allocate indx blocks ****
      do nb=0,1
        call D3dB_nfft3d(1,nfft3d)
        value = BA_alloc_get(mt_int,nfft3d,'i_indx0',
     >                     i_indx(2,nb),i_indx(1,nb))
        value = value.and.
     >          BA_alloc_get(mt_int,nfft3d,'j_indx0',
     >                     j_indx(2,nb),j_indx(1,nb))
        value = value.and.
     >          BA_alloc_get(mt_int,nfft3d,'k_indx0',
     >                     k_indx(2,nb),k_indx(1,nb))
        if (.not. value) 
     >  call errquit('strfac_init:out of heap memory',nb+1, MA_ERR)

        do k=0,nz-1
         !call D3dB_ktoqp(1,(k+1),q,p)
         !if (p.eq.taskid) then
         do j=0,ny-1
         do i=0,nxh
            call D3dB_ijktoindexp(1,i+1,j+1,k+1,index,p)
            if (p.eq.taskid) then
c             index = (q-1)*(nxh+1)*ny
c    >              + j*(nxh+1)
c    >              + i+1
              int_mb(i_indx(1,nb)+index-1) = i+1
              int_mb(j_indx(1,nb)+index-1) = j+1
              int_mb(k_indx(1,nb)+index-1) = k+1
            end if
         end do
         end do
        end do
        call Pack_i_pack(nb,int_mb(i_indx(1,nb)))
        call Pack_i_pack(nb,int_mb(j_indx(1,nb)))
        call Pack_i_pack(nb,int_mb(k_indx(1,nb)))
      end do

      !**** generate phase factors for charges ****
      if (ion_nion_q().gt.0) call phafac_q()

      return
      end


C     ***********************************
C     *                                 *
C     *         strfac_end              *
C     *                                 *
C     ***********************************

      subroutine strfac_end()
      implicit none
#include "errquit.fh"

#include "bafdecls.fh"
#include "strfac_common.fh"

      logical value
      integer nb

      value = BA_free_heap(ex1(2))
      value = value.and.BA_free_heap(ex2(2))
      value = value.and.BA_free_heap(ex3(2))
      if (.not. value) 
     >  call errquit('strfac_end:error freeing heap memory',0, MA_ERR)

      do nb=0,1
         value = BA_free_heap(i_indx(2,nb))
         value = value.and.BA_free_heap(j_indx(2,nb))
         value = value.and.BA_free_heap(k_indx(2,nb))
         if (.not. value) 
     >   call errquit('strfac_end:error freeing heap memory',nb, MA_ERR)
      end do

      return
      end


C     ***********************************
C     *                                 *
C     *         PHAFAC                  *
C     *                                 *
C     ***********************************
C
C    This file Sets up the phaze factors, that are used
C for creating structure factors.
C
C       Entry - 
C       Exit -
C
C       Uses - ion data structure
C
      subroutine phafac()
      implicit none

#include "bafdecls.fh"
#include "strfac_common.fh"

*     *** local variables
      integer    i,k
      integer    nx,ny,nz
      complex*16 cw1,cw2,cw3
      real*8     sw1,sw2,sw3
      real*8     pi

*     **** external functions ****
      integer  ion_nion
      real*8   lattice_unitg,ion_rion
      external ion_nion
      external lattice_unitg,ion_rion

      integer tid
      integer nthr
      integer  Parallel_threadid, Parallel_nthreads
      external Parallel_threadid, Parallel_nthreads

      call nwpw_timing_start(20)

      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()


      pi = 4.0d0*datan(1.0d0)

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)

*::::::::::::::::::  phase factor of ion positions  :::::::::::::::::::
!      do i=1,ion_nion()
      do i=tid+1,ion_nion(),nthr

        sw1= lattice_unitg(1,1)*ion_rion(1,i)
     >     + lattice_unitg(2,1)*ion_rion(2,i)
     >     + lattice_unitg(3,1)*ion_rion(3,i) + pi
        sw2= lattice_unitg(1,2)*ion_rion(1,i)
     >     + lattice_unitg(2,2)*ion_rion(2,i)
     >     + lattice_unitg(3,2)*ion_rion(3,i) + pi
        sw3= lattice_unitg(1,3)*ion_rion(1,i)
     >     + lattice_unitg(2,3)*ion_rion(2,i)
     >     + lattice_unitg(3,3)*ion_rion(3,i) + pi

c      !***DEBUG ***
c       write(*,*) "Warning: paw box center being used:"
c       sw1= lattice_unitg(1,1)*ion_rion(1,i)
c    >     + lattice_unitg(2,1)*ion_rion(2,i)
c    >     + lattice_unitg(3,1)*ion_rion(3,i) 
c       sw2= lattice_unitg(1,2)*ion_rion(1,i)
c    >     + lattice_unitg(2,2)*ion_rion(2,i)
c    >     + lattice_unitg(3,2)*ion_rion(3,i) 
c       sw3= lattice_unitg(1,3)*ion_rion(1,i)
c    >     + lattice_unitg(2,3)*ion_rion(2,i)
c    >     + lattice_unitg(3,3)*ion_rion(3,i) 
c      !***DEBUG ***


        cw1=dcmplx(dcos(sw1),-dsin(sw1))
        cw2=dcmplx(dcos(sw2),-dsin(sw2))
        cw3=dcmplx(dcos(sw3),-dsin(sw3))
c       ex1(0,i)=dcmplx(1.0d0,0.0d0)
c       ex2(0,i)=dcmplx(1.0d0,0.0d0)
c       ex3(0,i)=dcmplx(1.0d0,0.0d0)
        dcpl_mb(ex1(1)+(i-1)*nx)=dcmplx(1.0d0,0.0d0)
        dcpl_mb(ex2(1)+(i-1)*ny)=dcmplx(1.0d0,0.0d0)
        dcpl_mb(ex3(1)+(i-1)*nz)=dcmplx(1.0d0,0.0d0)

        do k=1,nx/2
c         ex1(k,I)=ex1(k-1,I)*cw1
c         ex1(nx-k,I)=dconjg(ex1(k,i))
          dcpl_mb(ex1(1)+k+(i-1)*nx)
     >         = dcpl_mb(ex1(1)+k-1+(i-1)*nx)*cw1
          dcpl_mb(ex1(1)+nx-k+(i-1)*nx)
     >         =dconjg(dcpl_mb(ex1(1)+k+(i-1)*nx))
        end do

        do k=1,ny/2
c         ex2(k,i)=ex2(k-1,i)*cw2
c         ex2(ny-k,i)=dconjg(ex2(k,i))
          dcpl_mb(ex2(1)+k+(i-1)*ny)
     >         = dcpl_mb(ex2(1)+k-1+(i-1)*ny)*cw2
          dcpl_mb(ex2(1)+ny-k+(i-1)*ny)
     >         =dconjg(dcpl_mb(ex2(1)+k+(i-1)*ny))
        end do

        do k=1,nz/2
c         ex3(k,I)=ex3(k-1,I)*cw3
c         ex3(nz-k,i)=dconjg(ex3(k,i))
          dcpl_mb(ex3(1)+k+(i-1)*nz)
     >         = dcpl_mb(ex3(1)+k-1+(i-1)*nz)*cw3
          dcpl_mb(ex3(1)+nz-k+(i-1)*nz)
     >         =dconjg(dcpl_mb(ex3(1)+k+(i-1)*nz))
        end do

c       ex1(nx/2,i)=dcmplx(0.0d0, 0.0d0)
c       ex2(ny/2,i)=dcmplx(0.0d0, 0.0d0)
c       ex3(nz/2,i)=dcmplx(0.0d0, 0.0d0)
        dcpl_mb(ex1(1)+nx/2+(i-1)*nx)=dcmplx(0.0d0, 0.0d0)
        dcpl_mb(ex2(1)+ny/2+(i-1)*ny)=dcmplx(0.0d0, 0.0d0)
        dcpl_mb(ex3(1)+nz/2+(i-1)*nz)=dcmplx(0.0d0, 0.0d0)
      end do

      call nwpw_timing_end(20)

      return
      end


C .. End PhaFac .... 


C     ***********************************
C     *                                 *
C     *         PHAFAC_Q                *
C     *                                 *
C     ***********************************
C
C    This file Sets up the phaze factors, that are used
C for creating structure factors for charges.
C
C       Entry - 
C       Exit -
C
C       Uses - ion data structure
C
      subroutine phafac_q()
      implicit none

#include "bafdecls.fh"
#include "strfac_common.fh"

*     *** local variables
      integer    i,k,nion,nion_q
      integer    nx,ny,nz
      complex*16 cw1,cw2,cw3
      real*8     sw1,sw2,sw3
      real*8     pi

*     **** external functions ****
      integer  ion_nion,ion_nion_q
      real*8   lattice_unitg,ion_rion
      external ion_nion,ion_nion_q
      external lattice_unitg,ion_rion

      integer tid
      integer nthr
      integer  Parallel_threadid, Parallel_nthreads
      external Parallel_threadid, Parallel_nthreads

      nion_q = ion_nion_q()

      if (nion_q.gt.0) then
         call nwpw_timing_start(20)

         tid  = Parallel_threadid()
         nthr = Parallel_nthreads()

         pi = 4.0d0*datan(1.0d0)

         call D3dB_nx(1,nx)
         call D3dB_ny(1,ny)
         call D3dB_nz(1,nz)

         nion   = ion_nion()

*::::::::::::::::::  phase factor of ion positions  :::::::::::::::::::
         do i=tid+nion+1,nion+nion_q,nthr

           sw1= lattice_unitg(1,1)*ion_rion(1,i)
     >        + lattice_unitg(2,1)*ion_rion(2,i)
     >        + lattice_unitg(3,1)*ion_rion(3,i) + pi
           sw2= lattice_unitg(1,2)*ion_rion(1,i)
     >        + lattice_unitg(2,2)*ion_rion(2,i)
     >        + lattice_unitg(3,2)*ion_rion(3,i) + pi
           sw3= lattice_unitg(1,3)*ion_rion(1,i)
     >        + lattice_unitg(2,3)*ion_rion(2,i)
     >        + lattice_unitg(3,3)*ion_rion(3,i) + pi

           cw1=dcmplx(dcos(sw1),-dsin(sw1))
           cw2=dcmplx(dcos(sw2),-dsin(sw2))
           cw3=dcmplx(dcos(sw3),-dsin(sw3))
           dcpl_mb(ex1(1)+(i-1)*nx)=dcmplx(1.0d0,0.0d0)
           dcpl_mb(ex2(1)+(i-1)*ny)=dcmplx(1.0d0,0.0d0)
           dcpl_mb(ex3(1)+(i-1)*nz)=dcmplx(1.0d0,0.0d0)

           do k=1,nx/2
             dcpl_mb(ex1(1)+k+(i-1)*nx)
     >            = dcpl_mb(ex1(1)+k-1+(i-1)*nx)*cw1
             dcpl_mb(ex1(1)+nx-k+(i-1)*nx)
     >            =dconjg(dcpl_mb(ex1(1)+k+(i-1)*nx))
           end do

           do k=1,ny/2
             dcpl_mb(ex2(1)+k+(i-1)*ny)
     >            = dcpl_mb(ex2(1)+k-1+(i-1)*ny)*cw2
             dcpl_mb(ex2(1)+ny-k+(i-1)*ny)
     >            =dconjg(dcpl_mb(ex2(1)+k+(i-1)*ny))
           end do

           do k=1,nz/2
             dcpl_mb(ex3(1)+k+(i-1)*nz)
     >            = dcpl_mb(ex3(1)+k-1+(i-1)*nz)*cw3
             dcpl_mb(ex3(1)+nz-k+(i-1)*nz)
     >            =dconjg(dcpl_mb(ex3(1)+k+(i-1)*nz))
           end do

           dcpl_mb(ex1(1)+nx/2+(i-1)*nx)=dcmplx(0.0d0, 0.0d0)
           dcpl_mb(ex2(1)+ny/2+(i-1)*ny)=dcmplx(0.0d0, 0.0d0)
           dcpl_mb(ex3(1)+nz/2+(i-1)*nz)=dcmplx(0.0d0, 0.0d0)
         end do

         call nwpw_timing_end(20)
      end if

      return
      end

C .. End PhaFac_Q .... 



C     ***********************************
C     *                                 *
C     *         phafac2                 *
C     *                                 *
C     ***********************************
      subroutine phafac2()
      implicit none
#include "bafdecls.fh"
#include "strfac_common.fh"
*     *** local variables
      integer    i,k
      integer    nx,ny,nz,tid,nthr
      complex*16 cw1,cw2,cw3
      real*8     sw1,sw2,sw3
      real*8     pi

*     **** external functions ****
      integer  ion_nion
      real*8   lattice_unitg,ion_rion2
      external ion_nion
      external lattice_unitg,ion_rion2
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      call nwpw_timing_start(20)
      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()
      pi = 4.0d0*datan(1.0d0)
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)

      do i=tid+1,ion_nion(),nthr

        sw1= lattice_unitg(1,1)*ion_rion2(1,i)
     >     + lattice_unitg(2,1)*ion_rion2(2,i)
     >     + lattice_unitg(3,1)*ion_rion2(3,i) + pi
        sw2= lattice_unitg(1,2)*ion_rion2(1,i)
     >     + lattice_unitg(2,2)*ion_rion2(2,i)
     >     + lattice_unitg(3,2)*ion_rion2(3,i) + pi
        sw3= lattice_unitg(1,3)*ion_rion2(1,i)
     >     + lattice_unitg(2,3)*ion_rion2(2,i)
     >     + lattice_unitg(3,3)*ion_rion2(3,i) + pi

c      !***DEBUG ***
c       write(*,*) "Warning: paw box center being used:"
c       sw1= lattice_unitg(1,1)*ion_rion2(1,i)
c    >     + lattice_unitg(2,1)*ion_rion2(2,i)
c    >     + lattice_unitg(3,1)*ion_rion2(3,i) 
c       sw2= lattice_unitg(1,2)*ion_rion2(1,i)
c    >     + lattice_unitg(2,2)*ion_rion2(2,i)
c    >     + lattice_unitg(3,2)*ion_rion2(3,i) 
c       sw3= lattice_unitg(1,3)*ion_rion2(1,i)
c    >     + lattice_unitg(2,3)*ion_rion2(2,i)
c    >     + lattice_unitg(3,3)*ion_rion2(3,i) 
c      !***DEBUG ***

        cw1=dcmplx(dcos(sw1),-dsin(sw1))
        cw2=dcmplx(dcos(sw2),-dsin(sw2))
        cw3=dcmplx(dcos(sw3),-dsin(sw3))
        dcpl_mb(ex1(1)+(i-1)*nx)=dcmplx(1.0d0,0.0d0)
        dcpl_mb(ex2(1)+(i-1)*ny)=dcmplx(1.0d0,0.0d0)
        dcpl_mb(ex3(1)+(i-1)*nz)=dcmplx(1.0d0,0.0d0)
        do k=1,nx/2
          dcpl_mb(ex1(1)+k+(i-1)*nx)
     >         = dcpl_mb(ex1(1)+k-1+(i-1)*nx)*cw1
          dcpl_mb(ex1(1)+nx-k+(i-1)*nx)
     >         =dconjg(dcpl_mb(ex1(1)+k+(i-1)*nx))
        end do
        do k=1,ny/2
          dcpl_mb(ex2(1)+k+(i-1)*ny)
     >         = dcpl_mb(ex2(1)+k-1+(i-1)*ny)*cw2
          dcpl_mb(ex2(1)+ny-k+(i-1)*ny)
     >         =dconjg(dcpl_mb(ex2(1)+k+(i-1)*ny))
        end do
        do k=1,nz/2
          dcpl_mb(ex3(1)+k+(i-1)*nz)
     >         = dcpl_mb(ex3(1)+k-1+(i-1)*nz)*cw3
          dcpl_mb(ex3(1)+nz-k+(i-1)*nz)
     >         =dconjg(dcpl_mb(ex3(1)+k+(i-1)*nz))
        end do
        dcpl_mb(ex1(1)+nx/2+(i-1)*nx)=dcmplx(0.0d0, 0.0d0)
        dcpl_mb(ex2(1)+ny/2+(i-1)*ny)=dcmplx(0.0d0, 0.0d0)
        dcpl_mb(ex3(1)+nz/2+(i-1)*nz)=dcmplx(0.0d0, 0.0d0)
      end do
      call nwpw_timing_end(20)
      return
      end

C     ***********************************
C     *                                 *
C     *         phafac_rion             *
C     *                                 *
C     ***********************************
C
C    This file Sets up the phaze factors, that are used
C for creating structure factors.
C
C       Entry - 
C       Exit -
C
C       Uses - ion data structure
C
      subroutine phafac_rion(ion_rion)
      implicit none
      real*8 ion_rion(3,*)

#include        "bafdecls.fh"
#include "strfac_common.fh"

*     *** local variables
      integer    i,k,tid,nthr
      integer    nx,ny,nz
      complex*16 cw1,cw2,cw3
      real*8     sw1,sw2,sw3
      real*8     pi

*     **** external functions ****
      integer  ion_nion
      real*8   lattice_unitg
      external ion_nion
      external lattice_unitg
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      call nwpw_timing_start(20)
      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

      pi = 4.0d0*datan(1.0d0)

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)

*::::::::::::::::::  phase factor of ion positions  :::::::::::::::::::
      do i=tid+1,ion_nion(),nthr
        sw1= lattice_unitg(1,1)*ion_rion(1,i)
     >     + lattice_unitg(2,1)*ion_rion(2,i)
     >     + lattice_unitg(3,1)*ion_rion(3,i) + pi
        sw2= lattice_unitg(1,2)*ion_rion(1,i)
     >     + lattice_unitg(2,2)*ion_rion(2,i)
     >     + lattice_unitg(3,2)*ion_rion(3,i) + pi
        sw3= lattice_unitg(1,3)*ion_rion(1,i)
     >     + lattice_unitg(2,3)*ion_rion(2,i)
     >     + lattice_unitg(3,3)*ion_rion(3,i) + pi
        cw1=dcmplx(dcos(sw1),-dsin(sw1))
        cw2=dcmplx(dcos(sw2),-dsin(sw2))
        cw3=dcmplx(dcos(sw3),-dsin(sw3))
c       ex1(0,i)=dcmplx(1.0d0,0.0d0)
c       ex2(0,i)=dcmplx(1.0d0,0.0d0)
c       ex3(0,i)=dcmplx(1.0d0,0.0d0)
        dcpl_mb(ex1(1)+(i-1)*nx)=dcmplx(1.0d0,0.0d0)
        dcpl_mb(ex2(1)+(i-1)*ny)=dcmplx(1.0d0,0.0d0)
        dcpl_mb(ex3(1)+(i-1)*nz)=dcmplx(1.0d0,0.0d0)

        do k=1,nx/2
c         ex1(k,I)=ex1(k-1,I)*cw1
c         ex1(nx-k,I)=dconjg(ex1(k,i))
          dcpl_mb(ex1(1)+k+(i-1)*nx)
     >         = dcpl_mb(ex1(1)+k-1+(i-1)*nx)*cw1
          dcpl_mb(ex1(1)+nx-k+(i-1)*nx)
     >         =dconjg(dcpl_mb(ex1(1)+k+(i-1)*nx))
        end do

        do k=1,ny/2
c         ex2(k,i)=ex2(k-1,i)*cw2
c         ex2(ny-k,i)=dconjg(ex2(k,i))
          dcpl_mb(ex2(1)+k+(i-1)*ny)
     >         = dcpl_mb(ex2(1)+k-1+(i-1)*ny)*cw2
          dcpl_mb(ex2(1)+ny-k+(i-1)*ny)
     >         =dconjg(dcpl_mb(ex2(1)+k+(i-1)*ny))
        end do

        do k=1,nz/2
c         ex3(k,I)=ex3(k-1,I)*cw3
c         ex3(nz-k,i)=dconjg(ex3(k,i))
          dcpl_mb(ex3(1)+k+(i-1)*nz)
     >         = dcpl_mb(ex3(1)+k-1+(i-1)*nz)*cw3
          dcpl_mb(ex3(1)+nz-k+(i-1)*nz)
     >         =dconjg(dcpl_mb(ex3(1)+k+(i-1)*nz))
        end do

c       ex1(nx/2,i)=dcmplx(0.0d0, 0.0d0)
c       ex2(ny/2,i)=dcmplx(0.0d0, 0.0d0)
c       ex3(nz/2,i)=dcmplx(0.0d0, 0.0d0)
        dcpl_mb(ex1(1)+nx/2+(i-1)*nx)=dcmplx(0.0d0, 0.0d0)
        dcpl_mb(ex2(1)+ny/2+(i-1)*ny)=dcmplx(0.0d0, 0.0d0)
        dcpl_mb(ex3(1)+nz/2+(i-1)*nz)=dcmplx(0.0d0, 0.0d0)
      end do

      call nwpw_timing_end(20)

      return
      end

C .. End PhaFac_rion .... 






C      **********************************
C      *                                *
C      *        STRFAC                  *
C      *                                *
C      **********************************
C
C       This routine returns the structure factor, EXI, for
C  a specified ion I.
C
C       Entry - II: the Ion index 
C       Exit  - EXI: the structure factor
C
C
      subroutine strfac(II,exi)
      implicit none
      integer    II
      complex*16 exi(*)

#include "bafdecls.fh"
#include "strfac_common.fh"
   
*     **** local variables ****
      integer i,j,k,index
      integer nx,ny,nz
      integer nxh
      integer p,q,taskid,tid,nthr

*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      call nwpw_timing_start(8)
      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()
      call Parallel2d_taskid_i(taskid)
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      nxh=nx/2

*     **** structure factor ****
      do k=0,nz-1
         !call D3dB_ktoqp(1,(k+1),q,p)
         !if (p.eq.taskid) then
         do j=0,ny-1
         do i=0,nxh
            call D3dB_ijktoindexp(1,i+1,j+1,k+1,index,p)
            if ((p.eq.taskid).and.(mod(index,nthr).eq.tid)) then
c               index = (q-1)*(nxh+1)*ny
c    >                + j*(nxh+1)
c    >                + i+1
                exi(index) = dcpl_mb(ex1(1)+i+(II-1)*nx)
     >                     * dcpl_mb(ex2(1)+j+(II-1)*ny)
     >                     * dcpl_mb(ex3(1)+k+(II-1)*nz)
            end if
         end do
         end do
      end do
!$OMP BARRIER

      call nwpw_timing_end(8)


      return
      end




C      **********************************
C      *                                *
C      *          strfac_pack           *
C      *                                *
C      **********************************
C
C       This routine returns the structure factor, EXI, for
C  a specified ion I.
C
C       Entry - nb:
c               II: the Ion index 
C       Exit  - EXI: the structure factor
C
C
      subroutine strfac_pack(nb,II,exi)
      implicit none
      integer    nb
      integer    II
      complex*16 exi(*)

#include "bafdecls.fh"
#include "strfac_common.fh"
   
*     **** local variables ****
      integer i,j,k,index,npack
      integer nx,ny,nz

      call nwpw_timing_start(8)

      call Pack_npack(nb,npack)      
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)

     
c*     **** structure factor in packed space ****
c      do index=1,npack
c        i = int_mb(i_indx(1,nb)+index-1)
c        j = int_mb(j_indx(1,nb)+index-1)
c        k = int_mb(k_indx(1,nb)+index-1)
c      
c        exi(index) = dcpl_mb(ex1(1)+i+(II-1)*nx)
c     >             * dcpl_mb(ex2(1)+j+(II-1)*ny)
c     >             * dcpl_mb(ex3(1)+k+(II-1)*nz)
c      end do
      call strfac_pack_sub(npack,
     >                     int_mb(i_indx(1,nb)),
     >                     int_mb(j_indx(1,nb)),
     >                     int_mb(k_indx(1,nb)),
     >                     dcpl_mb(ex1(1)+(II-1)*nx),
     >                     dcpl_mb(ex2(1)+(II-1)*ny),
     >                     dcpl_mb(ex3(1)+(II-1)*nz),
     >                     exi)

      call nwpw_timing_end(8)
      return
      end

      subroutine strfac_pack_sub(npack,
     >                           indxi,indxj,indxk,
     >                           ex1,ex2,ex3,
     >                           exi)
      implicit none
      integer npack,indxi(*),indxj(*),indxk(*)
      complex*16 ex1(*),ex2(*),ex3(*),exi(*)

      integer i,j,k,indx
!$OMP DO
      do indx=1,npack
        exi(indx) = ex1(indxi(indx))*ex2(indxj(indx))*ex3(indxk(indx))
      end do
!$OMP END DO
      return
      end

C ... End STRFAC ....



*     ***********************************
*     *                                 *
*     *         qstrfac_init            *
*     *                                 *
*     ***********************************

      subroutine qstrfac_init(nion)
      implicit none
      integer nion

#include "bafdecls.fh"
#include "strfac_common.fh"
#include "errquit.fh"

      integer nx,ny,nz
      logical value


*     **** allocate phase factors ****
      qnion = nion
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
    
      value = BA_alloc_get(mt_dcpl,(nx*qnion),'qex1',qex1(2),qex1(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,(ny*qnion),'qex2',qex2(2),qex2(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,(nz*qnion),'qex3',qex3(2),qex3(1))
      if (.not. value) 
     >  call errquit('qstrfac_init:out of heap memory',0, MA_ERR)

      return
      end


C     ***********************************
C     *                                 *
C     *         qstrfac_end             *
C     *                                 *
C     ***********************************

      subroutine qstrfac_end()
      implicit none

#include "bafdecls.fh"
#include "strfac_common.fh"
#include "errquit.fh"

      logical value

      value =           BA_free_heap(qex1(2))
      value = value.and.BA_free_heap(qex2(2))
      value = value.and.BA_free_heap(qex3(2))
      if (.not. value) 
     >  call errquit('qstrfac_end:error freeing heap memory',0,MA_ERR)

      return
      end


C     ***********************************
C     *                                 *
C     *         qphafac_rion            *
C     *                                 *
C     ***********************************
C
C    This file Sets up the phaze factors, that are used
C for creating structure factors.
C
C       Entry - 
C       Exit -
C
C       Uses - ion data structure
C
      subroutine qphafac_rion(ion_rion)
      implicit none
      real*8 ion_rion(3,*)

#include        "bafdecls.fh"
#include "strfac_common.fh"

*     *** local variables
      integer    i,k,tid,nthr
      integer    nx,ny,nz
      complex*16 cw1,cw2,cw3
      real*8     sw1,sw2,sw3
      real*8     pi

*     **** external functions ****
      integer  ion_nion
      real*8   lattice_unitg
      external ion_nion
      external lattice_unitg
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      call nwpw_timing_start(20)
      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

      pi = 4.0d0*datan(1.0d0)

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)

*::::::::::::::::::  phase factor of ion positions  :::::::::::::::::::
      do i=tid+1,qnion,nthr
        sw1= lattice_unitg(1,1)*ion_rion(1,i)
     >     + lattice_unitg(2,1)*ion_rion(2,i)
     >     + lattice_unitg(3,1)*ion_rion(3,i) + pi
        sw2= lattice_unitg(1,2)*ion_rion(1,i)
     >     + lattice_unitg(2,2)*ion_rion(2,i)
     >     + lattice_unitg(3,2)*ion_rion(3,i) + pi
        sw3= lattice_unitg(1,3)*ion_rion(1,i)
     >     + lattice_unitg(2,3)*ion_rion(2,i)
     >     + lattice_unitg(3,3)*ion_rion(3,i) + pi
        cw1=dcmplx(dcos(sw1),-dsin(sw1))
        cw2=dcmplx(dcos(sw2),-dsin(sw2))
        cw3=dcmplx(dcos(sw3),-dsin(sw3))
        dcpl_mb(qex1(1)+(i-1)*nx)=dcmplx(1.0d0,0.0d0)
        dcpl_mb(qex2(1)+(i-1)*ny)=dcmplx(1.0d0,0.0d0)
        dcpl_mb(qex3(1)+(i-1)*nz)=dcmplx(1.0d0,0.0d0)

        do k=1,nx/2
          dcpl_mb(qex1(1)+k+(i-1)*nx)
     >         = dcpl_mb(qex1(1)+k-1+(i-1)*nx)*cw1
          dcpl_mb(qex1(1)+nx-k+(i-1)*nx)
     >         =dconjg(dcpl_mb(qex1(1)+k+(i-1)*nx))
        end do

        do k=1,ny/2
          dcpl_mb(qex2(1)+k+(i-1)*ny)
     >         = dcpl_mb(qex2(1)+k-1+(i-1)*ny)*cw2
          dcpl_mb(qex2(1)+ny-k+(i-1)*ny)
     >         =dconjg(dcpl_mb(qex2(1)+k+(i-1)*ny))
        end do

        do k=1,nz/2
          dcpl_mb(qex3(1)+k+(i-1)*nz)
     >         = dcpl_mb(qex3(1)+k-1+(i-1)*nz)*cw3
          dcpl_mb(qex3(1)+nz-k+(i-1)*nz)
     >         =dconjg(dcpl_mb(qex3(1)+k+(i-1)*nz))
        end do

        dcpl_mb(qex1(1)+nx/2+(i-1)*nx)=dcmplx(0.0d0, 0.0d0)
        dcpl_mb(qex2(1)+ny/2+(i-1)*ny)=dcmplx(0.0d0, 0.0d0)
        dcpl_mb(qex3(1)+nz/2+(i-1)*nz)=dcmplx(0.0d0, 0.0d0)
      end do

      call nwpw_timing_end(20)

      return
      end

C .. End PhaFac_rion .... 






C      **********************************
C      *                                *
C      *        qstrfac                 *
C      *                                *
C      **********************************
C
C       This routine returns the structure factor, EXI, for
C  a specified ion I.
C
C       Entry - II: the Ion index 
C       Exit  - EXI: the structure factor
C
C
      subroutine qstrfac(II,exi)
      implicit none
      integer    II
      complex*16 exi(*)

#include "bafdecls.fh"
#include "strfac_common.fh"
   
*     **** local variables ****
      integer i,j,k,index
      integer nx,ny,nz
      integer nxh
      integer p,q,taskid,tid,nthr

*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      call nwpw_timing_start(8)
      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

      call Parallel2d_taskid_i(taskid)
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      nxh=nx/2

*     **** structure factor ****
      do k=0,nz-1
         do j=0,ny-1
         do i=0,nxh
            call D3dB_ijktoindexp(1,i+1,j+1,k+1,index,p)
            if ((p.eq.taskid).and.(mod(index,nthr).eq.tid)) then
                exi(index) = dcpl_mb(qex1(1)+i+(II-1)*nx)
     >                     * dcpl_mb(qex2(1)+j+(II-1)*ny)
     >                     * dcpl_mb(qex3(1)+k+(II-1)*nz)
            end if
         end do
         end do
      end do
!$OMP BARRIER

      call nwpw_timing_end(8)


      return
      end




C      **********************************
C      *                                *
C      *          qstrfac_pack          *
C      *                                *
C      **********************************
C
C       This routine returns the structure factor, EXI, for
C  a specified ion I.
C
C       Entry - nb:
c               II: the Ion index 
C       Exit  - EXI: the structure factor
C
C
      subroutine qstrfac_pack(nb,II,exi)
      implicit none
      integer    nb
      integer    II
      complex*16 exi(*)

#include "bafdecls.fh"
#include "strfac_common.fh"
   
*     **** local variables ****
      integer i,j,k,index,npack
      integer nx,ny,nz

      call nwpw_timing_start(8)

      call Pack_npack(nb,npack)      
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)

     
*     **** structure factor in packed space ****
c      do index=1,npack
c        i = int_mb(i_indx(1,nb)+index-1)
c        j = int_mb(j_indx(1,nb)+index-1)
c        k = int_mb(k_indx(1,nb)+index-1)
c      
c        exi(index) = dcpl_mb(qex1(1)+i+(II-1)*nx)
c     >             * dcpl_mb(qex2(1)+j+(II-1)*ny)
c     >             * dcpl_mb(qex3(1)+k+(II-1)*nz)
c      end do
      call strfac_pack_sub(npack,
     >                     int_mb(i_indx(1,nb)),
     >                     int_mb(j_indx(1,nb)),
     >                     int_mb(k_indx(1,nb)),
     >                     dcpl_mb(qex1(1)+(II-1)*nx),
     >                     dcpl_mb(qex2(1)+(II-1)*ny),
     >                     dcpl_mb(qex3(1)+(II-1)*nz),
     >                     exi)


      call nwpw_timing_end(8)

      return
      end

C ... End qSTRFAC ....



*     ***********************************
*     *                                 *
*     *         strfac_translate        *
*     *                                 *
*     ***********************************
      subroutine strfac_translate(trans,exi)
      implicit none
      real*8 trans(3)
      complex*16 exi(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      integer tid,taskid
      integer nthr

      logical    value
      integer    i,j,k,index,p
      integer    nx,ny,nz,nxh
      complex*16 cw1,cw2,cw3
      real*8     sw1,sw2,sw3
      real*8     pi,ion_tran(3)
      integer ex1(2),ex2(2),ex3(2)

*     **** external functions ****
      real*8   lattice_unitg,lattice_unita
      external lattice_unitg,lattice_unita
      integer  Parallel_threadid, Parallel_nthreads
      external Parallel_threadid, Parallel_nthreads

      call nwpw_timing_start(20)

      ion_tran(1) = trans(1)
     >            - (   lattice_unita(1,1)
     >                + lattice_unita(1,2)
     >                + lattice_unita(1,3) )/2.0d0
      ion_tran(2) = trans(2)
     >            - (   lattice_unita(2,1)
     >                + lattice_unita(2,2)
     >                + lattice_unita(2,3) )/2.0d0
      ion_tran(3) = trans(3)
     >            - (   lattice_unita(3,1)
     >                + lattice_unita(3,2)
     >                + lattice_unita(3,3) )/2.0d0


      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()
      call Parallel2d_taskid_i(taskid)

      pi = 4.0d0*datan(1.0d0)

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      nxh = nx/2

      value = BA_push_get(mt_dcpl,(nx),'ex1',ex1(2),ex1(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(ny),'ex2',ex2(2),ex2(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(nz),'ex3',ex3(2),ex3(1))
      if (.not. value)
     >  call errquit('strfac_translate:out of stack',0, MA_ERR)

*::::::::::::::::::  phase factor of translation  :::::::::::::::::::
!$OMP SINGLE

        sw1 = lattice_unitg(1,1)*ion_tran(1)
     >      + lattice_unitg(2,1)*ion_tran(2)
     >      + lattice_unitg(3,1)*ion_tran(3) + pi
        sw2 = lattice_unitg(1,2)*ion_tran(1)
     >      + lattice_unitg(2,2)*ion_tran(2)
     >      + lattice_unitg(3,2)*ion_tran(3) + pi
        sw3 = lattice_unitg(1,3)*ion_tran(1)
     >      + lattice_unitg(2,3)*ion_tran(2)
     >      + lattice_unitg(3,3)*ion_tran(3) + pi

        cw1 = dcmplx(dcos(sw1),-dsin(sw1))
        cw2 = dcmplx(dcos(sw2),-dsin(sw2))
        cw3 = dcmplx(dcos(sw3),-dsin(sw3))
        dcpl_mb(ex1(1)) = dcmplx(1.0d0,0.0d0)
        dcpl_mb(ex2(1)) = dcmplx(1.0d0,0.0d0)
        dcpl_mb(ex3(1)) = dcmplx(1.0d0,0.0d0)

        do k=1,nx/2
          dcpl_mb(ex1(1)+k) = dcpl_mb(ex1(1)+k-1)*cw1
          dcpl_mb(ex1(1)+nx-k) = dconjg(dcpl_mb(ex1(1)+k))
        end do

        do k=1,ny/2
          dcpl_mb(ex2(1)+k) = dcpl_mb(ex2(1)+k-1)*cw2
          dcpl_mb(ex2(1)+ny-k) = dconjg(dcpl_mb(ex2(1)+k))
        end do

        do k=1,nz/2
          dcpl_mb(ex3(1)+k) = dcpl_mb(ex3(1)+k-1)*cw3
          dcpl_mb(ex3(1)+nz-k) = dconjg(dcpl_mb(ex3(1)+k))
        end do

        dcpl_mb(ex1(1)+nx/2) = dcmplx(0.0d0, 0.0d0)
        dcpl_mb(ex2(1)+ny/2) = dcmplx(0.0d0, 0.0d0)
        dcpl_mb(ex3(1)+nz/2) = dcmplx(0.0d0, 0.0d0)

!$OMP END SINGLE

      call nwpw_timing_end(20)
      call nwpw_timing_start(8)

*     **** structure factor ****
      do k=0,nz-1
         do j=0,ny-1
         do i=0,nxh
            call D3dB_ijktoindexp(1,i+1,j+1,k+1,index,p)
            if ((p.eq.taskid).and.(mod(index,nthr).eq.tid)) then
                exi(index) = dcpl_mb(ex1(1)+i)
     >                     * dcpl_mb(ex2(1)+j)
     >                     * dcpl_mb(ex3(1)+k)
            end if
         end do
         end do
      end do
!$OMP BARRIER

      value = BA_pop_stack(ex3(2))
      value = value.and.BA_pop_stack(ex2(2))
      value = value.and.BA_pop_stack(ex1(2))
      if (.not. value)
     >  call errquit('strfac_translate:error popping stack',0, MA_ERR)

      call nwpw_timing_end(8)

      return
      end


