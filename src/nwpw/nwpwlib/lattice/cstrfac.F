*
* $Id$
*

*     ***********************************
*     *					*
*     *		cstrfac_init		*
*     *					*
*     ***********************************

      subroutine cstrfac_init()
      implicit none
#include "errquit.fh"

#include "bafdecls.fh"
#include "cstrfac_common.fh"

      integer nion,nx,ny,nz,nbrillq
      integer nb,i,j,k,index,taskid,p,q
      integer nfft3d
      integer i_indx(2),j_indx(2),k_indx(2)
      logical value

*     ***** external functions ****
      integer  ion_nion,brillioun_nbrillq
      external ion_nion,brillioun_nbrillq

      nion       = ion_nion()
      nbrillq    = brillioun_nbrillq()
      call Parallel3d_taskid_i(taskid)
      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
       

      value = BA_alloc_get(mt_dcpl,(nx*nion),'ex1',ex1(2),ex1(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,(ny*nion),'ex2',ex2(2),ex2(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,(nz*nion),'ex3',ex3(2),ex3(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,(nbrillq*nion),
     >                    'cxr',cxr(2),cxr(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(2*nbrillq+2),
     >                    'i_indx_list',i_indx_list(2),i_indx_list(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(2*nbrillq+2),
     >                    'j_indx_list',j_indx_list(2),j_indx_list(1))
      value = value.and.
     >        BA_alloc_get(mt_int,(2*nbrillq+2),
     >                    'k_indx_list',k_indx_list(2),k_indx_list(1))
      if (.not. value) call errquit('cstrfac_init:out of heap',0,
     &       MA_ERR)


*     **** allocate indx blocks ****
      do nb=0,nbrillq
        call C3dB_nfft3d(1,nfft3d)
        value = BA_alloc_get(mt_int,nfft3d,'i_indx0',
     >                     i_indx(2),i_indx(1))
        value = value.and.
     >          BA_alloc_get(mt_int,nfft3d,'j_indx0',
     >                     j_indx(2),j_indx(1))
        value = value.and.
     >          BA_alloc_get(mt_int,nfft3d,'k_indx0',
     >                     k_indx(2),k_indx(1))
        if (.not. value) 
     >  call errquit('cstrfac_init:out of heap memory',nb+1, MA_ERR)

        int_mb(i_indx_list(1)+2*nb)  =i_indx(1)
        int_mb(i_indx_list(1)+2*nb+1)=i_indx(2)

        int_mb(j_indx_list(1)+2*nb)  =j_indx(1)
        int_mb(j_indx_list(1)+2*nb+1)=j_indx(2)

        int_mb(k_indx_list(1)+2*nb)  =k_indx(1) 
        int_mb(k_indx_list(1)+2*nb+1)=k_indx(2) 
      
        do k=0,nz-1
         !call C3dB_ktoqp(1,(k+1),q,p)
         !if (p.eq.taskid) then
            do j=0,ny-1
            do i=0,nx-1
              call C3dB_ijktoindexp(1,i+1,j+1,k+1,index,p)
              if (p.eq.taskid) then
c               index = (q-1)*(nx)*ny
c    >                + j*(nx)
c    >                + i+1
                int_mb(i_indx(1)+index-1) = i
                int_mb(j_indx(1)+index-1) = j
                int_mb(k_indx(1)+index-1) = k
              end if
            end do
            end do
         !end if
        end do
        call Cram_i_pack(nb,int_mb(i_indx(1)))
        call Cram_i_pack(nb,int_mb(j_indx(1)))
        call Cram_i_pack(nb,int_mb(k_indx(1)))
      end do



      return
      end

C     ***********************************
C     *					*
C     *		cstrfac_end 		*
C     *					*
C     ***********************************

      subroutine cstrfac_end()
      implicit none
#include "errquit.fh"

#include "bafdecls.fh"
#include "cstrfac_common.fh"

*     **** local variables ****
      logical value
      integer i_indx,j_indx,k_indx
      integer nb,nbrillq

*     **** external functions ****
      integer  brillioun_nbrillq
      external brillioun_nbrillq

      nbrillq = brillioun_nbrillq()

*     **** deallocate index blocks ****
      do nb=0,nbrillq
         i_indx=int_mb(i_indx_list(1)+2*nb+1)
         j_indx=int_mb(j_indx_list(1)+2*nb+1)
         k_indx=int_mb(k_indx_list(1)+2*nb+1)
         value =           BA_free_heap(i_indx)
         value = value.and.BA_free_heap(j_indx)
         value = value.and.BA_free_heap(k_indx)
         if (.not.value) call errquit('cstrfac_end:freeing heap',nb,
     &       MA_ERR)
      end do

      value =           BA_free_heap(i_indx_list(2))
      value = value.and.BA_free_heap(j_indx_list(2))
      value = value.and.BA_free_heap(k_indx_list(2))
      value = value.and.BA_free_heap(ex1(2))
      value = value.and.BA_free_heap(ex2(2))
      value = value.and.BA_free_heap(ex3(2))
      value = value.and.BA_free_heap(cxr(2))
      if (.not. value) call errquit('cstrfac_end:freeing heap',990,
     &       MA_ERR)
 
      return
      end

C     ***********************************
C     *					*
C     *		CPHAFAC			*
C     *					*
C     ***********************************
C
C    This file Sets up the phaze factors, that are used
C for creating structure factors.
C
C	Entry - 
C	Exit -
C
C  	Uses - ion data structure
C
      subroutine cphafac()
      implicit none

#include	"bafdecls.fh"
#include "cstrfac_common.fh"

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

      call nwpw_timing_start(8)

      pi = 4.0d0*datan(1.0d0)

      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)

*::::::::::::::::::  phase factor of ion positions  :::::::::::::::::::
      do i=1,ion_nion()
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

      call nwpw_timing_end(8)

      return
      end

C .. End PhaFac .... 




C     ***********************************
C     *					*
C     *		cphafac_rion
C     *					*
C     ***********************************
C
C    This file Sets up the phaze factors, that are used
C for creating structure factors.
C
C	Entry - 
C	Exit -
C
C  	Uses - ion data structure
C
      subroutine cphafac_rion(ion_rion)
      implicit none
      real*8 ion_rion(3,*)

#include	"bafdecls.fh"
#include "cstrfac_common.fh"

*     *** local variables
      integer    i,k
      integer    nx,ny,nz
      complex*16 cw1,cw2,cw3
      real*8     sw1,sw2,sw3
      real*8     pi

*     **** external functions ****
      integer  ion_nion
      real*8   lattice_unitg
      external ion_nion
      external lattice_unitg

      call nwpw_timing_start(8)

      pi = 4.0d0*datan(1.0d0)

      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)

*::::::::::::::::::  phase factor of ion positions  :::::::::::::::::::
      do i=1,ion_nion()
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

      call nwpw_timing_end(8)

      return
      end

C .. End cPhaFac_rion .... 








C      **********************************
C      *				*
C      *	CSTRFAC 		*
C      *				*
C      **********************************
C
C 	This routine returns the structure factor, EXI, for
C  a specified ion I.
C
C	Entry - II: the Ion index 
C	Exit  - EXI: the structure factor
C
C
      subroutine cstrfac(II,exi)
      implicit none
      integer    II
      complex*16 exi(*)

#include "bafdecls.fh"
#include "cstrfac_common.fh"
   
*     **** local variables ****
      integer i,j,k,index
      integer nx,ny,nz
      integer nxh
      integer p,q,taskid

      call nwpw_timing_start(8)

      call Parallel3d_taskid_i(taskid)
      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      nxh=nx/2

*::::::::::::::::  structure factor ::::::::::::::::::::::::::::::::
      do k=0,nz-1
         !call C3dB_ktoqp(1,(k+1),q,p)
         !if (p.eq.taskid) then
            do j=0,ny-1
            do i=0,nx-1
              call C3dB_ijktoindexp(1,i+1,j+1,k+1,index,p)
              if (p.eq.taskid) then
c               index = (q-1)*(nx)*ny
c    >                + j*(nx)
c    >                + i+1
                exi(index) = dcpl_mb(ex1(1)+i+(II-1)*nx)
     >                     * dcpl_mb(ex2(1)+j+(II-1)*ny)
     >                     * dcpl_mb(ex3(1)+k+(II-1)*nz)
              end if
            end do
            end do
         !end if
      end do

      call nwpw_timing_end(8)



      return
      end

C ... End STRFAC ....



C      **********************************
C      *				*
C      *	CSTRFAC_pack 		*
C      *				*
C      **********************************
C
C 	This routine returns the structure factor, EXI, for
C  a specified ion I.
C
C	Entry - II: the Ion index 
C	Exit  - EXI: the structure factor
C
C
      subroutine cstrfac_pack(nb,II,exi)
      implicit none
      integer    nb
      integer    II
      complex*16 exi(*)

#include "bafdecls.fh"
#include "cstrfac_common.fh"
   
*     **** local variables ****
      integer i,j,k,index
      integer nx,ny,nz,npack
      integer i_indx,j_indx,k_indx

      call nwpw_timing_start(8)

      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      call Cram_npack(nb,npack)
      i_indx = int_mb(i_indx_list(1)+2*nb)
      j_indx = int_mb(j_indx_list(1)+2*nb)
      k_indx = int_mb(k_indx_list(1)+2*nb)

*     **** crammed structure factor ****
      do index=1,npack
         i = int_mb(i_indx+index-1)
         j = int_mb(j_indx+index-1)
         k = int_mb(k_indx+index-1)

         exi(index) = dcpl_mb(ex1(1)+i+(II-1)*nx)
     >              * dcpl_mb(ex2(1)+j+(II-1)*ny)
     >              * dcpl_mb(ex3(1)+k+(II-1)*nz)
      end do

      call nwpw_timing_end(8)

      return
      end


*      **********************************
*      *				*
*      *	CSTRFAC_k 	        *
*      *				*
*      **********************************
*
* 	This routine returns the exp(ik.R_I) for
*  a specified ion I.
*
*	Entry - II: the Ion index 
*	        nb: the band index
*	Exit  - ck: the structure factor
*
*
      subroutine cstrfac_k(II,nb,ck)
      implicit none
      integer    II,nb
      complex*16 ck

#include "bafdecls.fh"
#include "cstrfac_common.fh"

      integer  nion
      integer  ion_nion
      external ion_nion

      nion = ion_nion()
      ck = dcpl_mb(cxr(1)+ii-1+(nb-1)*nion)
      return
      end

*      **********************************
*      *				*
*      *	   cphafac_k 		*
*      *				*
*      **********************************
      subroutine cphafac_k()
      implicit none

#include "bafdecls.fh"
#include "cstrfac_common.fh"

*     *** local variables
      integer    i,nb,nion

*     **** external functions ****
      integer  brillioun_nbrillq,ion_nion
      real*8   brillioun_k,ion_rion
      external brillioun_nbrillq,ion_nion
      external brillioun_k,ion_rion

      nion = ion_nion()
      do nb=1,brillioun_nbrillq()
      do i =1,nion
         dcpl_mb(cxr(1)+i-1+(nb-1)*nion) =
     >        dcmplx(cos( (brillioun_k(1,nb))*ion_rion(1,i)
     >                  + (brillioun_k(2,nb))*ion_rion(2,i)
     >                  + (brillioun_k(3,nb))*ion_rion(3,i)),
     >               sin( (brillioun_k(1,nb))*ion_rion(1,i)
     >                  + (brillioun_k(2,nb))*ion_rion(2,i)
     >                  + (brillioun_k(3,nb))*ion_rion(3,i)))
      end do
      end do

      return
      end

      subroutine cstrfac_k_print()
      implicit none

#include "bafdecls.fh"
#include "cstrfac_common.fh"

*     *** local variables
      integer    i,nb,nion

*     **** external functions ****
      integer  brillioun_nbrillq,ion_nion
      external brillioun_nbrillq,ion_nion

      write(*,*)
      write(*,*) "cstrfac_k_print:"
      nion = ion_nion()
      do nb=1,brillioun_nbrillq()
      do i =1,nion
        write(*,*) "cstrfac_k nb,ii,cxr:",nb,i,
     >             dcpl_mb(cxr(1)+i-1+(nb-1)*nion)
      end do
      end do
      write(*,*)
      write(*,*)

      return
      end



*     ***********************************
*     *                                 *
*     *        cstrfac_translate        *
*     *                                 *
*     ***********************************
      subroutine cstrfac_translate(trans,exi)
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
      call Parallel3d_taskid_i(taskid)

      pi = 4.0d0*datan(1.0d0)

      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      nxh = nx/2

      value = BA_push_get(mt_dcpl,(nx),'ex1',ex1(2),ex1(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(ny),'ex2',ex2(2),ex2(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,(nz),'ex3',ex3(2),ex3(1))
      if (.not. value)
     >  call errquit('cstrfac_translate:out of stack',0, MA_ERR)

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

*::::::::::::::::  structure factor ::::::::::::::::::::::::::::::::
      do k=0,nz-1
            do j=0,ny-1
            do i=0,nx-1
              call C3dB_ijktoindexp(1,i+1,j+1,k+1,index,p)
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
     >  call errquit('cstrfac_translate:error popping stack',0, MA_ERR)

      call nwpw_timing_end(8)

      return
      end


