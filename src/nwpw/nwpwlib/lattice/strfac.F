C     ***********************************
C     *					*
C     *		strfac_init		*
C     *					*
C     ***********************************

      subroutine strfac_init()
*
* $Id: strfac.F,v 1.1 2001-08-30 18:31:38 bylaska Exp $
*
      implicit none

#include	"mafdecls.fh"

*     **** local common block ****
c     complex*16 ex1(0:nfft1-1,natmx)
c     complex*16 ex2(0:nfft2-1,natmx)
c     complex*16 ex3(0:nfft3-1,natmx)
      integer ex1(2),ex2(2),ex3(2)
      common / strfac_block / ex1,ex2,ex3

      integer nion,nx,ny,nz
      logical value

*     ***** external functions ****
      integer  ion_nion
      external ion_nion

      nion = ion_nion()
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
       

      value = MA_alloc_get(mt_dcpl,(nx*nion),'ex1',ex1(2),ex1(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,(ny*nion),'ex2',ex2(2),ex2(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,(nz*nion),'ex3',ex3(2),ex3(1))
      if (.not. value) call errquit('out of heap memory',0)

      return
      end

C     ***********************************
C     *					*
C     *		strfac_end 		*
C     *					*
C     ***********************************

      subroutine strfac_end()
      implicit none

#include "mafdecls.fh"

*     **** local common block ****
c     complex*16 ex1(0:nfft1-1,natmx)
c     complex*16 ex2(0:nfft2-1,natmx)
c     complex*16 ex3(0:nfft3-1,natmx)
      integer ex1(2),ex2(2),ex3(2)
      common / strfac_block / ex1,ex2,ex3

      logical value

      value = MA_free_heap(ex1(2))
      value = MA_free_heap(ex2(2))
      value = MA_free_heap(ex3(2))
       
      return
      end

C     ***********************************
C     *					*
C     *		PHAFAC			*
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
      subroutine phafac()
      implicit none

#include	"mafdecls.fh"

*     **** local common block ****
c     complex*16 ex1(0:nfft1-1,natmx)
c     complex*16 ex2(0:nfft2-1,natmx)
c     complex*16 ex3(0:nfft3-1,natmx)
      integer ex1(2),ex2(2),ex3(2)
      common / strfac_block / ex1,ex2,ex3



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

*     1 - fft_time
*     2 - dot_time
*     3 - ABC_time
*     4 - L_time
*     5 - vlocal_time
*     6 - vnl_time
*     7 - vh_time
*     8 - strfac_time
*     9 - maskC_time
      real*8 tim1,tim2
      real*8 times(50)
      common / timing / times

      call current_second(tim1)

      pi = 4.0d0*datan(1.0d0)

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)

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

      call current_second(tim2)
      times(8) = times(8) + (tim2-tim1)

      return
      end

C .. End PhaFac .... 






C      **********************************
C      *				*
C      *	STRFAC 			*
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
      subroutine strfac(II,exi)
      implicit none
      integer    II
      complex*16 exi(*)

#include "mafdecls.fh"
   
*     **** local common block ****
c     complex*16 ex1(0:nfft1-1,natmx)
c     complex*16 ex2(0:nfft2-1,natmx)
c     complex*16 ex3(0:nfft3-1,natmx)
      integer ex1(2),ex2(2),ex3(2)
      common / strfac_block / ex1,ex2,ex3


*     1 - fft_time
*     2 - dot_time
*     3 - ABC_time
*     4 - L_time
*     5 - vlocal_time
*     6 - vnl_time
*     7 - vh_time
*     8 - strfac_time
*     9 - maskC_time
      real*8 tim1,tim2
      real*8 times(50)
      common / timing / times



*     **** local variables ****
      integer i,j,k,index
      integer nx,ny,nz
      integer nxh
      integer p,q,taskid

      call current_second(tim1)

      call Parallel_taskid(taskid)
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      nxh=nx/2

*::::::::::::::::  structure factor ::::::::::::::::::::::::::::::::
      do k=0,nz-1
         call D3dB_ktoqp(1,(k+1),q,p)
         if (p.eq.taskid) then
            do j=0,ny-1
            do i=0,nxh
                index = (q-1)*(nxh+1)*ny
     >                + j*(nxh+1)
     >                + i+1
                exi(index) = dcpl_mb(ex1(1)+i+(II-1)*nx)
     >                     * dcpl_mb(ex2(1)+j+(II-1)*ny)
     >                     * dcpl_mb(ex3(1)+k+(II-1)*nz)
            end do
            end do
         end if
      end do

      call current_second(tim2)
      times(8) = times(8) + (tim2-tim1)


      return
      end

C ... End STRFAC ....

