*
* $Id: cstrfac.F,v 1.5 2002-01-12 00:50:16 bylaska Exp $
*

*     ***********************************
*     *					*
*     *		cstrfac_init		*
*     *					*
*     ***********************************

      subroutine cstrfac_init()
      implicit none

#include "mafdecls.fh"
#include "cstrfac_common.fh"

      integer nion,nx,ny,nz,nbrillioun
      logical value

*     ***** external functions ****
      integer  ion_nion,brillioun_nbrillioun
      external ion_nion,brillioun_nbrillioun

      nion       = ion_nion()
      nbrillioun = brillioun_nbrillioun()
      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
       

      value = MA_alloc_get(mt_dcpl,(nx*nion),'ex1',ex1(2),ex1(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,(ny*nion),'ex2',ex2(2),ex2(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,(nz*nion),'ex3',ex3(2),ex3(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,(nbrillioun*nion),
     >                    'cxr',cxr(2),cxr(1))
      if (.not. value) call errquit('out of heap memory',0)

      return
      end

C     ***********************************
C     *					*
C     *		cstrfac_end 		*
C     *					*
C     ***********************************

      subroutine cstrfac_end()
      implicit none

#include "mafdecls.fh"
#include "cstrfac_common.fh"

      logical value

      value = MA_free_heap(ex1(2))
      value = MA_free_heap(ex2(2))
      value = MA_free_heap(ex3(2))
      value = MA_free_heap(cxr(2))
       
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

#include	"mafdecls.fh"
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






C      **********************************
C      *				*
C      *	CSTRFAC 			*
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

#include "mafdecls.fh"
#include "cstrfac_common.fh"
   
*     **** local variables ****
      integer i,j,k,index
      integer nx,ny,nz
      integer nxh
      integer p,q,taskid

      call nwpw_timing_start(8)

      call Parallel_taskid(taskid)
      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      nxh=nx/2

*::::::::::::::::  structure factor ::::::::::::::::::::::::::::::::
      do k=0,nz-1
         call C3dB_ktoqp(1,(k+1),q,p)
         if (p.eq.taskid) then
            do j=0,ny-1
            do i=0,nx-1
                index = (q-1)*(nx)*ny
     >                + j*(nx)
     >                + i+1
                exi(index) = dcpl_mb(ex1(1)+i+(II-1)*nx)
     >                     * dcpl_mb(ex2(1)+j+(II-1)*ny)
     >                     * dcpl_mb(ex3(1)+k+(II-1)*nz)
            end do
            end do
         end if
      end do

      call nwpw_timing_end(8)



      return
      end

C ... End STRFAC ....


*      **********************************
*      *								*
*      *			CSTRFAC_k 			*
*      *								*
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

#include "mafdecls.fh"
#include "cstrfac_common.fh"

      integer  nion
      integer  ion_nion
      external ion_nion

      nion = ion_nion()
      ck = dcpl_mb(cxr(1)+ii-1+(nb-1)*nion)
      return
      end

*      **********************************
*      *								*
*      *			cphafac_k 			*
*      *								*
*      **********************************
      subroutine cphafac_k()
      implicit none

#include "mafdecls.fh"
#include "cstrfac_common.fh"

*     *** local variables
      integer    i,nb,nion

*     **** external functions ****
      integer  brillioun_nbrillioun,ion_nion
      real*8   brillioun_k,ion_rion
      external brillioun_nbrillioun,ion_nion
      external brillioun_k,ion_rion

      nion = ion_nion()
      do nb=1,brillioun_nbrillioun()
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

#include "mafdecls.fh"
#include "cstrfac_common.fh"

*     *** local variables
      integer    i,nb,nion

*     **** external functions ****
      integer  brillioun_nbrillioun,ion_nion
      external brillioun_nbrillioun,ion_nion

      write(*,*)
      write(*,*) "cstrfac_k_print:"
      nion = ion_nion()
      do nb=1,brillioun_nbrillioun()
      do i =1,nion
        write(*,*) "cstrfac_k nb,ii,cxr:",nb,i,
     >             dcpl_mb(cxr(1)+i-1+(nb-1)*nion)
      end do
      end do
      write(*,*)
      write(*,*)

      return
      end

