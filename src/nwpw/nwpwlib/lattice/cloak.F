*
* $Id: cloak.F,v 1.3 2001-12-08 01:06:38 bylaska Exp $
*

*     *******************************
*     *								*
*     *	        cloak_init          *
*     *								*
*     *******************************
      subroutine cloak_init()
      implicit none

#include "mafdecls.fh"
#include "cloak_common.fh"

      logical value
      integer nfft3d

*     **** allocate masker memory ****
      call C3dB_nfft3d(1,nfft3d)

      value = MA_alloc_get(mt_log,nfft3d,
     >                     'masker',masker(2),masker(1))
      if (.not. value) call errquit('out of heap memory',0)

      return
      end

*     *******************************
*     *								*
*     *	        cloak_end           *
*     *								*
*     *******************************
      subroutine cloak_end()
      implicit none

#include "mafdecls.fh"
#include "cloak_common.fh"

      logical value

      value = MA_free_heap(masker(2))
      if (.not. value) call errquit('out of heap memory',0)

      return
      end

*     *******************************
*     *								*
*     *	        cloak_set           *
*     *								*
*     *******************************
      subroutine cloak_set(kvector,ecut)
      implicit none
      real*8 kvector(3),ecut

#include "mafdecls.fh"
#include "cloak_common.fh"

*     **** local variables ****
      integer nfft3d
      integer i,j,k,p,q,index
      integer k1,k2,k3
      integer nx,ny,nz
      integer nxh,nyh,nzh
      real*8  ggcut,g1,g2,g3,gg
      integer taskid

*     **** external functions ***
      real*8   lattice_unitg
      external lattice_unitg

      call Parallel_taskid(taskid)
      call C3dB_nfft3d(1,nfft3d)
      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      nxh = nx/2
      nyh = ny/2
      nzh = nz/2

*     **** set all masker on ****
      do i=1,nfft3d
            log_mb(masker(1)+i-1) = .true.
      end do

*     **** get fermi sphere cut-off ****
      nwave = 0
      ggcut = 2.0d0*ecut


*     **** undo masker in sphere defined by ggcut ****
      do k3 = -nzh+1, nzh-1
         do k2 = -nyh+1, nyh-1
            do k1 = -nxh+1, nxh-1
               g1 = k1*lattice_unitg(1,1) 
     >            + k2*lattice_unitg(1,2) 
     >            + k3*lattice_unitg(1,3)
     >            + kvector(1)
               g2 = k1*lattice_unitg(2,1) 
     >            + k2*lattice_unitg(2,2) 
     >            + k3*lattice_unitg(2,3)
     >            + kvector(2)
               g3 = k1*lattice_unitg(3,1) 
     >            + k2*lattice_unitg(3,2) 
     >            + k3*lattice_unitg(3,3)
     >            + kvector(3)
               i=k1
               j=k2
               k=k3
               if (i .lt. 0) i = i + nx
               if (j .lt. 0) j = j + ny
               if (k .lt. 0) k = k + nz
                call C3dB_ktoqp(1,k+1,q,p)
                if (p .eq. taskid) then
                  gg = g1*g1 + g2*g2 + g3*g3
                  if (gg.lt.ggcut) then
                     index = (q-1)*(nx)*ny
     >                     + j*(nx)
     >                     + i+1
                     log_mb(masker(1)+index-1) = .false.
                     nwave = nwave + 1
                  end if
               end if
            end do
         end do
      end do
      nwave_all = nwave
      call C3dB_ISumAll(nwave_all)

      return
      end

*     *******************************
*     *								*
*     *	        cloak_C             *
*     *								*
*     *******************************
      subroutine cloak_C(A)
      implicit none
      complex*16 A(*)

#include "mafdecls.fh"
#include "cloak_common.fh"

*     **** local variables ****
      integer i,nfft3d

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


      call C3dB_nfft3d(1,nfft3d)
      do i=1,nfft3d
         if (log_mb(masker(1)+i-1)) A(i) = dcmplx(0.0d0,0.0d0)
      end do

      call current_second(tim2)
      times(9) = times(9) + (tim2-tim1)

      return 
      end

*     *******************************
*     *								*
*     *	        cloak_R             *
*     *								*
*     *******************************
      subroutine cloak_R(A)
      implicit none
      real*8    A(*)

#include "mafdecls.fh"
#include "cloak_common.fh"

*     **** local variables ****
      integer i,nfft3d

      call C3dB_nfft3d(1,nfft3d)
      do i=1,nfft3d
         if (log_mb(masker(1)+i-1)) A(i) = 0.0d0
      end do

      return 
      end


*     *******************************
*     *								*
*     *	        cloak_nwave         *
*     *								*
*     *******************************
      integer function cloak_nwave()
      implicit none

#include "cloak_common.fh"

      cloak_nwave = nwave
      return
      end

*     *******************************
*     *								*
*     *	        cloak_nwave_all     *
*     *								*
*     *******************************
      integer function cloak_nwave_all()
      implicit none

#include "cloak_common.fh"

      cloak_nwave_all = nwave_all
      return
      end
