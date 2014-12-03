*
* $Id$
*

*     *******************************
*     *			  	    *
*     *	        cloak_init          *
*     *				    *
*     *******************************
      subroutine cloak_init()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "cloak_common.fh"

      logical value
      integer nfft3d

*     **** allocate masker memory ****
      call C3dB_nfft3d(1,nfft3d)

      value = BA_alloc_get(mt_log,nfft3d,
     >                     'masker',masker(2),masker(1))
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)

      return
      end

*     *******************************
*     *				    *
*     *	        cloak_end           *
*     *				    *
*     *******************************
      subroutine cloak_end()
      implicit none
#include "errquit.fh"

#include "bafdecls.fh"
#include "cloak_common.fh"

      logical value

      value = BA_free_heap(masker(2))
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)

      return
      end

*     *******************************
*     *				    *
*     *	        cloak_set           *
*     *				    *
*     *******************************
      subroutine cloak_set(kvector,ecut)
      implicit none
      real*8 kvector(3),ecut

#include "bafdecls.fh"
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
      real*8   lattice_unitg_frozen
      external lattice_unitg_frozen

      call Parallel3d_taskid_i(taskid)
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
               g1 = k1*lattice_unitg_frozen(1,1) 
     >            + k2*lattice_unitg_frozen(1,2) 
     >            + k3*lattice_unitg_frozen(1,3)
     >            + kvector(1)
               g2 = k1*lattice_unitg_frozen(2,1) 
     >            + k2*lattice_unitg_frozen(2,2) 
     >            + k3*lattice_unitg_frozen(2,3)
     >            + kvector(2)
               g3 = k1*lattice_unitg_frozen(3,1) 
     >            + k2*lattice_unitg_frozen(3,2) 
     >            + k3*lattice_unitg_frozen(3,3)
     >            + kvector(3)
               i=k1
               j=k2
               k=k3
               if (i .lt. 0) i = i + nx
               if (j .lt. 0) j = j + ny
               if (k .lt. 0) k = k + nz

                !call C3dB_ktoqp(1,k+1,q,p)
                call C3dB_ijktoindexp(1,i+1,j+1,k+1,index,p)
                if (p .eq. taskid) then
                  gg = g1*g1 + g2*g2 + g3*g3
                  if (gg.lt.ggcut) then
c                    index = (q-1)*(nx)*ny
c    >                     + j*(nx)
c    >                     + i+1
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
*     *		 		    *
*     *	        cloak_C             *
*     *				    *
*     *******************************
      subroutine cloak_C(A)
      implicit none
      complex*16 A(*)

#include "bafdecls.fh"
#include "cloak_common.fh"

*     **** local variables ****
      integer i,nfft3d

      call nwpw_timing_start(9)

      call C3dB_nfft3d(1,nfft3d)
      do i=1,nfft3d
         if (log_mb(masker(1)+i-1)) A(i) = dcmplx(0.0d0,0.0d0)
      end do

      call nwpw_timing_end(9)
      return 
      end

*     *******************************
*     *				    *
*     *	        cloak_R             *
*     *				    *
*     *******************************
      subroutine cloak_R(A)
      implicit none
      real*8    A(*)

#include "bafdecls.fh"
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
*     *			  	    *
*     *	        cloak_nwave         *
*     *				    *
*     *******************************
      integer function cloak_nwave()
      implicit none

#include "cloak_common.fh"

      cloak_nwave = nwave
      return
      end

*     *******************************
*     *				    *
*     *	        cloak_nwave_all     *
*     *				    *
*     *******************************
      integer function cloak_nwave_all()
      implicit none

#include "cloak_common.fh"

      cloak_nwave_all = nwave_all
      return
      end

*     *******************************
*     *				    *
*     *	        cloak_masker        *
*     *				    *
*     *******************************
      logical function cloak_masker(i)
      implicit none
      integer i

#include "bafdecls.fh"
#include "cloak_common.fh"

      cloak_masker = log_mb(masker(1)+i-1) 
      return
      end

