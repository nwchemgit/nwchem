*
* $Id: mask_init.F,v 1.2 2001-12-12 02:42:27 bylaska Exp $
*

      subroutine mask_end()
      implicit none

#include "mafdecls.fh"

*     **** common block for masker ****
      integer masker(2,0:1),nwave(0:1),nwave_all(0:1)
      common / mask_block / masker,nwave,nwave_all

      logical value

      value = MA_free_heap(masker(2,0))
      value = value.and.
     >        MA_free_heap(masker(2,1))
      if (.not. value) call errquit('out of heap memory',0)

      return
      end



      subroutine mask_init()
      implicit none

#include "mafdecls.fh"

*     **** common block for masker ****
*     logical masker(nfft3d,0:1)
      integer masker(2,0:1),nwave(0:1),nwave_all(0:1)
      common / mask_block / masker,nwave,nwave_all

*     **** local variables ****
      integer nfft3d
      integer nb
      integer i,j,k,p,q,index
      integer k1,k2,k3
      integer nx,ny,nz
      integer nxh,nyh,nzh
      real*8  ggcut,g1,g2,g3,gg
      integer taskid
      logical value

*     **** external functions ***
      real*8   lattice_ggcut,lattice_wggcut,lattice_unitg
      external lattice_ggcut,lattice_wggcut,lattice_unitg

*     **** allocate masker memory ****
      call D3dB_nfft3d(1,nfft3d)

      value = MA_alloc_get(mt_log,nfft3d,
     >                     'masker',masker(2,0),masker(1,0))
      value = value.and.
     >        MA_alloc_get(mt_log,nfft3d,
     >                     'masker',masker(2,1),masker(1,1))
      if (.not. value) call errquit('out of heap memory',0)


      call Parallel_taskid(taskid)
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      nxh = nx/2
      nyh = ny/2
      nzh = nz/2

*     **** set all masker on ****
      do nb=0,1
      do i=1,nfft3d
            log_mb(masker(1,nb)+i-1) = .true.
      end do
      end do

*     **** get fermi sphere cut-off ****
      do nb=0,1
         nwave(nb) = 0
         if (nb.eq.0) then
            ggcut = lattice_ggcut()
         else
            ggcut = lattice_wggcut()
         end if


*        **** undo masker in sphere defined by ggcut ****
         do k3 = -nzh+1, nzh-1
            do k2 = -nyh+1, nyh-1
               do k1 = 0,nxh-1
                  g1 = k1*lattice_unitg(1,1) 
     >               + k2*lattice_unitg(1,2) 
     >               + k3*lattice_unitg(1,3)
                  g2 = k1*lattice_unitg(2,1) 
     >               + k2*lattice_unitg(2,2) 
     >               + k3*lattice_unitg(2,3)
                  g3 = k1*lattice_unitg(3,1) 
     >               + k2*lattice_unitg(3,2) 
     >               + k3*lattice_unitg(3,3)
                  i=k1
                  j=k2
                  k=k3
                  if (i .lt. 0) i = i + nx
                  if (j .lt. 0) j = j + ny
                  if (k .lt. 0) k = k + nz

                  call D3dB_ktoqp(1,k+1,q,p)

                  if (p .eq. taskid) then
                     gg = g1*g1 + g2*g2 + g3*g3
                     if (gg.lt.ggcut) then
                        index = (q-1)*(nxh+1)*ny
     >                        + j*(nxh+1)
     >                        + i+1
                        log_mb(masker(1,nb)+index-1) = .false.
                        nwave(nb) = nwave(nb) + 1
                     end if
                  end if
               end do
            end do
         end do
         nwave_all(nb) = nwave(nb)
         call D3dB_ISumAll(nwave_all(nb))
      end do

      return
      end

      subroutine mask_C(nb,A)
      implicit none
      integer    nb
      complex*16 A(*)

#include "mafdecls.fh"

*     **** common block for masker ****
      integer masker(2,0:1),nwave(0:1),nwave_all(0:1)
      common / mask_block / masker,nwave,nwave_all

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


      call D3dB_nfft3d(1,nfft3d)
      do i=1,nfft3d
         if (log_mb(masker(1,nb)+i-1)) A(i) = dcmplx(0.0d0,0.0d0)
      end do

      call current_second(tim2)
      times(9) = times(9) + (tim2-tim1)

      return 
      end





      subroutine mask_R(nb,A)
      implicit none
      integer   nb
      real*8    A(*)

#include "mafdecls.fh"

*     **** common block for masker ****
      integer masker(2,0:1),nwave(0:1),nwave_all(0:1)
      common / mask_block / masker,nwave,nwave_all

*     **** local variables ****
      integer i,nfft3d

      call D3dB_nfft3d(1,nfft3d)
      do i=1,nfft3d
         if (log_mb(masker(1,nb)+i-1)) A(i) = 0.0d0
      end do

      return 
      end


      integer function mask_nwave(nb)
      implicit none
      integer nb

*     **** common block for masker ****
      integer masker(2,0:1),nwave(0:1),nwave_all(0:1)
      common / mask_block / masker,nwave,nwave_all

      mask_nwave = nwave(nb)
      return
      end

      integer function mask_nwave_all(nb)
      implicit none
      integer nb

*     **** common block for masker ****
      integer masker(2,0:1),nwave(0:1),nwave_all(0:1)
      common / mask_block / masker,nwave,nwave_all

      mask_nwave_all = nwave_all(nb)
      return
      end
