*
* $Id$
*

*     *********************************************
*     *                                           *
*     *                 mask_end                  *
*     *                                           *
*     *********************************************
      subroutine mask_end()
      implicit none

#include "bafdecls.fh"
#include "mask_common.fh"
#include "errquit.fh"

      logical value

      value =           BA_free_heap(masker(2,0))
      value = value.and.BA_free_heap(masker(2,1))
      if (has_small) then
         value = value.and.BA_free_heap(masker(2,2))
         value = value.and.BA_free_heap(masker(2,3))
      end if
      if (.not.value) call errquit('mask_end: freeing heap',0, MA_ERR)

      return
      end


*     *********************************************
*     *                                           *
*     *                 mask_init                 *
*     *                                           *
*     *********************************************
      subroutine mask_init()
      implicit none

#include "bafdecls.fh"
#include "mask_common.fh"
#include "errquit.fh"


*     **** local variables ****
      integer nfft3d
      integer nb
      integer i,j,k,p,q,index
      integer k1,k2,k3
      integer nx,ny,nz
      integer nxh,nyh,nzh
      real*8  ggcut,g1,g2,g3,gg
      double precision eps
      parameter (eps=1d-12)
      integer taskid
      logical value

*     **** external functions ***
      real*8   lattice_ggcut,lattice_wggcut,lattice_unitg_frozen
      external lattice_ggcut,lattice_wggcut,lattice_unitg_frozen
      logical  lattice_has_small
      external lattice_has_small
      real*8   lattice_unitg_frozen_small
      external lattice_unitg_frozen_small

*     **** allocate masker memory ****
      call D3dB_nfft3d(1,nfft3d)

      value = BA_alloc_get(mt_log,nfft3d,
     >                     'masker',masker(2,0),masker(1,0))
      value = value.and.
     >        BA_alloc_get(mt_log,nfft3d,
     >                     'masker',masker(2,1),masker(1,1))
      if (.not.value) call errquit('out of heap memory',0, MA_ERR)


      call Parallel2d_taskid_i(taskid)
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
                  g1 = k1*lattice_unitg_frozen(1,1) 
     >               + k2*lattice_unitg_frozen(1,2) 
     >               + k3*lattice_unitg_frozen(1,3)
                  g2 = k1*lattice_unitg_frozen(2,1) 
     >               + k2*lattice_unitg_frozen(2,2) 
     >               + k3*lattice_unitg_frozen(2,3)
                  g3 = k1*lattice_unitg_frozen(3,1) 
     >               + k2*lattice_unitg_frozen(3,2) 
     >               + k3*lattice_unitg_frozen(3,3)
                  i=k1
                  j=k2
                  k=k3
                  if (i .lt. 0) i = i + nx
                  if (j .lt. 0) j = j + ny
                  if (k .lt. 0) k = k + nz

c                 call D3dB_ktoqp(1,k+1,q,p)
                  call D3dB_ijktoindexp(1,i+1,j+1,k+1,index,p)
                  if (p .eq. taskid) then
                     gg = g1*g1 + g2*g2 + g3*g3
                     gg= gg-ggcut
                     if (gg.lt.-eps) then
c                       index = (q-1)*(nxh+1)*ny
c    >                        + j*(nxh+1)
c    >                        + i+1
                        log_mb(masker(1,nb)+index-1) = .false.
                        nwave(nb) = nwave(nb) + 1
                     end if
                  end if
               end do
            end do
         end do
         nwave_entire(nb) = nwave(nb)
         call D3dB_ISumAll(nwave_entire(nb))

      end do

*     ************************************
*     ***** initialize small lattice *****
*     ************************************
      has_small = lattice_has_small()
      if (has_small) then

*        **** allocate masker memory ****
         call D3dB_nfft3d(3,nfft3d)

         value = BA_alloc_get(mt_log,nfft3d,
     >                        'masker',masker(2,2),masker(1,2))
         value = value.and.
     >           BA_alloc_get(mt_log,nfft3d,
     >                        'masker',masker(2,3),masker(1,3))
         if (.not.value) call errquit('mask_init:out of heap',1, MA_ERR)


         call Parallel2d_taskid_i(taskid)
         call D3dB_nx(3,nx)
         call D3dB_ny(3,ny)
         call D3dB_nz(3,nz)
         nxh = nx/2
         nyh = ny/2
         nzh = nz/2

*        **** set all masker on ****
         do nb=2,3
         do i=1,nfft3d
               log_mb(masker(1,nb)+i-1) = .true.
         end do
         end do

*        **** get fermi sphere cut-off ****
         do nb=2,3
            nwave(nb) = 0
            if (nb.eq.2) then
               ggcut = lattice_ggcut()
            else
               ggcut = lattice_wggcut()
            end if

*           **** undo masker in sphere defined by ggcut ****
            do k3 = -nzh+1, nzh-1
               do k2 = -nyh+1, nyh-1
                  do k1 = 0,nxh-1
                     g1 = k1*lattice_unitg_frozen_small(1,1)
     >                  + k2*lattice_unitg_frozen_small(1,2)
     >                  + k3*lattice_unitg_frozen_small(1,3)
                     g2 = k1*lattice_unitg_frozen_small(2,1)
     >                  + k2*lattice_unitg_frozen_small(2,2)
     >                  + k3*lattice_unitg_frozen_small(2,3)
                     g3 = k1*lattice_unitg_frozen_small(3,1)
     >                  + k2*lattice_unitg_frozen_small(3,2)
     >                  + k3*lattice_unitg_frozen_small(3,3)
                     i=k1
                     j=k2
                     k=k3
                     if (i .lt. 0) i = i + nx
                     if (j .lt. 0) j = j + ny
                     if (k .lt. 0) k = k + nz

c                    call D3dB_ktoqp(3,k+1,q,p)
                     call D3dB_ijktoindexp(3,i+1,j+1,k+1,index,p)
                     if (p .eq. taskid) then
                        gg = g1*g1 + g2*g2 + g3*g3
                        gg= gg-ggcut
                        if (gg.lt.-eps) then
c                          index = (q-1)*(nxh+1)*ny
c    >                           + j*(nxh+1)
c    >                           + i+1
                           log_mb(masker(1,nb)+index-1) = .false.
                           nwave(nb) = nwave(nb) + 1
                        end if
                     end if
                  end do
               end do
            end do
            nwave_entire(nb) = nwave(nb)
            call D3dB_ISumAll(nwave_entire(nb))

         end do

      end if

      return
      end

      subroutine mask_C(nb,A)
      implicit none
      integer    nb
      complex*16 A(*)

#include "bafdecls.fh"
#include "mask_common.fh"

*     **** local variables ****
      integer i,nfft3d

      call nwpw_timing_start(9)

      if (nb.lt.2) then
         call D3dB_nfft3d(1,nfft3d)
      else
         call D3dB_nfft3d(3,nfft3d)
      end if
!$OMP DO
      do i=1,nfft3d
         if (log_mb(masker(1,nb)+i-1)) A(i) = dcmplx(0.0d0,0.0d0)
      end do
!$OMP END DO

      call nwpw_timing_end(9)

      return 
      end





      subroutine mask_R(nb,A)
      implicit none
      integer   nb
      real*8    A(*)

#include "bafdecls.fh"
#include "mask_common.fh"

*     **** local variables ****
      integer i,nfft3d

      if (nb.lt.2) then
         call D3dB_nfft3d(1,nfft3d)
      else
         call D3dB_nfft3d(3,nfft3d)
      end if
      do i=1,nfft3d
         if (log_mb(masker(1,nb)+i-1)) A(i) = 0.0d0
      end do

      return 
      end


      integer function mask_nwave(nb)
      implicit none
      integer nb

#include "mask_common.fh"

      mask_nwave = nwave(nb)
      return
      end

      integer function mask_nwave_all(nb)
      implicit none
      integer nb

#include "mask_common.fh"

      mask_nwave_all = nwave_entire(nb)
      return
      end
