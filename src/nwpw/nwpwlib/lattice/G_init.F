c
c $Id$
c

      integer function G_indx(i)
      implicit none 
      integer i

#include "bafdecls.fh"

*     **** G_block common block ****
      integer G(3),g_hndl
      common / G_block / G,g_hndl

      G_indx = G(i)
      return
      end


      integer function Gsmall_indx(i)
      implicit none
      integer i

#include "bafdecls.fh"

*     **** G_small_block common block ****
      logical has_small
      integer Gsmall(3),gsmall_hndl
      common / G_small_block / Gsmall,gsmall_hndl,has_small

      Gsmall_indx = Gsmall(i)
      return
      end


      subroutine G_init()
      implicit none 


#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      integer taskid,p
      integer i,j,k
      integer k1,k2,k3,index
      integer nx,ny,nz,nfft3d
      integer nxh,nyh,nzh
      real*8  g1,g2,g3,dk1,dk2,dk3
      logical value

*     **** G_block common block ****
      integer G(3),g_hndl
      common / G_block / G,g_hndl

*     **** G_small_block common block ****
      logical has_small
      integer Gsmall(3),gsmall_hndl
      common / G_small_block / Gsmall,gsmall_hndl,has_small


*     **** external functions ***
      logical  lattice_has_small
      real*8   lattice_unitg,lattice_unitg_small
      external lattice_has_small
      external lattice_unitg,lattice_unitg_small


*     **************************
*     **** gamma point code ****
*     **************************

*       **** allocate G_block common block ****
        call D3dB_nfft3d(1,nfft3d)
        value = BA_alloc_get(mt_dbl,(3*nfft3d),'G',g_hndl,G(1))
        if (.not. value) call errquit('out of heap memory',0, MA_ERR)

        G(2) = G(1) + nfft3d
        G(3) = G(2) + nfft3d

        call Parallel2d_taskid_i(taskid)
        call D3dB_nx(1,nx)
        call D3dB_ny(1,ny)
        call D3dB_nz(1,nz)
        nxh = nx/2
        nyh = ny/2
        nzh = nz/2


*       **** inefficient but less confusing algorithm ****
        call dcopy(3*nfft3d,0.0d0,0,dbl_mb(G(1)),1)
        do k3 = -nzh+1, nzh
           do k2 = -nyh+1, nyh
              do k1 = 0,nxh
               dk1 = dble(k1)
               dk2 = dble(k2)
               dk3 = dble(k3)
               g1 = dk1*lattice_unitg(1,1) 
     >            + dk2*lattice_unitg(1,2) 
     >            + dk3*lattice_unitg(1,3)
               g2 = dk1*lattice_unitg(2,1) 
     >            + dk2*lattice_unitg(2,2) 
     >            + dk3*lattice_unitg(2,3)
               g3 = dk1*lattice_unitg(3,1) 
     >            + dk2*lattice_unitg(3,2) 
     >            + dk3*lattice_unitg(3,3)
               i=k1
               j=k2
               k=k3
               if (i .lt. 0) i = i + nx
               if (j .lt. 0) j = j + ny
               if (k .lt. 0) k = k + nz

               !call D3dB_ktoqp(1,k+1,q,p)
               call D3dB_ijktoindexp(1,i+1,j+1,k+1,index,p)
               if (p .eq. taskid) then
c                 index = (q-1)*(nxh+1)*ny
c    >                  + j*(nxh+1) 
c    >                  + i+1
                  dbl_mb(G(1)+index-1) = g1
                  dbl_mb(G(2)+index-1) = g2
                  dbl_mb(G(3)+index-1) = g3
               end if

              end do  
           end do
        end do


      has_small = lattice_has_small()
      if (has_small) then
*       **** allocate G_block common block ****
        call D3dB_nfft3d(3,nfft3d)
        value = BA_alloc_get(mt_dbl,(3*nfft3d),'Gsmall',
     >                               gsmall_hndl,Gsmall(1))
        if (.not. value) call errquit('out of heap memory',0, MA_ERR)

        Gsmall(2) = Gsmall(1) + nfft3d
        Gsmall(3) = Gsmall(2) + nfft3d

        call Parallel2d_taskid_i(taskid)
        call D3dB_nx(3,nx)
        call D3dB_ny(3,ny)
        call D3dB_nz(3,nz)
        nxh = nx/2
        nyh = ny/2
        nzh = nz/2

*       **** inefficient but less confusing algorithm ****
        call dcopy(3*nfft3d,0.0d0,0,dbl_mb(Gsmall(1)),1)
        do k3 = -nzh+1, nzh
           do k2 = -nyh+1, nyh
              do k1 = 0,nxh
               dk1 = dble(k1)
               dk2 = dble(k2)
               dk3 = dble(k3)
               g1 = dk1*lattice_unitg_small(1,1)
     >            + dk2*lattice_unitg_small(1,2)
     >            + dk3*lattice_unitg_small(1,3)
               g2 = dk1*lattice_unitg_small(2,1)
     >            + dk2*lattice_unitg_small(2,2)
     >            + dk3*lattice_unitg_small(2,3)
               g3 = dk1*lattice_unitg_small(3,1)
     >            + dk2*lattice_unitg_small(3,2)
     >            + dk3*lattice_unitg_small(3,3)
               i=k1
               j=k2
               k=k3
               if (i .lt. 0) i = i + nx
               if (j .lt. 0) j = j + ny
               if (k .lt. 0) k = k + nz

               !call D3dB_ktoqp(1,k+1,q,p)
               call D3dB_ijktoindexp(3,i+1,j+1,k+1,index,p)
               if (p .eq. taskid) then
c                 index = (q-1)*(nxh+1)*ny
c    >                  + j*(nxh+1) 
c    >                  + i+1
                  dbl_mb(Gsmall(1)+index-1) = g1
                  dbl_mb(Gsmall(2)+index-1) = g2
                  dbl_mb(Gsmall(3)+index-1) = g3
               end if

              end do
           end do
        end do

      end if


      return
      end

      subroutine G_end()
      implicit none 
#include "bafdecls.fh"
#include "errquit.fh"

*     **** G_block common block ****
      integer G(3),g_hndl
      common / G_block / G,g_hndl

*     **** G_small_block common block ****
      logical has_small
      integer Gsmall(3),gsmall_hndl
      common / G_small_block / Gsmall,gsmall_hndl,has_small

      logical value

      value = BA_free_heap(g_hndl)
      if (has_small)  value = value.and.BA_free_heap(gsmall_hndl)
      if (.not.value) call errquit('G_end:freeing heap',0,MA_ERR)

      return
      end

