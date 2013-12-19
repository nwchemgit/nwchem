*     **************************************************
*     *                                                *
*     *            nwpw_meta_gga_bd                    *
*     *                                                *
*     **************************************************
      block data nwpw_meta_gga_bd
      implicit none
#include "nwpw_meta_gga.fh"
      data meta_on/.false./
      end

*     **************************************************
*     *                                                *
*     *            nwpw_meta_gga_init                  *
*     *                                                *
*     **************************************************
      subroutine nwpw_meta_gga_init()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_meta_gga.fh"
#include "nwxc.fh"

*     **** local variables ****
      logical ok
      integer n2ft3d

      meta_on = nwxc_is_mgga()

      if (meta_on) then
         call D3dB_n2ft3d(1,n2ft3d)
         ok =        MA_alloc_get(mt_dbl, 2*n2ft3d,
     >                     'dfdtau',dfdtau(2),dfdtau(1))
         ok = ok.and.MA_alloc_get(mt_dbl, 2*n2ft3d,
     >                     'tau',tau(2),tau(1))
         if (.not.ok)
     >      call errquit('nwpw_meta_gga_init:out of heap',0,MA_ERR)
         call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(dfdtau(1)),1)
         call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(tau(1)),1)
      end if
      return
      end

*     **************************************************
*     *                                                *
*     *            nwpw_meta_gga_end                   *
*     *                                                *
*     **************************************************
      subroutine nwpw_meta_gga_end()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_meta_gga.fh"

      logical ok

      if (meta_on) then
         ok =        MA_free_heap(dfdtau(2))
         ok = ok.and.MA_free_heap(tau(2))
         if (.not.ok)
     >      call errquit('nwpw_meta_gga_end:freeing heap',0,MA_ERR)
      end if
      meta_on = .false.
      return
      end

*     **************************************************
*     *                                                *
*     *            nwpw_meta_gga_on                    *
*     *                                                *
*     **************************************************
      logical function nwpw_meta_gga_on()
      implicit none

#include "nwpw_meta_gga.fh"

      nwpw_meta_gga_on = meta_on
      return
      end



*     **************************************************
*     *                                                *
*     *            nwpw_meta_gga_gen_tau               *
*     *                                                *
*     **************************************************

*   This routine computes
*   tau = |grad psi|
*
      subroutine nwpw_meta_gga_gen_tau(ispin,neq,psi)
      implicit none
      integer ispin,neq(2)
      complex*16 psi(*)


#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_meta_gga.fh"

*     **** local variables ****
      integer ms,n,n2ft3d,nfft3d,npack1,xyz,nx,ny,nz
      integer G(3),n1(2),n2(2),dpsi(2)
      real*8 scal2

*     **** external functions ****
      integer  Pack_G_indx
      external Pack_G_indx
      real*8   lattice_omega
      external lattice_omega


      if (meta_on) then 
      n1(1) = 1
      n2(1) = neq(1)
      n1(2) = neq(1)+1
      n2(2) = neq(1)+neq(2)

      scal2 = 0.5d0/lattice_omega()

      call D3dB_nfft3d(1,nfft3d)
      call D3dB_n2ft3d(1,n2ft3d)
      call Pack_npack(1,npack1)

*     **** allocate stack memory ****
      if (.not.MA_push_get(mt_dcpl, nfft3d,'dpsi',dpsi(2),dpsi(1)))
     >   call errquit('nwpw_meta_gga_Hpsik:out of stack',0,MA_ERR)
      call dcopy(2*nfft3d,0.0d0,0,dcpl_mb(dpsi(1)),1)

      G(1)  = Pack_G_indx(1,1)
      G(2)  = Pack_G_indx(1,2)
      G(3)  = Pack_G_indx(1,3)

      call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(tau(1)),1)
      do ms=1,ispin
         do n=n1(ms),n2(ms)
            do xyz=1,3

*              **** compute dpsi(G,1:3) ****
               call Pack_tc_iMul(1,dbl_mb(G(xyz)),
     >                             psi(1+(n-1)*npack1),
     >                             dcpl_mb(dpsi(1)))

*              **** compute dpsi(r,1:3) ****
               call Pack_c_unpack(1,dcpl_mb(dpsi(1)))
               call D3dB_cr_fft3b(1,dcpl_mb(dpsi(1)))

*              **** tau(r,ms) += dpsi(r,ms)**2 ****
               call D3dB_rr_SqrAdd(1,dcpl_mb(dpsi(1)),
     >                               dbl_mb(tau(1)+(ms-1)*n2ft3d))
            end do
         end do
         call D3dB_r_SMul1(1,(scal2),dbl_mb(tau(1)+(ms-1)*n2ft3d))
         call D3dB_r_Zero_Ends(1,dbl_mb(tau(1)+(ms-1)*n2ft3d))
      end do

*     **** deallocate stack ****
      if (.not.MA_pop_stack(dpsi(2)))
     >   call errquit('nwpw_meta_gga_Hpsik:popping stack',0,MA_ERR)


      end if
      return
      end


*     **************************************************
*     *                                                *
*     *            nwpw_meta_gga_Hpsik                 *
*     *                                                *
*     **************************************************

*   This routine computes
*   Hpsik = del dot (dfdtau * del psi)
*
      subroutine nwpw_meta_gga_Hpsik(ispin,neq,psi,Hpsi)
      implicit none
      integer ispin,neq(2)
      complex*16 psi(*)
      complex*16 Hpsi(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_meta_gga.fh"

*     **** local variables ****
      integer ms,n,n2ft3d,nfft3d,npack1,nx,ny,nz,xyz
      integer G(3),n1(2),n2(2),dpsi(2)
      real*8 scal1,scal2

*     **** external functions ****
      integer  Pack_G_indx
      external Pack_G_indx
      real*8   lattice_omega
      external lattice_omega

      if (meta_on) then
      n1(1) = 1
      n2(1) = neq(1)
      n1(2) = neq(1)+1
      n2(2) = neq(1)+neq(2)

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 0.5d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()


      call D3dB_nfft3d(1,nfft3d)
      call D3dB_n2ft3d(1,n2ft3d)
      call Pack_npack(1,npack1)

*     **** allocate stack memory ****
      if (.not.MA_push_get(mt_dcpl, nfft3d,'dpsi',dpsi(2),dpsi(1)))
     >   call errquit('nwpw_meta_gga_Hpsik:out of stack',0,MA_ERR)
      call dcopy(2*nfft3d,0.0d0,0,dcpl_mb(dpsi(1)),1)

      G(1)  = Pack_G_indx(1,1)
      G(2)  = Pack_G_indx(1,2)
      G(3)  = Pack_G_indx(1,3)

      do ms=1,ispin
      do n=n1(ms),n2(ms)
         do xyz=1,3

*           **** compute dpsi(G,1:3) ****
            call Pack_tc_iMul(1,dbl_mb(G(xyz)),
     >                          psi(1+(n-1)*npack1),
     >                          dcpl_mb(dpsi(1)))

*           **** compute dpsi(r,1:3) ****
            !call Pack_c_SMul1(1,scal2,dcpl_mb(dpsi(1)))
            call Pack_c_unpack(1,dcpl_mb(dpsi(1)))
            call D3dB_cr_fft3b(1,dcpl_mb(dpsi(1)))

*           **** multiply dpsi(r,1:3) = dfdtau(r,ms)*dpsi(r,1:3) ****
            call D3dB_rr_Mul2(1,dbl_mb(dfdtau(1)+(ms-1)*n2ft3d),
     >                          dcpl_mb(dpsi(1)))

*           **** compute dpsi(G,1:3) ****
            call D3dB_r_Zero_Ends(1,dcpl_mb(dpsi(1)))
            call D3dB_rc_fft3f(1,dcpl_mb(dpsi(1)))
            call Pack_c_pack(1,dcpl_mb(dpsi(1)))

*           **** compute Hpsi = grad dot dpsi(G) ****
            call Pack_c_SMul1(1,(scal1),dcpl_mb(dpsi(1)))
            call Pack_tc_iMul2(1,dbl_mb(G(xyz)),
     >                          dcpl_mb(dpsi(1)))
            call Pack_cc_Sum2(1,dcpl_mb(dpsi(1)),Hpsi(1+(n-1)*npack1))

         end do
      end do
      end do

*     **** deallocate stack ****
      if (.not.MA_pop_stack(dpsi(2)))
     >   call errquit('nwpw_meta_gga_Hpsik:popping stack',0,MA_ERR)

      end if
      return
      end


*     **************************************************
*     *                                                *
*     *            nwpw_meta_gga_pxc                   *
*     *                                                *
*     **************************************************
      real*8 function nwpw_meta_gga_pxc(ispin,neq,psi)
      implicit none
      integer ispin,neq(2)
      complex*16 psi(*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_meta_gga.fh"

*     **** local variables ****
      integer ms,n,n2ft3d,nfft3d,npack1,nx,ny,nz,xyz
      integer G(3),n1(2),n2(2),dpsi(2)
      real*8 scal1,scal2,ee
      real*8 pmeta

*     **** external functions ****
      integer  Pack_G_indx
      external Pack_G_indx
      real*8   lattice_omega
      external lattice_omega

      pmeta = 0.0d0

      if (meta_on) then
      n1(1) = 1
      n2(1) = neq(1)
      n1(2) = neq(1)+1
      n2(2) = neq(1)+neq(2)

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 0.5d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()


      call D3dB_nfft3d(1,nfft3d)
      call D3dB_n2ft3d(1,n2ft3d)
      call Pack_npack(1,npack1)

*     **** allocate stack memory ****
      if (.not.MA_push_get(mt_dcpl, nfft3d,'dpsi',dpsi(2),dpsi(1)))
     >   call errquit('nwpw_meta_gga_Hpsik:out of stack',0,MA_ERR)
      call dcopy(2*nfft3d,0.0d0,0,dcpl_mb(dpsi(1)),1)

      G(1)  = Pack_G_indx(1,1)
      G(2)  = Pack_G_indx(1,2)
      G(3)  = Pack_G_indx(1,3)

      do ms=1,ispin
      do n=n1(ms),n2(ms)
         do xyz=1,3

*           **** compute dpsi(G,1:3) ****
            call Pack_tc_iMul(1,dbl_mb(G(xyz)),
     >                          psi(1+(n-1)*npack1),
     >                          dcpl_mb(dpsi(1)))

*           **** compute dpsi(r,1:3) ****
            !call Pack_c_SMul1(1,scal2,dcpl_mb(dpsi(1)))
            call Pack_c_unpack(1,dcpl_mb(dpsi(1)))
            call D3dB_cr_fft3b(1,dcpl_mb(dpsi(1)))

*           **** multiply dpsi(r,1:3) = dfdtau(r,ms)*dpsi(r,1:3) ****
            call D3dB_rr_Mul2(1,dbl_mb(dfdtau(1)+(ms-1)*n2ft3d),
     >                          dcpl_mb(dpsi(1)))

*           **** compute dpsi(G,1:3) ****
            call D3dB_r_Zero_Ends(1,dcpl_mb(dpsi(1)))
            call D3dB_rc_fft3f(1,dcpl_mb(dpsi(1)))
            call Pack_c_pack(1,dcpl_mb(dpsi(1)))

*           **** compute Hpsi = grad dot dpsi(G) ****
            call Pack_c_SMul1(1,(scal1),dcpl_mb(dpsi(1)))
            call Pack_tc_iMul2(1,dbl_mb(G(xyz)),
     >                          dcpl_mb(dpsi(1)))
            call Pack_cc_idot(1,psi(1+(n-1)*npack1),dcpl_mb(dpsi(1)),ee)
            pmeta = pmeta - ee

         end do
      end do
      end do

*     **** deallocate stack ****
      if (.not.MA_pop_stack(dpsi(2)))
     >   call errquit('nwpw_meta_gga_Hpsik:popping stack',0,MA_ERR)

      call Parallel_SumAll(pmeta)
      if (ispin.eq.1) pmeta = pmeta + pmeta

      end if

      nwpw_meta_gga_pxc = pmeta
      return
      end
c $Id$
