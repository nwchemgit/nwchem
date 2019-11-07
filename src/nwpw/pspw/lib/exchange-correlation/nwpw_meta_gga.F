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
      subroutine nwpw_meta_gga_init(gga)
      implicit none
      integer gga

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_meta_gga.fh"
#include "nwpwxc.fh"
#include "util.fh"

*     **** local variables ****
      logical ok
      integer n2ft3d

      meta_on = .false.
      meta_on = nwpwxc_is_mgga().or.(gga.ge.300)

      if (meta_on) then
         call D3dB_n2ft3d(1,n2ft3d)
         ok =        BA_alloc_get(mt_dbl, 2*n2ft3d,
     >                     'dfdtau',dfdtau(2),dfdtau(1))
         ok = ok.and.BA_alloc_get(mt_dbl, 2*n2ft3d,
     >                     'tau',tau(2),tau(1))
         if (.not.ok)
     >      call errquit('nwpw_meta_gga_init:out of heap',0,MA_ERR)
         !call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(dfdtau(1)),1)
         !call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(tau(1)),1)
         call Parallel_shared_vector_zero(.false.,2*n2ft3d,
     >                                    dbl_mb(dfdtau(1)))
         call Parallel_shared_vector_zero(.true.,2*n2ft3d,
     >                                    dbl_mb(tau(1)))
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

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_meta_gga.fh"

      logical ok

      if (meta_on) then
         ok =        BA_free_heap(dfdtau(2))
         ok = ok.and.BA_free_heap(tau(2))
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
*     *            nwpw_meta_tau_ptr                   *
*     *                                                *
*     **************************************************
      integer function nwpw_meta_tau_ptr()
      implicit none
#include "nwpw_meta_gga.fh"
      nwpw_meta_tau_ptr = tau(1)
      return 
      end

*     **************************************************
*     *                                                *
*     *            nwpw_meta_dfdtau_ptr                *
*     *                                                *
*     **************************************************
      integer function nwpw_meta_dfdtau_ptr()
      implicit none
#include "nwpw_meta_gga.fh"
      nwpw_meta_dfdtau_ptr = dfdtau(1)
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


#include "bafdecls.fh"
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
      if (.not.BA_push_get(mt_dcpl, nfft3d,'dpsi',dpsi(2),dpsi(1)))
     >   call errquit('nwpw_meta_gga_Hpsik:out of stack',0,MA_ERR)
      !call dcopy(2*nfft3d,0.0d0,0,dcpl_mb(dpsi(1)),1)
      call Parallel_shared_vector_zero(.false.,2*nfft3d,
     >                                 dcpl_mb(dpsi(1)))

      G(1)  = Pack_G_indx(1,1)
      G(2)  = Pack_G_indx(1,2)
      G(3)  = Pack_G_indx(1,3)

      !call dcopy(2*n2ft3d,0.0d0,0,dbl_mb(tau(1)),1)
      call Parallel_shared_vector_zero(.true.,2*n2ft3d,dbl_mb(tau(1)))
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
      call D1dB_Vector_SumAll(ispin*n2ft3d,dbl_mb(tau(1)))

*     **** deallocate stack ****
      if (.not.BA_pop_stack(dpsi(2)))
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
*   Hpsik = 0.5 * del dot (dfdtau * del psi)
*
*           /
*   Emeta = | f(tau) dr
*           /
*
*  dEmeta/dpsi*(r) = dEmeta/dtau(r') * dtau(r')/dpsi*(r)
*
      subroutine nwpw_meta_gga_Hpsik(ispin,neq,psi,Hpsi)
      implicit none
      integer ispin,neq(2)
      complex*16 psi(*)
      complex*16 Hpsi(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_meta_gga.fh"

*     **** local variables ****
      integer ms,n,n2ft3d,nfft3d,npack1,nx,ny,nz,xyz
      integer G(3),n1(2),n2(2),dpsi(2)
      real*8 scal1,scal2,scal
      !real*8 scal1,scal2,scal,pmeta0,ee

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
      scal  = scal1


      call D3dB_nfft3d(1,nfft3d)
      call D3dB_n2ft3d(1,n2ft3d)
      call Pack_npack(1,npack1)

*     **** allocate stack memory ****
      if (.not.BA_push_get(mt_dcpl, nfft3d,'dpsi',dpsi(2),dpsi(1)))
     >   call errquit('nwpw_meta_gga_Hpsik:out of stack',0,MA_ERR)
      !call dcopy(2*nfft3d,0.0d0,0,dcpl_mb(dpsi(1)),1)
      call Parallel_shared_vector_zero(.true.,2*nfft3d,dcpl_mb(dpsi(1)))

      G(1)  = Pack_G_indx(1,1)
      G(2)  = Pack_G_indx(1,2)
      G(3)  = Pack_G_indx(1,3)

c      pmeta0 = 0.0d0
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

*           **** compute Hpsi = -0.5 * grad dot dpsi(G) ****
            call Pack_c_SMul1(1,(scal),dcpl_mb(dpsi(1)))
            call Pack_tc_iMul2(1,dbl_mb(G(xyz)),
     >                          dcpl_mb(dpsi(1)))
            call Pack_cc_Sum2(1,dcpl_mb(dpsi(1)),Hpsi(1+(n-1)*npack1))
c            call Pack_cc_dot(1,psi(1+(n-1)*npack1),
c     >                        dcpl_mb(dpsi(1)),ee)
c            pmeta0 = pmeta0 + ee

         end do
      end do
      end do
c      if (ispin.eq.1) pmeta0 = pmeta0 + pmeta0
      !write(*,*) "PMETA0=",pmeta0

*     **** deallocate stack ****
      if (.not.BA_pop_stack(dpsi(2)))
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

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_meta_gga.fh"

*     **** local variables ****
      integer ms,n2ft3d,nx,ny,nz
c      integer n,nfft3d,npack1,xyz
c      integer G(3),n1(2),n2(2),dpsi(2)
c      real*8 scal1,scal2
      real*8 pmeta,dV,ee,pmeta2

*     **** external functions ****
c      integer  Pack_G_indx
c      external Pack_G_indx
      real*8   lattice_omega
      external lattice_omega

      pmeta = 0.0d0

      if (meta_on) then
c         n1(1) = 1
c         n2(1) = neq(1)
c         n1(2) = neq(1)+1
c         n2(2) = neq(1)+neq(2)

         call D3dB_nx(1,nx)
         call D3dB_ny(1,ny)
         call D3dB_nz(1,nz)
c         scal1 = 0.5d0/dble(nx*ny*nz)
c         scal2 = 1.0d0/lattice_omega()


c         call D3dB_nfft3d(1,nfft3d)
         call D3dB_n2ft3d(1,n2ft3d)
c         call Pack_npack(1,npack1)

c*        **** allocate stack memory ****
c         if (.not.BA_push_get(mt_dcpl, nfft3d,'dpsi',dpsi(2),dpsi(1)))
c     >      call errquit('nwpw_meta_gga_Hpsik:out of stack',0,MA_ERR)
c         call dcopy(2*nfft3d,0.0d0,0,dcpl_mb(dpsi(1)),1)
c
c         G(1)  = Pack_G_indx(1,1)
c         G(2)  = Pack_G_indx(1,2)
c         G(3)  = Pack_G_indx(1,3)
c
c         do ms=1,ispin
c         do n=n1(ms),n2(ms)
c            do xyz=1,3
c
c*              **** compute dpsi(G,1:3) ****
c               call Pack_tc_iMul(1,dbl_mb(G(xyz)),
c     >                             psi(1+(n-1)*npack1),
c     >                             dcpl_mb(dpsi(1)))
c
c*              **** compute dpsi(r,1:3) ****
c               !call Pack_c_SMul1(1,scal2,dcpl_mb(dpsi(1)))
c               call Pack_c_unpack(1,dcpl_mb(dpsi(1)))
c               call D3dB_cr_fft3b(1,dcpl_mb(dpsi(1)))
c
c*              **** multiply dpsi(r,1:3) = dfdtau(r,ms)*dpsi(r,1:3) ****
c               call D3dB_rr_Mul2(1,dbl_mb(dfdtau(1)+(ms-1)*n2ft3d),
c     >                             dcpl_mb(dpsi(1)))
c
c*              **** compute dpsi(G,1:3) ****
c               call D3dB_r_Zero_Ends(1,dcpl_mb(dpsi(1)))
c               call D3dB_rc_fft3f(1,dcpl_mb(dpsi(1)))
c               call Pack_c_pack(1,dcpl_mb(dpsi(1)))
c
c*              **** compute Hpsi = -0.5 *  grad dot dpsi(G) ****
c               call Pack_c_SMul1(1,(scal1),dcpl_mb(dpsi(1)))
c               call Pack_tc_iMul2(1,dbl_mb(G(xyz)),
c     >                             dcpl_mb(dpsi(1)))
c               call Pack_cc_idot(1,psi(1+(n-1)*npack1),
c     >                             dcpl_mb(dpsi(1)),ee)
c               pmeta = pmeta - ee
c
c            end do
c         end do
c         end do
c
c*        **** deallocate stack ****
c         if (.not.BA_pop_stack(dpsi(2)))
c     >      call errquit('nwpw_meta_gga_Hpsik:popping stack',0,MA_ERR)
c
c         call Parallel_SumAll(pmeta)

         dV    = lattice_omega()/dble(nx*ny*nz)
         do ms=1,ispin
            call D3dB_rr_dot(1,dbl_mb(dfdtau(1)+(ms-1)*n2ft3d),
     >                         dbl_mb(tau(1)+(ms-1)*n2ft3d),ee)
            pmeta = pmeta + ee*dV
         end do

         if (ispin.eq.1) pmeta = pmeta + pmeta

      end if

      !write(*,*) "PMETA=",pmeta

      nwpw_meta_gga_pxc = (pmeta)
      return
      end



*     **************************************************
*     *                                                *
*     *            nwpw_meta_gga_print                 *
*     *                                                *
*     **************************************************
      subroutine nwpw_meta_gga_print(unit)
      implicit none
      integer unit

#include "nwpw_meta_gga.fh"

      if (meta_on) then
         write(unit,1001)
      end if

      return
1001  format(6x,"- Meta GGA Functional")
      end



c $Id$
