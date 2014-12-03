*
* $Id$
*
*     *************************************************
*     *                                               *
*     *           c_coulomb_screened_init             *
*     *                                               *
*     *************************************************
*
      subroutine c_coulomb_screened_init(flag0,rcut0,pp0)
      implicit none 
      integer flag0
      real*8 rcut0,pp0

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block ****
      real*8 rcut,pp
      integer flag
      integer vc_indx,vc_hndl
      common / c_vc_screened_block / rcut,pp,flag,vc_indx,vc_hndl

*     **** local variables ****
      integer npack0

      call nwpw_timing_start(7)
      flag = flag0
      rcut = rcut0
      pp   = pp0

*     **** allocate vc memory ****
      call Cram_npack(0,npack0)
      if (.not.BA_alloc_get(mt_dbl,npack0,'vc',vc_hndl,vc_indx))
     > call errquit('c_coulomb_screened_init:out of heap memory',
     >              0,MA_ERR)


      call nwpw_timing_end(7)
      return
      end


*
*     *************************************************
*     *                                               *
*     *           c_coulomb_screened_end              *
*     *                                               *
*     *************************************************
      subroutine c_coulomb_screened_end()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block ****
      real*8 rcut,pp
      integer flag
      integer vc_indx,vc_hndl
      common / c_vc_screened_block / rcut,pp,flag,vc_indx,vc_hndl

      if (.not.BA_free_heap(vc_hndl))
     > call errquit('c_coulomb_screened_end:errror freeing heap',
     >              0,MA_ERR)

      return
      end


*     *************************************************
*     *                                               *
*     *           c_coulomb_screened_set              *
*     *                                               *
*     *************************************************
      subroutine c_coulomb_screened_set(kv1,kv2)
      implicit none 
      real*8 kv1(*),kv2(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block used for coulomb.f ****
      real*8 rcut,pp
      integer flag
      integer vc_indx,vc_hndl
      common / c_vc_screened_block / rcut,pp,flag,vc_indx,vc_hndl

*     **** local variables ****
      real*8 EPSILON
      logical value
      integer npack0,nfft3d,G(3),n2ft3d
      integer i,j,k,k1,k2,k3,p,q,index
      integer zero,qzero,pzero,taskid
      integer nx,ny,nz,nxh,nyh,nzh
      real*8  pi,twopi,fourpi,gg,x,y,z,temp,dV,pp2,scal1
      real*8  a(3,3),gx,gy,gz
      integer gk(2),glr(2),gk2(2)

*     **** external functions ****
      integer  c_G_indx
      real*8   util_erf,lattice_omega,lattice_unita
      external c_G_indx
      external util_erf,lattice_omega,lattice_unita


      call nwpw_timing_start(7)
      call Parallel2d_taskid_i(taskid)

      call C3dB_nfft3d(1,nfft3d)
      call C3dB_n2ft3d(1,n2ft3d)
      call Cram_npack(0,npack0)
      G(1) = c_G_indx(1)
      G(2) = c_G_indx(2)
      G(3) = c_G_indx(3)
      pp2 = pp + 2.0d0

      EPSILON=1.0d0/rcut

*     **** allocate vc memory ****
      value = BA_push_get(mt_dcpl,nfft3d,'gk',gk(2),gk(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nfft3d,'glr',glr(2),glr(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nfft3d,'gk2',gk2(2),gk2(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)


      call C3dB_nz(1,nz)
      call C3dB_ny(1,ny)
      call C3dB_nx(1,nx)
      nxh = nx/2
      nyh = ny/2
      nzh = nz/2

      pi     = 4.0d0*datan(1.0d0)
      twopi  = 2.0d0*pi
      fourpi = 4.0d0*pi

*     **** elemental vectors ****
      do i=1,3
         a(i,1) = lattice_unita(i,1)/nx
         a(i,2) = lattice_unita(i,2)/ny
         a(i,3) = lattice_unita(i,3)/nz
      end do
      dV     = lattice_omega()/dble(nx*ny*nz)
      scal1  = 1.0d0/dble(nx*ny*nz)


*     ***********************************************
*     ***** use aperiodic definitions of kernel ******
*     ***********************************************
      if ((flag.eq.0).or.(flag.eq.2)) then

*     ***********************************
*     *** Set up the Green's function ***
*     ***********************************

*     ***** find the G==0 point in the lattice *****
      i=0
      j=0
      k=0
      call C3dB_ijktoindexp(1,i+1,j+1,k+1,zero,pzero)


      call dcopy(2*nfft3d,0.0d0,0,dcpl_mb(gk(1)),1)

*     **** short-range part of Greens function set only for short-range ****
      if (flag.eq.0) then
         do i=1,nfft3d
            gx  =  dbl_mb(G(1)+i-1) + kv1(1) - kv2(1)
            gy  =  dbl_mb(G(2)+i-1) + kv1(2) - kv2(2)
            gz  =  dbl_mb(G(3)+i-1) + kv1(3) - kv2(3)
            gg  = gx*gx + gy*gy + gz*gz

            if (gg<1.0d-6) then
               temp = pi/EPSILON**2
            else
               temp = (fourpi/gg)
     >              * (1.0d0 - dexp(-gg/(4.0d0*EPSILON**2)))
            end if
           
            dcpl_mb(gk(1)+i-1) = dcmplx(temp,0.0d0)
         end do
         call C3dB_cr_fft3b(1,dcpl_mb(gk(1)))
         call C3dB_c_SMul1(1,scal1,dcpl_mb(gk(1)))
      end if


*     **** long-range part of Greens function ****
      call dcopy(nfft3d,0.0d0,0,dbl_mb(glr(1)),1)
      do k3=-nzh,nzh-1
      do k2=-nyh,nyh-1
      do k1=-nxh,nxh-1
         i = k1
         j = k2
         k = k3
         if (i.lt.0) i=i+nx
         if (j.lt.0) j=j+ny
         if (k.lt.0) k=k+nz

         call C3dB_ijktoindex2p(1,i+1,j+1,k+1,index,p)
         if (p .eq. taskid) then
            x = a(1,1)*k1 + a(1,2)*k2 + a(1,3)*k3
            y = a(2,1)*k1 + a(2,2)*k2 + a(2,3)*k3
            z = a(3,1)*k1 + a(3,2)*k2 + a(3,3)*k3

            temp = dsqrt(x*x + y*y + z*z)

            if (temp .gt. 1.0d-10) then
              temp = (util_erf(EPSILON*temp)/temp) 
            else
              temp = 2.0d0*EPSILON/dsqrt(pi)
            end if
 
            dbl_mb(glr(1)+index-1) = temp*dV
 
         end if

      end do
      end do
      end do


*     **** add long-range part to short-range part  ****
      call C3dB_rc_Sum2(1,dbl_mb(glr(1)),dcpl_mb(gk(1)))


*     **** Multiply by the screening function ****
      do k3=-nzh,nzh-1
      do k2=-nyh,nyh-1
      do k1=-nxh,nxh-1
         i = k1
         j = k2
         k = k3
         if (i.lt.0) i=i+nx 
         if (j.lt.0) j=j+ny 
         if (k.lt.0) k=k+nz 

         call C3dB_ijktoindex2p(1,i+1,j+1,k+1,index,p)
         if (p .eq. taskid) then

            x = a(1,1)*k1 + a(1,2)*k2 + a(1,3)*k3
            y = a(2,1)*k1 + a(2,2)*k2 + a(2,3)*k3
            z = a(3,1)*k1 + a(3,2)*k2 + a(3,3)*k3
            temp = dsqrt(x*x + y*y + z*z)

            dbl_mb(glr(1)+index-1) 
     >        = dbl_mb(glr(1)+index-1)
     >       * (1.0d0-(1.0d0-dexp(-(temp/rcut)**pp2))**pp)

         end if

      end do
      end do
      end do


*     **** Put glr in k-space ****
      call C3dB_rc_fft3f(1,dcpl_mb(gk(1)))
        

*     *** note that only real parts of tranformed grl are used ****
      do k=1,nfft3d
         dbl_mb(gk2(1)+k-1) =  dbl_mb(glr(1)+2*(k-1)) 
      end do




*     **************************************************************
*     ***** flag==4 use erfc definitions of cutoff-kernel     ******
*     **************************************************************
      else if (flag.eq.4) then

      EPSILON=1.0d0/rcut

*     ***** find the G==0 point in the lattice *****
c      i=0
c      j=0
c      k=0
c      call C3dB_ijktoindexp(1,i+1,j+1,k+1,zero,pzero)
      call dcopy(nfft3d,0.0d0,0,dbl_mb(gk2(1)),1)

*     **** short-range part of Greens function set only for short-range ****
      do i=1,nfft3d
         gx  =  dbl_mb(G(1)+i-1) + kv1(1) - kv2(1)
         gy  =  dbl_mb(G(2)+i-1) + kv1(2) - kv2(2)
         gz  =  dbl_mb(G(3)+i-1) + kv1(3) - kv2(3)
         gg  = gx*gx + gy*gy + gz*gz

         if (gg .lt. 1.0d-6) then
            temp = pi/EPSILON**2
         else
            temp = (fourpi/gg)
     >           * (1.0d0 - dexp(-gg/(4.0d0*EPSILON**2)))
         end if
           
         dbl_mb(gk2(1)+i-1) = temp
      end do 

      end if

*     *** Cram the kernel ****
      call Cram_r_pack(0,dbl_mb(gk2(1)))
      call Cram_r_Copy(0,dbl_mb(gk2(1)),dbl_mb(vc_indx))

      value = BA_pop_stack(gk2(2))
      value = value.and.BA_pop_stack(glr(2))
      value = value.and.BA_pop_stack(gk(2))
      if (.not. value) 
     >  call errquit('c_coulomb_screened_set: poppingstack memory',0,
     &       MA_ERR)


      call nwpw_timing_end(7)

      return
      end



*     *************************************************
*     *                                               *
*     *           c_coulomb_screened_v                *
*     *                                               *
*     *************************************************
      subroutine c_coulomb_screened_v(dng,vc_out)
      implicit none
      complex*16 dng(*)
      complex*16 vc_out(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block ****
      real*8 rcut,pp
      integer flag
      integer vc_indx,vc_hndl
      common / c_vc_screened_block / rcut,pp,flag,vc_indx,vc_hndl

      call nwpw_timing_start(7)
      call Cram_rc_Mul(0,dbl_mb(vc_indx),dng,vc_out)
      call nwpw_timing_end(7)
      return 
      end

*     *************************************************
*     *                                               *
*     *           c_coulomb_screened_e                *
*     *                                               *
*     *************************************************
      real*8 function c_coulomb_screened_e(dng)
      implicit none
      complex*16 dng(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block ****
      real*8 rcut,pp
      integer flag
      integer vc_indx,vc_hndl
      common / c_vc_screened_block / rcut,pp,flag,vc_indx,vc_hndl

*     **** local variables ****
      integer npack0
      real*8 ec
      integer tmp1(2)

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

      write(*,*) "fera,npack0=",npack0
      call util_flush(6)
      call nwpw_timing_start(7)

      call Cram_npack(0,npack0)
      write(*,*) "fera,npack0=",npack0
      if (.not.BA_push_get(mt_dbl,npack0,'tmp1',tmp1(2),tmp1(1)))
     > call errquit('c_coulomb_screened_e:out of stack memory',0,MA_ERR)
      
      call Cram_cr_Sqr(0,dng,dbl_mb(tmp1(1)))
      call Cram_rr_dot(0,dbl_mb(tmp1(1)),dbl_mb(vc_indx),ec)
      ec = 0.5d0*ec*lattice_omega()

      if (.not.BA_pop_stack(tmp1(2)))
     > call errquit('c_coulomb_screened_e:popping stack',0,MA_ERR)
      
      call nwpw_timing_end(7)
      c_coulomb_screened_e = ec
      return 
      end

*     *************************************************
*     *                                               *
*     *           c_icoulomb_screened_e               *
*     *                                               *
*     *************************************************
      real*8 function c_icoulomb_screened_e(dng)
      implicit none
      complex*16 dng(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block ****
      real*8 rcut,pp
      integer flag
      integer vc_indx,vc_hndl
      common / c_vc_screened_block / rcut,pp,flag,vc_indx,vc_hndl

*     **** local variables ****
      integer tmp1(2)
      integer npack0
      real*8 ec

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

      call nwpw_timing_start(7)
      call Cram_npack(0,npack0)
      if (.not.BA_push_get(mt_dbl,npack0,'tmp1',tmp1(2),tmp1(1)))
     > call errquit('out of stack',0, MA_ERR)

      call Cram_cr_Sqr(0,dng,dbl_mb(tmp1(1)))
      call Cram_rr_idot(0,dbl_mb(tmp1(1)),dbl_mb(vc_indx),ec)

      ec = 0.5d0*ec*lattice_omega()

      if (.not.BA_pop_stack(tmp1(2)))
     >  call errquit('popping stack memory',0, MA_ERR)
      call nwpw_timing_end(7)

      c_icoulomb_screened_e = ec
      return
      end




*     **********************************
*     *                                *
*     *   c_coulomb_screened_euv_init  *
*     *                                *
*     **********************************

      subroutine c_coulomb_screened_euv_init(flag,rcut,pp)
      implicit none 
      integer flag
      real*8 rcut,pp

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block used for coulomb_screened_euv ****
      integer dvscr_loc(2),vloc(2)
      common / vc_screened_euv_block / dvscr_loc,vloc

*     **** local variables ****
      real*8 EPSILON
      parameter (EPSILON=1.0d0)   !** should be inputed by user! **

      logical value
      integer npack0,nfft3d,G(3),n2ft3d,u,s
      integer i,j,k,k1,k2,k3,p,q,index
      integer zero,qzero,pzero,taskid
      integer nx,ny,nz,nxh,nyh,nzh
      real*8  pi,twopi,fourpi,gg,x,y,z,temp,temp2,dV,pp2,scal1
      real*8  a(3,3),xyz(3),r
      integer gk(2),dgk(2),gr(2),dgr(2),gk2(2),dvscr(3,3)

*     **** external functions ****
      integer  G_indx
      real*8   util_erf,lattice_omega,lattice_unita
      real*8   kiril_coulomb_transform, kiril_coulomb_transform0
      external G_indx
      external util_erf,lattice_omega,lattice_unita
      external kiril_coulomb_transform, kiril_coulomb_transform0

      call nwpw_timing_start(7)
      call Parallel2d_taskid_i(taskid)


      call C3dB_nfft3d(1,nfft3d)
      call C3dB_n2ft3d(1,n2ft3d)
      call Cram_npack(0,npack0)
      G(1) = G_indx(1)
      G(2) = G_indx(2)
      G(3) = G_indx(3)
      pp2 = pp + 2.0d0

*     **** allocate vc_screened_euv_block memory ****
      value = BA_alloc_get(mt_dbl,9*npack0,'dvscr_loc',
     >                     dvscr_loc(2),dvscr_loc(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,npack0,'vloc',
     >                     vloc(2),vloc(1))
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)
      j = 0
      do s=1,3
      do u=1,3
         dvscr(u,s)  = dvscr_loc(1)  + j
         j = j + npack0
      end do
      end do

      value = BA_push_get(mt_dcpl,nfft3d,'gk',gk(2),gk(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nfft3d,'dgk',dgk(2),dgk(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nfft3d,'gk2',gk2(2),gk2(1))
      value = value.and.
     >        BA_push_get(mt_dbl,n2ft3d,'gr',gr(2),gr(1))
      value = value.and.
     >        BA_push_get(mt_dbl,n2ft3d,'dgr',dgr(2),dgr(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)


      call C3dB_nz(1,nz)
      call C3dB_ny(1,ny)
      call C3dB_nx(1,nx)
      nxh = nx/2
      nyh = ny/2
      nzh = nz/2

      pi     = 4.0d0*datan(1.0d0)
      twopi  = 2.0d0*pi
      fourpi = 4.0d0*pi

*     **** elemental vectors ****
      do i=1,3
         a(i,1) = lattice_unita(i,1)/nx
         a(i,2) = lattice_unita(i,2)/ny
         a(i,3) = lattice_unita(i,3)/nz
      end do
      dV     = lattice_omega()/dble(nx*ny*nz)
      scal1  = 1.0d0/dble(nx*ny*nz)


*     ***********************************************
*     ***** use aperiodic definitions of kernel ******
*     ***********************************************
      if (flag.eq.0) then


*     ***********************************************
*     *** Set up the short range Green's function ***
*     ***********************************************

*     ***** find the G==0 point in the lattice *****
      i=0
      j=0
      k=0
      call C3dB_ijktoindexp(1,i+1,j+1,k+1,zero,pzero)


*     **** short-range part of Greens function ****
      call dcopy(2*nfft3d,0.0d0,0,dcpl_mb(gk(1)),1)
      do i=1,nfft3d

         gg  = ( dbl_mb(G(1)+i-1)*dbl_mb(G(1)+i-1)
     >         + dbl_mb(G(2)+i-1)*dbl_mb(G(2)+i-1)
     >         + dbl_mb(G(3)+i-1)*dbl_mb(G(3)+i-1) )

         if ((pzero.eq.taskid) .and. (i.eq.zero)) then
            temp = pi/EPSILON**2
         else
            temp = (fourpi/gg)
     >           * (1.0d0 - dexp(-gg/(4.0d0*EPSILON**2)))
         end if

         dcpl_mb(gk(1)+i-1)  = dcmplx(temp,0.0d0)
      end do
      call C3dB_cr_fft3b(1,dcpl_mb(gk(1)))
      call dcopy(n2ft3d,dcpl_mb(gk(1)), 1,dbl_mb(gr(1)), 1)


*     **** Multiply by the screening function ****
      do k3=-nzh,nzh-1
      do k2=-nyh,nyh-1
      do k1=-nxh,nxh-1
         i = k1
         j = k2
         k = k3
         if (i.lt.0) i=i+nx
         if (j.lt.0) j=j+ny
         if (k.lt.0) k=k+nz

         call C3dB_ijktoindex2p(1,i+1,j+1,k+1,index,p)
         if (p .eq. taskid) then
            x = a(1,1)*k1 + a(1,2)*k2 + a(1,3)*k3
            y = a(2,1)*k1 + a(2,2)*k2 + a(2,3)*k3
            z = a(3,1)*k1 + a(3,2)*k2 + a(3,3)*k3
            r = dsqrt(x*x + y*y + z*z)
            dbl_mb(gr(1)+index-1)
     >        = dbl_mb(gr(1)+index-1)
     >         *(1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
         end if

      end do
      end do
      end do
      call C3dB_r_SMul1(1,scal1,dbl_mb(gr(1)))


*     **** Put gk in k-space ****
      call C3dB_rc_fft3f(1,dbl_mb(gr(1)))


*     *** note that only real parts of tranformed gls are used ****
      do k=1,nfft3d
         dbl_mb(gk2(1)+k-1) =  dbl_mb(gr(1)+2*(k-1))
      end do

*     *** Cram the kernel ****
      call Cram_r_pack(0,dbl_mb(gk2(1)))
      call Cram_r_Copy(0,dbl_mb(gk2(1)),dbl_mb(vloc(1)))




*     *********************************************
*     *** Set up the Green's function gradients ***
*     *********************************************
      do s=1,3
      do u=1,3


*     **** short-range part of Greens function ****
      call dcopy(2*nfft3d,0.0d0,0,dcpl_mb(gk(1)),1)
      call dcopy(2*nfft3d,0.0d0,0,dcpl_mb(dgk(1)),1)
      do i=1,nfft3d

         gg  = ( dbl_mb(G(1)+i-1)*dbl_mb(G(1)+i-1)
     >         + dbl_mb(G(2)+i-1)*dbl_mb(G(2)+i-1)
     >         + dbl_mb(G(3)+i-1)*dbl_mb(G(3)+i-1) )

         if ((pzero.eq.taskid) .and. (i.eq.zero)) then
            temp = pi/EPSILON**2
            temp2 = 0.0d0
         else
            temp = (fourpi/gg)
     >           * (1.0d0 - dexp(-gg/(4.0d0*EPSILON**2)))
            temp2 = (-2.0d0*dbl_mb(G(u)+i-1)*dbl_mb(G(s)+i-1))
     >            * ( (-fourpi/(gg*gg))
     >               *(1.0d0 - dexp(-gg/(4.0d0*EPSILON**2)))
     >               + (pi/(gg*EPSILON**2))
     >                *dexp(-gg/(4.0d0*EPSILON**2)))
         end if

         dcpl_mb(gk(1)+i-1)  = dcmplx(temp,0.0d0)
         dcpl_mb(dgk(1)+i-1) = dcmplx(temp2,0.0d0)
      end do
      call C3dB_cr_fft3b(1,dcpl_mb(gk(1)))
      call C3dB_cr_fft3b(1,dcpl_mb(dgk(1)))

      call dcopy(n2ft3d,dcpl_mb(gk(1)), 1,dbl_mb(gr(1)), 1)
      call dcopy(n2ft3d,dcpl_mb(dgk(1)),1,dbl_mb(dgr(1)),1)


*     **** Multiply by the screening function ****
      do k3=-nzh,nzh-1
      do k2=-nyh,nyh-1
      do k1=-nxh,nxh-1
         i = k1
         j = k2
         k = k3
         if (i.lt.0) i=i+nx 
         if (j.lt.0) j=j+ny 
         if (k.lt.0) k=k+nz 

         call C3dB_ijktoindex2p(1,i+1,j+1,k+1,index,p)
         if (p .eq. taskid) then
            x = a(1,1)*k1 + a(1,2)*k2 + a(1,3)*k3
            y = a(2,1)*k1 + a(2,2)*k2 + a(2,3)*k3
            z = a(3,1)*k1 + a(3,2)*k2 + a(3,3)*k3
            r = dsqrt(x*x + y*y + z*z)
            xyz(1) = x
            xyz(2) = y
            xyz(3) = z

            if (r.gt.1.0d-10) then
              gg = -dexp(-(r/rcut)**pp2)
              gg = gg*(1.0d0-dexp(-(r/rcut)**pp2))**(pp-1)
              gg = gg*(pp*pp2)*(r/rcut)**(pp+1)
              gg = gg/rcut
              dbl_mb(gr(1)+index-1) = 
     >           dbl_mb(gr(1)+index-1)*gg*(xyz(u)*xyz(s))/r
            else
              dbl_mb(gr(1)+index-1) = 0.0d0
            end if

            dbl_mb(gr(1)+index-1) 
     >        = dbl_mb(gr(1)+index-1)
     >        + dbl_mb(dgr(1)+index-1)
     >         *(1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)

         end if

      end do
      end do
      end do
      call C3dB_r_SMul1(1,scal1,dbl_mb(gr(1)))


*     **** Put gr in k-space ****
      call C3dB_rc_fft3f(1,dbl_mb(gr(1)))
        

*     *** note that only real parts of tranformed gls are used ****
      do k=1,nfft3d
         dbl_mb(gk2(1)+k-1) =  dbl_mb(gr(1)+2*(k-1)) 
      end do

*     *** Cram the kernel ****
      call Cram_r_pack(0,dbl_mb(gk2(1)))
      call Cram_r_Copy(0,dbl_mb(gk2(1)),dbl_mb(dvscr(u,s)))


*     **** long-range part of Greens function ****
      call dcopy(n2ft3d,0.0d0,0,dbl_mb(gr(1)),1)
      call dcopy(n2ft3d,0.0d0,0,dbl_mb(dgr(1)),1)
      do k3=-nzh,nzh-1
      do k2=-nyh,nyh-1
      do k1=-nxh,nxh-1
         i = k1
         j = k2
         k = k3
         if (i.lt.0) i=i+nx
         if (j.lt.0) j=j+ny
         if (k.lt.0) k=k+nz

         call C3dB_ijktoindex2p(1,i+1,j+1,k+1,index,p)
         if (p .eq. taskid) then
            x = a(1,1)*k1 + a(1,2)*k2 + a(1,3)*k3
            y = a(2,1)*k1 + a(2,2)*k2 + a(2,3)*k3
            z = a(3,1)*k1 + a(3,2)*k2 + a(3,3)*k3
            r = dsqrt(x*x + y*y + z*z)
            xyz(1) = x
            xyz(2) = y
            xyz(3) = z

            if (r .gt. 1.0d-10) then
              temp  = (util_erf(EPSILON*r)/r)
              temp2 = ((2.0d0*EPSILON/dsqrt(pi))
     >              *dexp(-EPSILON**2*r**2)/r
     >             - (util_erf(EPSILON*r)/r**2))

            else
              temp  = 2.0d0*EPSILON/dsqrt(pi)
              temp2 = 0.0d0
            end if

            dbl_mb(gr(1)+index-1)  = temp
            dbl_mb(dgr(1)+index-1) = temp2
         end if

      end do
      end do
      end do


*     **** Multiply by the screening function ****
      do k3=-nzh,nzh-1
      do k2=-nyh,nyh-1
      do k1=-nxh,nxh-1
         i = k1
         j = k2
         k = k3
         if (i.lt.0) i=i+nx
         if (j.lt.0) j=j+ny
         if (k.lt.0) k=k+nz

         call C3dB_ijktoindex2p(1,i+1,j+1,k+1,index,p)
         if (p .eq. taskid) then
            x = a(1,1)*k1 + a(1,2)*k2 + a(1,3)*k3
            y = a(2,1)*k1 + a(2,2)*k2 + a(2,3)*k3
            z = a(3,1)*k1 + a(3,2)*k2 + a(3,3)*k3
            r = dsqrt(x*x + y*y + z*z)
            xyz(1) = x
            xyz(2) = y
            xyz(3) = z

            if (r.gt.1.0d-10) then
              gg = -dexp(-(r/rcut)**pp2)
              gg = gg*(1.0d0-dexp(-(r/rcut)**pp2))**(pp-1)
              gg = gg*(pp*pp2)*(r/rcut)**(pp+1)
              gg = gg/rcut
              dbl_mb(gr(1)+index-1) =
     >           dbl_mb(gr(1)+index-1)*gg
            else
              dbl_mb(gr(1)+index-1) = 0.0d0
            end if

            dbl_mb(gr(1)+index-1)
     >        = dbl_mb(gr(1)+index-1)
     >        + dbl_mb(dgr(1)+index-1)
     >         *(1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)

            if (r.gt.1.0d-10) then
              dbl_mb(gr(1)+index-1)
     >          = dbl_mb(gr(1)+index-1)*xyz(u)*xyz(s)/r
            else
              dbl_mb(gr(1)+index-1) = 0.0d0
            end if

         end if

      end do
      end do
      end do
      call C3dB_r_SMul1(1,dV,dbl_mb(gr(1)))


*     **** Put gk in k-space ****
      call C3dB_rc_fft3f(1,dbl_mb(gr(1)))


*     *** note that only real parts of tranformed gls are used ****
      do k=1,nfft3d
         dbl_mb(gk2(1)+k-1) =  dbl_mb(gr(1)+2*(k-1))
      end do

*     *** Cram the kernel ****
      call Cram_r_pack(0,dbl_mb(gk2(1)))
      call Cram_rr_Sum2(0,dbl_mb(gk2(1)),dbl_mb(dvscr(u,s)))

      end do 
      end do 

*     *******************************************
*     ***** use erfc definitions of kernel ******
*     ******************************************* 
      else if (flag.eq.4) then
*     ***********************************************
*     *** Set up the short range Green's function ***
*     ***********************************************
*     ***** find the G==0 point in the lattice *****
      i=0
      j=0
      k=0
      call C3dB_ijktoindexp(1,i+1,j+1,k+1,zero,pzero)

*     **** short-range part of Greens function ****
      call dcopy(2*nfft3d,0.0d0,0,dcpl_mb(gk(1)),1)
      do i=1,nfft3d

         gg  = ( dbl_mb(G(1)+i-1)*dbl_mb(G(1)+i-1)
     >         + dbl_mb(G(2)+i-1)*dbl_mb(G(2)+i-1)
     >         + dbl_mb(G(3)+i-1)*dbl_mb(G(3)+i-1) )

         if ((pzero.eq.taskid) .and. (i.eq.zero)) then
            temp = pi/EPSILON**2
         else
            temp = (fourpi/gg)
     >           * (1.0d0 - dexp(-gg/(4.0d0*EPSILON**2)))
         end if

         dbl_mb(gk2(1)+i-1)  = temp
      end do

*     *** Cram the kernel ****
      call Cram_r_pack(0,dbl_mb(gk2(1)))
      call Cram_r_Copy(0,dbl_mb(gk2(1)),dbl_mb(vloc(1)))






*     ***********************************************
*     ***** use periodic definitions of kernel ******
*     ***********************************************
      else
        write(*,*) "periodic screening stresses not implemented!"
        call errquit('coulomb_screened_euv_init:',0,0)
      end if



*     *** free up stack memory ***
      value =           BA_pop_stack(dgr(2))
      value =           BA_pop_stack(gr(2))
      value = value.and.BA_pop_stack(gk2(2))
      value = value.and.BA_pop_stack(dgk(2))
      value = value.and.BA_pop_stack(gk(2))
      if (.not. value) 
     >  call errquit('coulomb_screened_euv_init: popping stack',0,
     &       MA_ERR)

      call nwpw_timing_end(7)
      return
      end


*     **********************************
*     *                                *
*     *   c_coulomb_screened_euv_end   *
*     *                                *
*     **********************************

      subroutine c_coulomb_screened_euv_end()
      implicit none
#include "bafdecls.fh"

*     **** common block used for coulomb_screened_euv ****
      integer dvscr_loc(2),vloc(2)
      common / vc_screened_euv_block / dvscr_loc,vloc

      logical value

      value = BA_free_heap(dvscr_loc(2))
      value = value.and.BA_free_heap(vloc(2))
      if (.not. value) 
     >  call errquit('coulomb_screened_euv_end: freeing heap',0,0)
      return
      end



*     **********************************
*     *                                *
*     *    c_coulomb_screened_euv      *
*     *                                *
*     **********************************
      subroutine c_coulomb_screened_euv(dng,euv)
      implicit none 
      complex*16 dng(*)
      real*8 euv(3,3)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block used for coulomb_screened_euv ****
      integer dvscr_loc(2),vloc(2)
      common / vc_screened_euv_block / dvscr_loc,vloc


*     **** local variables ****
      integer npack0
      integer i,j
      integer u,v,s
      logical value

      real*8 pi,scal,ss,sum
      real*8 hm(3,3),Bus(3,3),ec
      integer tmp1(2),dvscr(3,3)

*     **** external functions ****
      real*8   lattice_unitg,lattice_omega
      external lattice_unitg,lattice_omega


      call Cram_npack(0,npack0)
      value = BA_push_get(mt_dbl,npack0,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     > call errquit('coulomb_screened_euv:out of stack memory',0,MA_ERR)
      j = 0
      do s=1,3
      do u=1,3
         dvscr(u,s)  = dvscr_loc(1)  + j
         j = j + npack0
      end do
      end do

*     *** compute local screened coulomb energy ***
      call Cram_cr_Sqr(0,dng,dbl_mb(tmp1(1)))
      call Cram_rr_dot(0,dbl_mb(tmp1(1)),dbl_mb(vloc(1)),ec)
      ec = 0.5d0*ec*lattice_omega()

*     **** Bus = Sum(G) (omega/2)*dvscr_loc(G,u,v) ****
      call dcopy(9,0.0d0,0,Bus,1)
      ss     = 0.5d0*lattice_omega()
      do u=1,3
      do s=u,3
        call Cram_rr_dot(0,dbl_mb(tmp1(1)),dbl_mb(dvscr(u,s)),sum)
        Bus(u,s) = ss*sum
      end do
      end do
      do u=1,3
      do s=u+1,3
         Bus(s,u) = Bus(u,s)
      end do
      end do

      pi     = 4.0d0*datan(1.0d0)
      scal   = 1.0d0/(2.0d0*pi)

*     *** define hm ****
      do j=1,3
      do i=1,3
         hm(i,j) = scal*lattice_unitg(i,j)
      end do
      end do      

      do v=1,3
      do u=1,3
         euv(u,v) = -ec*hm(u,v)
         do s=1,3
            euv(u,v) = euv(u,v) + Bus(u,s)*hm(s,v)
         end do
      end do
      end do

      if (.not.BA_pop_stack(tmp1(2)))
     > call errquit('c_coulomb_screened_euv:popping stack',0,MA_ERR)

      return
      end



