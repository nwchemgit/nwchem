*
* $Id: coulomb_screened_small.F 27911 2015-12-11 01:40:19Z bylaska $

*     ****************************************************
*     *                                                  *
*     *             coulomb_screened_small_init          *
*     *                                                  *
*     ****************************************************
*
      subroutine coulomb_screened_small_init(flag,rcut,pp)
      implicit none 
      integer flag
      real*8 rcut,pp

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block used for coulomb.f ****
      integer vc_indx,vc_hndl
      common / vc_screened_small_block / vc_indx,vc_hndl

*     **** local variables ****
      real*8 EPSILON
c      parameter (EPSILON=1.0d0)

      logical value
      integer npack0,nfft3d,G(3),n2ft3d
      integer i,j,k,k1,k2,k3,p,q,index
      integer zero,qzero,pzero,taskid
      integer nx,ny,nz,nxh,nyh,nzh
      real*8  pi,twopi,fourpi,gg,x,y,z,temp,dV,pp2,scal1
      real*8  a(3,3)
      integer gk(2),glr(2),gk2(2)

*     **** external functions ****
      integer  Gsmall_indx
      real*8   util_erf,lattice_omega_small,lattice_unita_small
      real*8   kiril_coulomb_transform, kiril_coulomb_transform0
      real*8  control_attenuation
      external Gsmall_indx
      external util_erf,lattice_omega_small,lattice_unita_small
      external kiril_coulomb_transform, kiril_coulomb_transform0
      external control_attenuation

      if (flag.eq.2) then
         EPSILON = control_attenuation()
      else
         EPSILON = 1.0d0
      end if

      call nwpw_timing_start(7)
      call Parallel2d_taskid_i(taskid)

      call D3dB_nfft3d(3,nfft3d)
      call D3dB_n2ft3d(3,n2ft3d)
      call Pack_npack(2,npack0)
      G(1) = Gsmall_indx(1)
      G(2) = Gsmall_indx(2)
      G(3) = Gsmall_indx(3)
      pp2 = pp + 2.0d0

*     **** allocate vc memory ****
      value = BA_alloc_get(mt_dbl,npack0,'vc',vc_hndl,vc_indx)
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)

      value = BA_push_get(mt_dcpl,nfft3d,'gk',gk(2),gk(1))
      value = value.and.
     >        BA_push_get(mt_dbl,nfft3d,'gk2',gk2(2),gk2(1))
      value = value.and.
     >        BA_push_get(mt_dbl,n2ft3d,'glr',glr(2),glr(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
      call dcopy(2*nfft3d,0.0d0,0,dcpl_mb(gk(1)),1)
      call dcopy(nfft3d,0.0d0,0,dbl_mb(gk2(1)),1)
      call dcopy(n2ft3d,0.0d0,0,dbl_mb(glr(1)),1)


      call D3dB_nz(3,nz)
      call D3dB_ny(3,ny)
      call D3dB_nx(3,nx)
      nxh = nx/2
      nyh = ny/2
      nzh = nz/2

      pi     = 4.0d0*datan(1.0d0)
      twopi  = 2.0d0*pi
      fourpi = 4.0d0*pi

*     **** elemental vectors ****
      do i=1,3
         a(i,1) = lattice_unita_small(i,1)/nx
         a(i,2) = lattice_unita_small(i,2)/ny
         a(i,3) = lattice_unita_small(i,3)/nz
      end do
      dV     = lattice_omega_small()/dble(nx*ny*nz)
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
      call D3dB_ijktoindexp(3,i+1,j+1,k+1,zero,pzero)


      call dcopy(2*nfft3d,0.0d0,0,dcpl_mb(gk(1)),1)

*     **** short-range part of Greens function set only for short-range ****
      if (flag.eq.0) then
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
           
         dcpl_mb(gk(1)+i-1) = dcmplx(temp,0.0d0)
      end do
      call D3dB_cr_fft3b(3,dcpl_mb(gk(1)))
      call D3dB_r_SMul1(3,scal1,dcpl_mb(gk(1)))
      end if



*     **** long-range part of Greens function ****
      call dcopy(n2ft3d,0.0d0,0,dbl_mb(glr(1)),1)
      do k3=-nzh,nzh-1
      do k2=-nyh,nyh-1
      do k1=-nxh,nxh-1
         i = k1
         j = k2
         k = k3
         if (i.lt.0) i=i+nx
         if (j.lt.0) j=j+ny
         if (k.lt.0) k=k+nz

         call D3dB_ijktoindex2p(3,i+1,j+1,k+1,index,p)
         if (p .eq. taskid) then

            x = a(1,1)*k1 + a(1,2)*k2 + a(1,3)*k3
            y = a(2,1)*k1 + a(2,2)*k2 + a(2,3)*k3
            z = a(3,1)*k1 + a(3,2)*k2 + a(3,3)*k3
c*           **** reverse y and z ****
c            x = a(1,1)*k1 + a(1,2)*k3 + a(1,3)*k2
c            y = a(2,1)*k1 + a(2,2)*k3 + a(2,3)*k2
c            z = a(3,1)*k1 + a(3,2)*k3 + a(3,3)*k2
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
      call D3dB_rr_Sum2(3,dcpl_mb(gk(1)),dbl_mb(glr(1)))


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

         call D3dB_ijktoindex2p(3,i+1,j+1,k+1,index,p)
         if (p .eq. taskid) then

            x = a(1,1)*k1 + a(1,2)*k2 + a(1,3)*k3
            y = a(2,1)*k1 + a(2,2)*k2 + a(2,3)*k3
            z = a(3,1)*k1 + a(3,2)*k2 + a(3,3)*k3
c*           **** reverse y and z ****
c            x = a(1,1)*k1 + a(1,2)*k3 + a(1,3)*k2
c            y = a(2,1)*k1 + a(2,2)*k3 + a(2,3)*k2
c            z = a(3,1)*k1 + a(3,2)*k3 + a(3,3)*k2
            temp = dsqrt(x*x + y*y + z*z)

            dbl_mb(glr(1)+index-1) 
     >        = dbl_mb(glr(1)+index-1)
     >       * (1.0d0-(1.0d0-dexp(-(temp/rcut)**pp2))**pp)

         end if

      end do
      end do
      end do


*     **** Put glr in k-space ****
      call D3dB_r_Zero_Ends(3,dbl_mb(glr(1)))
      call D3dB_rc_fft3f(3,dbl_mb(glr(1)))
        

*     *** note that only real parts of tranformed grl are used ****
      do k=1,nfft3d
         dbl_mb(gk2(1)+k-1) =  dbl_mb(glr(1)+2*(k-1)) 
      end do




*     *********************************************************
*     ***** flag == 1 use periodic definitions of kernel ******
*     *********************************************************
      else if (flag.eq.1) then
        i=0
        j=0
        k=0
        call D3dB_ijktoindexp(3,i+1,j+1,k+1,zero,pzero)

*       **** short-range part of Greens function ****
        call dcopy(nfft3d,0.0d0,0,dbl_mb(gk2(1)),1)
        do i=1,nfft3d

           gg  = ( dbl_mb(G(1)+i-1)*dbl_mb(G(1)+i-1)
     >           + dbl_mb(G(2)+i-1)*dbl_mb(G(2)+i-1)
     >           + dbl_mb(G(3)+i-1)*dbl_mb(G(3)+i-1) )

           if ((pzero.eq.taskid) .and. (i.eq.zero)) then
              !temp = 8.0d0*pi*pi*dsqrt(pi*rcut)
              !temp = 0.25d0*dsqrt(pi)*lattice_omega_small()/(pi*pi)
              !temp = (5.633714987781071d0)*0.5d0*lattice_omega_small()/(pi*pi)
              !temp =(5.633714987781071d0)*0.25d0*lattice_omega_small()/(pi*pi)
              !temp = (5.633714987781071d0)*fourpi
              temp = (5.633714987781071d0)*
     >              (lattice_omega_small()**(2.0d0/3.0d0))/(pi)
              !temp = 0.0d0
           else
              temp = fourpi/gg
           end if
           dbl_mb(gk2(1)+i-1) = temp
        end do

*     **************************************************************
*     ***** flag==3 use periodic definitions of cutoff-kernel ******
*     **************************************************************
      else if (flag.eq.3) then
        i=0
        j=0
        k=0
        call D3dB_ijktoindexp(3,i+1,j+1,k+1,zero,pzero)

*       **** short-range part of Greens function ****
        call dcopy(nfft3d,0.0d0,0,dbl_mb(gk2(1)),1)
        do i=1,nfft3d

           gg  = ( dbl_mb(G(1)+i-1)*dbl_mb(G(1)+i-1)
     >           + dbl_mb(G(2)+i-1)*dbl_mb(G(2)+i-1)
     >           + dbl_mb(G(3)+i-1)*dbl_mb(G(3)+i-1) )

           if ((pzero.eq.taskid) .and. (i.eq.zero)) then
              temp = kiril_coulomb_transform0(flag,rcut,pp)
           else
              temp = kiril_coulomb_transform(flag,gg,rcut,pp)
           end if
           dbl_mb(gk2(1)+i-1) = temp
        end do


*     **************************************************************
*     ***** flag==4 use erfc definitions of kernel            ******
*     **************************************************************
      else
*     ***********************************
*     *** Set up the Green's function ***
*     ***********************************
      EPSILON=1.0d0/rcut

*     ***** find the G==0 point in the lattice *****
      i=0
      j=0
      k=0
      call D3dB_ijktoindexp(3,i+1,j+1,k+1,zero,pzero)
      call dcopy(nfft3d,0.0d0,0,dbl_mb(gk2(1)),1)

*     **** short-range part of Greens function , fourier transform of erfc(e*r)/r ****
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

         dbl_mb(gk2(1)+i-1) = temp
      end do
      
      end if



*     *** Pack the kernel ****
      call Pack_t_pack(2,dbl_mb(gk2(1)))
      call Pack_t_Copy(2,dbl_mb(gk2(1)),dbl_mb(vc_indx))

      value = BA_pop_stack(glr(2))
      value = value.and.BA_pop_stack(gk2(2))
      value = value.and.BA_pop_stack(gk(2))
      if (.not. value) 
     >   call errquit('coulomb_screened_small_init: poppingstack',0,
     >                MA_ERR)


      call nwpw_timing_end(7)

      return
      end



*     ****************************************************
*     *                                                  *
*     *             coulomb_screened_small_end           *
*     *                                                  *
*     ****************************************************
      subroutine coulomb_screened_small_end()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block used for coulomb.f ****
      integer vc_indx,vc_hndl
      common / vc_screened_small_block / vc_indx,vc_hndl
      logical value

      value = BA_free_heap(vc_hndl)
      if (.not. value) 
     >   call errquit('coulomb_screened_small_end: free heap',0,
     >                MA_ERR)
      return
      end


*     ****************************************************
*     *                                                  *
*     *             coulomb_screened_small_v             *
*     *                                                  *
*     ****************************************************
      subroutine coulomb_screened_small_v(dng,vc_out)
      implicit none
      complex*16 dng(*)
      complex*16 vc_out(*)

#include "bafdecls.fh"

*     **** common block used for coulomb.f ****
      integer vc_indx,vc_hndl
      common / vc_screened_small_block / vc_indx,vc_hndl

      call nwpw_timing_start(7)
      
      call Pack_tc_Mul(2,dbl_mb(vc_indx),dng,vc_out)
     
      call nwpw_timing_end(7)

      return 
      end


*     ****************************************************
*     *                                                  *
*     *             coulomb_screened_small_e             *
*     *                                                  *
*     ****************************************************
      real*8 function coulomb_screened_small_e(dng)
      implicit none
      complex*16 dng(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block used for coulomb.f ****
      integer vc_indx,vc_hndl
      common / vc_screened_small_block / vc_indx,vc_hndl

*     **** local variables ****
      integer npack0
      real*8 ec

      integer tmp1(2)
      logical value

*     **** external functions ****
      real*8   lattice_omega_small
      external lattice_omega_small

      call nwpw_timing_start(7)

      call Pack_npack(2,npack0)
      value = BA_push_get(mt_dbl,npack0,'tmp1',tmp1(2),tmp1(1))
      if (.not.value) 
     >   call errquit('coulomb_screened_small_e:out of stack',0, MA_ERR)
      
      call Pack_ct_Sqr(2,dng,dbl_mb(tmp1(1)))
      call Pack_tt_dot(2,dbl_mb(tmp1(1)),dbl_mb(vc_indx),ec)
       
      ec = 0.5d0*ec*lattice_omega_small()

      value = BA_pop_stack(tmp1(2))
      if (.not.value) 
     >   call errquit('coulomb_screened_small_e:popping stack',1,MA_ERR)
      
      call nwpw_timing_end(7)
      coulomb_screened_small_e = ec
      return 
      end


*     ****************************************************
*     *                                                  *
*     *            icoulomb_screened_small_e             *
*     *                                                  *
*     ****************************************************
      real*8 function icoulomb_screened_small_e(dng)
      implicit none
      complex*16 dng(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block used for coulomb.f ****
      integer vc_indx,vc_hndl
      common / vc_screened_small_block / vc_indx,vc_hndl

*     **** local variables ****
      integer tmp1(2)
      integer npack0
      real*8 ec

*     **** external functions ****
      real*8   lattice_omega_small
      external lattice_omega_small

      call nwpw_timing_start(7)
      call Pack_npack(2,npack0)
      if (.not.BA_push_get(mt_dbl,npack0,'tmp1',tmp1(2),tmp1(1)))
     >   call errquit('icoulomb_screened_small_e:out stack',0, MA_ERR)

      call Pack_ct_Sqr(2,dng,dbl_mb(tmp1(1)))
      call Pack_tt_idot(2,dbl_mb(tmp1(1)),dbl_mb(vc_indx),ec)

      ec = 0.5d0*ec*lattice_omega_small()

      if (.not.BA_pop_stack(tmp1(2)))
     >   call errquit('icoulomb_screened_small_e:pop stack',0,MA_ERR)
      call nwpw_timing_end(7)

      icoulomb_screened_small_e = ec
      return
      end







