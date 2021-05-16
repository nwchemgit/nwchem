*
* $Id$

*     ****************************************************
*     *                                                  *
*     *             coulomb_screened_init                *
*     *                                                  *
*     ****************************************************
*
      subroutine coulomb_screened_init(flag,rcut,pp)
      implicit none 
      integer flag
      real*8 rcut,pp

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block used for coulomb.f ****
      integer vc_indx,vc_hndl
      common / vc_screened_block / vc_indx,vc_hndl

*     **** local variables ****
      real*8 EPSILON
c      parameter (EPSILON=1.0d0)

      logical value
      integer npack0,nfft3d,G(3),n2ft3d
      integer i,j,k,k1,k2,k3,p,q,index
      integer zero,qzero,pzero,taskid
      integer nx,ny,nz,nxh,nyh,nzh
      real*8  pi,twopi,fourpi,gg,x,y,z,temp,dV,pp2,scal1
      real*8  a(3,3),unita(3,3),unitg(3,3)
      integer gk(2),glr(2),gk2(2)

*     **** external functions ****
      integer  G_indx
      real*8   util_erf,lattice_omega,lattice_unita,lattice_unitg
      real*8   kiril_coulomb_transform, kiril_coulomb_transform0
      real*8  control_attenuation
      external G_indx
      external util_erf,lattice_omega,lattice_unita,lattice_unitg
      external kiril_coulomb_transform, kiril_coulomb_transform0
      external control_attenuation

      if (flag.eq.2) then
         EPSILON = control_attenuation()
      else
         EPSILON = 1.0d0
      end if

      call nwpw_timing_start(7)
      call Parallel2d_taskid_i(taskid)

      call D3dB_nfft3d(1,nfft3d)
      call D3dB_n2ft3d(1,n2ft3d)
      call Pack_npack(0,npack0)
      G(1) = G_indx(1)
      G(2) = G_indx(2)
      G(3) = G_indx(3)
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


      call D3dB_nz(1,nz)
      call D3dB_ny(1,ny)
      call D3dB_nx(1,nx)
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
c     call D3dB_ktoqp(1,k+1,qzero,pzero)
c     zero = (qzero-1)*(nxh+1)*ny
c    >     + j*(nxh+1)
c    >     + i+1
      call D3dB_ijktoindexp(1,i+1,j+1,k+1,zero,pzero)


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
      call D3dB_cr_fft3b(1,dcpl_mb(gk(1)))
c      call D3dB_r_SMul(1,scal1,dcpl_mb(gk(1)),dcpl_mb(gk(1)))
      call D3dB_r_SMul1(1,scal1,dcpl_mb(gk(1)))
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

         !call D3dB_ktoqp(1,k+1,q,p)
         call D3dB_ijktoindex2p(1,i+1,j+1,k+1,index,p)
         if (p .eq. taskid) then
c           index = (q-1)*(nx+2)*ny
c    >            + j    *(nx+2)
c    >            + i+1

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
      call D3dB_rr_Sum2(1,dcpl_mb(gk(1)),dbl_mb(glr(1)))


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

         !call D3dB_ktoqp(1,k+1,q,p)
         call D3dB_ijktoindex2p(1,i+1,j+1,k+1,index,p)
         if (p .eq. taskid) then
c           index = (q-1)*(nx+2)*ny
c    >            + j    *(nx+2)
c    >            + i+1

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
      call D3dB_r_Zero_Ends(1,dbl_mb(glr(1)))
      call D3dB_rc_fft3f(1,dbl_mb(glr(1)))
        

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
        call D3dB_ijktoindexp(1,i+1,j+1,k+1,zero,pzero)

*       **** short-range part of Greens function ****
        call dcopy(nfft3d,0.0d0,0,dbl_mb(gk2(1)),1)
        do i=1,nfft3d

           gg  = ( dbl_mb(G(1)+i-1)*dbl_mb(G(1)+i-1)
     >           + dbl_mb(G(2)+i-1)*dbl_mb(G(2)+i-1)
     >           + dbl_mb(G(3)+i-1)*dbl_mb(G(3)+i-1) )

           if ((pzero.eq.taskid) .and. (i.eq.zero)) then
              !temp = 8.0d0*pi*pi*dsqrt(pi*rcut)
              !temp = 0.25d0*dsqrt(pi)*lattice_omega()/(pi*pi)
              !temp = (5.633714987781071d0)*0.5d0*lattice_omega()/(pi*pi)
              !temp =(5.633714987781071d0)*0.25d0*lattice_omega()/(pi*pi)
              !temp = (5.633714987781071d0)*fourpi
              temp = (5.633714987781071d0)*
     >              (lattice_omega()**(2.0d0/3.0d0))/(pi)
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
        call D3dB_ijktoindexp(1,i+1,j+1,k+1,zero,pzero)

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
      else if (flag.eq.4) then
*     ***********************************
*     *** Set up the Green's function ***
*     ***********************************
      EPSILON=1.0d0/rcut

*     ***** find the G==0 point in the lattice *****
      i=0
      j=0
      k=0
      call D3dB_ijktoindexp(1,i+1,j+1,k+1,zero,pzero)
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

*     **************************************************************
*     ***** flag==5 Filtered exchange definitions of kernel   ******
*     **************************************************************
      else 
         !write(*,*) "Filtered Exchange not implemenented!"
         do j=1,3
         do i=1,3
            unita(i,j) = lattice_unita(i,j)
            unitg(i,j) = lattice_unitg(i,j)
         end do
         end do
         call coulomb_filter(unita,unitg,nfft3d,dbl_mb(G(1)),
     >                       dbl_mb(gk2(1)))
      end if



*     *** Pack the kernel ****
      call Pack_t_pack(0,dbl_mb(gk2(1)))
      call Pack_t_Copy(0,dbl_mb(gk2(1)),dbl_mb(vc_indx))

      value = BA_pop_stack(glr(2))
      value = value.and.BA_pop_stack(gk2(2))
      value = value.and.BA_pop_stack(gk(2))
      if (.not. value) 
     >  call errquit('coulomb_screened_init: poppingstack memory',0,
     &       MA_ERR)


      call nwpw_timing_end(7)

*     ***** DEBUG printing *****
c      call coulomb_screened_init_print(rcut,pp)



      return
      end




      subroutine coulomb_screened_end()
      implicit none
#include "bafdecls.fh"

*     **** common block used for coulomb.f ****
      integer vc_indx,vc_hndl
      common / vc_screened_block / vc_indx,vc_hndl
      logical value

      value = BA_free_heap(vc_hndl)
      return
      end


      subroutine coulomb_screened_v(dng,vc_out)
      implicit none
      complex*16 dng(*)
      complex*16 vc_out(*)

#include "bafdecls.fh"


*     **** common block used for coulomb.f ****
      integer vc_indx,vc_hndl
      common / vc_screened_block / vc_indx,vc_hndl

      call nwpw_timing_start(7)
      
      call Pack_tc_Mul(0,dbl_mb(vc_indx),dng,vc_out)
     
      call nwpw_timing_end(7)

      return 
      end


      real*8 function coulomb_screened_e(dng)
      implicit none
      complex*16 dng(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block used for coulomb.f ****
      integer vc_indx,vc_hndl
      common / vc_screened_block / vc_indx,vc_hndl

*     **** local variables ****
      integer npack0
      real*8 ec

      integer tmp1(2)
      logical value

*     **** external functions ****
      real*8  lattice_omega
      external lattice_omega

      call nwpw_timing_start(7)

      call Pack_npack(0,npack0)
      value = BA_push_get(mt_dbl,npack0,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
      
      call Pack_ct_Sqr(0,dng,dbl_mb(tmp1(1)))
      call Pack_tt_dot(0,dbl_mb(tmp1(1)),dbl_mb(vc_indx),ec)
       
      ec = 0.5d0*ec*lattice_omega()

      value = BA_pop_stack(tmp1(2))
      
      call nwpw_timing_end(7)
      coulomb_screened_e = ec
      return 
      end


      real*8 function icoulomb_screened_e(dng)
      implicit none
      complex*16 dng(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block used for coulomb.f ****
      integer vc_indx,vc_hndl
      common / vc_screened_block / vc_indx,vc_hndl

*     **** local variables ****
      integer tmp1(2)
      integer npack0
      real*8 ec

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

      call nwpw_timing_start(7)
      call Pack_npack(0,npack0)
      if (.not.BA_push_get(mt_dbl,npack0,'tmp1',tmp1(2),tmp1(1)))
     > call errquit('out of stack',0, MA_ERR)

      call Pack_ct_Sqr(0,dng,dbl_mb(tmp1(1)))
      call Pack_tt_idot(0,dbl_mb(tmp1(1)),dbl_mb(vc_indx),ec)

      ec = 0.5d0*ec*lattice_omega()

      if (.not.BA_pop_stack(tmp1(2)))
     >  call errquit('popping stack memory',0, MA_ERR)
      call nwpw_timing_end(7)

      icoulomb_screened_e = ec
      return
      end




*     *******************************
*     *                             *
*     *   kiril_coulomb_transform   *
*     *                             *
*     *******************************

*     This function returns the fourier transform of 
*
*             if flag==1  v_kiril = exp(-(r/rcut)**pp)/r
*     or      if flag==2  v_kiril = (1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)/r
*
*     Entry - gg: g squared
*             rcut: 
*             pp:
*                                           
*     Exit - returns 
*                              /infty
*                             | 
*      v_kiril(g) = (4*pi)  * | r**2 * v_kiril(r)* j0(gr) dr
*                             |
*                            / 0

      real*8 function kiril_coulomb_transform(flag,gg,rcut,pp)
      implicit none
      integer flag
      real*8 gg,rcut,pp

*     **** local variables ****
      integer nrho
      parameter (nrho=15000)

      integer i
      real*8 q,fourpi,r,sum,drho,pp2

      pp2    = pp + 2.0d0
      drho   = 2.0d0*rcut/dble(nrho)
      q      = dsqrt(gg)
      fourpi = 4.0d0*(4.0d0*datan(1.0d0))

      sum = 0.0d0
      r = 0.0d0
      if (flag.eq.1) then
       do i=2,nrho-1
         r = r + drho
         sum = sum + sin(q*r)*dexp(-(r/rcut)**pp)
        end do
        r = drho*(nrho-1)
        sum = sum + 0.5d0*sin(q*r)*dexp(-(r/rcut)**pp)
      else
       do i=2,nrho-1
c        r = (i-1)*drho
         r = r + drho
         sum = sum + sin(q*r)
     >        *(1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
        end do
        r = drho*(nrho-1)
        sum = sum + 0.5d0*sin(q*r)
     >     *(1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
      end if

      kiril_coulomb_transform = (fourpi/q)*sum*drho
      return
      end

*     *******************************
*     *                             *
*     *   kiril_coulomb_transform0  *
*     *                             *
*     *******************************

*     This function returns the fourier transform of 
*
*           if flag==1   v_kiril = exp(-(r/rcut)**pp)/r
*     or    if flag==2   v_kiril = (1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)/r
*
*     Entry - 
*             rcut: 
*             pp:
*                                           
*     Exit - returns 
*                              /infty
*                             | 
*      v_kiril(g=0) = (4*pi)* | r**2 * v_kiril(r) dr
*                             |
*                            / 0

      real*8 function kiril_coulomb_transform0(flag,rcut,pp)
      implicit none
      integer flag
      real*8 rcut,pp

*     **** local variables ****
      integer nrho
      parameter (nrho=15000)

      integer i
      real*8  fourpi,r,sum,drho,pp2

      pp2    = pp + 2.0d0
      drho   = 2.0d0*rcut/dble(nrho)
      fourpi = 4.0d0*(4.0d0*datan(1.0d0))

      sum = 0.0d0
      r   = 0.0d0
      if (flag.eq.1) then
        do i=2,nrho-1
           r = r + drho
           sum = sum + r*dexp(-(r/rcut)**pp)
        end do
        r = drho*(nrho-1)
        sum = sum + 0.5d0*r*dexp(-(r/rcut)**pp)
      else
        do i=2,nrho-1
c          r = (i-1)*drho
           r = r + drho
           sum = sum + r*
     >          *(1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
        end do
        r = drho*(nrho-1)
        sum = sum + 0.5d0*r*
     >       *(1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
      end if

      kiril_coulomb_transform0 = (fourpi)*sum*drho
      return
      end




*
      subroutine coulomb_screened_init_print(rcut,pp)
      implicit none 
#include "errquit.fh"
      real*8 rcut,pp

#include "bafdecls.fh"

*     **** common block used for coulomb.f ****
      integer vc_indx,vc_hndl
      common / vc_screened_block / vc_indx,vc_hndl

*     
      logical value
      integer i,npack0,tmp1(2),nfft3d


      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(0,npack0)

*     **** allocate memory ****
      value = BA_push_get(mt_dcpl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     >  call errquit('coulomb_screened_init_print:push stack',0, MA_ERR)


      call dcopy(2*nfft3d,0.0d0,0,dcpl_mb(tmp1(1)),1)
      do i=1,npack0
         !dcpl_mb(tmp1(1)+i-1) = dcmplx(dbl_mb(vc_indx+i-1),0.0d0)
         dcpl_mb(tmp1(1)+i-1) = dbl_mb(vc_indx+i-1)
      end do


      call Pack_c_unpack(0,dcpl_mb(tmp1(1)))

      call D3dB_cr_fft3b(1,dcpl_mb(tmp1(1)))

      call sub_coulomb_screened_init_print(dcpl_mb(tmp1(1)),rcut,pp)

      value = BA_pop_stack(tmp1(2))
      if (.not. value) 
     >  call errquit('coulomb_screened_init_print:pop stack',1, MA_ERR)


      return
      end



      subroutine sub_coulomb_screened_init_print(vc,rcut,pp)
      implicit none
      real*8 vc(*)
      real*8 rcut,pp


*     **** local variables ****
      integer nfft3d,n2ft3d
      integer i,j,k,p,q,taskid
      integer index,k1,k2,k3
      integer np1,np2,np3
      integer nph1,nph2,nph3
      real*8  a(3,3),r,pp2,dV,scal1
      real*8  x,y,z

      integer rtdb
      integer  control_rtdb
      external control_rtdb

*     **** external functions ****
      real*8   lattice_unita,lattice_omega
      external lattice_unita,lattice_omega


!$OMP MASTER
*     **** constants ****
      call Parallel_taskid(taskid)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d
      call D3dB_nx(1,np1)
      call D3dB_ny(1,np2)
      call D3dB_nz(1,np3)
      pp2    = pp + 2.0d0

      nph1 = np1/2
      nph2 = np2/2
      nph3 = np3/2

*     **** elemental vectors ****
      do i=1,3
         a(i,1) = lattice_unita(i,1)/np1
         a(i,2) = lattice_unita(i,2)/np2
         a(i,3) = lattice_unita(i,3)/np3
      end do
      scal1 = 1.0d0/dble(np1*np2*np3)
      dV    = lattice_omega()/dble(np1*np2*np3)


      !call D3dB_r_SMul(1,scal1,vc,vc)

*     **** grid points in coordination space ****
      k3=0
      k2=0
      do k1 = -nph1, nph1-1
         i = k1
         j = k2
         k = k3
         if (i.lt.0) i=i+np1
         if (j.lt.0) j=j+np2
         if (k.lt.0) k=k+np3

         !call D3dB_ktoqp(1,k+1,q,p)
         call D3dB_ijktoindex2p(1,i+1,j+1,k+1,index,p)
         if (p .eq. taskid) then
c           index = (q-1)*(np1+2)*np2
c    >            + j    *(np1+2)
c    >            + i+1

*         **** reverse y and z ****
          r = a(1,1)*k1
          !r = a(1,1)*k1 + a(1,2)*k3 + a(1,3)*k2
          !r(2) = a(2,1)*k1 + a(2,2)*k3 + a(2,3)*k2
          !r(3) = a(3,1)*k1 + a(3,2)*k3 + a(3,3)*k2

          if (dabs(r).gt.1.0d-4) then
c          write(56,*) r,vc(index),
c     >     (1.0d0-(1.0d0-dexp(-(dabs(r)/rcut)**pp2))**pp)/dabs(r) * dV
          write(56,'(4E14.4)') r,vc(index),pp/dabs(r),pp
          else
          write(56,'(4E14.4)') r,vc(index),0.0d0,pp
          end if

         end if
      end do
      call util_flush(56)

      rtdb = control_rtdb()
      call dplot_gcube_write(rtdb,"kernel.cube",0,
     >          "Screened Coulomb KERNEL",vc)

!$OMP END MASTER
!$OMP BARRIER
      return
      end





*     **********************************
*     *                                *
*     *    coulomb_screened_euv_init   *
*     *                                *
*     **********************************

      subroutine coulomb_screened_euv_init(flag,rcut,pp)
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


      call D3dB_nfft3d(1,nfft3d)
      call D3dB_n2ft3d(1,n2ft3d)
      call Pack_npack(0,npack0)
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
      call dcopy(2*nfft3d,0.0d0,0,dcpl_mb(gk(1)),1)
      call dcopy(2*nfft3d,0.0d0,0,dcpl_mb(dgk(1)),1)
      call dcopy(nfft3d,0.0d0,0,dbl_mb(gk2(1)),1)
      call dcopy(n2ft3d,0.0d0,0,dbl_mb(gr(1)),1)
      call dcopy(n2ft3d,0.0d0,0,dbl_mb(dgr(1)),1)


      call D3dB_nz(1,nz)
      call D3dB_ny(1,ny)
      call D3dB_nx(1,nx)
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
      call D3dB_ijktoindexp(1,i+1,j+1,k+1,zero,pzero)


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
      call D3dB_cr_fft3b(1,dcpl_mb(gk(1)))
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

         call D3dB_ijktoindex2p(1,i+1,j+1,k+1,index,p)
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
c      call D3dB_r_SMul(1,scal1,dbl_mb(gr(1)),dbl_mb(gr(1)))
      call D3dB_r_SMul1(1,scal1,dbl_mb(gr(1)))


*     **** Put gk in k-space ****
      call D3dB_r_Zero_Ends(1,dbl_mb(gr(1)))
      call D3dB_rc_fft3f(1,dbl_mb(gr(1)))


*     *** note that only real parts of tranformed gls are used ****
      do k=1,nfft3d
         dbl_mb(gk2(1)+k-1) =  dbl_mb(gr(1)+2*(k-1))
      end do

*     *** Pack the kernel ****
      call Pack_t_pack(0,dbl_mb(gk2(1)))
      call Pack_t_Copy(0,dbl_mb(gk2(1)),dbl_mb(vloc(1)))




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
      call D3dB_cr_fft3b(1,dcpl_mb(gk(1)))
      call D3dB_cr_fft3b(1,dcpl_mb(dgk(1)))

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

         call D3dB_ijktoindex2p(1,i+1,j+1,k+1,index,p)
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
c      call D3dB_r_SMul(1,scal1,dbl_mb(gr(1)),dbl_mb(gr(1)))
      call D3dB_r_SMul1(1,scal1,dbl_mb(gr(1)))


*     **** Put gr in k-space ****
      call D3dB_r_Zero_Ends(1,dbl_mb(gr(1)))
      call D3dB_rc_fft3f(1,dbl_mb(gr(1)))
        

*     *** note that only real parts of tranformed gls are used ****
      do k=1,nfft3d
         dbl_mb(gk2(1)+k-1) =  dbl_mb(gr(1)+2*(k-1)) 
      end do

*     *** Pack the kernel ****
      call Pack_t_pack(0,dbl_mb(gk2(1)))
      call Pack_t_Copy(0,dbl_mb(gk2(1)),dbl_mb(dvscr(u,s)))


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

         call D3dB_ijktoindex2p(1,i+1,j+1,k+1,index,p)
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

         call D3dB_ijktoindex2p(1,i+1,j+1,k+1,index,p)
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
c      call D3dB_r_SMul(1,dV,dbl_mb(gr(1)),dbl_mb(gr(1)))  !** note extra omega ***
      call D3dB_r_SMul1(1,dV,dbl_mb(gr(1)))


*     **** Put gk in k-space ****
      call D3dB_r_Zero_Ends(1,dbl_mb(gr(1)))
      call D3dB_rc_fft3f(1,dbl_mb(gr(1)))



*     *** note that only real parts of tranformed gls are used ****
      do k=1,nfft3d
         dbl_mb(gk2(1)+k-1) =  dbl_mb(gr(1)+2*(k-1))
      end do

*     *** Pack the kernel ****
      call Pack_t_pack(0,dbl_mb(gk2(1)))
c      call Pack_tt_Sum(0,dbl_mb(gk2(1)),
c     >                 dbl_mb(dvscr(u,s)),
c     >                 dbl_mb(dvscr(u,s)))
      call Pack_tt_Sum2(0,dbl_mb(gk2(1)),dbl_mb(dvscr(u,s)))


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
      call D3dB_ijktoindexp(1,i+1,j+1,k+1,zero,pzero)

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

*     *** Pack the kernel ****
      call Pack_t_pack(0,dbl_mb(gk2(1)))
      call Pack_t_Copy(0,dbl_mb(gk2(1)),dbl_mb(vloc(1)))

*     *********************************************
*     *** Set up the Green's function gradients ***
*     *********************************************
      do s=1,3
      do u=1,3

*     **** short-range part of Greens function ****
      call dcopy(nfft3d,0.0d0,0,dbl_mb(gk2(1)),1)
      do i=1,nfft3d
         gg  = ( dbl_mb(G(1)+i-1)*dbl_mb(G(1)+i-1)
     >         + dbl_mb(G(2)+i-1)*dbl_mb(G(2)+i-1)
     >         + dbl_mb(G(3)+i-1)*dbl_mb(G(3)+i-1) )
         if ((pzero.eq.taskid) .and. (i.eq.zero)) then
            temp2 = 0.0d0
         else
            temp2 = (-2.0d0*dbl_mb(G(u)+i-1)*dbl_mb(G(s)+i-1))
     >            * ( (-fourpi/(gg*gg))
     >               *(1.0d0 - dexp(-gg/(4.0d0*EPSILON**2)))
     >               + (pi/(gg*EPSILON**2))
     >                *dexp(-gg/(4.0d0*EPSILON**2)))
         end if

         dbl_mb(gk2(1)+i-1)  = temp2
      end do
      call Pack_t_pack(0,dbl_mb(gk2(1)))
      call Pack_t_Copy(0,dbl_mb(gk2(1)),dbl_mb(dvscr(u,s)))
        
      end do 
      end do 

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
*     *    coulomb_screened_euv_end    *
*     *                                *
*     **********************************

      subroutine coulomb_screened_euv_end()
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
*     *    coulomb_screened_euv        *
*     *                                *
*     **********************************

      subroutine coulomb_screened_euv(dng,euv)
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



      call Pack_npack(0,npack0)
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
      call Pack_ct_Sqr(0,dng,dbl_mb(tmp1(1)))
      call Pack_tt_dot(0,dbl_mb(tmp1(1)),dbl_mb(vloc(1)),ec)
      ec = 0.5d0*ec*lattice_omega()

*     **** Bus = Sum(G) (omega/2)*dvscr_loc(G,u,v) ****
      call dcopy(9,0.0d0,0,Bus,1)
      ss     = 0.5d0*lattice_omega()
      do u=1,3
      do s=u,3
        call Pack_tt_dot(0,dbl_mb(tmp1(1)),dbl_mb(dvscr(u,s)),sum)
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

      value = BA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('coulomb_screened_euv:popping stack',0,MA_ERR)

      return
      end



