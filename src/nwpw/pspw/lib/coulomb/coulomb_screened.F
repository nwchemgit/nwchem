*
* $Id: coulomb_screened.F,v 1.12 2003-10-17 22:57:27 carlfahl Exp $


*
      subroutine coulomb_screened_init(rcut,pp)
      implicit none 
#include "errquit.fh"
      real*8 rcut,pp

#include "mafdecls.fh"

*     **** common block used for coulomb.f ****
      integer vc_indx,vc_hndl
      common / vc_screened_block / vc_indx,vc_hndl

*     **** local variables ****
      real*8 EPSILON
      parameter (EPSILON=1.0d0)

      logical value
      integer npack0,nfft3d,G(3),n2ft3d
      integer i,j,k,k1,k2,k3,p,q,index
      integer zero,qzero,pzero,taskid
      integer nx,ny,nz,nxh,nyh,nzh
      real*8  pi,twopi,fourpi,gg,x,y,z,temp,dV,pp2,scal1
      real*8  a(3,3)
      integer gk(2),glr(2),gk2(2)

*     **** external functions ****
      integer  G_indx
      real*8   util_erf,lattice_omega,lattice_unita
      external G_indx
      external util_erf,lattice_omega,lattice_unita

      call nwpw_timing_start(7)
      call Parallel_taskid(taskid)

      call D3dB_nfft3d(1,nfft3d)
      call D3dB_n2ft3d(1,n2ft3d)
      call Pack_npack(0,npack0)
      G(1) = G_indx(1)
      G(2) = G_indx(2)
      G(3) = G_indx(3)
      pp2 = pp + 2.0d0

*     **** allocate vc memory ****
      value = MA_alloc_get(mt_dbl,npack0,'vc',vc_hndl,vc_indx)
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)

      value = MA_push_get(mt_dcpl,nfft3d,'gk',gk(2),gk(1))
      value = value.and.
     >        MA_push_get(mt_dbl,nfft3d,'gk2',gk2(2),gk2(1))
      value = value.and.
     >        MA_push_get(mt_dbl,n2ft3d,'glr',glr(2),glr(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)


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

*     ***********************************
*     *** Set up the Green's function ***
*     ***********************************

*     ***** find the G==0 point in the lattice *****
      i=0
      j=0
      k=0
      call D3dB_ktoqp(1,k+1,qzero,pzero)
      zero = (qzero-1)*(nxh+1)*ny
     >     + j*(nxh+1)
     >     + i+1

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
           
         dcpl_mb(gk(1)+i-1) = dcmplx(temp,0.0d0)
      end do
      call D3dB_cr_fft3b(1,dcpl_mb(gk(1)))
      call D3dB_r_SMul(1,scal1,dcpl_mb(gk(1)),dcpl_mb(gk(1)))


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

         call D3dB_ktoqp(1,k+1,q,p)
         if (p .eq. taskid) then
            index = (q-1)*(nx+2)*ny
     >            + j    *(nx+2)
     >            + i+1

*           **** reverse y and z ****
            x = a(1,1)*k1 + a(1,2)*k3 + a(1,3)*k2
            y = a(2,1)*k1 + a(2,2)*k3 + a(2,3)*k2
            z = a(3,1)*k1 + a(3,2)*k3 + a(3,3)*k2
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
      call D3dB_rr_Sum(1,dcpl_mb(gk(1)),
     >                 dbl_mb(glr(1)),
     >                 dbl_mb(glr(1)))

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

         call D3dB_ktoqp(1,k+1,q,p)
         if (p .eq. taskid) then
            index = (q-1)*(nx+2)*ny
     >            + j    *(nx+2)
     >            + i+1

*           **** reverse y and z ****
            x = a(1,1)*k1 + a(1,2)*k3 + a(1,3)*k2
            y = a(2,1)*k1 + a(2,2)*k3 + a(2,3)*k2
            z = a(3,1)*k1 + a(3,2)*k3 + a(3,3)*k2
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

*     *** Pack the kernel ****
      call Pack_t_pack(0,dbl_mb(gk2(1)))
      call Pack_t_Copy(0,dbl_mb(gk2(1)),dbl_mb(vc_indx))

      value = MA_pop_stack(glr(2))
      value = value.and.MA_pop_stack(gk2(2))
      value = value.and.MA_pop_stack(gk(2))
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
#include "mafdecls.fh"

*     **** common block used for coulomb.f ****
      integer vc_indx,vc_hndl
      common / vc_screened_block / vc_indx,vc_hndl
      logical value

      value = MA_free_heap(vc_hndl)
      return
      end


      subroutine coulomb_screened_v(dng,vc_out)
      implicit none
      complex*16 dng(*)
      complex*16 vc_out(*)

#include "mafdecls.fh"


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
#include "errquit.fh"
      complex*16 dng(*)

#include "mafdecls.fh"

*     **** common block used for coulomb.f ****
      integer vc_indx,vc_hndl
      common / vc_screened_block / vc_indx,vc_hndl


*     **** local variables ****
      integer npack0
      real*8 ec

c     real*8  tmp1(*)
      integer tmp1(2)
      logical value

*     **** external functions ****
      real*8 lattice_omega
      external lattice_omega

      
      call nwpw_timing_start(7)

      call Pack_npack(0,npack0)
      value = MA_push_get(mt_dbl,npack0,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
      
      call Pack_ct_Sqr(0,dng,dbl_mb(tmp1(1)))
      call Pack_tt_dot(0,dbl_mb(tmp1(1)),dbl_mb(vc_indx),ec)
       
      ec = 0.5d0*ec*lattice_omega()

      value = MA_pop_stack(tmp1(2))
      
      call nwpw_timing_end(7)

      coulomb_screened_e = ec
      return 
      end

*     *******************************
*     *                             *
*     *   kiril_coulomb_transform   *
*     *                             *
*     *******************************

*     This function returns the fourier transform of 
*
*            v_kiril = exp(-(r/rcut)**pp)/r
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

      real*8 function kiril_coulomb_transform(gg,rcut,pp)
      implicit none
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
      do i=2,nrho-1
         r = (i-1)*drho
c        sum = sum + sin(q*r)*dexp(-(r/rcut)**pp)
         sum = sum + sin(q*r)
     >        *(1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
      end do
      r = drho*(nrho-1)
c      sum = sum + 0.5d0*sin(q*r)*dexp(-(r/rcut)**pp)
      sum = sum + 0.5d0*sin(q*r)
     >     *(1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)

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
*            v_kiril = exp(-(r/rcut)**pp)/r
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

      real*8 function kiril_coulomb_transform0(rcut,pp)
      implicit none
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
      do i=2,nrho-1
         r = (i-1)*drho
c        sum = sum + r*dexp(-(r/rcut)**pp)
         sum = sum + r*
     >        *(1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
      end do
      r = drho*(nrho-1)
c      sum = sum + 0.5d0*r*dexp(-(r/rcut)**pp)
      sum = sum + 0.5d0*r*
     >     *(1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)

      kiril_coulomb_transform0 = (fourpi)*sum*drho
      return
      end




*
      subroutine coulomb_screened_init_print(rcut,pp)
      implicit none 
#include "errquit.fh"
      real*8 rcut,pp

#include "mafdecls.fh"

*     **** common block used for coulomb.f ****
      integer vc_indx,vc_hndl
      common / vc_screened_block / vc_indx,vc_hndl

*     
      logical value
      integer i,npack0,tmp1(2),nfft3d


      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(0,npack0)

*     **** allocate memory ****
      value = MA_push_get(mt_dcpl,nfft3d,'tmp1',tmp1(2),tmp1(1))
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

      value = MA_pop_stack(tmp1(2))
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

      integer rtdb
      integer  control_rtdb
      external control_rtdb

*     **** external functions ****
      real*8   lattice_unita,lattice_omega
      external lattice_unita,lattice_omega


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


      call D3dB_r_SMul(1,scal1,vc,vc)

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

         call D3dB_ktoqp(1,k+1,q,p)
         if (p .eq. taskid) then
            index = (q-1)*(np1+2)*np2
     >            + j    *(np1+2)
     >            + i+1

*         **** reverse y and z ****
          r = a(1,1)*k1
          !r = a(1,1)*k1 + a(1,2)*k3 + a(1,3)*k2
          !r(2) = a(2,1)*k1 + a(2,2)*k3 + a(2,3)*k2
          !r(3) = a(3,1)*k1 + a(3,2)*k3 + a(3,3)*k2

          if (dabs(r).gt.1.0d-4) then
          write(56,*) r,vc(index),
     >     (1.0d0-(1.0d0-dexp(-(dabs(r)/rcut)**pp2))**pp)/dabs(r) * dV
          else
          write(56,*) r,vc(index),0.0d0
          end if

         end if
      end do
      call util_flush(56)

      rtdb = control_rtdb()
      call dplot_gcube_write(rtdb,"kernel.cube",0,
     >          "Screened Coulomb KERNEL",vc)

      return
      end

