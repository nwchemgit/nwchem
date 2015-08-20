
*
* $Id$
*

*     ***************************
*     *                         *
*     *      generate_ELF       *
*     *                         *
*     ***************************

*   This routine returns the electron localization function (ELF) defined 
* by Becke and Edgecombe, J. Chem. Phys., vol 92, page 5397-5403, 1990.
*
*    ELF = (1+ (D/D0)^2)^-1
*
*    where
*
*       D0 = (3/5)*(6*pi^2)^(2/3) * dn^(5/3)
*       D  = Sum(i) |del*psi_i|^2  - 0.25*|del*dn|^2/dn
*
      subroutine generate_ELF(npack1,ne,psi,dn,elf)
      implicit none
#include "errquit.fh"
      integer    npack1,ne
      complex*16 psi(npack1,ne)
      real*8     dn(*)
      real*8     elf(*)

#include "bafdecls.fh"

*     **** local variables ****
      logical value
      integer i,k,npack0,nfft3d,n2ft3d
      integer dng(2),dng2(2),dns(2),dnsp(2)
      integer nx,ny,nz
      real*8 tt,vv,uu,scaluu
      real*8 scal1,scal2

*     **** external functions ****
      real*8  lattice_omega
      external lattice_omega

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()
      
      call Pack_npack(0,npack0)
      call D3dB_nfft3d(1,nfft3d)
      call D3dB_n2ft3d(1,n2ft3d)

      call dcopy(n2ft3d,0.0d0,0,elf,1)

      value = BA_push_get(mt_dbl,n2ft3d,'dns',dns(2),dns(1))
      value = value.and.
     >        BA_push_get(mt_dbl,n2ft3d,'dnsp',dnsp(2),dnsp(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nfft3d,'dng',dng(2),dng(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nfft3d,'dng2',dng2(2),dng2(1))
      if (.not.value) call errquit('generate_ELF:pushing stack',0,
     &       MA_ERR)
      call Ursenbach_smoother(n2ft3d,dn,
     >                        dbl_mb(dns(1)),
     >                        dbl_mb(dnsp(1)))
      do i=1,n2ft3d
       dbl_mb(dns(1)+i-1) = log(dn(i))
      end do

*     **** elf <-- elf + |del*rho|^2/rho ****
      do i=1,n2ft3d
       dbl_mb(dns(1)+i-1) = log(dn(i))
      end do
      !call D3dB_r_SMul(1,scal1,dbl_mb(dns(1)),dcpl_mb(dng(1)))
      call D3dB_r_SMul(1,scal1,dbl_mb(int(dn(1))),dcpl_mb(dng(1)))
      call D3dB_r_SMul(1,scal1,dbl_mb(dns(1)),dcpl_mb(dng2(1)))
      call D3dB_r_Zero_Ends(1,dcpl_mb(dng(1)))
      call D3dB_r_Zero_Ends(1,dcpl_mb(dng2(1)))
      call D3dB_rc_fft3f(1,dcpl_mb(dng(1)))
      call D3dB_rc_fft3f(1,dcpl_mb(dng2(1)))
      call Pack_c_pack(0,dcpl_mb(dng(1)))
      call Pack_c_pack(0,dcpl_mb(dng2(1)))
      call gradient_energy_density2(0,npack0,
     >                             dcpl_mb(dng(1)),
     >                             dcpl_mb(dng2(2)),
     >                             n2ft3d,elf)
      !call D3dB_rr_Mul(1,dbl_mb(dnsp(1)),elf,elf)
      !call D3dB_rr_Mul(1,dbl_mb(dnsp(1)),elf,elf)

*     **** elf <-- -0.25*elf ****
      !call D3dB_rr_Divide(1,elf,dn,elf)
c      call D3dB_r_SMul(1,(-0.25d0),elf,elf)
      call D3dB_r_SMul1(1,(-0.25d0),elf)

*     **** elf <-- elf + Sum(i) |del*psi_i|^2 ****
      do i=1,ne
         call dcopy(n2ft3d,0.0d0,0,dcpl_mb(dng(1)),1)
         call Pack_c_SMul(1,dsqrt(scal2),psi(1,i),dcpl_mb(dng(1)))
         call gradient_energy_density(1,npack1,dcpl_mb(dng(1)),
     >                                n2ft3d,elf)
      end do      

*     **** elf <-- (1+(elf/D0)^2)^(-1) ****
      scaluu = 4.0d0*datan(1.0d0)
      scaluu = (3.0d0/5.0d0)*(6.0d0*scaluu**2)**(2.0d0/3.0d0)
      do k=1,n2ft3d
         tt = elf(k)
         vv = dn(k)
         vv = scaluu*(vv**(5.0d0/3.0d0))
         !uu = 1.0d0+(tt/vv)**2
         !elf(k) = 1.0d0/uu
         uu = vv*vv + tt*tt
         uu = 1.0d0/uu
         elf(k) = vv*vv*uu
      end do

 10   continue

      value =           BA_pop_stack(dng2(2))
      value = value.and.BA_pop_stack(dng(2))
      value = value.and.BA_pop_stack(dnsp(2))
      value = value.and.BA_pop_stack(dns(2))
      if (.not.value) call errquit('generate_ELF:popping stack',1,
     &       MA_ERR)

      return
      end


*     ***********************************
*     *                                 *
*     *     gradient_energy_density     *
*     *                                 *
*     ***********************************
* 
*   This routine calculates
*
*     df = df + |del*f|^2
*
      subroutine gradient_energy_density(nb,npack,f,
     >                                   n2ft3d,df)
      implicit none
#include "errquit.fh"
      integer nb,npack
      complex*16 f(*)
      integer n2ft3d
      real*8 df(*)

#include "bafdecls.fh"

*     **** local variables ****
      logical value
      integer nfft3d
      integer fx(2),fy(2),fz(2)
      integer Gx(2),Gy(2),Gz(2)

*     **** external functions ****
      integer  G_indx
      external G_indx


      nfft3d = n2ft3d/2
      value =           BA_push_get(mt_dbl,n2ft3d,
     >                               'fx',fx(2),fx(1))
      value = value.and.BA_push_get(mt_dbl,n2ft3d,
     >                                'fy',fy(2),fy(1))
      value = value.and.BA_push_get(mt_dbl,n2ft3d,
     >                                'fz',fz(2),fz(1))
      value = value.and.BA_push_get(mt_dbl,nfft3d,
     >                                'Gx',Gx(2),Gx(1))
      value = value.and.BA_push_get(mt_dbl,nfft3d,
     >                                'Gy',Gy(2),Gy(1))
      value = value.and.BA_push_get(mt_dbl,nfft3d,
     >                                'Gz',Gz(2),Gz(1))
      if (.not.value) 
     >  call errquit('gradient_energy_density:push stack',0, MA_ERR)

*     **** define Gx,Gy,Gz in packed space ****
      call dcopy(nfft3d,dbl_mb(G_indx(1)),1,dbl_mb(Gx(1)),1)
      call dcopy(nfft3d,dbl_mb(G_indx(2)),1,dbl_mb(Gy(1)),1)
      call dcopy(nfft3d,dbl_mb(G_indx(3)),1,dbl_mb(Gz(1)),1)
      call Pack_t_pack(nb,dbl_mb(Gx(1)))
      call Pack_t_pack(nb,dbl_mb(Gy(1)))
      call Pack_t_pack(nb,dbl_mb(Gz(1)))

*     ***** calculate:        ****
*     **** fx(G)=i*Gx(G)*f(G) ****
*     **** fy(G)=i*Gy(G)*f(G) ****
*     **** fz(G)=i*Gz(G)*f(G) ****
      call Pack_itc_Mul(nb,dbl_mb(Gx(1)),f,dbl_mb(fx(1)))
      call Pack_itc_Mul(nb,dbl_mb(Gy(1)),f,dbl_mb(fy(1)))
      call Pack_itc_Mul(nb,dbl_mb(Gz(1)),f,dbl_mb(fz(1)))

      value =           BA_pop_stack(Gz(2))
      value = value.and.BA_pop_stack(Gy(2))
      value = value.and.BA_pop_stack(Gx(2))
      if (.not.value) 
     >  call errquit('gradient_energy_density:pop stack',1, MA_ERR)

*     **** Fourier transform fx,fy,fz to real space ****
      call Pack_c_unpack(nb,dbl_mb(fx(1)))
      call Pack_c_unpack(nb,dbl_mb(fy(1)))
      call Pack_c_unpack(nb,dbl_mb(fz(1)))
      call D3dB_cr_fft3b(1,dbl_mb(fx(1)))
      call D3dB_cr_fft3b(1,dbl_mb(fy(1)))
      call D3dB_cr_fft3b(1,dbl_mb(fz(1)))

*     **** calculate df = df + fx*fx + fy*fy + fz*fz ****
c      call D3dB_rr_Sqr(1,dbl_mb(fx(1)),dbl_mb(fx(1)))
c      call D3dB_rr_Sqr(1,dbl_mb(fy(1)),dbl_mb(fy(1)))
c      call D3dB_rr_Sqr(1,dbl_mb(fz(1)),dbl_mb(fz(1)))
c      call D3dB_rr_Sum(1,dbl_mb(fx(1)),df,df)
c      call D3dB_rr_Sum(1,dbl_mb(fy(1)),df,df)
c      call D3dB_rr_Sum(1,dbl_mb(fz(1)),df,df)
      call D3dB_rr_Sqr1(1,dbl_mb(fx(1)))
      call D3dB_rr_Sqr1(1,dbl_mb(fy(1)))
      call D3dB_rr_Sqr1(1,dbl_mb(fz(1)))
      call D3dB_rr_Sum2(1,dbl_mb(fx(1)),df)
      call D3dB_rr_Sum2(1,dbl_mb(fy(1)),df)
      call D3dB_rr_Sum2(1,dbl_mb(fz(1)),df)

      value =           BA_pop_stack(fz(2))
      value = value.and.BA_pop_stack(fy(2))
      value = value.and.BA_pop_stack(fx(2))
      if (.not.value) 
     >  call errquit('gradient_energy_density:pop stack',2, MA_ERR)

      return
      end

*     ***********************************
*     *                                 *
*     *     gradient_energy_density2    *
*     *                                 *
*     ***********************************
* 
*   This routine calculates
*
*     df = df + del*f * del*h
*
      subroutine gradient_energy_density2(nb,npack,f,h,
     >                                   n2ft3d,df)
      implicit none
#include "errquit.fh"
      integer nb,npack
      complex*16 f(*),h(*)
      integer n2ft3d
      real*8 df(*)

#include "bafdecls.fh"

*     **** local variables ****
      logical value
      integer nfft3d
      integer fx(2),fy(2),fz(2)
      integer hx(2),hy(2),hz(2)
      integer Gx(2),Gy(2),Gz(2)

*     **** external functions ****
      integer  G_indx
      external G_indx


      nfft3d = n2ft3d/2
      value =           BA_push_get(mt_dbl,n2ft3d,
     >                                'hx',hx(2),hx(1))
      value = value.and.BA_push_get(mt_dbl,n2ft3d,
     >                                'hy',hy(2),hy(1))
      value = value.and.BA_push_get(mt_dbl,n2ft3d,
     >                                'hz',hz(2),hz(1))
      value = value.and.BA_push_get(mt_dbl,n2ft3d,
     >                               'fx',fx(2),fx(1))
      value = value.and.BA_push_get(mt_dbl,n2ft3d,
     >                                'fy',fy(2),fy(1))
      value = value.and.BA_push_get(mt_dbl,n2ft3d,
     >                                'fz',fz(2),fz(1))

      value = value.and.BA_push_get(mt_dbl,nfft3d,
     >                                'Gx',Gx(2),Gx(1))
      value = value.and.BA_push_get(mt_dbl,nfft3d,
     >                                'Gy',Gy(2),Gy(1))
      value = value.and.BA_push_get(mt_dbl,nfft3d,
     >                                'Gz',Gz(2),Gz(1))
      if (.not.value) 
     >  call errquit('gradient_energy_density:push stack',0, MA_ERR)

*     **** define Gx,Gy,Gz in packed space ****
      call dcopy(nfft3d,dbl_mb(G_indx(1)),1,dbl_mb(Gx(1)),1)
      call dcopy(nfft3d,dbl_mb(G_indx(2)),1,dbl_mb(Gy(1)),1)
      call dcopy(nfft3d,dbl_mb(G_indx(3)),1,dbl_mb(Gz(1)),1)
      call Pack_t_pack(nb,dbl_mb(Gx(1)))
      call Pack_t_pack(nb,dbl_mb(Gy(1)))
      call Pack_t_pack(nb,dbl_mb(Gz(1)))

*     ***** calculate:        ****
*     **** fx(G)=i*Gx(G)*f(G) ****
*     **** fy(G)=i*Gy(G)*f(G) ****
*     **** fz(G)=i*Gz(G)*f(G) ****
      call Pack_itc_Mul(nb,dbl_mb(Gx(1)),f,dbl_mb(fx(1)))
      call Pack_itc_Mul(nb,dbl_mb(Gy(1)),f,dbl_mb(fy(1)))
      call Pack_itc_Mul(nb,dbl_mb(Gz(1)),f,dbl_mb(fz(1)))
      call Pack_itc_Mul(nb,dbl_mb(Gx(1)),h,dbl_mb(hx(1)))
      call Pack_itc_Mul(nb,dbl_mb(Gy(1)),h,dbl_mb(hy(1)))
      call Pack_itc_Mul(nb,dbl_mb(Gz(1)),h,dbl_mb(hz(1)))

      value =           BA_pop_stack(Gz(2))
      value = value.and.BA_pop_stack(Gy(2))
      value = value.and.BA_pop_stack(Gx(2))
      if (.not.value) 
     >  call errquit('gradient_energy_density2:pop stack',1, MA_ERR)

*     **** Fourier transform fx,fy,fz to real space ****
      call Pack_c_unpack(nb,dbl_mb(fx(1)))
      call Pack_c_unpack(nb,dbl_mb(fy(1)))
      call Pack_c_unpack(nb,dbl_mb(fz(1)))
      call D3dB_cr_fft3b(1,dbl_mb(fx(1)))
      call D3dB_cr_fft3b(1,dbl_mb(fy(1)))
      call D3dB_cr_fft3b(1,dbl_mb(fz(1)))

      call Pack_c_unpack(nb,dbl_mb(hx(1)))
      call Pack_c_unpack(nb,dbl_mb(hy(1)))
      call Pack_c_unpack(nb,dbl_mb(hz(1)))
      call D3dB_cr_fft3b(1,dbl_mb(hx(1)))
      call D3dB_cr_fft3b(1,dbl_mb(hy(1)))
      call D3dB_cr_fft3b(1,dbl_mb(hz(1)))

*     **** calculate df = df + hx*fx + hy*fy + hz*fz ****
c      call D3dB_rr_Mul(1,dbl_mb(hx(1)),dbl_mb(fx(1)),dbl_mb(fx(1)))
c      call D3dB_rr_Mul(1,dbl_mb(hy(1)),dbl_mb(fy(1)),dbl_mb(fy(1)))
c      call D3dB_rr_Mul(1,dbl_mb(hz(1)),dbl_mb(fz(1)),dbl_mb(fz(1)))
c      call D3dB_rr_Sum(1,dbl_mb(fx(1)),df,df)
c      call D3dB_rr_Sum(1,dbl_mb(fy(1)),df,df)
c      call D3dB_rr_Sum(1,dbl_mb(fz(1)),df,df)
      call D3dB_rr_Mul2(1,dbl_mb(hx(1)),dbl_mb(fx(1)))
      call D3dB_rr_Mul2(1,dbl_mb(hy(1)),dbl_mb(fy(1)))
      call D3dB_rr_Mul2(1,dbl_mb(hz(1)),dbl_mb(fz(1)))
      call D3dB_rr_Sum2(1,dbl_mb(fx(1)),df)
      call D3dB_rr_Sum2(1,dbl_mb(fy(1)),df)
      call D3dB_rr_Sum2(1,dbl_mb(fz(1)),df)

      value =           BA_pop_stack(fz(2))
      value = value.and.BA_pop_stack(fy(2))
      value = value.and.BA_pop_stack(fx(2))
      value = value.and.BA_pop_stack(hz(2))
      value = value.and.BA_pop_stack(hy(2))
      value = value.and.BA_pop_stack(hx(2))
      if (.not.value) 
     >  call errquit('gradient_energy_density2:pop stack',2, MA_ERR)

      return
      end




*     **********************************
*     *                                *
*     *      generate_dmatrix_column   *
*     *                                *
*     **********************************

      subroutine generate_dmatrix_column(ms,xyz,
     >                                   ispin,ne,n2ft3d,psi_r,rho)
      implicit none
      integer    ms
      real*8     xyz(3)
      integer    ispin,ne(2),n2ft3d
      real*8     psi_r(n2ft3d,ne(1)+ne(2))
      real*8     rho(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer i,k,jj,jjindex(100),rgrid(2),nx,ny,nz
      real*8 summ(100),tsum,msum,x,y,z,dv,dx,dy,dz,rr,rrmin

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

      if (.not.BA_push_get(mt_dbl,3*n2ft3d,'rgrd',rgrid(2),rgrid(1)))
     > call errquit('pspw_dplot_loop:push stack',0,MA_ERR)
      call lattice_r_grid(dbl_mb(rgrid(1)))

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      dv = lattice_omega()/dble(nx*ny*nz)

c     *** find nearest i ***
      rrmin = 9.9d99
      jj = 0
      do i=1,n2ft3d
         x = dbl_mb(rgrid(1)+3*(i-1))
         y = dbl_mb(rgrid(1)+3*(i-1)+1)
         z = dbl_mb(rgrid(1)+3*(i-1)+2)
         dx = x-xyz(1)
         dy = y-xyz(2)
         dz = z-xyz(3)
         rr = dx**2 + dy**2 + dz**2
         if (rr.lt.rrmin) then
            rrmin = rr
            jj = i
         end if
      end do

c      jj = j
c      msum = 0.0d0
c      do i=1,n2ft3d
c         tsum = 0.0d0
c         do k=1,ne(ms)
c            tsum = tsum + psi_r(i,k+(ms-1)*ne(1))**2
c         end do
c         if (tsum.gt.msum) then
c            jj = i
c            msum = tsum
c         end if
c      end do
      write(*,*) "jj=",jj,rrmin
      x = dbl_mb(rgrid(1)+3*(jj-1))
      y = dbl_mb(rgrid(1)+3*(jj-1)+1)
      z = dbl_mb(rgrid(1)+3*(jj-1)+2)

      do k=1,ne(ms)
         summ(k) = psi_r(jj,k+(ms-1)*ne(1))
      end do
      call dcopy(n2ft3d,0.0d0,0,rho,1)
      do i=1,n2ft3d
         do k=1,ne(ms)
            rho(i) = rho(i) + psi_r(i,k+(ms-1)*ne(1))*summ(k)
         end do
      end do
      do i=1,n2ft3d
         tsum = tsum + rho(i)
      end do
      tsum = tsum * dv
      write(*,*) "<rho(:,j)> = ",jj,x,y,z,tsum

      if (.not.BA_pop_stack(rgrid(2)))
     > call errquit('pspw_dplot_loop:pop stack',0,MA_ERR)
      return
      end


