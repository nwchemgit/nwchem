c
c $Id: switching_charge.F,v 1.5 2002-03-03 17:42:42 bylaska Exp $
c

***********************************************************
*                                                         *
*               Switching Charge module                   *
*                                                         *
*          Interfaced to nwchem-PSPW code                 *
*                                                         *
*   -- developed by Eric J. Bylaska on February 27, 2002  *
*                                                         *
***********************************************************
*
*
*

*     **********************************
*     *                                *
*     *     Switching_Charge_init      *
*     *                                *
*     **********************************

      subroutine Switching_Charge_init(rtdb)
      implicit none
      integer rtdb

#include "mafdecls.fh"
#include "rtdb.fh"
#include "switching_charge.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      integer npack0,nfft3d,G(3)
      integer i,j,k
      integer zero,qzero,pzero,taskid
      integer nx,ny,nxh,nyh
      real*8  gg
      logical value
      integer tmp1(2)


*     **** external functions ****
      integer  G_indx
      real*8   Truncating_Function_transform
      external G_indx
      external Truncating_Function_transform


      call Parallel_taskid(taskid)

      if (.not.rtdb_get(rtdb,
     >      'nwpw:charge_dielectric_d',mt_dbl,1,s_d))
     >    s_d = 0.5d0/0.529177d0

      if (.not.rtdb_get(rtdb,
     >      'nwpw:charge_dielectric_rho',mt_dbl,1,s_rho))
     >    s_rho = 2.5d0/0.529177d0

      if (.not.rtdb_get(rtdb,
     >      'nwpw:charge_rcut',mt_dbl,1,s_sigma))
     >    s_sigma = 0.8d0/0.529177d0

      if (.not.rtdb_get(rtdb,'nwpw:truncating_rcut',mt_dbl,1,rcut))
     >   rcut = 2.0d0/0.529177d0

      if (.not.rtdb_get(rtdb,'nwpw:truncating_power',mt_dbl,1,pp))
     >   pp = 8.0d0


      if (taskid.eq.MASTER) then
          write(*,1160)
          write(*,1161) rcut
          write(*,1162 pp
          write(*,1170)
          write(*,1171) s_d,s_rho,s_sigma
      end if
 1160 FORMAT(/' Truncating Function definition:')
 1161 FORMAT(5X, '- Truncating rcut  = (', E11,3,' )  ')
 1162 FORMAT(5X, '- Truncating power = (', E11,3,' )  ')
 1170 FORMAT(/' Switching charges used:')
 1171 FORMAT(5X, '- (d,rho,sigma)    = (',3E11.3,' )  ')


*     **** define Truncating Function ****      
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(0,npack0)
      G(1) = G_indx(1)
      G(2) = G_indx(2)
      G(3) = G_indx(3)

      value = MA_alloc_get(mt_dbl,(npack0),'Trnc',Trnc(2),Trnc(1))
      if (.not. value) 
     >  call errquit('Switching_Charge_init: out of heap memory',0)
      call dcopy(npack0,0.0d0,0,dbl_mb(Trnc(1)), 1)


      value = MA_push_get(mt_dbl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     >   call errquit('Switching_Charge_init: pushing stack',0)


*     ***** Generate Fourier transform of Trunctaing Function *****
      do i = 1,nfft3d
         gg  = ( dbl_mb(G(1)+i-1)*dbl_mb(G(1)+i-1)
     >         + dbl_mb(G(2)+i-1)*dbl_mb(G(2)+i-1)
     >         + dbl_mb(G(3)+i-1)*dbl_mb(G(3)+i-1) )
         dbl_mb(tmp1(1)+i-1) = 
     >       Truncating_Function_transform(gg,rcut,pp)
      end do
      call Pack_t_pack(0,dbl_mb(tmp1(1)))
      call Pack_t_Copy(0,dbl_mb(tmp1(1)),dbl_mb(Trnc(1)))
      value = MA_pop_stack(tmp1(2))

      return
      end

*     *************************************
*     *                                   *
*     *   Truncating_Function_transform   *
*     *                                   *
*     *************************************

*     This function returns the fourier transform of
*
*            T(r) = (1.0d0-(1.0d0-dexp(-(r/rcut)**(pp+2)))**pp)
*
*      May want to use:
*            T(r) = exp(-(r/rcut)**pp)
*
*     Entry - gg: g squared
*             rcut:
*             pp:
*
*     Exit - returns
*                      /infty
*                     |
*      T(g) = (4*pi)* | r**2 * T(r)* j0(gr) dr
*                     |
*                     / 0

      real*8 function Truncating_Function_transform(gg,rcut,pp)
      implicit none
      real*8 gg,rcut,pp

*     **** local variables ****
      integer nrho
      parameter (nrho=5000)

      integer i
      real*8 q,fourpi,r,sum,drho,pp2,transform

      pp2    = pp + 2.0d0
      drho   = 2.0d0*rcut/dble(nrho)
      fourpi = 4.0d0*(4.0d0*datan(1.0d0))

*     **** gg > 0 ****
      if (gg.gt.1.0d-9) then
         q  = dsqrt(gg)
         sum = 0.0d0
         do i=2,nrho-1
            r = (i-1)*drho
            sum = sum + r*sin(q*r)
     >           *(1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
         end do
         r = drho*(nrho-1)
         sum = sum + r*sin(q*r)
     >        *(1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
         transform = (fourpi/q)*sum*drho

*     **** gg == 0 *****
      else
         sum = 0.0d0
         do i=2,nrho-1
            r = (i-1)*drho
            sum = sum + r*r*
     >           *(1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
         end do
         r = drho*(nrho-1)
         sum = sum + r*r*
     >        *(1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
         transform = (fourpi)*sum*drho
      end if

      Truncating_Function_transform = transform
      return
      end


*     **********************************
*     *                                *
*     *     Switching_Charge_end       *
*     *                                *
*     **********************************

      subroutine Switching_Charge_end()
      implicit none

#include "mafdecls.fh"
#include "switching_charge.fh"

      if (.not. MA_free_heap(Trnc(2))) 
     >  call errquit('Switching_Charge_end: freeing heap',0)
      return
      end


*     *************************************
*     *           	                  *
*     *   Switching_Charge_Generate_V     *
*     *                                   *
*     *************************************
      subroutine Switching_Charge_Generate_V(n2ft3d,rgrid,Vqm)
      implicit none
      integer n2ft3d
      real*8  rgrid(3,*)
      real*8  Vqm(*)

#include "switching_charge.fh"

*     **** local variables ****
      integer ii,k,nion
      real*8 x1,y1,z1,q1,r,epsilon

*     **** external functions ****
      integer  pspw_qmmm_nion
      real*8   pspw_qmmm_rion,util_erf
      external pspw_qmmm_nion
      external pspw_qmmm_rion,util_erf

      nion = pspw_qmmm_nion()

*     **** switching charge potential ****
      do ii=1,nion
        x1 = pspw_qmmm_rion(1,ii)
        y1 = pspw_qmmm_rion(2,ii)
        z1 = pspw_qmmm_rion(3,ii)
        call pspw_Q_param_MM(ii,q1)

        do k=1,n2ft3d
          r = (rgrid(1,k)-x1)**2
     >      + (rgrid(2,k)-y1)**2
     >      + (rgrid(3,k)-z1)**2
          r = dsqrt(r)

*         **** define dielectric switching function ****
          if (r.le.s_d) then
             epsilon = 0.0d0
          else if (r.lt.(s_d+s_rho)) then
             epsilon = 1.0d0-(1.0d0-(r-s_d)**2/s_rho**2)**2
          else
             epsilon = 1.0d0
          end if

          Vqm(k) = Vqm(k) - epsilon*q1*util_erf(r/s_sigma)/r
        end do
      end do

      return
      end


*     *************************************
*     *                                   *
*     *  Switching_Charge_Potential_Fmm   *
*     *                                   *
*     *************************************

      subroutine Switching_Charge_Potential_Fmm(n2ft3d,rgrid,rho,dv,fmm)
      implicit none
      integer n2ft3d
      real*8 rgrid(3,*)
      real*8 rho(*)
      real*8 dv
      real*8 fmm(3,*)

#include "mafdecls.fh"
#include "switching_charge.fh"

*     ***** local variables ****
      integer i,j,nion
      real*8 rx,ry,rz,c
      real*8 fx,fy,fz,verf,yerf,v
      real*8 x,y,z,q,r,epsilon,depsilon,sqrt_pi

*     **** external functions ****
      integer  pspw_qmmm_nion
      real*8   pspw_qmmm_rion,util_erf
      external pspw_qmmm_nion
      external pspw_qmmm_rion,util_erf

      nion = pspw_qmmm_nion()

      sqrt_pi = dsqrt(4.0d0*datan(1.0d0))
      c = 1.0d0/s_sigma

      do j=1,nion
         x = pspw_qmmm_rion(1,j)
         y = pspw_qmmm_rion(2,j)
         z = pspw_qmmm_rion(3,j)
         call pspw_Q_param_MM(j,q)
         q = -q
         fx = 0.0d0
         fy = 0.0d0
         fz = 0.0d0
         do i=1,n2ft3d
            rx = x - rgrid(1,i)
            ry = y - rgrid(2,i)
            rz = z - rgrid(3,i)
            r  = dsqrt( rx**2 + ry**2 + rz**2)

            if (r .gt. 1.0d-8) then

*             **** define dielectric switching function ****
              if (r.le.s_d) then
                 epsilon  = 0.0d0
                 depsilon = 0.0d0
              else if (r.lt.(s_d+s_rho)) then
                 epsilon = 1.0d0-(1.0d0-(r-s_d)**2/s_rho**2)**2
                 depsilon = 4.0d0*((r-s_d)/s_rho**2)
     >                           *(1.0d0-(r-s_d)**2/s_rho**2)
              else
                 epsilon  = 1.0d0
                 depsilon = 0.0d0
              end if
              yerf=r*c
              verf = util_erf(yerf)
              v    = epsilon*q*( (2.0d0/sqrt_pi)*(r*c)*exp(-(r*c)**2)
     >                          - verf)/r**3
     >             + depsilon*q*verf/r**2
            else
              v = 0.0d0
            end if

            fx = fx + rho(i)*rx*v
            fy = fy + rho(i)*ry*v
            fz = fz + rho(i)*rz*v
         end do

         call D3dB_SumAll(fx)
         call D3dB_SumAll(fy)
         call D3dB_SumAll(fz)
         fmm(1,j) = fmm(1,j) - fx*dv
         fmm(2,j) = fmm(2,j) - fy*dv
         fmm(3,j) = fmm(3,j) - fz*dv
      end do


      return
      end


*     *************************************
*     *           	                  *
*     *   Generate_Truncating_Function    *
*     *                                   *
*     *************************************
      subroutine Generate_Truncating_Function(T)

C
      data (vander(i),i=1,36)
     1                  /1.20D+00,1.20D+00,1.37D+00,1.45D+00,
     2 1.45D+00,1.50D+00,1.50D+00,1.40D+00,1.35D+00,1.30D+00,
     3 1.57D+00,1.36D+00,1.24D+00,1.17D+00,1.80D+00,1.75D+00,
     4 1.70D+00,19*0.0D+00/

