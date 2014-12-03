c
c $Id$
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

*     **********************************
*     *                                *
*     *     Switching_Charge_init      *
*     *                                *
*     **********************************

      subroutine Switching_Charge_init(rtdb)
      implicit none
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"
#include "stdio.fh"
#include "switching_charge.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)


      call Parallel_taskid(taskid)

      if (.not.btdb_get(rtdb,
     >      'nwpw:charge_rcut',mt_dbl,1,s_sigma))
     >    s_sigma = 0.8d0/0.529177d0


      if (taskid.eq.MASTER) then
          write(luout,1170)
          write(luout,1171) s_sigma
      end if
 1170 FORMAT(/' Electrostatic point charges used:')
 1171 FORMAT(5X, '- sigma = (',E11.3,' )  ')


      call Truncating_Function_init(rtdb)

      return
      end


*     **********************************
*     *                                *
*     *     Switching_Charge_end       *
*     *                                *
*     **********************************

      subroutine Switching_Charge_end()
      implicit none

      call Truncating_Function_end()
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

#include "errquit.fh"
#include "stdio.fh"
#include "bafdecls.fh"
#include "switching_charge.fh"

*     **** local variables ****
      logical value
      real*8  scal2
      integer xtmp(2),i,count

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

*     scal2   = 1.0d0/dsqrt(lattice_omega())
      scal2   = 1.0d0

*     **** switching charge potential ****
      call Switching_Charge_Vmm(n2ft3d,rgrid,Vqm)

*     *** generate Tall(r) ****
      value = BA_push_get(mt_dbl,n2ft3d,'xtmp',xtmp(2),xtmp(1))
      if (.not. value) 
     >  call errquit(' Switching_Charge_Generate_V:pushing stack',0,
     &       MA_ERR)
      call Truncating_Function(rgrid,dbl_mb(xtmp(1)))

c*     **** debug ****
c       call dcopy(n2ft3d,1.0d0,0,dbl_mb(xtmp(1)),1)
       scal2 = 0.0d0
       count = 0
       do i=1,n2ft3d
           
          if ( dabs(dbl_mb(xtmp(1)+i-1)) .gt.0.1d0) then
          write(luout,*) "T:",i,dbl_mb(xtmp(1)+i-1)
          end if
          if ( dbl_mb(xtmp(1)+i-1) .lt.0.0d0) then
             scal2 = scal2 + dbl_mb(xtmp(1)+i-1)
             count = count + 1
          end if
       end do
       write(luout,*) "ERROR:",count,scal2

*     **** Tall*Vqm ****
c      call D3dB_rr_Mul(1,Vqm,dbl_mb(xtmp(1)),Vqm)
      call D3dB_rr_Mul2(1,dbl_mb(xtmp(1)),Vqm)

      value = BA_pop_stack(xtmp(2))
      if (.not. value)
     >  call errquit(' Switching_Charge_Generate_V:popping stack',1,
     &       MA_ERR)

      return
      end


*     *************************************
*     *                                   *
*     *  Switching_Charge_Potential_Fion  *
*     *                                   *
*     *************************************

      subroutine Switching_Charge_Potential_Fion(
     >           n2ft3d,rgrid,rho,dv,fion)
      implicit none
#include "errquit.fh"
      integer n2ft3d
      real*8 rgrid(3,*)
      real*8 rho(*)
      real*8 dv
      real*8 fion(3,*)

#include "bafdecls.fh"
#include "switching_charge.fh"

*     ***** local variables *****
      logical value
      integer nx,ny,nz,Amm(2)
      real*8  scal1

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)

      value = BA_push_get(mt_dbl,n2ft3d,'Amm',Amm(2),Amm(1))
      if (.not. value)
     >  call errquit(' Switching_Charge_Generate_Fion:pushing stack',0,
     &       MA_ERR)

*     **** switching charge potential ****
      call D3dB_r_Zero(1,dbl_mb(Amm(1)))
      call Switching_Charge_Vmm(n2ft3d,rgrid,dbl_mb(Amm(1)))

*     **** Amm(r) = rho(r)*Vmm(r) ****
c      call D3dB_rr_Mul(1,dbl_mb(Amm(1)),rho,dbl_mb(Amm(1)))
      call D3dB_rr_Mul2(1,rho,dbl_mb(Amm(1)))


*     **** Amm(G) = FFT(Amm(r)) ****
      call D3dB_rc_fft3f(1,dbl_mb(Amm(1)))
      call Pack_c_pack(0,dbl_mb(Amm(1)))
c      call Pack_c_SMul(0,scal1,dbl_mb(Amm(1)),dbl_mb(Amm(1)))
      call Pack_c_SMul1(0,scal1,dbl_mb(Amm(1)))
      

      call grad_Truncating_Function(dbl_mb(Amm(1)),fion)

      value = BA_pop_stack(Amm(2))
      if (.not. value)
     >  call errquit(' Switching_Charge_Generate_Fion:popping stack',1,
     &       MA_ERR)

      return
      end


*     *************************************
*     *                                   *
*     *  Switching_Charge_Potential_Fmm   *
*     *                                   *
*     *************************************

      subroutine Switching_Charge_Potential_Fmm(n2ft3d,rgrid,rho,dv,fmm)
      implicit none
#include "errquit.fh"
      integer n2ft3d
      real*8 rgrid(3,*)
      real*8 rho(*)
      real*8 dv
      real*8 fmm(3,*)

#include "bafdecls.fh"
#include "switching_charge.fh"

*     ***** local variables ****
      logical value
      integer i,j,nion
      real*8 rx,ry,rz,c
      real*8 fx,fy,fz,verf,yerf,v
      real*8 x,y,z,q,r,sqrt_pi,scal2
      integer xtmp(2)

*     **** external functions ****
      integer  pspw_qmmm_nion
      real*8   pspw_qmmm_rion,util_erf,lattice_omega
      external pspw_qmmm_nion
      external pspw_qmmm_rion,util_erf,lattice_omega

      nion = pspw_qmmm_nion()

      sqrt_pi = dsqrt(4.0d0*datan(1.0d0))
      c       = 1.0d0/s_sigma
      scal2   = 1.0d0/lattice_omega()

*     *** generate Tall ****
      value = BA_push_get(mt_dbl,n2ft3d,'xtmp',xtmp(2),xtmp(1))
      if (.not. value) 
     >  call errquit(' Switching_Charge_Potentail_Fmm:pushing stack',0,
     &       MA_ERR)
      call Truncating_Function(rgrid,dbl_mb(xtmp(1)))

c*     **** debug ****
c      call dcopy(n2ft3d,1.0d0,0,dbl_mb(xtmp(1)),1)

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

              yerf=r*c
              verf = util_erf(yerf)
              v    = q*( (2.0d0/sqrt_pi)*(r*c)*exp(-(r*c)**2)
     >                          - verf)/r**3
            else
              v = 0.0d0
            end if

            fx = fx + dbl_mb(xtmp(1)+i-1)*rho(i)*rx*v
            fy = fy + dbl_mb(xtmp(1)+i-1)*rho(i)*ry*v
            fz = fz + dbl_mb(xtmp(1)+i-1)*rho(i)*rz*v
         end do

         call D3dB_SumAll(fx)
         call D3dB_SumAll(fy)
         call D3dB_SumAll(fz)
         fmm(1,j) = fmm(1,j) - fx*dv
         fmm(2,j) = fmm(2,j) - fy*dv
         fmm(3,j) = fmm(3,j) - fz*dv
      end do

      value = BA_pop_stack(xtmp(2))
      if (.not. value)
     >  call errquit(' Switching_Charge_Potentail_Fmm:popping stack',1,
     &       MA_ERR)

      return
      end

*     *************************************
*     *           	                      *
*     *   Switching_Charge_Vmm            *
*     *                                   *
*     *************************************
      subroutine Switching_Charge_Vmm(n2ft3d,rgrid,Vmm)
      implicit none
      integer n2ft3d
      real*8  rgrid(3,*)
      real*8  Vmm(*)

#include "switching_charge.fh"

*     **** local variables ****
      integer ii,k,nion
      real*8 x1,y1,z1,q1,r,c,sqrt_pi

*     **** external functions ****
      integer  pspw_qmmm_nion
      real*8   pspw_qmmm_rion,util_erf
      external pspw_qmmm_nion
      external pspw_qmmm_rion,util_erf

      nion    = pspw_qmmm_nion()
      sqrt_pi = dsqrt(4.0d0*datan(1.0d0))
      c       = 1.0d0/s_sigma

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

          if (r .gt. 1.0d-15) then
            Vmm(k) = Vmm(k) - q1*util_erf(c*r)/r
          else
            Vmm(k) = Vmm(k) - 2.0d0*q1*c/sqrt_pi
          end if

        end do
      end do

      return
      end


