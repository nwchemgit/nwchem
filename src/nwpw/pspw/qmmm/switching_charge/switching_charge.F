
*     *************************************
*     *           			              *
*     *   Switching_Charge_Generate_V     *
*     *                         			*
*     *************************************
      subroutine Switching_Charge_Generate_V()
      implicit none

*     **** external functions ****
      integer  pspw_qmmm_nion
      real*8   pspw_qmmm_rion
      external pspw_qmmm_nion
      external pspw_qmmm_rion

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

*     **********************************
*     *                                *
*     *     pspw_charge_rho_Fcharge        *
*     *                                *
*     **********************************

      subroutine pspw_charge_rho_Fcharge(n2ft3d,rgrid,rho,dv,fcharge)
      implicit none
      integer n2ft3d
      real*8 rgrid(3,*)
      real*8 rho(*)
      real*8 dv
      real*8 fcharge(3,*)

#include "mafdecls.fh"
#include "pspw_charge.fh"

*     ***** local variables ****
      integer i,j
      real*8 rx,ry,rz,c
      real*8 fx,fy,fz,verf,yerf,v
      real*8 x,y,z,q,r,epsilon,depsilon,sqrt_pi

*     **** external functions ****
      real*8   util_erf
      external util_erf


      if (nion.gt.0) then
      sqrt_pi = dsqrt(4.0d0*datan(1.0d0))

      do j=1,nion
         x = dbl_mb(r1(1)+3*(j-1))
         y = dbl_mb(r1(1)+3*(j-1)+1)
         z = dbl_mb(r1(1)+3*(j-1)+2)
         q = -dbl_mb(charge(1)+j-1)
         c = 1.0d0/s_sigma
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
     >             + depsilon*q*verf/r
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
         fcharge(1,j) = fcharge(1,j) - fx*dv
         fcharge(2,j) = fcharge(2,j) - fy*dv
         fcharge(3,j) = fcharge(3,j) - fz*dv
      end do


      end if

      return
      end



