c $Id$

*    ************************************
*    *                                  *
*    *          nwpw_gga                *
*    *                                  *
*    ************************************
      subroutine nwpw_gga(gga,ic,ispin,
     >                   rho,agr,
     >                   xce,fn,fdn,tmp)
      implicit none
      integer gga,ic,ispin
      double precision rho(ic,ispin)
      double precision agr(ic,*)

      double precision xce(ic)
      double precision fn(ic,ispin)
      double precision fdn(ic,*)
      double precision tmp(ic)

#include "nwpwxc.fh"
#include "util.fh"

*     **** local variables ****
      double precision grad(ispin+1)
      double precision dtmp
      double precision dncut
      parameter(dncut = 1.0d-30)
      
      integer i
      logical use_nwpwxc

      
      use_nwpwxc = nwpwxc_is_on()

c     **** restricted calculations ****
      if (ispin.eq.1) then

!$OMP DO
         do i=1,ic
           tmp(i) = (rho(i,1)+rho(i,1))
         end do
!$OMP END DO

cccccc USE_NWXC needs to be checked cccccccc
         if (use_nwpwxc) then
           !call D3db_rr_Sqr1(1,agr)
!OMP DO
           do i=1,ic
              agr(i,1) = agr(i,1)**2
           end do
!OMP END DO
!$OMP MASTER
           call nwpwxc_eval_df(1,ic,tmp,agr,
     >                       dtmp,xce,
     >                       fn,fdn,dtmp)
!$OMP END MASTER
!$OMP BARRIER
c
c          Combine (df/d|grad a|) with (df/d(grad a|grad b))
c
           !call D3dB_r_SMul(1,0.5d0,fdn(1,2),fdn)
c
c          Calculate energy density from energy
c
!$OMP DO
           do i = 1, ic
              fdn(i,1) = fdn(i,1) + 0.5d0*fdn(i,2)
              xce(i) = xce(i)/(tmp(i)+dncut)
              agr(i,1) = dsqrt(agr(i,1))
              fdn(i,1) = fdn(i,1)*agr(i,1)
           enddo
!$OMP END DO
           !call D3dB_rr_Sqrt1(1,agr)
           !call D3dB_rr_Mul2(1,agr,fdn)
cccccc USE_NWXC needs to be checked cccccccc

         else
         if (gga.eq.10) then
         call gen_PBE96_BW_restricted(ic,
     >                                tmp,
     >                                agr,
     >                                1.0d0,1.0d0,
     >                                xce,
     >                                fn,
     >                                fdn)
         else if (gga.eq.11) then
         call gen_BLYP_BW_restricted(ic,
     >                                tmp,
     >                                agr,
     >                                1.0d0,1.0d0,
     >                                xce,
     >                                fn,
     >                                fdn)

         else if (gga.eq.12) then
         call gen_revPBE_BW_restricted(ic,
     >                                tmp,
     >                                agr,
     >                                1.0d0,1.0d0,
     >                                xce,
     >                                fn,
     >                                fdn)
         else
          call errquit('bad gga',0,0)
         end if
         endif ! nwpwxc_is_on


c     **** unrestricted calculations ****
      else if (ispin.eq.2) then

         if (use_nwpwxc) then
!$OMP DO
           do i = 1, ic
             grad(1) = agr(i,1)**2
             grad(2) = agr(i,3)**2
             grad(3) = agr(i,2)**2
             grad(2) = 0.5d0*(grad(2) - grad(1) - grad(3))
             agr(i,1) = grad(1)
             agr(i,2) = grad(2)
             agr(i,3) = grad(3)
           enddo
!$OMP END DO
!$OMP MASTER
           call nwpwxc_eval_df(2,ic,rho,agr,
     >                       dtmp,xce,
     >                       fn,fdn,dtmp)
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO
           do i = 1, ic
             xce(i)   = xce(i)/(rho(i,1)+rho(i,2)+dncut)
             dtmp     = fdn(i,2)
             fdn(i,2) = fdn(i,3)
             fdn(i,3) = dtmp
             grad(1)  = dsqrt(agr(i,1))
             grad(2)  = dsqrt(agr(i,3))
             grad(3)  = dsqrt(2.0d0*agr(i,2)+agr(i,1)+agr(i,3))
             agr(i,1) = grad(1)
             agr(i,2) = grad(2)
             agr(i,3) = grad(3)
             fdn(i,1) = (2.0d0*fdn(i,1)-fdn(i,3))*agr(i,1)
             fdn(i,2) = (2.0d0*fdn(i,2)-fdn(i,3))*agr(i,2)
             fdn(i,3) = fdn(i,3)*agr(i,3)  
           enddo
!$OMP END DO
         else
         if (gga.eq.10) then
         call gen_PBE96_BW_unrestricted(ic,
     >                                rho,
     >                                agr,
     >                                1.0d0,1.0d0,
     >                                xce,
     >                                fn,
     >                                fdn)
         else if (gga.eq.11) then
         call gen_BLYP_BW_unrestricted(ic,
     >                                rho,
     >                                agr,
     >                                1.0d0,1.0d0,
     >                                xce,
     >                                fn,
     >                                fdn)

         else if (gga.eq.12) then
         call gen_revPBE_BW_unrestricted(ic,
     >                                rho,
     >                                agr,
     >                                1.0d0,1.0d0,
     >                                xce,
     >                                fn,
     >                                fdn)
         else
          call errquit('bad gga',0,0)
         end if
         end if ! nwpwxc_is_on

      else
        call errquit('bad spin',0,ispin)
      end if
      return
      end
