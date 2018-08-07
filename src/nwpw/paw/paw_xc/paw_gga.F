
*    ************************************
*    *                                  *
*    *          paw_gga                 *
*    *                                  *
*    ************************************
      subroutine paw_gga(gga,ic,ispin,
     >                   rho,agr,
     >                   xce,fn,fdn,tmp)
      implicit none
#include "nwpwxc.fh"
#include "util.fh"
      integer gga,ic,ispin
      logical use_nwpwxc
      double precision rho(ic,ispin)
      double precision agr(ic,*)
      double precision grad(ispin+1)

      double precision xce(ic)
      double precision fn(ic,ispin)
      double precision fdn(ic,*)
      double precision tmp(ic)

      double precision dtmp, dncut, dumtau
      parameter(dncut = 1.0d-30)
      
*     **** local variables ****
      integer i

c     **** restricted calculations ****
      if (ispin.eq.1) then

         do i=1,ic
           tmp(i) = (rho(i,1)+rho(i,1))
         end do
         use_nwpwxc = .false.
         use_nwpwxc = nwpwxc_is_on()
        
         if (use_nwpwxc) then
           call D3db_rr_Sqr1(1,agr)
           call nwpwxc_eval_df(1,ic,tmp,agr,
     >                       dumtau,xce,
     >                       fn,fdn,dumtau)
c
c          Combine (df/d|grad a|) with (df/d(grad a|grad b))
c
           call D3dB_rr_daxpy(1,0.5d0,fdn(1,2),fdn)
c
c          Calculate energy density from energy
c
           do i = 1, ic
             xce(i) = xce(i)/(tmp(i)+dncut)
           enddo
           call D3dB_rr_Sqrt1(1,agr)
           call D3dB_rr_Mul2(1,agr,fdn)
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

         use_nwpwxc = .false.
         use_nwpwxc = nwpwxc_is_on()
        
         if (use_nwpwxc) then
           do i = 1, ic
             grad(1) = agr(i,1)**2
             grad(2) = agr(i,3)**2
             grad(3) = agr(i,2)**2
             grad(2) = 0.5d0*(grad(2) - grad(1) - grad(3))
             agr(i,1) = grad(1)
             agr(i,2) = grad(2)
             agr(i,3) = grad(3)
           enddo
           call nwpwxc_eval_df(2,ic,rho,grad,
     >                       dumtau,xce,
     >                       fn,fdn,dumtau)
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
         endif ! nwpwxc_is_on

      else
        call errquit('bad spin',0,ispin)
      end if
      return
      end
c $Id$
