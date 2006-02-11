
*    ************************************
*    *                                  *
*    *          paw_gga                 *
*    *                                  *
*    ************************************
      subroutine paw_gga(gga,ic,ispin,
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
      
*     **** local variables ****
      integer i

c     **** restricted calculations ****
      if (ispin.eq.1) then

         do i=1,ic
           tmp(i) = (rho(i,1)+rho(i,1))
         end do
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


c     **** unrestricted calculations ****
      else if (ispin.eq.2) then

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

      else
        call errquit('bad spin',0,ispin)
      end if
      return
      end
