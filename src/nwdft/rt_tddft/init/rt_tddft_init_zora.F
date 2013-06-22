C
C     rt_tddft_init_zora.F
C
C     Initialization for ZORA, one routine for scalar and one for spin orbit.
C
      subroutine rt_tddft_init_zora_scalar (params)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "cdft.fh"
#include "zora.fh"
#include "rt_tddft.fh"


C     == In/out ==
      type(rt_params_t) params  !zora params/handles to be stored in here

      
C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_init_zora: "


C     == Variables ==
      integer g_zora_sf(2)
      integer g_zora_scale_sf(2)
      character*256 fname
      integer me
      integer g_ztmp
      integer nbf_array(2)
      integer i


C     == External ==
      logical, external ::  dft_zora_read

      
      me = ga_nodeid()


C
C     Checks
C     
      if (.not. do_zora)
     $     call errquit (pname//"called without ZORA active",0,0)

      if ((ipol.ne.1).and.(ipol.ne.2))
     $    call errquit (pname//"invalid ipol", ipol, 0) 


      params%lzorasf = .true.
      params%lzoraso = .false.
      

      
c$$$      if (.not. ga_create (mt_dbl, params%nbf_ao, params%nbf_ao,
c$$$     $     "SF", 0, 0, g_zora_sf))
c$$$     $     call errquit (pname//"create failed",0,GA_ERR)
c$$$
c$$$      if (.not. ga_create (mt_dbl, params%nbf_ao, params%nbf_ao,
c$$$     $     "scale SF", 0, 0, g_zora_scale_sf))
c$$$     $     call errquit (pname//"create failed",0,GA_ERR)
c$$$
c$$$      call ga_zero (g_zora_sf)
c$$$      call ga_zero (g_zora_scale_sf)
      
C
C     Allocation
C
      do i = 1, ipol
         if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $        "sf", 0, 0, g_zora_sf(i)))
     $        call errquit (pname//"couldnt make g_zora_sf", 0, GA_ERR)

         if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $        "scale sf", 0, 0, g_zora_scale_sf(i)))
     $        call errquit (pname//"couldnt make g_zora_sf", 0, GA_ERR)

         call ga_zero (g_zora_sf(i))
         call ga_zero (g_zora_scale_sf(i))

         nbf_array(i) = params%nbf_ao
      
      enddo


C
C     Get ZORA data from file.  If we can't read it, quit since this is
C     safer than recoding the reconstruction here (in case it chages).
C
      call util_file_name ("zora_sf", .false., .false., fname)

      if (.not.dft_zora_read(fname, params%nbf_ao, ipol,
     $     nbf_array, mult, g_zora_sf, g_zora_scale_sf))
     $     call errquit (pname//
     $     "failed to read scalar ZORA data from: "//trim(fname),0,0)


c$$$      if (ipol .eq. 1) then
c$$$         if (.not.dft_zora_read(fname, params%nbf_ao, 1,
c$$$     $        nbf_array, mult, g_zora_sf, g_zora_scale_sf))
c$$$     $        call errquit (pname//
c$$$     $        "failed to read scalar ZORA data from: "//trim(fname),0,0)
c$$$      elseif (ipol .eq. 2) then
c$$$         call errquit (pname//"ipol = 2 not implemented", ipol, 0)
c$$$      else
c$$$         call errquit (pname//"invalid ipol", ipol, 0)
c$$$      endif

      do i = 1, ipol
         params%g_zora_sf(i) = g_zora_sf(i)
         params%g_zora_scale_sf(i) = g_zora_scale_sf(i)
      enddo


C
C     Check that matrices we just read in are consistent with those in
C     the ground state code.
C
      if (params%matrix_checks) then
         if (.not. ga_create(mt_dcpl, params%nbf_ao, params%nbf_ao,
     $        "ztmp", 0, 0, g_ztmp))
     $        call errquit (pname//"couldnt make g_ztmp", 0, GA_ERR)
         
         if (ipol.eq.1) then
            call ga_zero (g_ztmp)
            call convert_d2z (1d0, g_zora_sf(1), 0d0, g_zora_sf(1),
     $           g_ztmp)
            call zmat_compare_and_print (params%nbf_ao, g_ztmp,
     $           "zora_cs_sf", params%tol_zero,
     $           "GS and RT-TDDFT CS ZORA matrices")
         elseif (ipol.eq.2) then
            call ga_zero (g_ztmp)
            call convert_d2z (1d0, g_zora_sf(1), 0d0, g_zora_sf(1),
     $           g_ztmp)
            call zmat_compare_and_print (params%nbf_ao, g_ztmp,
     $           "zora_alpha_sf", params%tol_zero,
     $           "GS and RT-TDDFT alpha ZORA matrices")
            
            call ga_zero (g_ztmp)
            call convert_d2z (1d0, g_zora_sf(2), 0d0, g_zora_sf(2),
     $           g_ztmp)
            call zmat_compare_and_print (params%nbf_ao, g_ztmp,
     $           "zora_beta_sf", params%tol_zero,
     $           "GS and RT-TDDFT beta ZORA matrices")
         else
            call errquit (pname//"invalid ipol", ipol, 0)
         endif
         
         if (.not. ga_destroy (g_ztmp))
     $        call errquit (pname//"destroy failed", 0,0)
      endif


C     XXX CLEAN UP ZORA ARRAYS AFTER PROPAGATION
      end subroutine

c $Id$
