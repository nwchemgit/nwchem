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
      integer, parameter   :: npol = 1
      character(*), parameter :: pname = "rt_tddft_init_zora: "


C     == Variables ==
      integer g_zora_sf
      integer g_zora_scale_sf
      character*256 fname
      integer me
      integer g_ztmp
      integer nbf_array


C     == External ==
      logical, external ::  dft_zora_read
      logical, external  :: zmat_compare

      
      me = ga_nodeid()

      if (me.eq.0) write (luout, *) pname//"Before alloc"

      if (.not. do_zora)
     $     call errquit (pname//"called without ZORA active",0,0)

      params%lzorasf = .true.
      params%lzoraso = .false.

      
      if (.not. ga_create (mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "SF", 0, 0, g_zora_sf))
     $     call errquit (pname//"create failed",0,GA_ERR)

      if (.not. ga_create (mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "scale SF", 0, 0, g_zora_scale_sf))
     $     call errquit (pname//"create failed",0,GA_ERR)

      call ga_zero (g_zora_sf)
      call ga_zero (g_zora_scale_sf)
      call halt ()
      
C
C     Allocation
C
      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "sf(1)", 0, 0, g_zora_sf))
     $     call errquit (pname//"couldnt make g_zora_sf(1)", 0, GA_ERR)

      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "scale sf(1)", 0, 0, g_zora_scale_sf))
     $     call errquit (pname//"couldnt make g_zora_sf(1)", 0, GA_ERR)


      if (me.eq.0) write (luout, *) pname//"Before file_name"

C
C     Get ZORA data from file.  If we can't read it, quit since this is
C     safer than recoding the reconstruction here (in case it chages).
C
      call util_file_name ("zora_sf", .false., .false., fname)

      if (me.eq.0) write (luout, *) pname//"Before zero"

      call ga_zero (g_zora_sf)
      call ga_zero (g_zora_scale_sf)

      if (me.eq.0) write (luout, *) pname//" XXXX DID NOT ZERO"

      if (me.eq.0) write (luout, *) pname//"Before params"

      nbf_array = params%nbf_ao

      if (me.eq.0) write (luout, *) pname//"Before read"

      if (ipol .eq. 1) then
         if (.not.dft_zora_read(fname, params%nbf_ao, 1,
     $        nbf_array, mult, g_zora_sf, g_zora_scale_sf))
     $        call errquit (pname//
     $        "failed to read scalar ZORA data from: "//trim(fname),0,0)
      elseif (ipol .eq. 2) then
         call errquit (pname//"ipol = 2 not implemented", ipol, 0)
      else
         call errquit (pname//"invalid ipol", ipol, 0)
      endif

      params%g_zora_sf(1) = g_zora_sf
      params%g_zora_scale_sf(1) = g_zora_scale_sf


      if (me.eq.0) write (luout, *) pname//"After read"

C
C     Check that matrices we just read in are consistent with those in
C     the ground state code.
C
      if (.not. ga_create(mt_dcpl, params%nbf_ao, params%nbf_ao,
     $     "ztmp", 0, 0, g_ztmp))
     $     call errquit (pname//"couldnt make g_ztmp", 0, GA_ERR)
      
      if (ipol.eq.1) then
         call convert_d2z (1d0, g_zora_sf, 0d0, g_zora_sf, g_ztmp)
         call zmat_compare_and_print (params%nbf_ao, g_ztmp,
     $        "g_zora_cs_sf", params%tol_zero,
     $        "GS and RT-TDDFT CS ZORA matrices")
      else
         call errquit (pname//"ipol = 2 not implemented", ipol, 0)
      endif

      if (.not. ga_destroy (g_ztmp))
     $     call errquit (pname//"destroy failed", 0,0)


C     XXX CLEAN UP ZORA ARRAYS AFTER PROPAGATION
      end subroutine


#if 0
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
      integer, parameter   :: npol = 1
      character(*), parameter :: pname = "rt_tddft_init_zora: "


C     == Variables ==
      integer g_zora_sf(2)           ! ZORA scalar arrays; real nbf_ao x nbf_ao
      integer g_zora_scale_sf(2)     ! ZORA scalar energy scaling arrays; real nbf_ao x nbf_ao
      character*256 fname
      integer me
      integer g_ztmp
      integer nbf_array(2)


C     == External ==
      logical, external ::  dft_zora_read
      logical, external  :: zmat_compare

      
      me = ga_nodeid()

      if (me.eq.0) write (luout, *) pname//"Before alloc"

      if (.not. do_zora)
     $     call errquit (pname//"called without ZORA active",0,0)

      params%lzorasf = .true.
      params%lzoraso = .false.

C
C     Allocation
C
      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "sf(1)", 0, 0, g_zora_sf(1)))
     $     call errquit (pname//"couldnt make g_zora_sf(1)", 0, GA_ERR)
      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "sf(2)", 0, 0, g_zora_sf(2)))
     $     call errquit (pname//"couldnt make g_zora_sf(2)", 0, GA_ERR)

      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "scale sf(1)", 0, 0, g_zora_scale_sf(1)))
     $     call errquit (pname//"couldnt make g_zora_sf(1)", 0, GA_ERR)
      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "scale sf(2)", 0, 0, g_zora_scale_sf(2)))
     $     call errquit (pname//"couldnt make g_zora_sf(2)", 0, GA_ERR)


      if (me.eq.0) write (luout, *) pname//"Before file_name"

C
C     Get ZORA data from file.  If we can't read it, quit since this is
C     safer than recoding the reconstruction here (in case it chages).
C
      call util_file_name ("zora_sf", .false., .false., fname)

      if (me.eq.0) write (luout, *) pname//"Before zero"

C      call ga_zero (g_zora_sf(1))
C      call ga_zero (g_zora_sf(2))
C      call ga_zero (g_zora_scale_sf(1))
C      call ga_zero (g_zora_scale_sf(2))

      if (me.eq.0) write (luout, *) pname//" XXXX DID NOT ZERO"

      if (me.eq.0) write (luout, *) pname//"Before params"

      nbf_array(1) = params%nbf_ao
      nbf_array(2) = params%nbf_ao

      if (me.eq.0) write (luout, *) pname//"Before read"

      if (ipol .eq. 1) then
         if (.not.dft_zora_read(fname, params%nbf_ao, 1,
     $        nbf_array, mult, g_zora_sf(1), g_zora_scale_sf(1)))
     $        call errquit (pname//
     $        "failed to read scalar ZORA data from: "//trim(fname),0,0)
      elseif (ipol .eq. 2) then
         if (.not.dft_zora_read(fname, params%nbf_ao, 2,
     $        nbf_array, mult, g_zora_sf, g_zora_scale_sf))
     $        call errquit (pname//
     $        "failed to read scalar ZORA data from: "//trim(fname),0,0)
      else
         call errquit (pname//"invalid ipol", ipol, 0)
      endif

      params%g_zora_sf(1) = g_zora_sf(1)
      params%g_zora_sf(2) = g_zora_sf(2)

      params%g_zora_scale_sf(1) = g_zora_scale_sf(1)
      params%g_zora_scale_sf(2) = g_zora_scale_sf(2)


      if (me.eq.0) write (luout, *) pname//"After read"

C
C     Check that matrices we just read in are consistent with those in
C     the ground state code.
C
      if (.not. ga_create(mt_dcpl, params%nbf_ao, params%nbf_ao,
     $     "ztmp", 0, 0, g_ztmp))
     $     call errquit (pname//"couldnt make g_ztmp", 0, GA_ERR)
      
      if (ipol.eq.1) then
         call convert_d2z (1d0, g_zora_sf(1), 0d0, g_zora_sf(1), g_ztmp)
         call zmat_compare_and_print (params%nbf_ao, g_ztmp,
     $        "g_zora_cs_sf", params%tol_zero,
     $        "GS and RT-TDDFT CS ZORA matrices")
      else

         call ga_zero (g_ztmp)
         call convert_d2z (1d0, g_zora_sf(1), 0d0, g_zora_sf(1), g_ztmp)
         call zmat_compare_and_print (params%nbf_ao, g_ztmp,
     $        "g_zora_alpha_sf", params%tol_zero,
     $        "GS and RT-TDDFT alpha ZORA SF matrices")

         call ga_zero (g_ztmp)
         call convert_d2z (1d0, g_zora_sf(2), 0d0, g_zora_sf(2), g_ztmp)
         call zmat_compare_and_print (params%nbf_ao, g_ztmp,
     $        "g_zora_beta_sf", params%tol_zero,
     $        "GS and RT-TDDFT beta ZORA SF matrices")
      endif

      if (.not. ga_destroy (g_ztmp))
     $     call errquit (pname//"destroy failed", 0,0)


C     XXX CLEAN UP ZORA ARRAYS AFTER PROPAGATION
      end subroutine
#endif
