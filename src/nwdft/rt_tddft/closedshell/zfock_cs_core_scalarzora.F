C     
C     zfock_cs_core_scalarzora.F
C
C     Calculates 1 electron scalar ZORA part (Hcore) of a closed shell
C     complex Fock matrix.  Note that input is the real part of the
C     density matrix, and output is real part of the Fock matrix.
C


      subroutine zfock_cs_core_scalarzora (params, g_densre, g_fcore)
      implicit none

#include "rtdb.fh"
#include "bas.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "global.fh"
#include "zora.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_densre

      
C     == Outputs ==
      integer, intent(in)           :: g_fcore 

      
C     == Parameters ==
      character(*), parameter  :: pname = "zfock_cs_core_scalarzora: "

      

      if (.not. do_zora)
     $     call errquit (pname//"called without ZORA active",0,0)

      call rt_tddft_cs_confirm (params,'zfock_cs_core_scalarzora.F')

      call ga_zero (g_fcore)      
      call ga_add (1d0, g_fcore, 1d0, params%g_zora_sf(1), g_fcore)

      end subroutine




      
#if 0      
      subroutine zfock_cs_core_scalarzora (params, g_densre, g_fcore)
      implicit none

#include "rtdb.fh"
#include "bas.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "global.fh"
#include "zora.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_densre

      
C     == Outputs ==
      integer, intent(in)           :: g_fcore 

      
C     == Parameters ==
      character(*), parameter  :: pname = "zfock_cs_core_scalarzora: "

      
C     == Variables ==
      double precision elapsed
      integer g_sf(2)           ! ZORA scalar arrays; real nbf_ao x nbf_ao
      integer g_scale_sf(2)     ! ZORA scalar energy scaling arrays; real nbf_ao x nbf_ao
      character*255 fname


C     == External ==
      logical, external ::  dft_zora_read


      if (.not. do_zora)
     $     call errquit (pname//"called without ZORA active",0,0)

      call rt_tddft_cs_confirm (params,'zfock_cs_core_scalarzora.F')

      if (params%prof) call prof_start (elapsed)

C
C     Allocation
C
      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "sf(1)", 0, 0, g_sf(1)))
     $     call errquit (pname//"couldnt make g_sf(1)", 0, GA_ERR)
      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "sf(2)", 0, 0, g_sf(2)))
     $     call errquit (pname//"couldnt make g_sf(2)", 0, GA_ERR)

      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "scale sf(1)", 0, 0, g_scale_sf(1)))
     $     call errquit (pname//"couldnt make g_sf(1)", 0, GA_ERR)
      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "scale sf(2)", 0, 0, g_scale_sf(2)))
     $     call errquit (pname//"couldnt make g_sf(2)", 0, GA_ERR)


C
C     Get ZORA data from file.  If we can't read it, quit--rebuilding it
C     on the fly is far too slow.  We only use g_so and g_sf, as the
C     scaling arrays are only applied to the eigenvalues in the SCF.
C
      call util_file_name ("zora_sf", .false., .false., fname)

      call ga_zero (g_sf)
      call ga_zero (g_scale_sf)

      if (ga_nodeid().eq.0)
     $     write (luout, *) pname//"todo: read in once only"


C     1 for ipol = 1
      if (.not.dft_zora_read(fname, params%nbf_ao, 1,
     $     params%nbf_ao, mult, g_sf, g_scale_sf))
     $     call errquit (pname//
     $     "failed to read scalar ZORA data from: "//trim(fname), 0, 0)
      
C      call ga_print (g_sf(1))

      call ga_zero (g_fcore)      
      call ga_add (1d0, g_fcore, 1d0, g_sf(1), g_fcore)

      if (.not. ga_destroy (g_sf(1)))
     $     call errquit (pname//"failed to destroy g_sf(1)", 0, 0)
      if (.not. ga_destroy (g_sf(2)))
     $     call errquit (pname//"failed to destroy g_sf(1)", 0, 0)

      if (.not. ga_destroy (g_scale_sf(1)))
     $     call errquit (pname//"failed to destroy g_scale_sf(1)", 0, 0)
      if (.not. ga_destroy (g_scale_sf(2)))
     $     call errquit (pname//"failed to destroy g_scale_sf(1)", 0, 0)

      if (params%prof)
     $     call prof_end (elapsed, "Fock CS core scalar ZORA")
      
      end subroutine
#endif
c $Id$
