C 
C     Compute the scaled real part of the AO ZORA density matrix.  This
C     is used before computing time-dependent properties.  Note this is
C     constructed using the ground state movecs and the time-dependent
C     MO occupations.
C
C     All matrices are real, since properties only depend on real part
C     anyways.
C
      subroutine rt_tddft_zora_scaled_densre (params, g_movecs_gs,
     $     moocc, g_densre_ao_scaled)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "geom.fh"
#include "util.fh"
#include "cdft.fh"
#include "rtdb.fh"
#include "zora.fh"      
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_movecs_gs !ground state movecs in AO basis ns x ns (0 padded if lindep)
      double precision, intent(in)  :: moocc(*) !MO occupations

C     == Output  ==
      integer, intent(in)           :: g_densre_ao_scaled !in=unscaled dens mat; out=scaled dens mat
      
C     == Parameters ==
      character(*),parameter :: pname="rt_tddft_zora_scale: "

C     == Variables ==
      integer iorb
      integer g_tmp, g_tmp2
      double precision orb_scale
      double precision zora_eint, zora_scale


      if (.not. ga_duplicate (g_movecs_gs, g_tmp, "tmp"))
     $     call errquit (pname//"failed to duplicate g_tmp",0,0)

      if (.not. ga_duplicate (g_movecs_gs, g_tmp2, "tmp2"))
     $     call errquit (pname//"failed to duplicate g_tmp2",0,0)

C
C     First, g_tmp is 0 except iorb column with movecs scaled by
C     time-dependent occupation.  Then, it becomes the density matrix
C     due to this orbital.
C
      call ga_zero (g_densre_ao_scaled)

      do iorb = 1, params%ns_ao
         call ga_zero (g_tmp)
         call ga_copy_patch ("N", 
     $        g_movecs_gs, 1, params%ns_ao, iorb, iorb,
     $        g_tmp, 1, params%ns_ao, iorb, iorb)
         orb_scale = moocc(iorb)/2d0 !XXX FOR CS ONLY

         call ga_scale (g_tmp, orb_scale)

         call ga_zero (g_tmp2)
         call ga_dgemm ("N", "T",
     $        params%ns_ao, params%ns_ao, params%ns_ao,
     $        1d0, g_tmp, g_tmp, 0d0, g_tmp2)

C     (XXX TMP DISABLED)
C         zora_eint = ga_ddot (g_tmp2, params%g_zora_scale_sf(1)) !xxx 1 for ipol; XXX make general
C         zora_scale = 1d0 / (1d0 + zora_eint)
C         call ga_scale (g_tmp2, zora_scale)
         
C         if (ga_nodeid().eq.0)
C     $        write (luout, *) "XXXRT zora_eint = ", zora_eint
         

C     XXX NOTE: 2 for CS; MAKE THIS MORE ELEGANT
         call ga_add (1d0, g_densre_ao_scaled, 2d0, g_tmp2,
     $        g_densre_ao_scaled)
      enddo

      if (.not. ga_destroy (g_tmp))
     $     call errquit (pname//"failed to destroy g_tmp",0,0)

      if (.not. ga_destroy (g_tmp2))
     $     call errquit (pname//"failed to destroy g_tmp2",0,0)

      end subroutine
c $Id$
