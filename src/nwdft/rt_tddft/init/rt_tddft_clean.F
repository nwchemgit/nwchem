C
C     rt_tddft_params_clean.F
C
C     Deallocates all universal GAs, parameters, etc.
C
      subroutine rt_tddft_clean (params)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "cdft.fh"
#include "geom.fh"      
#include "rt_tddft.fh"
#include "bas.fh"


C     == Inputs ==
      type(rt_params_t), intent(in) :: params


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_clean: "


C     == Variables ==
      logical ok
      integer ik, ig


      
C
C     Clean up dipole and quadrupole integrals.
C
      ok = .true.

      ok = ok .and. ga_destroy(params%g_dip%x)
      ok = ok .and. ga_destroy(params%g_dip%y)
      ok = ok .and. ga_destroy(params%g_dip%z)

      ok = ok .and. ga_destroy (params%g_quad%xx)
      ok = ok .and. ga_destroy (params%g_quad%xy) 
      ok = ok .and. ga_destroy (params%g_quad%xz)
      ok = ok .and. ga_destroy (params%g_quad%yy)
      ok = ok .and. ga_destroy (params%g_quad%yz)
      ok = ok .and. ga_destroy (params%g_quad%zz)

      if (.not.ok)
     $     call errquit(pname//"g_dip, g_quad destroy failed", 0,GA_ERR)


C
C     Clean up 1e integrals
C
      if (.not. ga_destroy (params%g_fock_ao_core))
     $     call errquit (pname//"destroy failed", 0, GA_ERR)

      
C
C     Clean up 2e and CD fitting integrals.
C
      call fock_2e_tidy (params%rtdb)

      if (cdfit) then
         if (.not. ma_pop_stack (params%l_3cwhat)) then
            call ma_summarize_allocated_blocks()
             call errquit (pname//"failed to pop 3cwhat", 
     C           params%l_3cwhat, MA_ERR)
          endif
      
         if (.not. ma_pop_stack (params%l_3ceri))
     $        call errquit (pname//"failed to pop 3ceri", 0, MA_ERR)
      endif


C
C     Clean up geometry masks.
C
      do ig = 1, params%ngeoms
         if (.not. ga_destroy (params%geom(ig)%g_mask))
     $        call errquit (pname//"failed to destroy mask",0,GA_ERR)
      enddo


C
C     Clean up CAP.
C      
      if (params%cap_active) then
         if (.not. ga_destroy (params%g_zcap))
     $        call errquit (pname//"failed to destroy zcap",0,GA_ERR)
      endif

      
C
C     Clean up canonical orthogonalization transforms.
C
      call canorg_clean (params)

CXXX  [KAL]: CLEAN UP ZORA HERE
      

CXXX  [KAL]: Needed if using ECP?  FIX THIS (need clean slate for dplot)!
C      call int_ecp_terminate ()
         
      end subroutine rt_tddft_clean


c $Id$
