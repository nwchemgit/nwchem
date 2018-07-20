C
C     rt_tddft_cs_init.F
C
C     Initialize closed shell RT-TDDFT.
C
      subroutine rt_tddft_cs_init (params)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "matutils.fh"
#include "cdft.fh"
#include "rt_tddft.fh"


C     == In/out ==
      type(rt_params_t), intent(inout) ::  params

      
C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_cs_init: "


C     == Variables ==
      integer me
      double precision nuclear_charge

      
      me = ga_nodeid ()


      params%theory = "CS"

      call rt_tddft_init_common (params)

      
C
C     Set up number of electrons for closed shell case.
C
      if (params%ipol .ne. 1)
     $     call errquit (pname//"only valid for ipol = 1", 0, 0)

      if (.not. rtdb_get(params%rtdb, 'dft:noc', mt_int, 2, noc))
     &     call errquit(pname//"rtdb_get of noc failed", 0, RTDB_ERR)

      params%nalpha = noc(1)
      params%nbeta = params%nalpha
      params%ntotel = params%nalpha + params%nbeta


C
C     Nuclear charge
CXXX  [KAL]: this doesnt seem to work with ECPs ...
C
      if (.not. geom_nuc_charge(geom, nuclear_charge))
     $     call errquit(pname//"geom_nuc_charge failed", 0, GEOM_ERR)

      params%syscharge = dint (nuclear_charge - params%ntotel)
      


C
C     Rebuild overlap matrix and initialze canonical orthogonalization
C     transforms.  This was done in dft_main0d, but we do it here again
C     so that the rt_tddft module is completely self-contained.  This
C     might also be done during the propagation if the nuclei move.
C
      call rt_tddft_init_overlap_canorg (params)


C
C     Super matrix size is just nbf_ao or nbf_mo (since this is vanilla
C     closed shell).
C      
      params%ns_ao = params%nbf_ao
      params%ns_mo = params%nbf_mo


      call rt_tddft_print_initialization (params)
      call rt_tddft_cs_confirm (params,'rt_tddft_cs_init.F')
      
      end subroutine rt_tddft_cs_init


c $Id$
