C====================================================================
C
C     Calculates the time-dependent occupations of the molecular orbitals.
C
C     n_k(t) = C'_k^+ P'(t) C'_k
C
C     where C'_k is the k^th eigenvector of the ground state Fock
C     matrix, and P' is the density matrix in the MO basis.  Note that P
C     is complex, but C' is real since it is from the SCF.  This means
C     that we can just use the real part of the dens mat, and take
C     transpose of C' instead of conjg transpose.
C
C     Note, you can send this either the full dens mat, or just the
C     alpha or beta spin part.
C
C     Note, can also compute using matrix multiplications (XXX double check):
C     n_k(t) = [C'^+ P'(t) C']_kk
C
      subroutine rt_tddft_moocc_calc (params, g_densre_mo, g_movecs_gs,
     $     moocc)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "geom.fh"
#include "util.fh"
#include "rt_tddft.fh"
#include "matutils.fh"      

C     == Inputs ==
      type (rt_params_t), intent(in) :: params
      integer, intent(in)            :: g_densre_mo  !re part of density matrix in MO basis
      integer, intent(in)            :: g_movecs_gs  !ground state movecs


C     == Outputs ==
      double precision, intent(out) :: moocc(*) !MO occupations

      
C     == Parameters ==
      character(*),parameter :: pname = "rt_tddft_moocc: "


C     == Variables ==
      integer me
      integer i
      integer lveck, iveck      !handle and index for kth eigenvector
      double precision occk     !occupation of orbital k
      integer g_veck, g_tmp, g_veckt
      integer k
      double precision dval

      me = ga_nodeid()

      
      if (.not. ga_create(mt_dbl, params%ns_mo, 1,
     $     "k^th evec", 0, 0, g_veck))
     $     call errquit ("failed to create veck", 0, GA_ERR)

      if (.not. ga_duplicate(g_veck, g_tmp, "moocc tmp"))
     $     call errquit ("failed to create g_tmp", 0, GA_ERR)

      if (.not. ga_duplicate(g_veck, g_veckt, "col of transpose C"))
     $     call errquit ("failed to create g_tmp", 0, GA_ERR)

C
C     Load k^th evec in g_veck and k^th column of C in g_vectk.
C
      do k = 1, params%ns_mo
CXXX  [KAL]: its redundant to have two g_veck

         call ga_zero (g_veck)
         call ga_zero (g_veckt)
         call ga_zero (g_tmp)

         call ga_copy_patch ("N",
     $        g_movecs_gs, 1, params%ns_mo, k, k,
     $        g_veck, 1, params%ns_mo, 1, 1)
         call ga_copy_patch ("T",
     $        g_movecs_gs, 1, params%ns_mo, k, k,
     $        g_veckt, 1, params%ns_mo, 1, 1)


C     Compute P'(t) C'_k.
         call ga_dgemm ("N", "N", params%ns_mo, 1, params%ns_mo,
     $        1d0, g_densre_mo, g_veck, 0d0, g_tmp)

         
C     Compute n_k = C'_k^T P'(t) C'_k.
         
         moocc(k) = ga_ddot (g_veckt, g_tmp)
      enddo


C     == Clean up ==
      if (.not. ga_destroy(g_veck))
     $     call errquit ("failed to destroy g_veck", 0, GA_ERR)
      if (.not. ga_destroy(g_tmp))
     $     call errquit ("failed to destroy g_tmp", 0, GA_ERR)
      if (.not. ga_destroy(g_veckt))
     $     call errquit ("failed to destroy g_tmp", 0, GA_ERR)

      end subroutine


C====================================================================
      subroutine rt_tddft_moocc_print (params, tt, moocc, moocc_tag)
      implicit none

#include "stdio.fh"
#include "global.fh"
#include "rt_tddft.fh"


C     == Inputs ==
      type(rt_params_t), intent(in) :: params
C      integer, intent(in)           :: it
      double precision, intent(in)  :: tt
      double precision, intent(in)  :: moocc(params%ns_mo)
      character(*), intent(in)      :: moocc_tag


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_moocc_print: "

      
C     == Variables ==      
      integer me
      integer k

      me = ga_nodeid()

      if (params%nt < 1) call errquit (pname//"nt must be > 0", 0, 0)
      
      if (me.eq.0) then


C
C     Write rt-tddft tag and current time
C
         write (luout, "(a,2x,1f11.5)", advance="no")
     $        trim(params%tag), tt
         
C     $        trim(params%tag), it*100/params%nt, "% ", tt
         
         
C
C     Loop over all MOs and print
C     
         do k = 1, params%ns_mo
            write (luout, "(1es22.12e3)", advance="no") moocc(k)
         enddo

C
C     Print tag and finish line
C     
         write(luout, *) "    "//moocc_tag
         
      endif
      
      call util_flush (luout)

      end subroutine


c $Id$
