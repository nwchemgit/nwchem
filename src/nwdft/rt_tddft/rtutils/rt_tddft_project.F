C====================================================================
C
C     Calculates A~, the projection of the matrix A onto matrix B
C
C     A~ = C'^+ A C'
C
C     where C' is the eigenvector matrix of the matrix B.
C
C     XXX CHECK:  only valid when B is the MO Fock matrix?
C
      subroutine rt_tddft_project (params, dir, g_za, g_zb, g_zaproj)
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
      character, intent(in)          :: dir !="F" for forward, "R" for reverse
      integer, intent(in)            :: g_za
      integer, intent(in)            :: g_zb

C     == Outputs ==
      integer, intent(in)            :: g_zaproj

C     == Parameters ==
      character(*),parameter :: pname = "rt_tddft_project: "

C     == Variables ==
      integer n
      integer ivals, lvals
      integer g_zevecs, g_zevecs_hconj, g_ztmp

      
C     XXX ASSUME SIZE: nmo
      n = params%ns_mo

      if ((dir .ne. "F").and.(dir .ne. "R")
     $     .and.(dir .ne. "f").and.(dir .ne. "r"))
     $     call errquit (pname//"invalid flag: "//dir, 0,0)

      if (.not. ga_duplicate (g_za, g_zevecs, "zevecs"))
     $     call errquit (pname//"failed to create zevecs", 0, GA_ERR)

      if (.not. ga_duplicate (g_za, g_zevecs_hconj, "zvecs_hconj"))
     $     call errquit (pname//"failed to create ztmp1", 0, GA_ERR)

      if (.not. ga_duplicate (g_za, g_ztmp, "ztmp"))
     $     call errquit (pname//"failed to create ztmp", 0, GA_ERR)

      if (.not.ma_push_get(mt_dcpl, params%ns_ao,'vals',lvals,ivals))
     &     call errquit(pname//'cannot allocate vals',0, MA_ERR)


C     Diagonalize B matrix to get C matrix and eigenvalues
      call ga_zero (g_zevecs)
C      call zmat_diag (g_zb, g_zevecs, dcpl_mb(ivals), .true.) !true = reorder evals by energy
      call zmat_diag (g_zb, g_zevecs, dcpl_mb(ivals), .false., .false.) !true = reorder evals by energy
      call ga_zero (g_zevecs_hconj)
      call mat_hermconj (g_zevecs, g_zevecs_hconj)
      

      if ( (dir.eq."F").or.(dir.eq."f")) then  ! C^+ A C
         call ga_zero (g_ztmp)
         call ga_zgemm("N", "N", n, n, n,
     $        z1, g_za, g_zevecs, z0, g_ztmp)
         call ga_zgemm("N", "N", n, n, n,
     $        z1, g_zevecs_hconj, g_ztmp, z0, g_zaproj)
         
      elseif ( (dir.eq."R").or.(dir.eq."r")) then   ! C A C^+
         call ga_zero (g_ztmp)
         call ga_zgemm("N", "N", n, n, n,
     $        z1, g_za, g_zevecs_hconj, z0, g_ztmp)
         call ga_zgemm("N", "N", n, n, n,
     $        z1, g_zevecs, g_ztmp, z0, g_zaproj)
         
      else
         call errquit (pname//"invalid flag: "//dir, 0,0)
      endif
      

      if (.not. ma_chop_stack (lvals))
     &   call errquit(pname//'cannot chop stack',0, MA_ERR)

      if (.not. ga_destroy (g_zevecs))
     $     call errquit (pname//"failed to destrory zevecs", 0, GA_ERR)
      
      if (.not. ga_destroy (g_zevecs_hconj))
     $     call errquit (pname//"failed to destrory ztmp1", 0, GA_ERR)

      if (.not. ga_destroy (g_ztmp))
     $     call errquit (pname//"failed to destrory ztmp", 0, GA_ERR)

      end subroutine


C====================================================================
C
C     Check the projecting A->B and reversing yields A again.
C
      subroutine rt_tddft_project_check (params, g_za, g_zb)
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
      integer, intent(in)            :: g_za
      integer, intent(in)            :: g_zb

C     == Parameters ==
      character(*),parameter :: pname = "rt_tddft_project_check: "


C     == Variables ==
      integer g_zaproj, g_zaback

      
      if (.not. ga_duplicate (g_za, g_zaproj, "zaproj"))
     $     call errquit (pname//"failed to create zaproj", 0, GA_ERR)

      if (.not. ga_duplicate (g_za, g_zaback, "zaback"))
     $     call errquit (pname//"failed to create zaback", 0, GA_ERR)

      
      call ga_zero (g_zaproj)
      call rt_tddft_project (params, "F", g_za, g_zb, g_zaproj)
      call rt_tddft_project (params, "R", g_zaproj, g_zb, g_zaback)

      if (.not. mat_is_same (g_za, g_zaback, params%tol_zero))
     $     call errquit (pname//"failed",0,0)

C      call ga_print (g_zaproj)
      

      if (.not. ga_destroy (g_zaproj))
     $     call errquit (pname//"failed to destrory zaproj", 0, GA_ERR)

      if (.not. ga_destroy (g_zaback))
     $     call errquit (pname//"failed to destrory zaback", 0, GA_ERR)

      end subroutine



C====================================================================
C
C     Project complex TD matrix g_za (in MO basis) onto ground state
C     representation.
C
      subroutine rt_tddft_project_gs (params, dir, g_za,
     $     g_movecs_gs, g_zaproj)
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
      character, intent(in)          :: dir !="F" for forward, "R" for reverse
      integer, intent(in)            :: g_za
      integer, intent(in)            :: g_movecs_gs

C     == Outputs ==
      integer, intent(in)            :: g_zaproj

C     == Parameters ==
      character(*),parameter :: pname = "rt_tddft_project_gs: "

C     == Variables ==
      integer g_zmovecs_gs
      integer n
      integer g_ztmp

      

C     XXX HARDCODED SIZE
      n = params%ns_mo
      

      if ((dir .ne. "F").and.(dir .ne. "R")
     $     .and.(dir .ne. "f").and.(dir .ne. "r"))
     $     call errquit (pname//"invalid flag: "//dir, 0,0)

      if (.not. ga_duplicate (g_za, g_zmovecs_gs, "zmovecs"))
     $     call errquit (pname//"alloc failed", 0, GA_ERR)

      if (.not. ga_duplicate (g_za, g_ztmp, "ztmp"))
     $     call errquit (pname//"failed to create ztmp", 0, GA_ERR)

      
      call ga_zero (g_zmovecs_gs)
      call convert_d2z(1d0, g_movecs_gs, 0d0, g_movecs_gs, g_zmovecs_gs)

      
      if ( (dir.eq."F").or.(dir.eq."f")) then  ! C^T A C
         call ga_zero (g_ztmp)
         call ga_zgemm("N", "N", n, n, n,
     $        z1, g_za, g_zmovecs_gs, z0, g_ztmp)
         call ga_zgemm("T", "N", n, n, n,
     $        z1, g_zmovecs_gs, g_ztmp, z0, g_zaproj)
         
      elseif ( (dir.eq."R").or.(dir.eq."r")) then   ! C A C^T
         call ga_zero (g_ztmp)
         call ga_zgemm("N", "T", n, n, n,
     $        z1, g_za, g_zmovecs_gs, z0, g_ztmp)
         call ga_zgemm("N", "N", n, n, n,
     $        z1, g_zmovecs_gs, g_ztmp, z0, g_zaproj)
         
      else
         call errquit (pname//"invalid flag: "//dir, 0,0)
      endif


      if (.not. ga_destroy (g_zmovecs_gs))
     $     call errquit (pname//"destroy failed", 0, GA_ERR)

      if (.not. ga_destroy (g_ztmp))
     $     call errquit (pname//"failed to destrory ztmp", 0, GA_ERR)
      

      end subroutine




C====================================================================
C
C     Check the GS projections
C
      subroutine rt_tddft_project_gs_check (params, g_za, g_movecs_gs)
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
      integer, intent(in)            :: g_za
      integer, intent(in)            :: g_movecs_gs


C     == Parameters ==
      character(*),parameter :: pname="rt_tddft_project_gs_check: "


C     == Variables ==
      integer g_zaproj, g_zaback

      
      if (.not. ga_duplicate (g_za, g_zaproj, "zaproj"))
     $     call errquit (pname//"failed to create zaproj", 0, GA_ERR)

      if (.not. ga_duplicate (g_za, g_zaback, "zaback"))
     $     call errquit (pname//"failed to create zaback", 0, GA_ERR)

      call ga_zero (g_zaproj)
      call rt_tddft_project_gs (params, "F",
     $     g_za, g_movecs_gs, g_zaproj)
      call rt_tddft_project_gs (params, "R",
     $     g_zaproj, g_movecs_gs, g_zaback)

      if (.not. mat_is_same (g_za, g_zaback, params%tol_zero)) then
         call ga_print (g_za)
         call ga_print (g_zaback)
         call errquit (pname//"failed",0,0)
      endif

C      call ga_print (g_zaproj)


      if (.not. ga_destroy (g_zaproj))
     $     call errquit (pname//"failed to destrory zaproj", 0, GA_ERR)

      if (.not. ga_destroy (g_zaback))
     $     call errquit (pname//"failed to destrory zaback", 0, GA_ERR)

      end subroutine
