      subroutine dimqm_tidy()
      implicit none
#include "dimqm.fh"
#include "errquit.fh"
#include "global.fh"
      dimqm_seeded = .false.
      xyz_seeded = .false.
c     Destroy diagonal array
      if(.not. ga_destroy(g_DIM_diag))
     $ call errquit("destroy DIM diag failed", 1, GA_ERR)
      

      end subroutine dimqm_tidy
