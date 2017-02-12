      subroutine scf_dump_evp_input(g_f,g_s)
      implicit none
!#include "errquit.fh"
!#include "mafdecls.fh"
#include "global.fh"
!#include "rtdb.fh"
!#include "geom.fh"
!#include "util.fh"
!
! call ga_diag( g_fock, g_overlap, g_mos, eigenvalues)
!
      integer g_f ! Fock matrix
      integer g_s ! overlap matrix
      call ga_print(g_f)
      call ga_print(g_s)
      end subroutine
