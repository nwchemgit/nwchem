c
c     Create a global array to hold the TDDFT gradients.
c
      subroutine tddft_grad_create_g(nat,nroots,g_g)
c
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
c
c     Input:
c
      integer nat    ! the number of atoms
      integer nroots ! the number of roots
c
c     Output:
c
      integer g_g       ! global array handle for the gradients
c
c     Local:
c
      integer idim(3)   ! the dimensions of the global array
      integer ichnk(3)  ! the chunk sizes of the global array
c
      idim(1)  = nroots
      idim(2)  = 3
      idim(3)  = nat
      ichnk(1) = -1
      ichnk(2) = -1
      ichnk(3) = -1
      if (.not.nga_create(mt_dbl,3,idim,"g_g",ichnk,g_g)) call errquit
     +  ("tddft_grad_create_g: could create global array g_g",0, GA_ERR)
c
      call ga_zero(g_g)
c
      end
c
c     Destroy the global array to hold the TDDFT gradients.
c
      subroutine tddft_grad_destroy_g(g_g)
c
      implicit none
c
#include "errquit.fh"
#include "global.fh"
c
c     Input:
c
      integer g_g       ! global array handle for the gradients
c
      if (.not.ga_destroy(g_g)) call errquit
     +  ("tddft_grad_destroy_g: could not destroy g_g",0, GA_ERR)
c
      end
c $Id$
