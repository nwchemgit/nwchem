c
c     Computes the 1-electron density matrices in AO basis from the MOs. 
c
      subroutine tddft_grad_compute_dao(ipol,nao,nocc,g_mo,g_dao)
c
      implicit none
c
#include "global.fh"
c
c     Input:
c
      integer ipol    ! =1 (restricted), =2 (unrestricted)
      integer nao     ! the number of atomic basis functions
      integer nocc(2) ! the number of occupied orbitals
      integer g_mo(2) ! the global array handle for the MOs
c
c     Output:
c
      integer g_dao(2) ! the global array handle for the density 
                       ! matrices
c
c     Local:
c
      integer ip       ! counter for the spin components
c
      do ip = 1, ipol
        call ga_dgemm('n','t',nao,nao,nocc(ip),1.0d0,g_mo(ip),g_mo(ip),
     +                0.0d0,g_dao(ip))
      enddo
c
      end
c $Id$
