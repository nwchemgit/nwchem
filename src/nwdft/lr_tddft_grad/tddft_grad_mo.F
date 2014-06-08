c
c     Create and Destroy an atom blocked global array to hold the MO coefficients.
c     HvD 5/2006, NG 11/2012
c
      subroutine tddft_grad_create_mo(ipol,ihdl_geom,ihdl_bas,g_mo)
c
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
c
c     Input:
c
      integer ipol      ! =1 (restricted), =2 (unrestricted)
      integer ihdl_geom ! handle for the geometry
      integer ihdl_bas  ! handle for the AO basis
c
c     Output:
c
      integer g_mo(2)   ! global array handles for the MOs
c
c     Functions:
c
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
c     Local:
c
      integer i         ! counter
c
      do i = 1,ipol
        g_mo(i) = ga_create_atom_blocked(ihdl_geom,ihdl_bas,
     +            'MO eigenvectors')
      enddo
c
      end
c
      subroutine tddft_grad_destroy_mo(ipol,g_mo)
c
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
c
c     Input:
c
      integer ipol      ! =1 (restricted), =2 (unrestricted)
      integer g_mo(2)   ! global array handles for the MOs
c
c     Local:
c
      integer i         ! counter
      character*32 pname
      pname = "tddft_grad_destroy_mo :"
c
      do i = 1,ipol
        if (.not.ga_destroy(g_mo(i))) 
     +    call errquit(pname//'failed to destroy g_mo',0, GA_ERR)
      enddo
c
      end
c $Id$
