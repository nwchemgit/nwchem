c
c     Create and Destroy workspace: orbital energies
c     HVD 5/2006, NG 11/2012
c
      subroutine tddft_grad_create_e(ipol,nao,ihdl_e,iptr_e)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
c
      integer ipol    ! =1 (restricted), =2 (unrestricted)
      integer nao     ! the number of AOs
      integer ihdl_e  ! the memory handle for the MO energies
      integer iptr_e  ! the memory index  for the MO energies
c
      character*32 pname
      pname ="tddft_grad_create_e: "
c
      if (.not.ma_push_get(mt_dbl,nao*ipol,'e_mo',ihdl_e,iptr_e))
     +  call errquit(pname//'failed to create e_mo',0,MA_ERR)
c
      end
c
      subroutine tddft_grad_destroy_e(ihdl_e)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
c
c
      integer ihdl_e  ! the memory handle for the MO energies
      character*32 pname
      pname ="tddft_grad_destroy_e: "
c
      if (.not.ma_pop_stack(ihdl_e))
     +  call errquit(pname//'failed to release e_mo',0,MA_ERR)
c
      end
c $Id$
