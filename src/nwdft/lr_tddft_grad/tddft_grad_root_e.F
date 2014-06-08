c
c     Create and destroy workspace for the excitation energies
c     HvD 5/2006, NG 11/2012
c
      subroutine tddft_grad_create_root_e(ntotroots,ihdl_root_e,
     &    iptr_root_e)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
c
c     Input:
c
      integer ntotroots 
c
c     Output:
c
      integer ihdl_root_e
      integer iptr_root_e 
c
      character*32 pname
c
      pname="tddft_grad_create_root_e: "
c
      if (.not.ma_push_get(mt_dbl,ntotroots,'root_e',
     +   ihdl_root_e, iptr_root_e))
     + call errquit(pname//'failed to allocate root_e',ntotroots,MA_ERR)
c
      end
c
      subroutine tddft_grad_destroy_root_e(ihdl_root_e)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
c
      integer ihdl_root_e  ! the memory handle for the MO energies
c
      character*32 pname
      pname = "tddft_grad_destroy_root_e: "
c
      if (.not.ma_pop_stack(ihdl_root_e))
     +  call errquit(pname//'failed to release root_e',0, MA_ERR)
c
      end
c $Id$
