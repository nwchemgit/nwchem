c-----------------------------------------------------------------------
c
      subroutine wfn1_allocate_mat(wfn1_inst,wfn1_mat)
      implicit none
#include "wfn1_param.fh"
#include "wfn1_matcnst.fh"
#include "global.fh"
!>
!> \brief Allocate 1-electron integral matrices
!>
!> Allocate global arrays to accomodate the overlap matrix, the
!> 1-electron Hamiltonian, and the Cauchy-Schwarz integrals.
!>
      type(wfn1_prminst), intent(in)  :: wfn1_inst
      type(wfn1_matcnst), intent(out) :: wfn1_mat
      integer, external :: ga_create_atom_blocked
!>
      wfn1_mat%g_s   = ga_create_atom_blocked(
     &                   wfn1_inst%geom,wfn1_inst%basis,"overlap")
      wfn1_mat%g_s12 = ga_create_atom_blocked(
     &                   wfn1_inst%geom,wfn1_inst%basis,"S^{-1/2}")
      wfn1_mat%g_tv  = ga_create_atom_blocked(
     &                   wfn1_inst%geom,wfn1_inst%basis,"1-electron")
      wfn1_mat%g_cs  = ga_create_atom_blocked(
     &                   wfn1_inst%geom,wfn1_inst%basis,"schwarz")
!>
      end
c
c-----------------------------------------------------------------------
c
      subroutine wfn1_deallocate_mat(wfn1_mat)
      implicit none
#include "wfn1_param.fh"
#include "wfn1_matcnst.fh"
#include "global.fh"
#include "errquit.fh"
!>
!> \brief Allocate 1-electron integral matrices
!>
!> Allocate global arrays to accomodate the overlap matrix, the
!> 1-electron Hamiltonian, and the Cauchy-Schwarz integrals.
!>
      type(wfn1_matcnst), intent(inout) :: wfn1_mat
      character*20 pname
      pname = "wfn1_deallocate_mat:"
!>
      if (.not.ga_destroy(wfn1_mat%g_s)) then
        call errquit(pname//" failed to dealloc overlap",0,GA_ERR)
      endif
      if (.not.ga_destroy(wfn1_mat%g_s12)) then
        call errquit(pname//" failed to dealloc S^{-1/2}",0,GA_ERR)
      endif
      if (.not.ga_destroy(wfn1_mat%g_tv)) then
        call errquit(pname//" failed to dealloc 1-electron",0,GA_ERR)
      endif
      if (.not.ga_destroy(wfn1_mat%g_cs)) then
        call errquit(pname//" failed to dealloc schwarz",0,GA_ERR)
      endif
!>
      end
c
c-----------------------------------------------------------------------
