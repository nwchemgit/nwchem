C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C NAME
C     rhf_fock -- Compute and return a complete Fock matrix, along with
C     one- and two-electron energies.
C
C REVISION
C     $Id$
C
C SYNOPSIS
      Subroutine RHF_Fock( Geom, Basis, G_Dens, Tol2e, OSkel,
     $   G_Fock, EOne, ETwo)
      Implicit NONE
#include "errquit.fh"
      Integer Geom           ! [in]  Geometry handle
      Integer Basis          ! [in]  Basis set handle
      Integer G_Dens         ! [in]  GA handle for AO density
      Double Precision Tol2e ! [in]  selection threshold for 2e integs
      Logical OSkel          ! [in]  Use symmetry (petite-list skeleton)
      Integer G_Fock         ! [in]  GA handle for Fock matrix (GA overwritten)
      Double Precision EOne  ! [out] Energy of one-electron terms
      Double Precision ETwo  ! [out] Energy of two-electron terms
C
C DESCRIPTION
C     Given a density matrix construct a full RHF Fock matrix (both one-
C     and two-electron terms), returning the Fock matrix along with the
C     energies.
C
C     E(SCF) = EOne + ETwo + (Nuclear Repulsion Energy)
C
C INCLUDED FILES
#include "global.fh"
c$$$#include "cscfps.fh"
C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LOCAL VARIABLES
      Integer G_HCore ! Temporary GA, cloned from G_Fock
c     
c     The guess density must be symmetrized if sym is used
c     
      if (OSkel) then
c$$$         if (oscfps) call pstat_on(ps_sym_sym)
         call sym_symmetrize(geom, basis, .true., G_Dens)
c$$$         if (oscfps) call pstat_off(ps_sym_sym)
      endif
c     
      If ( .NOT. GA_Duplicate( G_Fock, G_HCore, 'Temp HCore') )
     $   Call ErrQuit('RHF_Fock: Can''t create GA for HCore', 0, GA_ERR)
c     
c     Build the fock_matrix
c     
      call ga_zero(G_Fock)
      call ga_zero(g_hcore)
      call int_1e_ga(basis, basis, g_hcore,'kinetic', OSkel)
      call int_1e_ga(basis, basis, g_hcore,'potential', OSkel)
c     
      call rhf_fock_2e(geom, basis, G_Dens, G_Fock, Tol2e,
     &   .true., .true., OSkel)
c     
c     Compute contributions to the energy and symmetrize
c     
      EOne = ga_ddot(g_hcore,G_Dens)
      ETwo = 0.5d0*ga_ddot(G_Fock,G_Dens)
      call ga_dadd(1.d0,g_hcore,1.d0,G_Fock,G_Fock)
      if (.not. ga_destroy(g_hcore)) call errquit
     $     ('rhf_fock: destroy of hcore?', 0, GA_ERR)
c     
      if (OSkel) then
c$$$         if (oscfps) call pstat_on(ps_sym_sym)
         call sym_symmetrize(geom, basis, .false., G_Fock)
c$$$         if (oscfps) call pstat_off(ps_sym_sym)
      endif
c     
      Return
      End
