c
c  Evaluate electronic energy expectation value
c  This follows closely "mcscf_fock"
c
c
       subroutine mcscf_etrace( geom, basis, nbf, nclosed, nact,
     $                          oskel, tol2e, dm1, dm2, g_movecs,
     $                          g_coul, eone, etwo )
       implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "util.fh"
#include "geom.fh"
       integer geom, basis                            ! [input] Geometry and basis handles
       integer nbf                                    ! [input] Number of basis functions
       integer nclosed                                ! [input] Number of closed shells
       integer nact                                   ! [input] Number of open shells
       logical oskel                                  ! [input] Symmetry toggle
       double precision tol2e                         ! [input] Integral tolerance
       double precision dm1(nact,nact)                ! [input] Active space 1PDM
       double precision dm2(nact,nact,nact,nact)      ! [input] Active space 2PDM
       integer g_movecs                               ! [input] MO coefficients
       integer g_coul                                 ! [input] Coulomb integrals
       double precision eone, etwo                    ! [output] Energy components
c
c
       integer nset
       parameter(nset=1)
       integer g_ifock, g_cdens, g_adens, g_tmp
       double precision e2ii, e2ia, e2ai, e2aa
       integer iv_dens(2), iv_fock(2)
       double precision jfac(2), kfac(2)
       integer ga_create_atom_blocked
       double precision ga_trace_diag, mcscf_act_ener
       external ga_create_atom_blocked, ga_trace_diag, mcscf_act_ener
       data jfac/1.0d0, 1.d0/, kfac/-0.5d0, -0.5d0/
c
c
c
       g_tmp = ga_create_atom_blocked(geom, basis, 'temp1')
       g_ifock = ga_create_atom_blocked(geom, basis, 'fock')
       g_cdens = ga_create_atom_blocked(geom, basis, 'closed dens')
       g_adens = ga_create_atom_blocked(geom, basis, 'active dens')
       call ga_zero(g_ifock)
       call mcscf_occ2dens(nbf, nclosed, nact, dm1, g_movecs,
     $                     g_cdens, g_adens )
c
c One-electron component
c
       call ga_zero(g_tmp)
       call int_1e_ga( basis, basis, g_tmp, 'kinetic', oskel)
       call int_1e_ga( basis, basis, g_tmp, 'potential', oskel)
       if (oskel) call sym_symmetrize(geom, basis, .false., g_tmp)
       eone = ga_ddot(g_tmp,g_cdens) + ga_ddot(g_tmp,g_adens)
c
c Two-electron component of the AO Fock matrices
c
       iv_dens(1) = g_cdens
       iv_fock(1) = g_ifock
       call fock_2e( geom, basis, nset, jfac, kfac, tol2e, oskel,
     $               iv_dens, iv_fock )
c
c  Symmetrize Fock AO components
c
       if (oskel)
     $    call sym_symmetrize(geom, basis, .false., g_ifock)
       e2ii = ga_ddot(g_ifock,g_cdens)
       e2ia = ga_ddot(g_ifock,g_adens)
       e2ai = e2ia
c
c  MO integral contribution
c
       e2aa = mcscf_act_ener( nbf, nclosed, nact, dm2, g_coul )
       etwo = (e2ii + e2ia + e2ai)*0.5d0 + e2aa

       if ((ga_nodeid().eq.0).and.
     $     (util_print('energy trace',print_debug))) then
         write(6,911) e2ii*0.5d0,e2ai*0.5d0,e2aa
 911     format(  'energy components'
     $          /,'    eii: ',f12.6,
     $          /,'    eia: ',f12.6,
     $          /,'    eaa: ',f12.6)
       endif
c
c  Clean up
c
       if (.not.ga_destroy(g_tmp))
     $   call errquit('mcscf_fock: cannot destroy',0)
       if (.not.ga_destroy(g_cdens))
     $   call errquit('mcscf_fock: cannot destroy',0)
       if (.not.ga_destroy(g_adens))
     $   call errquit('mcscf_fock: cannot destroy',0)
       if (.not.ga_destroy(g_ifock))
     $   call errquit('mcscf_fock: cannot destroy',0)
       return
       end






c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
c
c          DEBUGGING VERSION
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
       subroutine mcscf_etrace00( geom, basis, nbf, nclosed, nact,
     $                            oskel, tol2e, dm1, dm2, g_movecs,
     $                            g_coul, eone, e2ii, e2ia, e2aa )
       implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
       integer geom, basis                            ! [input] Geometry and basis handles
       integer nbf                                    ! [input] Number of basis functions
       integer nclosed                                ! [input] Number of closed shells
       integer nact                                   ! [input] Number of open shells
       logical oskel                                  ! [input] Symmetry toggle
       double precision tol2e                         ! [input] Integral tolerance
       double precision dm1(nact,nact)                ! [input] Active space 1PDM
       double precision dm2(nact,nact,nact,nact)      ! [input] Active space 2PDM
       integer g_movecs                               ! [input] MO coefficients
       integer g_coul                                 ! [input] Coulomb integrals
       double precision eone                          ! [output] Energy components
       double precision e2ii, e2ia, e2aa
c
c
       integer nset
       parameter(nset=1)
       integer g_ifock, g_cdens, g_adens, g_tmp
       integer iv_dens(2), iv_fock(2)
       double precision jfac(2), kfac(2)
       integer ga_create_atom_blocked
       double precision ga_trace_diag, mcscf_act_ener
       external ga_create_atom_blocked, ga_trace_diag, mcscf_act_ener
       data jfac/1.0d0, 1.d0/, kfac/-0.5d0, -0.5d0/
c
c
c
       g_tmp = ga_create_atom_blocked(geom, basis, 'temp1')
       g_ifock = ga_create_atom_blocked(geom, basis, 'fock')
       g_cdens = ga_create_atom_blocked(geom, basis, 'closed dens')
       g_adens = ga_create_atom_blocked(geom, basis, 'active dens')
       call ga_zero(g_ifock)
       call mcscf_occ2dens(nbf, nclosed, nact, dm1, g_movecs,
     $                     g_cdens, g_adens )
c
c One-electron component
c
       call ga_zero(g_tmp)
       call int_1e_ga( basis, basis, g_tmp, 'kinetic', oskel)
       call int_1e_ga( basis, basis, g_tmp, 'potential', oskel)
       if (oskel) call sym_symmetrize(geom, basis, .false., g_tmp)
       eone = ga_ddot(g_tmp,g_cdens) + ga_ddot(g_tmp,g_adens)
c
c Two-electron component of the AO Fock matrices
c
       iv_dens(1) = g_cdens
       iv_fock(1) = g_ifock
       call fock_2e( geom, basis, nset, jfac, kfac, tol2e, oskel,
     $               iv_dens, iv_fock )
c
c  Symmetrize Fock AO components
c
       if (oskel)
     $    call sym_symmetrize(geom, basis, .false., g_ifock)
       e2ii = ga_ddot(g_ifock,g_cdens)*0.5d0
       e2ia = ga_ddot(g_ifock,g_adens)
c
c  MO integral contribution
c
C       e2aa = mcscf_act_ener( nbf, nclosed, nact, dm2, g_coul )
       e2aa = 0.d0
c
c  Clean up
c
       if (.not.ga_destroy(g_tmp))
     $   call errquit('mcscf_fock: cannot destroy',0)
       if (.not.ga_destroy(g_cdens))
     $   call errquit('mcscf_fock: cannot destroy',0)
       if (.not.ga_destroy(g_adens))
     $   call errquit('mcscf_fock: cannot destroy',0)
       if (.not.ga_destroy(g_ifock))
     $   call errquit('mcscf_fock: cannot destroy',0)
       return
       end
       


