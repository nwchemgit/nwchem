c
c  Evaluate electronic energy expectation value
c  This follows closely "mcscf_fock"
c
c
       subroutine mcscf_etrace( geom, basis, nbf, nclosed, nact,
     $                          oact, oskel, tol2e, dm1, dm2,
     $                          g_movecs, g_coul, eone, etwo, ecore )
*
* $Id$
*
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "util.fh"
#include "geom.fh"
       integer geom, basis                            ! [input] Geometry and basis handles
       integer nbf                                    ! [input] Number of basis functions
       integer nclosed                                ! [input] Number of closed shells
       integer nact                                   ! [input] Number of open shells
       logical oact                                   ! [input] Toggle include active space 
       logical oskel                                  ! [input] Symmetry toggle
       double precision tol2e                         ! [input] Integral tolerance
       double precision dm1(nact,nact)                ! [input] Active space 1PDM
       double precision dm2(nact,nact,nact,nact)      ! [input] Active space 2PDM
       integer g_movecs                               ! [input] MO coefficients
       integer g_coul                                 ! [input] Coulomb integrals
       double precision eone, etwo, ecore             ! [output] Energy components
c
       integer nset
       parameter(nset=1)
c
       integer g_ifock, g_cdens, g_adens, g_tmp
       double precision e2ii, e2ia, e2ai, e2aa
       integer i
       double precision xx
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
       if (oact) then
         eone = ga_ddot(g_tmp,g_cdens) + ga_ddot(g_tmp,g_adens)
       else
         eone = ga_ddot(g_tmp,g_cdens)
       endif
c
c Two-electron component of the AO Fock matrices
c
       iv_dens(1) = g_cdens
       iv_fock(1) = g_ifock
       call fock_2e( geom, basis, nset, jfac, kfac, tol2e, oskel,
     $               iv_dens, iv_fock, .false. )
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
       if (oact) then
         e2aa = mcscf_act_ener( nbf, nclosed, nact, dm2, g_coul )
       else
         e2aa = 0.d0
       endif
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
c  Inactive core energy
c
       call ga_dadd( 1.d0, g_tmp, 1.d0, g_ifock, g_adens )
       call two_index_transf( g_adens, g_movecs, g_movecs,
     $                       g_tmp, g_ifock )
       ecore = 0.d0
       do i=ga_nodeid()+1,nclosed,ga_nnodes()
         call ga_get(g_ifock,i,i,i,i,xx,1)
         ecore = ecore + 2.d0*xx
       enddo
       call ga_sync()
       call ga_dgop(msg_mcscf_etrace,ecore,1,'+')
       ecore = ecore - 0.5d0*e2ii
c
c  Clean up
c
       if (.not.ga_destroy(g_tmp))
     $   call errquit('mcscf_etr: cannot destroy',0, GA_ERR)
       if (.not.ga_destroy(g_cdens))
     $   call errquit('mcscf_etr: cannot destroy',0, GA_ERR)
       if (.not.ga_destroy(g_adens))
     $   call errquit('mcscf_etr: cannot destroy',0, GA_ERR)
       if (.not.ga_destroy(g_ifock))
     $   call errquit('mcscf_etr: cannot destroy',0, GA_ERR)
       return
       end






c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
c
c          DEBUGGING VERSION
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
c$$$       subroutine mcscf_etrace00( geom, basis, nbf, nclosed, nact,
c$$$     $                            oskel, tol2e, dm1, dm2, g_movecs,
c$$$     $                            g_coul, eone, e2ii, e2ia, e2aa )
c$$$       implicit none
c$$$#include "errquit.fh"
c$$$#include "mafdecls.fh"
c$$$#include "global.fh"
c$$$#include "rtdb.fh"
c$$$#include "bas.fh"
c$$$#include "geom.fh"
c$$$       integer geom, basis                            ! [input] Geometry and basis handles
c$$$       integer nbf                                    ! [input] Number of basis functions
c$$$       integer nclosed                                ! [input] Number of closed shells
c$$$       integer nact                                   ! [input] Number of open shells
c$$$       logical oskel                                  ! [input] Symmetry toggle
c$$$       double precision tol2e                         ! [input] Integral tolerance
c$$$       double precision dm1(nact,nact)                ! [input] Active space 1PDM
c$$$       double precision dm2(nact,nact,nact,nact)      ! [input] Active space 2PDM
c$$$       integer g_movecs                               ! [input] MO coefficients
c$$$       integer g_coul                                 ! [input] Coulomb integrals
c$$$       double precision eone                          ! [output] Energy components
c$$$       double precision e2ii, e2ia, e2aa
c$$$c
c$$$c
c$$$       integer nset
c$$$       parameter(nset=1)
c$$$       integer g_ifock, g_cdens, g_adens, g_tmp
c$$$       integer iv_dens(2), iv_fock(2)
c$$$       double precision jfac(2), kfac(2)
c$$$       integer ga_create_atom_blocked
c$$$       double precision ga_trace_diag, mcscf_act_ener
c$$$       external ga_create_atom_blocked, ga_trace_diag, mcscf_act_ener
c$$$       data jfac/1.0d0, 1.d0/, kfac/-0.5d0, -0.5d0/
c$$$c
c$$$c
c$$$c
c$$$       g_tmp = ga_create_atom_blocked(geom, basis, 'temp1')
c$$$       g_ifock = ga_create_atom_blocked(geom, basis, 'fock')
c$$$       g_cdens = ga_create_atom_blocked(geom, basis, 'closed dens')
c$$$       g_adens = ga_create_atom_blocked(geom, basis, 'active dens')
c$$$       call ga_zero(g_ifock)
c$$$       call mcscf_occ2dens(nbf, nclosed, nact, dm1, g_movecs,
c$$$     $                     g_cdens, g_adens )
c$$$c
c$$$c One-electron component
c$$$c
c$$$       call ga_zero(g_tmp)
c$$$       call int_1e_ga( basis, basis, g_tmp, 'kinetic', oskel)
c$$$       call int_1e_ga( basis, basis, g_tmp, 'potential', oskel)
c$$$       if (oskel) call sym_symmetrize(geom, basis, .false., g_tmp)
c$$$       eone = ga_ddot(g_tmp,g_cdens) + ga_ddot(g_tmp,g_adens)
c$$$c
c$$$c Two-electron component of the AO Fock matrices
c$$$c
c$$$       iv_dens(1) = g_cdens
c$$$       iv_fock(1) = g_ifock
c$$$       call fock_2e( geom, basis, nset, jfac, kfac, tol2e, oskel,
c$$$     $               iv_dens, iv_fock, .false. )
c$$$c
c$$$c  Symmetrize Fock AO components
c$$$c
c$$$       if (oskel)
c$$$     $    call sym_symmetrize(geom, basis, .false., g_ifock)
c$$$       e2ii = ga_ddot(g_ifock,g_cdens)*0.5d0
c$$$       e2ia = ga_ddot(g_ifock,g_adens)
c$$$c
c$$$c  MO integral contribution
c$$$c
c$$$C       e2aa = mcscf_act_ener( nbf, nclosed, nact, dm2, g_coul )
c$$$       e2aa = 0.d0
c$$$c
c$$$c  Clean up
c$$$c
c$$$       if (.not.ga_destroy(g_tmp))
c$$$     $   call errquit('mcscf_etr: cannot destroy',0)
c$$$       if (.not.ga_destroy(g_cdens))
c$$$     $   call errquit('mcscf_etr: cannot destroy',0)
c$$$       if (.not.ga_destroy(g_adens))
c$$$     $   call errquit('mcscf_etr: cannot destroy',0)
c$$$       if (.not.ga_destroy(g_ifock))
c$$$     $   call errquit('mcscf_etr: cannot destroy',0)
c$$$       return
c$$$       end
c$$$       
c$$$
c$$$
