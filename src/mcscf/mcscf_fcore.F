c
c        
c  Form
c
c       F    = h                                core Hamiltonian
c        ij     ij
c
c                              rs    rs
c       F    =  h   +  sum ( 2J   - K   )       inactive Fock element
c        rs      rs       i    ii    ii
c
c
       subroutine mcscf_fcore( basis, nbf, nclosed, nact,
     $                         g_movecs, g_coul, g_exch, g_fcore )
*
* $Id$
*
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"       
#include "global.fh"
#include "util.fh"
       integer basis
       integer nbf
       integer nclosed
       integer nact
       integer g_movecs
       integer g_coul
       integer g_exch
       integer g_fcore
       integer nn, t, u, tu, jlo, jhi, k_j, k_k, ld1
       logical ga_check_JKblocked
       external ga_check_JKblocked
c
c  1e Hamiltononian
c
       if (util_print('fcore', print_never)) then
          write(6,*) ' MOVECS in fcore '
          call ga_print(g_movecs)
       endif
       call moints_1e( nbf, basis, g_movecs, g_fcore)
       if (util_print('fcore', print_never))  call ga_print(g_fcore)
c
c  Inactive-active Coulomb interation
c
       nn = nbf*nbf
       if (.not.ga_check_JKblocked(g_coul,nact,nbf,jlo,jhi))
     $      call errquit('mcscf_fcore: wrong distrib. for Coulomb',0,
     &       GA_ERR)
       do t=1,nact
         do u=1,t
           tu = ((t-1)*t)/2 + u
           if ((tu.ge.jlo).and.(tu.le.jhi)) then
             call ga_access(g_coul,1,nn,tu,tu,k_j,ld1)
             call mcscf_fcore01( nbf, nclosed, nact, t, u, 2.d0,
     $                           dbl_mb(k_j), g_fcore)
             call ga_release(g_coul,1,nn,tu,tu)
           endif
         enddo
       enddo
c
c  Inactive-active Exchange interaction
c
       if (.not.ga_check_JKblocked(g_exch,nact,nbf,jlo,jhi))
     $      call errquit('mcscf_fcore: wrong distrib. for exchange',0,
     &       GA_ERR)
       do t=1,nact
         do u=1,t
           tu = ((t-1)*t)/2 + u
           if ((tu.ge.jlo).and.(tu.le.jhi)) then
             call ga_access(g_exch,1,nn,tu,tu,k_k,ld1)
             call mcscf_fcore01( nbf, nclosed, nact, t, u, -1.d0,
     $                           dbl_mb(k_k), g_fcore)
             call ga_release(g_exch,1,nn,tu,tu)
           endif
         enddo
       enddo
c
c  Properly synchronize
c
       call ga_sync()
c
c  Done
c
       return
       end






       subroutine mcscf_fcore01( nbf, nclosed, nact, t, u, fact,
     $                           xmo, g_fcore )
       implicit none
#include "global.fh"
       integer nbf
       integer nclosed
       integer nact
       integer t, u
       double precision fact
       double precision xmo(nbf,nbf)
       integer g_fcore
       double precision xx
       integer i, tt, uu

       xx = 0.d0
       do i=1,nclosed
         xx = xx + xmo(i,i)
       enddo
       tt = nclosed + t
       uu = nclosed + u
       call ga_acc(g_fcore, tt, tt, uu, uu, xx, 1, fact )
       if (t.ne.u) call ga_acc(g_fcore, uu, uu, tt, tt, xx, 1, fact )
       return
       end

         
           
           






           


c
c  Evaluate inactive Fock piece separately 
c  (see mcscf_fock)
c
c
       subroutine mcscf_ifock( geom, basis, nbf, nclosed, nact,
     $                         oskel, tol2e, g_movecs, eone, etwo,
     $                         ecore, g_ifock )
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
       logical oskel                                  ! [input] Symmetry toggle
       double precision tol2e                         ! [input] Integral tolerance
       integer g_movecs                               ! [input] MO coefficients
       double precision eone, etwo, ecore             ! [output] Energy components
       integer g_ifock                                ! [output] Inactive Fock matrix
c
       integer nset
       parameter(nset=1)
c
       integer g_cdens, g_tmp
       double precision e2ii
       integer i
       double precision xx
       integer iv_dens(2), iv_fock(2)
       double precision jfac(2), kfac(2)
       integer ga_create_atom_blocked
       external ga_create_atom_blocked
       data jfac/1.0d0, 1.d0/, kfac/-0.5d0, -0.5d0/
c
c
c
       g_tmp = ga_create_atom_blocked(geom, basis, 'temp1')
       g_cdens = ga_create_atom_blocked(geom, basis, 'closed dens')
       call ga_dgemm( 'n', 't', nbf, nbf, nclosed, 2.d0,
     $                g_movecs, g_movecs, 0.d0, g_cdens )
c
c One-electron component
c
       call ga_zero(g_ifock)
       call int_1e_ga( basis, basis, g_ifock, 'kinetic', oskel)
       call int_1e_ga( basis, basis, g_ifock, 'potential', oskel)
       if (oskel) call sym_symmetrize(geom, basis, .false., g_ifock)
       eone = ga_ddot(g_ifock,g_cdens)
c
c Two-electron component of the AO Fock matrices
c
       call ga_zero(g_tmp)
       iv_dens(1) = g_cdens
       iv_fock(1) = g_tmp
       call fock_2e( geom, basis, nset, jfac, kfac, tol2e, oskel,
     $               iv_dens, iv_fock, .false. )
c
c  Symmetrize Fock AO components
c
       if (oskel)
     $    call sym_symmetrize(geom, basis, .false., g_tmp )
       e2ii = ga_ddot(g_tmp,g_cdens)
c
c  MO integral contribution
c
       etwo = e2ii*0.5d0
       if ((ga_nodeid().eq.0).and.
     $     (util_print('energy trace',print_debug))) then
         write(6,911) e2ii*0.5d0
 911     format(  'energy components  eii: ',f12.6)
       endif
c
c  Inactive core energy
c
       call ga_dadd( 1.d0, g_tmp, 1.d0, g_ifock, g_cdens )
       call two_index_transf( g_cdens, g_movecs, g_movecs,
     $                       g_tmp, g_ifock )
       ecore = 0.d0
       do i=ga_nodeid()+1,nclosed,ga_nnodes()
         call ga_get(g_ifock,i,i,i,i,xx,1)
         ecore = ecore + 2.d0*xx
       enddo
       call ga_sync()
       call ga_dgop(msg_mcscf_ifocktrace,ecore,1,'+')
       ecore = ecore - 0.5d0*e2ii
c
c  Clean up
c
       if (.not.ga_destroy(g_tmp))
     $   call errquit('mcscf_etr: cannot destroy',0, GA_ERR)
       if (.not.ga_destroy(g_cdens))
     $   call errquit('mcscf_etr: cannot destroy',0, GA_ERR)
       return
       end
