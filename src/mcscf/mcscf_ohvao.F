*
* $Id$
*
       subroutine mcscf_hessv_2e_ao( geom, basis, nbf, nclosed, nact,
     $                               tol2e, oskel, dm1, g_movecs,
     $                               g_x, g_tmp, g_ax )
#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "mcscfprof.fh"
c
c
       integer geom, basis                          ! [input] Handles
       integer nbf                                  ! [input] Basis functions
       integer nclosed                              ! [input] Closed shells
       integer nact                                 ! [input] Active shells
       double precision tol2e                       ! [input] Integral tolerance
       logical oskel                                ! [input] Symmetry selection
       double precision dm1(nact,nact)              ! [input] 1PDM
       integer g_movecs                             ! [input] MO coefficients
       integer g_tmp                                ! [input] Temporary  (nbf * nbf)
       integer g_x                                  ! [input] Argument parameter matrix 
       integer g_ax                                 ! [output] Hessian product (in matrix rep)
c
c
c
       integer g_f1, g_f2, g_f3, g_f4, g_f5, g_f6   ! Atom-blocked Fock matrices
       integer g_d1, g_d2, g_d3, g_d4, g_d5, g_d6   ! Atom-blocked densities
       integer nvir, vlen, voff, aoff, aend, nn
       integer g_tmp1, g_tmp2
c
c
       integer nsets
       parameter(nsets=6)
       integer iv_dens(nsets), iv_fock(nsets)
       double precision jfac(nsets),kfac(nsets)
c
c
       integer ga_create_atom_blocked
       external ga_create_atom_blocked
c
       data jfac/6*1.d0/
       data kfac/6*-0.5d0/
c
c
       if (omcscfprof) call pstat_on(ps_hv2eao)
       nvir = nbf - nclosed - nact
       vlen = (nclosed+nact)*nvir + nclosed*nact
       voff = nclosed + nact + 1
       aoff = nclosed + 1
       aend = nclosed + nact
c
c
       if (.not.ga_duplicate(g_tmp,g_tmp1,'temp1'))
     $      call errquit('mcscf_hessv_2e_ao: cannot duplicate',0,
     &       GA_ERR)
       if (.not.ga_duplicate(g_tmp,g_tmp2,'temp2'))
     $      call errquit('mcscf_hessv_2e_ao: cannot duplicate',0,
     &       GA_ERR)
       call ga_zero(g_tmp1)
       call ga_zero(g_tmp2)
       call ga_zero(g_tmp)
       if (ga_nodeid().eq.0)
     $      call ga_put(g_tmp,aoff,aend,aoff,aend,dm1,nact)
       call ga_sync()
       g_d1 = ga_create_atom_blocked( geom, basis, 'Density 1')
       g_d2 = ga_create_atom_blocked( geom, basis, 'Density 2')
       g_d3 = ga_create_atom_blocked( geom, basis, 'Density 3')
       g_d4 = ga_create_atom_blocked( geom, basis, 'Density 4')
       g_d5 = ga_create_atom_blocked( geom, basis, 'Density 5')
       g_d6 = ga_create_atom_blocked( geom, basis, 'Density 6')
       g_f1 = ga_create_atom_blocked( geom, basis, 'Fock 1')
       g_f2 = ga_create_atom_blocked( geom, basis, 'Fock 2')
       g_f3 = ga_create_atom_blocked( geom, basis, 'Fock 3')
       g_f4 = ga_create_atom_blocked( geom, basis, 'Fock 4')
       g_f5 = ga_create_atom_blocked( geom, basis, 'Fock 5')
       g_f6 = ga_create_atom_blocked( geom, basis, 'Fock 6')
c
c Inactive-Virtual "density"
c                                               iv     iv t
c                                              D  = C.X .C
c
       call ga_matmul_patch( 'n', 't', 1.d0, 0.d0,
     $                       g_x, voff, nbf, 1, nclosed,
     $                       g_movecs, 1, nclosed, 1, nbf,
     $                       g_tmp1, voff, nbf, 1, nbf )
       call ga_matmul_patch( 'n', 'n', 1.d0, 0.d0,
     $                       g_movecs, 1, nbf, voff, nbf,
     $                       g_tmp1, voff, nbf, 1, nbf,
     $                       g_tmp2, 1, nbf, 1, nbf )
       call ga_symmetrize(g_tmp2)
       call ga_copy( g_tmp2, g_d1 )
       call ga_copy( g_d1, g_d2 )
       call ga_copy( g_d1, g_d3 )
c
c Inactive-Virtual, Active-Virtual "density"    
c                                               ,av     av   t
c                                              D   = C.X .d.C
c
       call ga_matmul_patch( 'n', 'n', 1.d0, 0.d0,
     $                       g_x, voff, nbf, aoff, aend,
     $                       g_tmp, aoff, aend, aoff, aend,
     $                       g_tmp1, voff, nbf, aoff, aend )
       call ga_matmul_patch( 'n', 't', 1.d0, 0.d0,
     $                       g_tmp1, voff, nbf, aoff, aend,
     $                       g_movecs, aoff, aend, 1, nbf,
     $                       g_tmp2, voff, nbf, 1, nbf )
       call ga_matmul_patch( 'n', 'n', 1.d0, 0.d0,
     $                       g_movecs, 1, nbf, voff, nbf,
     $                       g_tmp2, voff, nbf, 1, nbf,
     $                       g_tmp1, 1, nbf, 1, nbf )
       call ga_symmetrize(g_tmp1)
       call ga_copy( g_tmp1, g_d4 )
       call ga_dadd( 1.d0, g_d4, 2.d0, g_d1, g_d1 )                  
c
c Inactive-Virtual, Inactive-Active "density"
c                                               ,ia   t      ia
c                                              D   = C.(2-d)X .C
c
       call ga_zero(g_tmp2)
       call diagfill_patch(g_tmp2, 2.d0, aoff, aend )
       call ga_dadd(1.d0, g_tmp2, -1.d0, g_tmp, g_tmp2 )
       call ga_zero(g_tmp1)
       call ga_matmul_patch('n', 'n', 1.d0, 0.d0,
     $                       g_tmp2, aoff, aend, aoff, aend,
     $                       g_x, aoff, aend, 1, nclosed,
     $                       g_tmp1, aoff, aend, 1, nclosed )
       call ga_matmul_patch('n', 't', 1.d0, 0.d0,
     $                       g_tmp1, aoff, aend, 1, nclosed,
     $                       g_movecs, 1, nclosed, 1, nbf,
     $                       g_tmp2, aoff, aend, 1, nbf )
       call ga_matmul_patch('n', 'n', 1.d0, 0.d0,
     $                       g_movecs, 1, nbf, aoff, aend,
     $                       g_tmp2, aoff, aend, 1, nbf,
     $                       g_tmp1, 1, nbf, 1, nbf )
       call ga_symmetrize(g_tmp1)
       call ga_copy( g_tmp1, g_d5 )
       call ga_dadd( 1.d0, g_d5, 1.d0, g_d1, g_d1 )
c
c Inactive-Active density
c                                               ia    t ia
c                                              D   = C.X .C
c
       call ga_zero(g_tmp1)
       call ga_zero(g_tmp2)
       call ga_matmul_patch('n', 't', 1.d0, 0.d0,
     $                      g_x, aoff, aend, 1, nclosed,
     $                      g_movecs, 1, nclosed, 1, nbf,
     $                      g_tmp1, aoff, aend, 1, nbf )
       call ga_matmul_patch('n', 'n', 1.d0, 0.d0,
     $                      g_movecs, 1, nbf, aoff, aend,
     $                      g_tmp1, aoff, aend, 1, nbf,
     $                      g_tmp2, 1, nbf, 1, nbf )
       call ga_symmetrize(g_tmp2)
       call ga_copy( g_tmp2, g_d5 )
       call ga_dadd( 1.d0, g_d5, 1.d0, g_d2, g_d2 )
c
c
c CA, CA density
c                                               ,,ia    t       ia
c                                              D    = -C.(1-d).X .C
c
       call ga_zero(g_tmp1)
       call diagfill_patch(g_tmp1, 1.d0, aoff, aend )                  
       call ga_dadd( 1.d0, g_tmp1, -1.d0, g_tmp, g_tmp1)               
       call ga_matmul_patch( 'n', 'n', 1.d0, 0.d0,
     $                       g_tmp1, aoff, aend, aoff, aend,
     $                       g_x, aoff, aend, 1, nclosed,
     $                       g_tmp2, aoff, aend, 1, nclosed )
       call ga_matmul_patch( 'n', 't', 1.d0, 0.d0,
     $                       g_tmp2, aoff, aend, 1, nclosed, 
     $                       g_movecs, 1, nclosed, 1, nbf,
     $                       g_tmp1, aoff, aend, 1, nbf )
       call ga_matmul_patch( 'n', 'n', 1.d0, 0.d0,
     $                       g_movecs, 1, nbf, aoff, aend,
     $                       g_tmp1, aoff, aend, 1, nbf,
     $                       g_tmp2, 1, nbf, 1, nbf )
       call ga_symmetrize(g_tmp2)
       call ga_copy(g_tmp2, g_d6)
c 
c
c
c  Density summary
c
c              iv   ,av   ,ia                 ,av
c       d1 : 2D  + D   + D              d4 : D
c
c             iv   ia                         ia
c       d2 : D  + D                     d5 : D
c
c             iv                              ,,ia
c       d3 : D                          d6 : D
c
c
       iv_dens(1) = g_d1
       iv_dens(2) = g_d2
       iv_dens(3) = g_d3
       iv_dens(4) = g_d4
       iv_dens(5) = g_d5
       iv_dens(6) = g_d6
c
c Fock build
c
       nn = 6
       call ga_zero(g_f1)
       call ga_zero(g_f2)
       call ga_zero(g_f3)
       call ga_zero(g_f4)
       call ga_zero(g_f5)
       call ga_zero(g_f6)
       iv_fock(1) = g_f1
       iv_fock(2) = g_f2
       iv_fock(3) = g_f3
       iv_fock(4) = g_f4
       iv_fock(5) = g_f5
       iv_fock(6) = g_f6
       call fock_2e( geom, basis, nn, jfac, kfac, tol2e, oskel,
     $               iv_dens, iv_fock, .false. )
c
c  Symmetrize AO Fock matrices
c
       if (oskel) then
         do i=1,nn
           call sym_symmetrize(geom, basis, .false., iv_fock(i))
         enddo
       endif
c
c Inactive-Virtual contribution
c 
c
       call ga_zero(g_tmp2)
       call ga_matmul_patch( 't', 'n', 1.d0, 0.d0,
     $                       g_movecs, voff, nbf, 1, nbf,
     $                       g_f1, 1, nbf, 1, nbf,
     $                       g_tmp1, voff, nbf, 1, nbf )
       call ga_matmul_patch( 'n', 'n', 8.d0, 0.d0,                      ! Where does this factor come from?
     $                       g_tmp1, voff, nbf, 1, nbf,
     $                       g_movecs, 1, nbf, 1, nclosed,
     $                       g_tmp2, voff, nbf, 1, nclosed )
       call ga_dadd_patch( 1.d0, g_tmp2, voff, nbf, 1, nclosed,
     $                     1.d0, g_ax, voff, nbf, 1, nclosed,
     $                           g_ax, voff, nbf, 1, nclosed )
c
c Active-Virtual contribution
c
c
       call ga_matmul_patch( 't', 'n', 1.d0, 0.d0,
     $                       g_movecs, voff, nbf, 1, nbf,
     $                       g_f2, 1, nbf, 1, nbf,
     $                       g_tmp1, voff, nbf, 1, nbf )
       call ga_zero(g_tmp2)
       call ga_matmul_patch( 'n', 'n', 1.d0, 0.d0,
     $                       g_tmp1, voff, nbf, 1, nbf,
     $                       g_movecs, 1, nbf, aoff, aend,
     $                       g_tmp2, voff, nbf, aoff, aend )
       call ga_matmul_patch( 'n', 'n', 8.d0, 0.d0,
     $                       g_tmp2, voff, nbf, aoff, aend,
     $                       g_tmp, aoff, aend, aoff, aend,
     $                       g_tmp1, voff, nbf, aoff, aend )
       call ga_dadd_patch( 1.d0, g_tmp1, voff, nbf, aoff, aend,
     $                     1.d0, g_ax, voff, nbf, aoff, aend,
     $                           g_ax, voff, nbf, aoff, aend )
c
c Inactive-Active contributions
c                                          t   iv
c                                 (2 - d) C.F[D  ].C
c
       call ga_zero(g_tmp1)
       call ga_zero(g_tmp2)
       call ga_matmul_patch( 't', 'n', 1.d0, 0.d0,
     $                       g_movecs, aoff, aend, 1, nbf,
     $                       g_f3, 1, nbf, 1, nbf,
     $                       g_tmp1, aoff, aend, 1, nbf )
       call ga_matmul_patch( 'n', 'n', 1.d0, 0.d0,
     $                       g_tmp1, aoff, aend, 1, nbf,
     $                       g_movecs, 1, nbf, 1, nclosed,
     $                       g_tmp2, aoff, aend, 1, nclosed )
       call ga_zero(g_tmp1)
       call diagfill_patch( g_tmp1, 2.d0, aoff, aend )
       call ga_dadd( 1.d0, g_tmp1, -1.d0, g_tmp, g_tmp1 )
       call ga_matmul_patch( 'n', 'n', 8.d0, 0.d0,
     $                       g_tmp1, aoff, aend, aoff, aend,
     $                       g_tmp2, aoff, aend, 1, nclosed,
     $                       g_tmp1, aoff, aend, 1, nclosed )
c
c                                     t   ,av  
c                                 -d C.F[D   ].C
c
       call ga_zero(g_tmp2)
       call ga_matmul_patch( 't', 'n', 1.d0, 0.d0,
     $                       g_movecs, aoff, aend, 1, nbf,
     $                       g_f4, 1, nbf, 1, nbf,
     $                       g_tmp2, aoff, aend, 1, nbf )
       call ga_matmul_patch( 'n', 'n', 8.d0, 1.d0,
     $                       g_tmp2, aoff, aend, 1, nbf,
     $                       g_movecs, 1, nbf, 1, nclosed,
     $                       g_tmp1, aoff, aend, 1, nclosed )
       call ga_dadd_patch( 1.d0, g_tmp1, aoff, aend, 1, nclosed,
     $                     1.d0, g_ax, aoff, aend, 1, nclosed,
     $                           g_ax, aoff, aend, 1, nclosed )
c
c
c Inactive-Active, Inactive-Active contribution
c    (note this last section has zero contribution
c    in ROHF theory and has not been debugged against the ROHF
c    Hessian product)
c
c
       call ga_matmul_patch( 't', 'n', 1.d0, 0.d0,                         
     $                       g_movecs, aoff, aend, 1, nbf,                
     $                       g_f6, 1, nbf, 1, nbf,                        
     $                       g_tmp1, aoff, aend, 1, nbf )                 
       call ga_matmul_patch( 'n', 'n', 8.d0, 0.d0,                        
     $                       g_tmp1, aoff, aend, 1, nbf,
     $                       g_movecs, 1, nbf, 1, nclosed,
     $                       g_tmp2, aoff, aend, 1, nclosed )
       call ga_dadd_patch( 1.d0, g_tmp2, aoff, aend, 1, nclosed,
     $                     1.d0, g_ax, aoff, aend, 1, nclosed,
     $                           g_ax, aoff, aend, 1, nclosed )
c
c
c
       call ga_matmul_patch( 't', 'n', 1.d0, 0.d0,
     $                       g_movecs, aoff, aend, 1, nbf,
     $                       g_f5, 1, nbf, 1, nbf,
     $                       g_tmp1, aoff, aend, 1, nbf )
       call ga_matmul_patch( 'n', 'n', 8.d0, 0.d0,
     $                       g_tmp1, aoff, aend, 1, nbf,
     $                       g_movecs, 1, nbf, 1, nclosed,
     $                       g_tmp2, aoff, aend, 1, nclosed )
       call diagfill_patch( g_tmp1, 1.d0, aoff, aend )
       call ga_dadd( 1.d0, g_tmp1, -1.d0, g_tmp, g_tmp1 )
       call ga_matmul_patch( 'n', 'n', 1.d0, 0.d0,
     $                       g_tmp1, aoff, aend, aoff, aend,
     $                       g_tmp2, aoff, aend, 1, nclosed,
     $                       g_tmp1, aoff, aend, 1, nclosed )
       call ga_dadd_patch( 1.d0, g_tmp1, aoff, aend, 1, nclosed,
     $                     1.d0, g_ax, aoff, aend, 1, nclosed,
     $                           g_ax, aoff, aend, 1, nclosed )
c
c
c
       if (.not.ga_destroy(g_tmp1))
     $      call errquit('mcscf_hessv_2e_ao: cannot destroy',0, GA_ERR)
       if (.not.ga_destroy(g_tmp2))
     $      call errquit('mcscf_hessv_2e_ao: cannot destroy',0, GA_ERR)
       if (.not.ga_destroy(g_d1))
     $      call errquit('mcscf_hessv_2e_ao: cannot destroy',0, GA_ERR)
       if (.not.ga_destroy(g_d2))
     $      call errquit('mcscf_hessv_2e_ao: cannot destroy',0, GA_ERR)
       if (.not.ga_destroy(g_d3))
     $      call errquit('mcscf_hessv_2e_ao: cannot destroy',0, GA_ERR)
       if (.not.ga_destroy(g_d4))
     $      call errquit('mcscf_hessv_2e_ao: cannot destroy',0, GA_ERR)
       if (.not.ga_destroy(g_d5))
     $      call errquit('mcscf_hessv_2e_ao: cannot destroy',0, GA_ERR)
       if (.not.ga_destroy(g_d6))
     $      call errquit('mcscf_hessv_2e_ao: cannot destroy',0, GA_ERR)
       if (.not.ga_destroy(g_f1))
     $      call errquit('mcscf_hessv_2e_ao: cannot destroy',0, GA_ERR)
       if (.not.ga_destroy(g_f2))
     $      call errquit('mcscf_hessv_2e_ao: cannot destroy',0, GA_ERR)
       if (.not.ga_destroy(g_f3))
     $      call errquit('mcscf_hessv_2e_ao: cannot destroy',0, GA_ERR)
       if (.not.ga_destroy(g_f4))
     $      call errquit('mcscf_hessv_2e_ao: cannot destroy',0, GA_ERR)
       if (.not.ga_destroy(g_f5))
     $      call errquit('mcscf_hessv_2e_ao: cannot destroy',0, GA_ERR)
       if (.not.ga_destroy(g_f6))
     $      call errquit('mcscf_hessv_2e_ao: cannot destroy',0, GA_ERR)
c
c
c
       if (omcscfprof) call pstat_off(ps_hv2eao)
       return
       end


