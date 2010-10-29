c
c  $Id$
c
c Return the MCSCF 2e Hessian-vector product
c
c
       subroutine mcscf_hessv_2e( geom, basis, nbf, nclosed, nact,
     $                            tol2e, oskel, dm1, dm2, g_movecs,
     $                            g_coul, g_exch, g_tmp1, g_tmp2,
     $                            g_tmp3, g_x, g_ax )
       implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "mcscfprof.fh"
c
c
       integer geom, basis                          ! [input] Handles
       integer nbf                                  ! [input] Basis functions
       integer nclosed                              ! [input] Closed shells
       integer nact                                 ! [input] Active shells
       double precision tol2e                       ! [input] Integral tolerance
       logical oskel                                ! [input] Symmetry toggle
       double precision dm1(nact,nact)              ! [input] 1PDM
       double precision dm2(nact,nact,nact,nact)    ! [input] 2PDM
       integer g_movecs                             ! [input] MO coefficients
       integer g_coul                               ! [input] Coulomb operator
       integer g_exch                               ! [input] Exchange operator
       integer g_tmp1, g_tmp2, g_tmp3               ! [input] Temporaries (nbf * nbf)
       integer g_x                                  ! [input] Argument vector
       integer g_ax                                 ! [output] Hessian vector product
c
       integer nvir, vlen, voff, aoff, aend
c
c
       if (omcscfprof) call pstat_on(ps_hv2e)
       nvir = nbf - nclosed - nact
       vlen = (nclosed+nact)*nvir + nclosed*nact
       voff = nclosed + nact + 1
       aoff = nclosed + 1
       aend = nclosed + nact
c
c
c
       call ga_zero( g_tmp3 )
       call ga_zero( g_tmp1 )
       call mcscf_vec2mat( nbf, nclosed, nact, 1.d0, g_x, 0.d0, g_tmp1 )
       call mcscf_hessv_2e_ao( geom, basis, nbf, nclosed, nact, tol2e,
     $                         oskel, dm1, g_movecs, g_tmp1, g_tmp2,
     $                         g_tmp3 )
c
c
c
       call mcscf_hessv_2e_mo( nbf, nclosed, nact, dm2,
     $                         g_coul, g_exch, g_tmp1, g_tmp3 )
c
c
c
       call mcscf_mat2vec(nbf, nclosed, nact, 1.d0, g_tmp3, 1.d0, g_ax )
c
c
c
       if (omcscfprof) call pstat_off(ps_hv2e)
       return
       end
