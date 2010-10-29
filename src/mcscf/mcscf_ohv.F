       subroutine mcscf_hessv( geom, basis, nbf, nclosed, nact,
     $                         tol2e, oskel, pflg, lshift, dm1, dm2,
     $                         g_movecs, g_ifock, g_afock, g_gfock,
     $                         g_coul, g_exch, g_x, g_ax )
*
* $Id$
*
       implicit none
#include "errquit.fh"
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
       integer pflg                                 ! [input] Hessian-type
       double precision lshift                      ! [input] Level-shift
       double precision dm1(nact,nact)              ! [input] 1PDM
       double precision dm2(nact,nact,nact,nact)    ! [input] 2PDM
       integer g_movecs                             ! [input] MO coefficients
       integer g_ifock                              ! [input] Inactive Fock matrix
       integer g_afock                              ! [input] Active Fock matrix
       integer g_gfock                              ! [input] Generalized Fock matrix
       integer g_coul                               ! [input] Coulomb operator
       integer g_exch                               ! [input] Exchange operator
       integer g_x                                  ! [input] Argument vector
       integer g_ax                                 ! [output] Hessian vector product
c
c
c
       integer g_tmp1, g_tmp2, g_tmp3               ! [input] Temporaries (nbf * nbf)
c
c
c
       if (omcscfprof) call pstat_on(ps_hv)
       if (.not.ga_duplicate(g_afock, g_tmp1, 'temp1'))
     $      call errquit('mcscf_hessv: cannot duplicate',0, GA_ERR)
       if (.not.ga_duplicate(g_afock, g_tmp2, 'temp2'))
     $      call errquit('mcscf_hessv: cannot duplicate',0, GA_ERR)
       if (.not.ga_duplicate(g_afock, g_tmp3, 'temp3'))
     $      call errquit('mcscf_hessv: cannot duplicate',0, GA_ERR)
c
c
c
       call ga_zero(g_ax)
       call mcscf_hessv_1e( geom, basis, nbf, nclosed, nact,
     $                      lshift, dm1, g_ifock, g_afock, g_gfock,
     $                      g_tmp1, g_tmp2, g_tmp3,
     $                      g_x, g_ax )
c
c
c
       if (pflg.gt.1) then
         call mcscf_hessv_2e( geom, basis, nbf, nclosed, nact, tol2e,
     $                        oskel, dm1, dm2, g_movecs,
     $                        g_coul, g_exch, g_tmp1, g_tmp2, g_tmp3,
     $                        g_x, g_ax )
       endif
c
c
c

       if (.not.ga_destroy(g_tmp1))
     $      call errquit('mcscf_hessv: cannot destroy temp1',0, GA_ERR)
       if (.not.ga_destroy(g_tmp2))
     $      call errquit('mcscf_hessv: cannot destroy temp2',0, GA_ERR)
       if (.not.ga_destroy(g_tmp3))
     $      call errquit('mcscf_hessv: cannot destroy temp3',0, GA_ERR)
c
c
c
       if (omcscfprof) call pstat_off(ps_hv)
       return
       end
