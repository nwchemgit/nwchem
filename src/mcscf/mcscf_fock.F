c
c  %A%
c  $Id: mcscf_fock.F,v 1.4 1996-03-13 22:54:31 rg240 Exp $
c
c Compute the MCSCF MO-Fock matrices;
c inactive, active and generalized
c
c
       subroutine mcscf_fock( geom, basis, nbf, nclosed, nact,
     $                        oskel, tol2e, dm1, dm2, g_movecs,
     $                        g_coul, eone, etwo,
     $                        g_ifock, g_afock, g_gfock )
       implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
#include "mcscfprof.fh"
c
c
c The definition for the *generalized* Fock
c elements is defined as twice the normal Fock elements 
c (see Siegbahn et al.) and therefore the gradient 
c is twice the generalized Fock element. This routine
c follows this convention, i.e.,
c
c                   gen    gen
c          g   = 2(F    - F   )  = 4F
c           pq      pq     qp        pq
c
c
c Note the asymmetry in the generalized Fock,
c
c
c          F   =  0      for p not closed-shell
c           pi
c
c          F   =  0      for p not open-shell
c           pt
c
c
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
       integer g_ifock                                ! [output] Inactive Fock
       integer g_afock                                ! [output] Active Fock
       integer g_gfock                                ! [output] Generalized Fock
c
c
       integer g_cdens, g_adens, g_tmp, g_tmp2, g_tmp3
       double precision e2ii, e2ia, e2ai, e2aa
       integer nset, aoff, aend
       integer iv_dens(2), iv_fock(2)
       double precision jfac(2), kfac(2)
       integer ga_create_atom_blocked
       external ga_create_atom_blocked
       double precision ga_trace_diag
       external ga_trace_diag
       data nset/2/
       data jfac/1.0d0, 1.d0/
       data kfac/-0.5d0, -0.5d0/
c
c
c
       if (omcscfprof) call pstat_on(ps_gfock)
       aoff = nclosed + 1
       aend = nclosed + nact
       call ga_zero(g_ifock)
       call ga_zero(g_afock)
       g_tmp = ga_create_atom_blocked(geom, basis, 'temp1')
       g_tmp2 = ga_create_atom_blocked(geom, basis, 'temp2')
       g_tmp3 = ga_create_atom_blocked(geom, basis, 'temp3')
       g_cdens = ga_create_atom_blocked(geom, basis, 'closed dens')
       g_adens = ga_create_atom_blocked(geom, basis, 'active dens')
       call mcscf_occ2dens(nbf, nclosed, nact, dm1, g_movecs,
     $                     g_cdens, g_adens )
c
c One-electron component
c
       call ga_zero(g_tmp)
       call ga_zero(g_tmp2)
       call ga_zero(g_tmp3)
       call int_1e_ga( basis, basis, g_tmp, 'kinetic', oskel)
       call int_1e_ga( basis, basis, g_tmp, 'potential', oskel)
       if (oskel) call sym_symmetrize(geom, basis, .false., g_tmp)
       eone = ga_ddot(g_tmp,g_cdens) + ga_ddot(g_tmp,g_adens)
c
c Two-electron component of the Fock matrices
c
       iv_dens(1) = g_cdens
       iv_dens(2) = g_adens
       iv_fock(1) = g_tmp2
       iv_fock(2) = g_tmp3
       call fock_2e( geom, basis, nset, jfac, kfac, tol2e, oskel,
     $               iv_dens, iv_fock )
c
c
c  Symmetrize Fock AO components
c
       if (oskel) then
         call sym_symmetrize(geom, basis, .false., g_tmp2)
         call sym_symmetrize(geom, basis, .false., g_tmp3)
       endif
       e2ii = ga_ddot(g_tmp2,g_cdens)
       e2ia = ga_ddot(g_tmp2,g_adens)
       e2ai = ga_ddot(g_tmp3,g_cdens)
       call ga_dadd(1.d0, g_tmp, 1.d0, g_tmp2, g_tmp2 )
c
c  Transform Fock to MO representation
c
       call two_index_transf(g_tmp2, g_movecs, g_movecs,
     $                       g_cdens, g_ifock )
       call two_index_transf(g_tmp3, g_movecs, g_movecs,
     $                       g_cdens, g_afock )
c
c Generalized Fock matrix elements
c               I      A
c    F   =  2 (F   +  F  )
c     iq        iq     iq
c
c
       call ga_zero(g_gfock)
       call ga_dadd_patch( 2.d0, g_ifock, 1, nclosed, 1, nbf,
     $                     2.d0, g_afock, 1, nclosed, 1, nbf,
     $                           g_gfock, 1, nclosed, 1, nbf )
c
c                     I 
c    F   =  sum[ D  .F   ]  +  2.sum[ P    (qu|vx) ]
c     tq          tu  qu               tuvx
c
c
       call ga_zero(g_tmp)
       call ga_zero(g_cdens)
       call mcscf_denfock_trace(nclosed, nact, nbf, dm1,
     $                          g_ifock, g_tmp)
       call mcscf_3index_trace( nbf, nclosed, nact, dm2,
     $                          g_coul, g_cdens )
       e2aa = ga_trace_diag(g_cdens)                            ! <--- NB. Fii and Faa are zero in this
       call ga_dadd_patch( 1.d0, g_tmp, aoff, aend, 1, nbf,     !          case, only Ftt is non-zero
     $                     1.d0, g_cdens, aoff, aend, 1, nbf,
     $                           g_gfock, aoff, aend, 1, nbf )
c
c  Sum 2e-energy components
c
       etwo = (e2ii + e2ia + e2ai + e2aa)*0.5d0

c$$$       WRITE(6,881) EONE,(E2II*0.5D0),E2IA,E2AA
c$$$ 881   FORMAT(//,5X,'EONE:',F20.14,
c$$$     $        /,5X,'EII:',1X,F20.14,
c$$$     $        /,5X,'EIA:',1X,F20.14,
c$$$     $        /,5X,'EAA:',1X,F20.14)
c
c  Clean up
c
       if (.not.ga_destroy(g_tmp))
     $   call errquit('mcscf_fock: cannot destroy',0)
       if (.not.ga_destroy(g_tmp2))
     $   call errquit('mcscf_fock: cannot destroy',0)
       if (.not.ga_destroy(g_tmp3))
     $   call errquit('mcscf_fock: cannot destroy',0)
       if (.not.ga_destroy(g_cdens))
     $   call errquit('mcscf_fock: cannot destroy',0)
       if (.not.ga_destroy(g_adens))
     $   call errquit('mcscf_fock: cannot destroy',0)
       if (omcscfprof) call pstat_off(ps_gfock)
       return
       end






       subroutine mcscf_denfock_trace(nclosed,nact,nbf,dm1,g_if,g_a)
       implicit none
#include "mafdecls.fh"
#include "global.fh"
       integer nclosed, nact, nbf
       double precision dm1(nact,nact)
       integer g_if, g_a
c
       integer l_ff, k_ff
       integer clo,chi,msize
       integer k_if,ld
       logical check_moblocked 
       external check_moblocked
c
       if (.not.check_moblocked(g_if,nbf,clo,chi))
     $      call errquit('mcscf_denfock_trace: wrong distrib.',0)
       if ((clo.gt.0).and.(chi.gt.0)) then
         msize = (chi-clo+1)*nact
         if (.not.ma_push_get(MT_DBL,msize,'ff',l_ff,k_ff))
     $     call errquit('mcscf_denfock_trace: cannot allocate',0)
         call dfill(msize,0.d0,dbl_mb(k_ff),1)
         call ga_access(g_if,1,nbf,clo,chi,k_if,ld)
         call mcscf_trace_xyz( nbf,nclosed,nact,clo,chi,dbl_mb(k_if),
     $                         dm1,dbl_mb(k_ff))
         call ga_release(g_if,1,nbf,clo,chi)
         call ga_acc(g_a,nclosed+1,nclosed+nact,clo,chi,dbl_mb(k_ff),
     $               nact,1.d0)
         if (.not.ma_pop_stack(l_ff))
     $     call errquit('mcscf_denfock_trace: cannot pop stack',0)
       endif
       call ga_sync()
       return
       end









       subroutine mcscf_trace_xyz(nbf,nclosed,nact,clo,chi,inf,
     $                            dm1,ff)
       implicit none
       integer nbf,nclosed,nact
       integer clo,chi
       double precision inf(nbf,clo:chi)
       double precision dm1(nact,nact)
       double precision ff(nact,clo:chi)
       integer clen,offset
       
       clen = chi - clo + 1
       offset = nclosed + 1
       call dgemm('n','n',nact,clen,nact,1.d0,dm1,nact,inf(offset,clo),
     $            nbf,0.d0,ff,nact)
       return
       end







c
c  Take this trace
c
c      F    =  sum     Gamma    (qu|vx)
c       tq        uvx       tuvx
c
c                                vx
c           =  sum  sum  Gamma  J        (1st summation is local,
c                 vx   u      tu qu       2nd, 3rd via ga_acc)
c
c
       subroutine mcscf_3index_trace( nbf, nclosed, nact, dm2,
     $                                g_coul, g_a )
       implicit none
#include "mafdecls.fh"
#include "global.fh"
       integer nbf,nclosed,nact
       double precision dm2(nact,nact,nact,nact)
       integer g_coul
       integer g_a
c
       integer msize, nn, clo, chi, v, x, vx
       integer k_ff, l_ff, ld, k_c
       logical ga_check_JKblocked
       external ga_check_JKblocked
c
c
c
       msize = nbf*nact
       nn = nbf*nbf
       if (.not.ma_push_get( MT_DBL, msize, 'ff', l_ff, k_ff ))
     $   call errquit('mcscf_3index_trace: cannot allocate',0)
       call dfill(msize,0.d0,dbl_mb(k_ff),1)
       if (.not.ga_check_JKblocked( g_coul, nact, nbf, clo, chi ))
     $   call errquit('mcscf_3index_trace: wrong distrib operator',0)
c
c   Loop over (vx), do summation if local
c
       do v=1,nact
         do x=1,v
           vx = (v*(v-1))/2 + x
           if ((vx.ge.clo).and.(vx.le.chi)) then
             call ga_access( g_coul, 1, nn, vx, vx, k_c, ld )
             call mcscf_trace_pqr( nbf, nclosed, nact, v, x,
     $                             dbl_mb(k_c), dm2, dbl_mb(k_ff) )
             call ga_release( g_coul, 1, nn, vx, vx )
           endif
         enddo
       enddo
       call ga_acc( g_a, nclosed+1, nclosed+nact, 1, nbf,
     $              dbl_mb(k_ff), nact, 1.d0 )
       if (.not.ma_pop_stack(l_ff))
     $   call errquit('mcscf_3index_trace: cannot pop stack',0)
       call ga_sync()
       return
       end







       subroutine mcscf_trace_pqr( nbf, nclosed, nact, v, x, c,
     $                             dm2, ff )
       implicit none
       integer nbf, nclosed, nact, v, x
       double precision c(nbf,nbf)
       double precision dm2(nact,nact,nact,nact)
       double precision ff(nact,nbf)
       integer offset
       double precision scale

       offset = nclosed + 1
       scale = 2.d0
       if (v.eq.x) scale = 1.d0
       call dgemm( 'n', 'n', nact, nbf, nact, scale, dm2(1,1,x,v),
     $             nact, c(offset,1), nbf, 1.d0, ff, nact )
       call dgemm( 'n', 't', nact, nbf, nact, scale, dm2(1,1,v,x),
     $             nact, c(1,offset), nbf, 1.d0, ff, nact )
       return
       end






c
c
c   Return  (sum over active diagonals of above quantity)
c   This is the active-active energy component
c
c      sum    Gamma    (tu|vx)
c         tuvx     tuvx
c
c

       double precision function mcscf_act_ener( nbf, nclosed, nact,
     $                                           dm2, g_coul )
       implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
       integer nbf,nclosed,nact
       double precision dm2(nact,nact,nact,nact)
       integer g_coul
c
       integer nn, clo, chi, v, x, vx, ld, k_c
       double precision xx, s
       logical ga_check_JKblocked
       double precision mcscf_trace_pq
       external ga_check_JKblocked
       external mcscf_trace_pq
c
c
c
       nn = nbf*nbf
       if (.not.ga_check_JKblocked( g_coul, nact, nbf, clo, chi ))
     $   call errquit('mcscf_3index_trace: wrong distrib operator',0)
c
c   Loop over (vx), do summation if local
c
       xx = 0.d0
       do v=1,nact
         do x=1,v
           vx = (v*(v-1))/2 + x
           s = 1.d0
           if (v.eq.x) s = 0.5d0
           if ((vx.ge.clo).and.(vx.le.chi)) then
             call ga_access( g_coul, 1, nn, vx, vx, k_c, ld )
             xx = xx + s*mcscf_trace_pq( nbf, nclosed, nact,
     $                                   dbl_mb(k_c), dm2(1,1,v,x))
             call ga_release( g_coul, 1, nn, vx, vx )
           endif
         enddo
       enddo
       call ga_sync()
       call ga_dgop(Msg_Act_Energy,xx,1,'+')
       mcscf_act_ener = xx
       return
       end












      double precision function mcscf_trace_pq( nbf, nclo, nact, c, d )
      implicit none
      integer nbf, nclo, nact
      double precision c(nbf,nbf), d(nact,nact)
      double precision xx
      integer i,j

      xx = 0.d0
      do i=1,nact
        do j=1,nact
          xx = xx + c(nclo+j,nclo+i)*(d(j,i)+d(i,j))
        enddo
      enddo
      mcscf_trace_pq = xx
      return
      end


      




