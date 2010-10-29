c
c  Given the unitrary transformation U
c  generate first order response of the gradient, B
c  (Werner's).
c
c  Involves 1-index transformation of J and K
c  and Fock build with first-order change in the
c  density
c
c
      subroutine mcscf_b( geom, basis, nbf, nclosed, nact, orbsym,
     $                    dm1, dm2, oskel, tol2e, bshift, g_vecs0,
     $                    g_ifock, g_afock, g_coul, g_exch, g_u, g_b )
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"      
      integer geom, basis                              ! [input] Handles
      integer nbf, nclosed, nact                       ! [input] Orbital parameters
      integer orbsym(nbf)                              ! [input] Orbital irreps
      double precision dm1(nact,nact)                  ! [input] 1-PDM
      double precision dm2(nact,nact,nact,nact)        ! [input] 2-PDM
      logical oskel                                    ! [input] Symmetry toggle
      double precision tol2e                           ! [input] AO integral tolerance
      double precision bshift                          ! [input] Shift for B-matrix
      integer g_vecs0                                  ! [input] MO vectors at alpha=0
      integer g_ifock                                  ! [input] Inactive Fock
      integer g_afock                                  ! [input] Active Fock
      integer g_coul                                   ! [input] Coulomb operator
      integer g_exch                                   ! [input] Exchange operator
      integer g_u                                      ! [input] Unitrary transformation
      integer g_b                                      ! [output] B matrix (gradient response)
c
      integer g_t, g_tmp, g_tmp2, g_tmp3
      integer i, nocc, aoff, aend, voff
c
c
c
      nocc = nclosed + nact
      aoff = nclosed + 1
      aend = nclosed + nact
      voff = aend + 1
      call ga_zero(g_b)
c
c
c
*ga:1:0
      if (.not.ga_create(MT_DBL, nbf, nbf, 'tmp', nbf, 0, g_tmp))
     $   call errquit('mcscf_b: cannot allocate tmp',0, GA_ERR)
*ga:1:0
      if (.not.ga_create(MT_DBL, nbf, nbf, 'tmp', nbf, 0, g_tmp2))
     $   call errquit('mcscf_b: cannot allocate tmp',0, GA_ERR)
*ga:1:0
      if (.not.ga_create(MT_DBL, nbf, nbf, 'tmp', nbf, 0, g_tmp3))
     $   call errquit('mcscf_b: cannot allocate tmp',0, GA_ERR)
*ga:1:0
      if (.not.ga_create(MT_DBL, nbf, nbf, 'T', nbf, 0, g_t))
     $   call errquit('mcscf_b: cannot allocate T',0, GA_ERR)
c
c  T = U - 1    
c      
      call ga_copy(g_u, g_t)
      do i=ga_nodeid()+1,nbf,ga_nnodes()
        call ga_acc(g_t, i, i, i, i, -1.d0, 1, 1.d0 )
      enddo
      call ga_sync()
c
c  AO contribution
c
      call mcscf_b_ao( geom, basis, nbf, nclosed, nact,
     $                 oskel, tol2e, dm1, g_vecs0, g_t,
     $                 g_b )
c
c  Transformed Fock contributions
c
c            I    A
c  B   = 2((F  + F ).U)
c   pi                 pi
c
      call ga_dadd( 1.d0, g_ifock, 1.d0, g_afock, g_tmp )
      call ga_dgemm( 'n', 'n', nbf, nbf, nbf, 1.d0, g_tmp, g_u,
     $               0.d0, g_tmp2 )
      call ga_dadd_patch( 2.d0, g_tmp2, 1, nbf, 1, nclosed,
     $                    1.d0, g_b,    1, nbf, 1, nclosed,
     $                          g_b,    1, nbf, 1, nclosed )
c
c          I
c  B   = (F .U.dm1)
c   pt             pt
c
c
      call ga_zero(g_tmp)
      if (ga_nodeid().eq.0)
     $  call ga_put(g_tmp, aoff, aend, aoff, aend, dm1, nact )
      call ga_sync()
      call ga_dgemm( 'n', 'n', nbf, nbf, nbf, 1.d0, g_ifock, g_u,
     $               0.d0, g_tmp2 )
      call ga_dgemm( 'n', 'n', nbf, nbf, nbf, 1.d0, g_tmp2, g_tmp,
     $               0.d0, g_tmp3 )
      call ga_dadd_patch( 1.d0, g_tmp3, 1, nbf, aoff, aend,
     $                    1.d0, g_b, 1, nbf, aoff, aend,
     $                          g_b, 1, nbf, aoff, aend )
c
c  Transform and trace J 
c
      call ga_zero( g_tmp )
      call mcscf_J_trace( nbf, nclosed, nact, dm2, g_u, g_coul, g_tmp )
      call ga_dadd(1.d0, g_tmp, 1.d0, g_b, g_b )
c
c  Transform and trace K
c
      call ga_zero( g_tmp )
      call mcscf_K_trace( nbf, nclosed, nact, dm2, g_t, g_exch, g_tmp )
      call ga_dadd(2.d0, g_tmp, 1.d0, g_b, g_b )
c
      if (.not.ga_destroy(g_t))
     $  call errquit('mcscf_b: cannot destroy T',0, GA_ERR)
      if (.not.ga_destroy(g_tmp3))
     $  call errquit('mcscf_b: cannot destroy tmp',0, GA_ERR)
      if (.not.ga_destroy(g_tmp2))
     $  call errquit('mcscf_b: cannot destroy tmp',0, GA_ERR)
      if (.not.ga_destroy(g_tmp))
     $  call errquit('mcscf_b: cannot destroy tmp',0, GA_ERR)
      return
      end










c
c  B matrix - AO integral contributions
c
      subroutine mcscf_b_ao( geom, basis, nbf, nclosed, nact,
     $                       oskel, tol2e, dm1, g_vecs0, g_t,
     $                       g_b )
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer geom
      integer basis
      integer nbf
      integer nclosed
      integer nact
      logical oskel
      double precision tol2e
      double precision dm1(nact,nact)
      integer g_vecs0
      integer g_t
      integer g_b
c
c
c
      integer nfock
      parameter(nfock=2)
      integer g_dens
      integer g_deltad, g_deltadc
      integer g_d1, g_d2, g_f1, g_f2
      integer i
      integer aoff, aend
      integer iv_fock(nfock), iv_dens(nfock)
      double precision jfac(nfock), kfac(nfock)
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c
      data jfac/2*1.d0/
      data kfac/2*-0.5d0/
c
c
      aoff = nclosed + 1
      aend = nclosed + nact
c
c
c
*ga:1:0
      if (.not.ga_create(MT_DBL, nbf, nbf, 'tmp', nbf, 0, g_dens))
     $   call errquit('mcscf_b: cannot allocate tmp',0, GA_ERR)
*ga:1:0
      if (.not.ga_create(MT_DBL, nbf, nbf, 'tmp', nbf, 0, g_deltadc))
     $   call errquit('mcscf_b: cannot allocate dD',0, GA_ERR)
*ga:1:0
      if (.not.ga_create(MT_DBL, nbf, nbf, 'tmp', nbf, 0, g_deltad))
     $   call errquit('mcscf_b: cannot allocate dD',0, GA_ERR)
      g_d1 = ga_create_atom_blocked( geom, basis, 'Density 1')
      g_d2 = ga_create_atom_blocked( geom, basis, 'Density 2')
      g_f1 = ga_create_atom_blocked( geom, basis, 'Fock 1')
      g_f2 = ga_create_atom_blocked( geom, basis, 'Fock 2')
      call ga_zero(g_dens)
c
c   D_ij = 2*delta_ij,   D_tu = dm1_tu
c
      do i=ga_nodeid()+1,nclosed,ga_nnodes()
        call ga_put(g_dens, i, i, i, i, 2.d0, 1 )
      enddo
      call ga_sync()
c
c  Inactive density diff (only double occupation)   There must be a trivial
c                                                   way to do this!
c    c          t
c  dD  = TD + DT
c
      call ga_dgemm( 'n', 'n', nbf, nbf, nbf, 1.d0, g_t, g_dens,
     $               0.d0, g_deltadc )
      call ga_dgemm( 'n', 't', nbf, nbf, nbf, 1.d0, g_dens, g_t,
     $               1.d0, g_deltadc )
c
c
c  Full density diff (with active occupation)
c               t
c  dD  = TD + DT
c
      if (ga_nodeid().eq.0)
     $  call ga_put(g_dens, aoff, aend, aoff, aend, dm1, nact )
      call ga_dgemm( 'n', 'n', nbf, nbf, nbf, 1.d0, g_t, g_dens,
     $               0.d0, g_deltad )
      call ga_dgemm( 'n', 't', nbf, nbf, nbf, 1.d0, g_dens, g_t,
     $               1.d0, g_deltad )
c
c                   t
c   dD(AO) =  C.dD.C
c
c     c           c  t
c   dD(AO) =  C.dD .C
c
c
      call ga_dgemm( 'n', 't', nbf, nbf, nbf, 1.d0, g_deltadc,
     $               g_vecs0, 0.d0, g_f1 )
      call ga_dgemm( 'n', 'n', nbf, nbf, nbf, 1.d0, g_vecs0,
     $               g_f1, 0.d0, g_d1 )
      call ga_dgemm( 'n', 't', nbf, nbf, nbf, 1.d0, g_deltad,
     $               g_vecs0, 0.d0, g_f1 )
      call ga_dgemm( 'n', 'n', nbf, nbf, nbf, 1.d0, g_vecs0,
     $               g_f1, 0.d0, g_d2 )
c
c  Compute Fock matrices
c
      call ga_zero(g_f1)
      call ga_zero(g_f2)
      iv_dens(1) = g_d1
      iv_dens(2) = g_d2
      iv_fock(1) = g_f1
      iv_fock(2) = g_f2
      call fock_2e( geom, basis, nfock, jfac, kfac, tol2e, oskel,
     $              iv_dens, iv_fock, .false. )
      if (oskel) then
        call sym_symmetrize(geom, basis, .false., g_f1)
        call sym_symmetrize(geom, basis, .false., g_f2)
      endif
c      
c  Transform to MO
c
      call two_index_transf(g_f1, g_vecs0, g_vecs0, g_d1,
     $                      g_deltadc )
      call two_index_transf(g_f2, g_vecs0, g_vecs0, g_d1,
     $                      g_deltad )
c
c  B   +=  F[dD]
c   ri          ri
c
      call ga_copy_patch('n', g_deltad, 1, nbf, 1, nclosed,
     $                        g_b,      1, nbf, 1, nclosed )
c
c               c
c  B   += ( F[dD ].D )
c   rt                rt
c
      call ga_dgemm( 'n', 'n', nbf, nbf, nbf, 1.d0, g_deltadc,
     $               g_dens, 0.d0, g_deltad )
      call ga_copy_patch('n', g_deltad, 1, nbf, aoff, aend,
     $                        g_b, 1, nbf, aoff, aend )

c
c  Clean up
c
      if (.not.ga_destroy(g_d1))
     $  call errquit('mcscf_b: cannot destroy AO dens1',0, GA_ERR)
      if (.not.ga_destroy(g_d2))
     $  call errquit('mcscf_b: cannot destroy AO dens2',0, GA_ERR)
      if (.not.ga_destroy(g_f1))
     $  call errquit('mcscf_b: cannot destroy AO fock1',0, GA_ERR)
      if (.not.ga_destroy(g_f2))
     $  call errquit('mcscf_b: cannot destroy AO fock2',0, GA_ERR)
      if (.not.ga_destroy(g_deltad))
     $  call errquit('mcscf_b: cannot destroy deltad',0, GA_ERR)
      if (.not.ga_destroy(g_deltadc))
     $  call errquit('mcscf_b: cannot destroy deltadc',0, GA_ERR)
      if (.not.ga_destroy(g_dens))
     $  call errquit('mcscf_b: cannot destroy dens',0, GA_ERR)
c
c
c
      return
      end







c
c  Do this trace:
c
c              vx          vx
c      B   =  J   U   Gamma  
c       pt     pq  qu      ut
c
c             ~vx      vx
c          =  J   Gamma  
c              pu      ut
c
c  Take advantage of symmetry
c
c          vx      xv
c         J    =  J
c          pq      pq
c
       subroutine mcscf_J_trace( nbf, nclosed, nact, dm2, g_u,
     $                           g_coul, g_b )
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
       integer nbf
       integer nclosed
       integer nact
       double precision dm2(nact,nact,nact,nact)
       integer g_u
       integer g_coul
       integer g_b
c
c
       integer v, x, vx, aoff, aend
       integer clo, chi
       integer l_b, k_b, l_t1, k_t1, l_u, k_u, k_xj, ld
       logical ga_check_JKblocked
       external ga_check_JKblocked
c
       aoff = nclosed + 1
       aend = nclosed + nact
       if (.not.ma_push_get(MT_DBL,(nbf*nact),'b',l_b,k_b))
     $   call errquit('mcscf_J_trace: cannot allocate b',0, MA_ERR)
       if (.not.ma_push_get(MT_DBL,(nbf*nbf),'a tmp',l_t1,k_t1))
     $   call errquit('mcscf_J_trace: cannot allocate t1',0, MA_ERR)
       if (.not.ma_push_get(MT_DBL,(nbf*nbf),'u',l_u,k_u))
     $   call errquit('mcscf_J_trace: cannot allocate u',0, MA_ERR)
       call ga_get(g_u,1,nbf,1,nbf,dbl_mb(k_u),nbf)
       if (.not.ga_check_JKblocked( g_coul, nact, nbf, clo, chi ))
     $   call errquit('mcscf_J_trace: wrong distrib operator',0, GA_ERR)
       call dfill((nact*nbf), 0.d0, dbl_mb(k_b), 1)
c
       do v=1,nact
         do x=1,v
           vx = (v*(v-1))/2 + x
           if ((vx.ge.clo).and.(vx.le.chi)) then
             call ga_access(g_coul, 1, (nbf*nbf), vx, vx, k_xj, ld )
             call dgemm('n', 'n', nbf, nbf, nbf, 1.d0, dbl_mb(k_xj),
     $                  nbf, dbl_mb(k_u), nbf, 0.d0, dbl_mb(k_t1), nbf )
             call mcscf_J_trace_a(nbf, nclosed, nact, v, x, dm2,
     $                            dbl_mb(k_t1), dbl_mb(k_b) )
             call ga_release(g_coul, 1, (nbf*nbf), vx, vx )
           endif
         enddo
       enddo
       call ga_acc(g_b, 1, nbf, aoff, aend, dbl_mb(k_b), nbf, 1.d0 )
       call ga_sync()
c
c
c
       if (.not.ma_pop_stack(l_u))
     $   call errquit('mcscf_J_trace: failed pop stack',0, MA_ERR)
       if (.not.ma_pop_stack(l_t1))
     $   call errquit('mcscf_J_trace: failed pop stack',0, MA_ERR)
       if (.not.ma_pop_stack(l_b))
     $   call errquit('mcscf_J_trace: failed pop stack',0, MA_ERR)
       return
       end
       


c
c          ~vx       vx        vx
c   B   =  J   (Gamma   + Gamma  )
c    tp     pu       ut        tu
c
c

      subroutine mcscf_J_trace_a( nbf, nclosed, nact, v, x, dm2,
     $                            xj1, b )
      implicit none
      integer nbf
      integer nclosed
      integer nact
      integer v, x
      double precision dm2(nact,nact,nact,nact)
      double precision xj1(nbf,nbf)
      double precision b(nbf,nact)
      integer p, t, u
      double precision scale

      scale = 1.d0
      if (v.eq.x) scale = 0.5d0
      do p=1,nbf
        do t=1,nact
          do u=1,nact
            b(p,t) = b(p,t) + scale*xj1(p,nclosed+u)*
     $                        (dm2(u,t,x,v)+dm2(t,u,x,v))
            b(p,t) = b(p,t) + scale*xj1(p,nclosed+u)*
     $                        (dm2(u,t,v,x)+dm2(t,u,v,x))
          enddo
        enddo
      enddo
      return
      end











c
c  Do this trace:
c
c              vx           vt        tv
c      B   =  K   T   (Gamma   + Gamma  )
c       pt     pq  qu       ux        ux
c
c             ~vx       vt        tv
c          =  K   (Gamma   + Gamma  )
c              pu       ux        ux
c
c  Note!
c
c          vx      xv
c         K   !=  K
c          pq      pq
c
      subroutine mcscf_K_trace( nbf, nclosed, nact, dm2, g_t,
     $                          g_exch, g_b )
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
       integer nbf
       integer nclosed
       integer nact
       double precision dm2(nact,nact,nact,nact)
       integer g_t
       integer g_exch
       integer g_b
c
c
       integer v, x, vx, aoff, aend
       integer clo, chi
       integer l_b, k_b, l_t1, k_t1, l_t, k_t, k_xj, ld
       logical ga_check_JKblocked
       external ga_check_JKblocked
c
       aoff = nclosed + 1
       aend = nclosed + nact
       if (.not.ma_push_get(MT_DBL,(nbf*nact),'b',l_b,k_b))
     $   call errquit('mcscf_K_trace: cannot allocate b',0, MA_ERR)
       if (.not.ma_push_get(MT_DBL,(nbf*nbf),'a tmp',l_t1,k_t1))
     $   call errquit('mcscf_K_trace: cannot allocate t1',0, MA_ERR)
       if (.not.ma_push_get(MT_DBL,(nbf*nbf),'u',l_t,k_t))
     $   call errquit('mcscf_K_trace: cannot allocate u',0, MA_ERR)
       call ga_get(g_t,1,nbf,1,nbf,dbl_mb(k_t),nbf)
       if (.not.ga_check_JKblocked( g_exch, nact, nbf, clo, chi ))
     $   call errquit('mcscf_K_trace: wrong distrib operator',0, GA_ERR)
       call dfill((nact*nbf), 0.d0, dbl_mb(k_b), 1)
c
       do v=1,nact
         do x=1,v
           vx = (v*(v-1))/2 + x
           if ((vx.ge.clo).and.(vx.le.chi)) then
             call ga_access(g_exch, 1, (nbf*nbf), vx, vx, k_xj, ld )
             call dgemm('n', 'n', nbf, nbf, nbf, 1.d0, dbl_mb(k_xj),
     $                  nbf, dbl_mb(k_t), nbf, 0.d0,
     $                  dbl_mb(k_t1), nbf )
             call mcscf_K_trace_a(nbf, nclosed, nact, v, x, dm2,
     $                            dbl_mb(k_t1), dbl_mb(k_b) )
             if (v.ne.x) then
               call dgemm('t', 'n', nbf, nbf, nbf, 1.d0, dbl_mb(k_xj),
     $                     nbf, dbl_mb(k_t), nbf, 0.d0,
     $                     dbl_mb(k_t1), nbf )
               call mcscf_K_trace_a(nbf, nclosed, nact, x, v, dm2,
     $                              dbl_mb(k_t1), dbl_mb(k_b) )
             endif
             call ga_release(g_exch, 1, (nbf*nbf), vx, vx )
           endif
         enddo
       enddo
       call ga_acc(g_b, 1, nbf, aoff, aend, dbl_mb(k_b), nbf, 1.d0 )
       call ga_sync()
c
c
c
       if (.not.ma_pop_stack(l_t))
     $   call errquit('mcscf_K_trace: failed pop stack',0, MA_ERR)
       if (.not.ma_pop_stack(l_t1))
     $   call errquit('mcscf_K_trace: failed pop stack',0, MA_ERR)
       if (.not.ma_pop_stack(l_b))
     $   call errquit('mcscf_K_trace: failed pop stack',0, MA_ERR)
       return
       end




c
c          ~vx       tv       tv
c   B   =  K   (Gamma  + Gamma  )
c    tp     pu       ux       xu
c
c
      subroutine mcscf_K_trace_a( nbf, nclosed, nact, v, x, dm2,
     $                            xk1, b )
      implicit none
      integer nbf
      integer nclosed
      integer nact
      integer v, x
      double precision dm2(nact,nact,nact,nact)
      double precision xk1(nbf,nbf)
      double precision b(nbf,nact)
      integer p, t, u

      do p=1,nbf
        do t=1,nact
          do u=1,nact
            b(p,t) = b(p,t) + xk1(p,nclosed+u)*
     $                        (dm2(u,x,t,v)+dm2(x,u,t,v))
          enddo
        enddo
      enddo
      return
      end









c
c   Return the new inactive Fock matrix
c   for a given unitrary rotation
c
c   To a good approximation, if the inactive
c   is frozen, just do a two-index
c   transformation, Fnew = U.Fold.U
c  
c   Must use the correct form if inactive orbitals 
c   have changed significantly, odelta_core = .true.
c   Requires a Fock build.
c
c
      subroutine mcscf_b_trf_ifock( geom, basis, nbf, nclosed, nact,
     $                              oskel, tol2e, g_vecs0, g_ifock,
     $                              g_u, odelta_core, g_ifock_new )
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"      
      integer geom, basis                              ! [input]  Handles
      integer nbf, nclosed, nact                       ! [input]  Orbital parameters
      logical oskel                                    ! [input]  Symmetry toggle
      double precision tol2e                           ! [input]  AO integral tolerance
      integer g_vecs0                                  ! [input]  MO vectors at alpha=0
      integer g_ifock                                  ! [input]  Inactive Fock
      integer g_u                                      ! [input]  Unitrary transformation
      logical odelta_core                              ! [input]  Toggle inactive core response (expensive!)
      integer g_ifock_new                              ! [output] New Inactive Fock
c
      integer nfock
      parameter(nfock=1)
      integer g_deltadc, g_d1, g_f1
      integer i, aoff, aend
      integer iv_fock(nfock), iv_dens(nfock)
      double precision jfac(nfock), kfac(nfock)
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c
      data jfac/1.d0/
      data kfac/-0.5d0/
c
c
      aoff = nclosed + 1
      aend = nclosed + nact
*ga:1:0
      if (.not.ga_create(MT_DBL, nbf, nbf, 'tmp', nbf, 0, g_deltadc))
     $   call errquit('mcscf_b: cannot allocate dD',0, GA_ERR)
      call ga_zero(g_deltadc)
      if (odelta_core) then
        g_d1 = ga_create_atom_blocked( geom, basis, 'Density 1')
        g_f1 = ga_create_atom_blocked( geom, basis, 'Fock 1')
c
c  T = U - 1    
c      
        call ga_copy(g_u, g_d1)
        do i=ga_nodeid()+1,nbf,ga_nnodes()
          call ga_acc(g_d1, i, i, i, i, -1.d0, 1, 1.d0 )
        enddo
        call ga_sync()
c
c  Inactive density diff (only double occupation)   
c
c    c          t
c  dD  = TD + DT 
c      
        call ga_zero(g_f1)
        do i=ga_nodeid()+1,nclosed,ga_nnodes()
          call ga_put(g_f1, i, i, i, i, 2.d0, 1 )
        enddo
        call ga_sync()
        call ga_dgemm( 'n', 'n', nbf, nbf, nbf, 1.d0, g_d1, g_f1,
     $                 0.d0, g_deltadc )
        call ga_dgemm( 'n', 't', nbf, nbf, nbf, 1.d0, g_f1, g_d1,
     $                 1.d0, g_deltadc )
c
c                   t
c   dD(AO) =  C.dD.C
c
        call ga_dgemm( 'n', 't', nbf, nbf, nbf, 1.d0, g_deltadc,
     $                 g_vecs0, 0.d0, g_f1 )
        call ga_dgemm( 'n', 'n', nbf, nbf, nbf, 1.d0, g_vecs0,
     $                 g_f1, 0.d0, g_d1 )
c
c  Compute Fock matrix
c
        call ga_zero(g_f1)
        iv_dens(1) = g_d1
        iv_fock(1) = g_f1
        call fock_2e( geom, basis, nfock, jfac, kfac, tol2e, oskel,
     $                iv_dens, iv_fock, .false. )
        if (oskel) call sym_symmetrize(geom, basis, .false., g_f1)
c     
c  Transform to MO
c
        call two_index_transf( g_f1, g_vecs0, g_vecs0, g_d1,
     $                         g_deltadc )
c
c  Clean up
c
        if (.not.ga_destroy(g_d1))
     $     call errquit('mcscf_b: cannot destroy AO dens1',0, GA_ERR)
        if (.not.ga_destroy(g_f1))
     $       call errquit('mcscf_b: cannot destroy AO fock1',0, GA_ERR)
      endif
c
c  Transform old inactive fock
c
c    I     I        c
c  dF  =  F.U + F[dD ]       second term computed already in AO section
c
c   I      t   I
c  F   =  U .dF
c
      call ga_dgemm( 'n', 'n', nbf, nbf, nbf, 1.d0, g_ifock, g_u,
     $               1.d0, g_deltadc )
      call ga_dgemm( 't', 'n', nbf, nbf, nbf, 1.d0, g_u, g_deltadc,
     $               0.d0, g_ifock_new )
c
c  Clean up
c
      if (.not.ga_destroy(g_deltadc))
     $  call errquit('mcscf_b: cannot destroy deltadc',0, GA_ERR)
      return
      end












#ifdef NOCOMPILE
c
c #######################################################
c                        ij
c  Routines to generate G   operator
c  Noddy version and not quite working...leave it out
c
c #######################################################
c
c
      subroutine mcscf_g( nbf, nclosed, nact, dm1, dm2,
     $                    g_ifock, g_afock )
      implicit none
      integer nbf
      integer nclosed
      integer nact
      double precision dm1(nact,nact)
      double precision dm2(nact,nact,nact,nact)
      integer g_ifock
      integer g_afock

      integer nocc
      double precision l(200000)
      double precision g(200000)

      nocc = nclosed + nact
      call dfill((nbf*nbf*nocc*nocc),0.d0,g,1)

      call mcscf_f2g( nbf, nclosed, nact, dm1, g_ifock,
     $                g_afock, g )

      call mcscf_lmat( nbf, nocc, l )
      call mcscf_l2g( nbf, nclosed, nact, dm1, l, g )
      call mcscf_la2g( nbf, nclosed, nact, dm2, g )

      call mcscf_g2file( nbf, nclosed, nact, g )


      return
      end


      subroutine mcscf_g2file( nbf, nclosed, nact, g )
      implicit none
      integer nbf
      integer nclosed
      integer nact
      double precision g(nbf*nbf,(nclosed+nact),(nclosed+nact))
      integer i,j,nn

      nn = nbf*nbf
      open(unit=89,file='Goper',status='unknown',form='unformatted')
      do i=1,nclosed+nact
        do j=1,nclosed+nact
          call hv_writev(89,nn,g(1,i,j))
          write(6,771) i,j
 771      format(' G oper: ',2i5)
          call moints_matprint(nbf,nbf,g(1,i,j))
        enddo
      enddo
      close(89)
      return
      end

      



c
c         tu               vx       tu       tu
c        G    +=    sum   J   (Gamma  + Gamma  )
c         pq           vx  pq       vx       xv
c
c
c         tu               vx       xu       xu
c        G    +=  2 sum   K   (Gamma  + Gamma  )
c         pq           vx  pq       vt       tv
c
c
      subroutine mcscf_la2g( nbf, nclosed, nact, dm2, g )
      implicit none
      integer nbf
      integer nclosed
      integer nact
      double precision dm2(nact,nact,nact,nact)
      double precision g(nbf*nbf,(nclosed+nact),(nclosed+nact))

      double precision hlp(1000), hlp2(1000), scale, dd
      integer t,u,tt,uu,nn,vv,xx,v,x

      nn = nbf*nbf
      open(unit=88,file='JKints',status='old',form='unformatted')
      do vv=1,nclosed+nact
        do xx=1,vv
          call hv_readv(88,nn,hlp)
          if ((vv.gt.nclosed).and.(xx.gt.nclosed)) then
            v = vv - nclosed
            x = xx - nclosed
            scale = 1.d0
            if (v.eq.x) scale = 0.5d0
            do t=1,nact
              tt = nclosed + t
              do u=1,nact
                uu = nclosed + u
                dd = (dm2(v,x,t,u) + dm2(x,v,t,u))*scale
                call daxpy(nn,dd,hlp,1,g(1,tt,uu),1)
                dd = (dm2(x,v,u,t) + dm2(v,x,u,t))*scale
                call daxpy(nn,dd,hlp,1,g(1,uu,tt),1)
              enddo
            enddo
          endif
        enddo
      enddo
      do vv=1,nclosed+nact
        do xx=1,vv
          call hv_readv(88,nn,hlp)
          if ((vv.gt.nclosed).and.(xx.gt.nclosed)) then
            v = vv - nclosed
            x = xx - nclosed
            scale = 2.d0
            if (v.eq.x) scale = 1.0d0
            do t=1,nact
              tt = nclosed + t
              do u=1,nact
                uu = nclosed + u
                dd = (dm2(v,t,x,u) + dm2(t,v,x,u))*scale
                call daxpy(nn,dd,hlp,1,g(1,tt,uu),1)
                call transpose_nw(nbf,nbf,hlp,hlp2)
                dd = (dm2(x,u,v,t) + dm2(u,x,v,t))*scale
                call daxpy(nn,dd,hlp2,1,g(1,uu,tt),1)
              enddo
            enddo
          endif
        enddo
      enddo
      close(88)

      return
      end






      subroutine mcscf_f2g( nbf, nclosed, nact, dm1, g_ifock,
     $                       g_afock, g )
      implicit none
      integer nbf
      integer nclosed
      integer nact
      double precision dm1(nact,nact)
      integer g_ifock
      integer g_afock
      double precision g(nbf*nbf,(nclosed+nact),(nclosed+nact))
      double precision x(1000),y(1000)

      integer i,t,u,tt,uu
      INTEGER AA,A,NVIR,J
      DOUBLE PRECISION XXX

      call ga_get(g_afock,1,nbf,1,nbf,x,nbf)
      call ga_get(g_ifock,1,nbf,1,nbf,y,nbf)
      do i=1,nclosed
        call daxpy((nbf*nbf),2.d0,x,1,g(1,i,i),1)
        call daxpy((nbf*nbf),2.d0,y,1,g(1,i,i),1)
      enddo

      do t=1,nact
        tt = nclosed + t
        do u=1,nact
          uu = nclosed + u
          call daxpy((nbf*nbf), dm1(t,u), y, 1, g(1,tt,uu), 1 )
        enddo
      enddo
      return
      end

        



      subroutine mcscf_l2g( nbf, nclosed, nact, dm1, l, g )
      implicit none
      integer nbf
      integer nclosed
      integer nact
      double precision dm1(nact,nact)
      double precision l(nbf*nbf,(nclosed+nact),(nclosed+nact))
      double precision g(nbf*nbf,(nclosed+nact),(nclosed+nact))
      
      integer i,j,t,u,tt,uu,nn
      
      nn = nbf*nbf
      do i=1,nclosed
        do j=1,nclosed
          call daxpy(nn, 2.d0, l(1,i,j), 1, g(1,i,j), 1 )
        enddo
      enddo
      do i=1,nclosed
        do t=1,nact
          tt = nclosed + t
          do u=1,nact
            uu =  nclosed + u
            call daxpy(nn, dm1(t,u), l(1,i,uu), 1, g(1,i,tt), 1)
            call daxpy(nn, dm1(t,u), l(1,uu,i), 1, g(1,tt,i), 1)
          enddo
        enddo
      enddo
      return
      end

          
        


c
c
c  Form MO supermatrix from J & K ints on disk
c
c          pq        pq     qp     pq
c         L    =  4 K   -  K   -  J
c          rs        rs     rs     rs
c
c
c
      subroutine mcscf_lmat( nbf, nocc, l )
      implicit none
      integer nbf
      integer nocc
      double precision l(nbf*nbf,nocc,nocc)

      integer i,j
      double precision hlp(10000), hlp2(10000)
      double precision scale,s4

      call dfill((nbf*nbf*nocc*nocc),0.d0,l,1)
      open(unit=88,file='JKints',status='old',form='unformatted')
      do i=1,nocc
        do j=1,i
          scale = -1.d0
          if (i.eq.j) scale = -0.5d0
          call hv_readv(88,(nbf*nbf),hlp)
          call daxpy((nbf*nbf), scale, hlp, 1, l(1,i,j), 1 )
          call daxpy((nbf*nbf), scale, hlp, 1, l(1,j,i), 1 )
        enddo
      enddo
      do i=1,nocc
        do j=1,i
          scale = -1.d0
          if (i.eq.j) scale = -0.5d0
          s4 = 4.d0
          if (i.eq.j) s4 = 2.d0
          call hv_readv(88,(nbf*nbf),hlp)
          call daxpy((nbf*nbf), s4, hlp, 1, l(1,i,j), 1 )
          call daxpy((nbf*nbf), scale, hlp, 1, l(1,j,i), 1 )
          call transpose_nw( nbf, nbf, hlp, hlp2 )
          call daxpy((nbf*nbf), scale, hlp2, 1, l(1,i,j), 1 )
          call daxpy((nbf*nbf), s4, hlp2, 1, l(1,j,i), 1 )
        enddo
      enddo
      close(88)
      return
      end







      subroutine mcscf_gt( nbf, nclosed, nact, tt, etgt )
      implicit none
      integer nbf
      integer nclosed
      integer nact
      double precision tt(nbf,nbf)
      double precision etgt

      double precision hlp(1000),vtmp(10000)
      integer i,j,nn
      double precision ddot

      nn = nbf*nbf
      etgt = 0.d0
      open(unit=89,file='Goper',status='unknown',form='unformatted')
      do i=1,nclosed+nact
        do j=1,nclosed+nact
          call hv_readv(89,nn,hlp)
          call dgemv( 'n', nbf, nbf, 1.d0, hlp, nbf, tt(1,j), 1,
     $                0.d0, vtmp, 1 )
          etgt = etgt + ddot(nbf, vtmp, 1, tt(1,i), 1 )
        enddo
      enddo
      close(89)
      return
      end

#endif




