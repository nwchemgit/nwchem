c
c Return the MCSCF 1e Hessian-vector product
c
c NB: some analogous 1e-terms from ROHF expression
c     cannot be computed from the Fock matrices
c     and are actually included in the 2e component
c
       subroutine mcscf_hessv_1e( geom, basis, nbf, nclosed, nact,
     $                            lshift, dm1, g_ifock, g_afock,
     $                            g_gfock, g_tmp1, g_tmp2, g_tmp3,
     $                            g_x, g_ax )
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
       double precision lshift                      ! [input] Level shift
       double precision dm1(nact,nact)              ! [input] 1PDM
       integer g_ifock                              ! [input] Inactive Fock matrix
       integer g_afock                              ! [input] Active Fock matrix
       integer g_gfock                              ! [input] Generalized Fock matrix
       integer g_tmp1, g_tmp2, g_tmp3               ! [input] Temporaries (nbf * nbf)
       integer g_x                                  ! [input] Argument vector
       integer g_ax                                 ! [output] Hessian vector product
c
c
       integer g_tmp4
       integer nvir, vlen, voff, aoff, aend
c
c
       if (omcscfprof) call pstat_on(ps_hv1e)
       nvir = nbf - nclosed - nact
       vlen = (nclosed+nact)*nvir + nclosed*nact
       voff = nclosed + nact + 1
       aoff = nclosed + 1
       aend = nclosed + nact
       if (.not.ga_duplicate(g_tmp3, g_tmp4, 'temp4'))
     $      call errquit('mcscf_hessv_1e: cannot duplicate',0, GA_ERR)
       call ga_zero(g_tmp4)
       call ga_zero(g_tmp3)
c
c Reshape argument vector by copying into matrix
c
       call ga_zero(g_tmp1)
       call mcscf_vec2mat( nbf, nclosed, nact, 1.d0, g_x, 0.d0, g_tmp1)
c
c 1-PDM density matrix is hidden in unused 
c active-active subblock
c
       call ga_put(g_tmp1, aoff, aend, aoff, aend, dm1, nact )
c
c                                          cv
c Useful intermediate (equivalent to ROHF F  )
c
c    Z     I     A
c   F  =  F  +  F
c
c
       call ga_dadd(1.d0, g_ifock, 1.d0, g_afock, g_tmp2 )
c
c Closed-virtual bit
c
c      iv             vv iv
c    Z1  = 4.(iF + aF)  X                    (eqn A4)
c       
       call ga_matmul_patch( 'n', 'n', 4.d0, 0.d0,
     $                       g_tmp2, voff, nbf, voff, nbf,
     $                       g_tmp1, voff, nbf, 1, nclosed,
     $                       g_tmp3, voff, nbf, 1, nclosed )
c
c      iv             ii iv
c    Z2  = 4.(iF + aF)  X                    (eqn A4)
c
       call ga_matmul_patch( 'n', 'n', -4.d0, 1.d0,
     $                       g_tmp1, voff, nbf, 1, nclosed,
     $                       g_tmp2, 1, nclosed, 1, nclosed,
     $                       g_tmp3, voff, nbf, 1, nclosed )
c
c      iv              ia av
c    Z3  = -2.(iF + aF)  X                   (eqn A5)
c
       call ga_matmul_patch( 'n', 'n', -2.d0, 1.d0,
     $                       g_tmp1, voff, nbf, aoff, aend,
     $                       g_tmp2, aoff, aend, 1, nclosed,
     $                       g_tmp3, voff, nbf, 1, nclosed )
c
c      iv    ia av
c    Z4  = -F  X                             (eqn A5)
c
       call ga_matmul_patch( 'n', 'n', -1.d0, 1.d0,             ! Possibly the transpose? 
     $                       g_tmp1, voff, nbf, aoff, aend,     ! Careful here ...        
     $                       g_gfock, aoff, aend, 1, nclosed,
     $                       g_tmp3, voff, nbf, 1, nclosed )
c
c      iv             av ia
c    Z5  = 4.(iF + aF)  X                    (eqn A2)
c
       call ga_matmul_patch( 'n', 'n', 4.d0, 1.d0,
     $                       g_tmp2, voff, nbf, aoff, aend,
     $                       g_tmp1, aoff, aend, 1, nclosed,
     $                       g_tmp3, voff, nbf, 1, nclosed )
c
c      iv    av ia
c    Z6  = -F  X                             (eqn A2)
c
       call ga_matmul_patch( 't', 'n', -1.d0, 1.d0,             ! Possibly the transpose? 
     $                       g_gfock, voff, nbf, aoff, aend,    ! Careful here ...        
     $                       g_tmp1, aoff, aend, 1, nclosed,
     $                       g_tmp3, voff, nbf, 1, nclosed )
c
c
c Active-virtual bit
c
c      av              ia iv
c    Z1  = -2.(iF + aF)  X                   (eqn A5)
c
       call ga_matmul_patch( 'n', 'n', -2.d0, 1.d0,
     $                       g_tmp1, voff, nbf, 1, nclosed,
     $                       g_tmp2, 1, nclosed, aoff, aend,
     $                       g_tmp3, voff, nbf, aoff, aend )
c
c      av    ia iv
c    Z2  = -F  X                             (eqn A5)
c
       call ga_matmul_patch( 'n', 't', -1.d0, 1.d0,
     $                       g_tmp1, voff, nbf, 1, nclosed,
     $                       g_gfock, 1, nclosed, aoff, aend,
     $                       g_tmp3, voff, nbf, aoff, aend )
c
c      av       av  vv
c    Z3  = 2.d.X  iF                         (eqn A6)
c
       call ga_matmul_patch( 'n', 'n', 1.d0, 0.d0,
     $                       g_tmp1, voff, nbf, aoff, aend,     
     $                       g_tmp1, aoff, aend, aoff, aend,    !                                     
     $                       g_tmp4, voff, nbf, aoff, aend )    ! X' = X.D    (multiply arg by 1PDM)  
       call ga_matmul_patch( 'n', 'n', 2.d0, 1.d0,              !                                     
     $                       g_ifock, voff, nbf, voff, nbf,     ! Z = F.X' = F.X.D                    
     $                       g_tmp4, voff, nbf, aoff, aend,
     $                       g_tmp3, voff, nbf, aoff, aend )
c
c      av      aa  av
c    Z4  = -2.F  X                           (eqn A6)
c
       call ga_matmul_patch( 'n', 'n', -2.d0, 1.d0,             ! Possibly the transpose? 
     $                       g_tmp1, voff, nbf, aoff, aend,     ! Careful here ...        
     $                       g_gfock, aoff, aend, aoff, aend,
     $                       g_tmp3, voff, nbf, aoff, aend )
c
c      av        ia   iv
c    Z5  = -2.d.X  iF                        (eqn A3)
c
       call ga_matmul_patch( 't', 'n', 1.d0, 0.d0,
     $                       g_tmp1, 1, nclosed, aoff, aend,
     $                       g_tmp1, aoff, aend, aoff, aend,    !                                         
     $                       g_tmp4, 1, nclosed, aoff, aend )   ! X' = tr(X).D    (multiply arg by 1PDM)  
       call ga_matmul_patch( 'n', 'n', -2.d0, 1.d0,             !                                         
     $                       g_ifock, voff, nbf, 1, nclosed,    ! Z = F.X' = F.tr(X).D                    
     $                       g_tmp4, 1, nclosed, aoff, aend,
     $                       g_tmp3, voff, nbf, aoff, aend )
c
c      av             iv ia
c    Z6  = 2.(iF + aF)  X                    (eqn A3)
c
       call ga_matmul_patch( 'n', 't', 2.d0, 1.d0,
     $                       g_tmp2, voff, nbf, 1, nclosed,
     $                       g_tmp1, 1, nclosed, aoff, aend,
     $                       g_tmp3, voff, nbf, aoff, aend )
c
c
c
c Closed-active bit
c
c      ia      iv         av
c    Y1   = 4.X  (iF + aF)                    (eqn A2)
c
       call ga_matmul_patch( 'n', 'n', 4.d0, 1.d0,
     $                       g_tmp2, aoff, aend, voff, nbf,
     $                       g_tmp1, voff, nbf, 1, nclosed,
     $                       g_tmp3, aoff, aend, 1, nclosed )
c
c      ia    iv av
c    Y2   = X  F                              (eqn A2)
c
       call ga_matmul_patch( 'n', 'n', -1.d0, 1.d0,
     $                       g_gfock, aoff, aend, voff, nbf,
     $                       g_tmp1, voff, nbf, 1, nclosed,
     $                       g_tmp3, aoff, aend, 1, nclosed )
c
c      ia         av   iv
c    Y1   = -2.d.X  .iF                       (eqn A3)
c               
c
       call ga_matmul_patch( 'n', 't', 1.d0, 0.d0,
     $                       g_tmp1, aoff, aend, aoff, aend,
     $                       g_tmp1, aoff, aend, voff, nbf,    !                                         
     $                       g_tmp4, aoff, aend, voff, nbf )   ! X' = D.tr(X)    (multiply arg by 1PDM)  
       call ga_matmul_patch( 'n', 'n', -2.d0, 1.d0,            !                                         
     $                       g_tmp4, aoff, aend, voff, nbf,    ! Z = X'.F = D.tr(X).F                  
     $                       g_ifock, voff, nbf, 1, nclosed,
     $                       g_tmp3, aoff, aend, 1, nclosed )
c
c      ia      av          iv
c    Y2   = 2.X  .(iF + aF)                  (eqn A3)
c
       call ga_matmul_patch( 't', 'n', 2.d0, 1.d0,
     $                       g_tmp1, aoff, aend, voff, nbf,
     $                       g_tmp2, voff, nbf, 1, nclosed,
     $                       g_tmp3, aoff, aend, 1, nclosed )
c
c      ia        ia   ii
c    Y3   = 2.d.X  .iF                       (eqn A1)
c
       call ga_matmul_patch( 'n', 'n', 1.d0, 0.d0,
     $                       g_tmp1, aoff, aend, aoff, aend,
     $                       g_tmp1, aoff, aend, 1, nclosed,   !                                       
     $                       g_tmp4, aoff, aend, 1, nclosed )  ! X' = D.X    (multiply arg by 1PDM)
       call ga_matmul_patch( 'n', 'n', 2.d0, 1.d0,             !                                       
     $                       g_tmp4, aoff, aend, 1, nclosed,   ! Z = X'.F = D.X.F                  
     $                       g_ifock, 1, nclosed, 1, nclosed,
     $                       g_tmp3, aoff, aend, 1, nclosed )
c
c      ia      ia               aa
c    Y4   = 2.X  (2iF + 2aF - F)             (eqn A1)
c             
                                                               !  x       z         Should really use patch routine
       call ga_dadd(2.d0, g_tmp2, -1.d0, g_gfock, g_tmp4 )     ! F  =  2.F  - F     here but only active-active
                                                               !                    block of g_tmp4 used.
       call ga_matmul_patch( 'n', 'n', 2.d0, 1.d0,
     $                       g_tmp4, aoff, aend, aoff, aend,
     $                       g_tmp1, aoff, aend, 1, nclosed,
     $                       g_tmp3, aoff, aend, 1, nclosed )
c
c      ia       ia         ii
c    Y5   = -4.X  (iF + aF)                  (eqn A1)
c
       call ga_matmul_patch( 'n', 'n', -4.d0, 1.d0,
     $                       g_tmp1, aoff, aend, 1, nclosed,
     $                       g_tmp2, 1, nclosed, 1, nclosed,
     $                       g_tmp3, aoff, aend, 1, nclosed )
c
c Copy result from matrix format to vector format
c
       call ga_dadd( lshift, g_x, 0.d0, g_ax, g_ax )
       call mcscf_mat2vec( nbf, nclosed, nact, 1.d0, g_tmp3, 1.d0, g_ax)
c
c
c
       if (.not.ga_destroy(g_tmp4))
     $      call errquit('mcscf_hessv_1e: cannot destroy',0, GA_ERR)
       if (omcscfprof) call pstat_off(ps_hv1e)
       return
       end

