       subroutine mcscf( rtdb )
       implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "context.fh"
#include "util.fh"
#include "sym.fh"
c
c
       integer rtdb
c
c
       integer geom, basis
       integer nbf, nclosed, nvir, nact
       integer nactel, nela, nelb, multip, orlen
       integer noper
       integer l_occ, k_occ, l_sym, k_sym
       integer l_evals, k_evals
       integer l_dm1, k_dm1, l_dm2, k_dm2
       integer l_tmp, k_tmp
       integer g_movecs, g_coul, g_exch
       integer g_hcore
       integer g_conj, g_gradold
       integer g_afock, g_ifock, g_gfock
       integer g_grad, g_prod, g_x
       integer mo_lo, mo_hi
       integer i, j, iter, maxiter, nmixed, clo, chi
       integer avstart, iastart
       integer icgreset                                   ! Number of iterations between CG resets
       double precision pfac
       double precision thresh
       double precision eone, etwo, energy, enrep
       double precision tol2e, gnorm, xx
       logical oconj
       logical oskel
       logical osym
       logical oorbsolve
       logical oprinfo
c
c  Debugging variables
c
       integer mclosed, mopen                             ! ROHF occupation
       integer g_fcv, g_fpv, g_fcp                        ! ROHF Fock matrices
       double precision crap(1000), hh(100000)
       integer ssiz
       parameter(ssiz=1000)
       double precision scr(ssiz), ev(ssiz)
       integer info
c
c
       integer ga_create_atom_blocked, ga_create_JKblocked
       external ga_create_atom_blocked, ga_create_JKblocked
       integer mcscf_rohf_den2occ
       external mcscf_rohf_den2occ
c
c
c
       data osym/.true./
       data oorbsolve/.true./
c
c
c
       if (.not. context_push('mcscf'))
     $   call errquit('mcscf: context_push failed', 0)
       call util_print_push
       call util_print_rtdb_load(rtdb, 'mcscf')
       oprinfo = util_print('information',print_default)
c
c  Get w.f. parameters
c
       call mcscf_getinfo(rtdb, geom, basis, nbf, nclosed, nact )
       nvir = nbf - nclosed - nact
       orlen = (nclosed*nvir) + (nact*nvir) + (nclosed*nact)
       avstart = nclosed*nvir + 1
       iastart = ((nclosed+nact)*nvir) + 1
       if (.not.rtdb_get(rtdb,'scf:skeleton',MT_LOG,1,oskel)) then
         oskel = sym_number_ops(geom).gt.0
       endif
       if (.not.geom_nuc_rep_energy( geom, enrep ))
     $   call errquit('mcscf: cannot retrieve nuclear repulsion',0)
       if (.not.rtdb_get(rtdb,'mcscf:tol2e',MT_DBL,1,tol2e))        
     $      tol2e = 0.d0                                                   ! Redundant recovered later
c
c  Get electron and spin multiplicity (NB: for info only)
c  Active elec and multiplicity must be set --- no defaults
c
      if (.not.rtdb_get(rtdb,'mcscf:nactelec',MT_INT,1,nactel))
     $   call errquit('number of active electrons not set',0)
      if (.not.rtdb_get(rtdb,'mcscf:multiplicity',MT_INT,1,multip))
     $   call errquit('spin multiplicity not set',0)
      nela = (nactel + multip - 1)/2
      nelb = nactel - nela
      if ((mod((nactel + multip - 1),2).ne.0).or.
     $    (nela.lt.0).or.(nelb.lt.0))
     $    call errquit('mcscf: incompatible elec and spin',0)
c
c  Get convergence info
c
      if (.not.rtdb_get(rtdb,'mcscf:maxiter',MT_INT,1,maxiter)) then
        if (.not.rtdb_get(rtdb,'scf:maxiter',MT_INT,1,maxiter))
     $     maxiter = 5
      endif
      if (.not.rtdb_get(rtdb,'mcscf:cgreset',MT_INT,1,icgreset))  
     $     icgreset = 20
      if (.not.rtdb_get(rtdb,'mcscf:thresh',MT_DBL,1,thresh)) then
        if (.not.rtdb_get(rtdb,'scf:thresh',MT_DBL,1,thresh))
     $     thresh = 1.d-6
      endif
c
c  Print info
c
       if (oprinfo) then
         write(6,900)
 900     format(/,27x,14('='),/
     $            28x,'Direct MCSCF',/
     $            27x,14('='),/)
         write(6,901) nclosed, nact, nactel, multip, orlen,
     $                1, avstart, iastart
 901     format(18x,35('-'),/
     $          20x,'Inactive shells:',10x,i5,/
     $          20x,'Active shells:',12x,i5,/,
     $          20x,'Active electrons:',9x,i5,/,
     $          20x,'Multiplicity:',13x,i5,/,
     $          20x,'Orbital rotations:',8x,i5,/
     $          25x,'Inact - Virt',9x,i5,/,
     $          25x,'Act - Virt',11x,i5,/,
     $          25x,'Inact - Act',10x,i5,/,
     $          18x,35('-'))
         write(6,902) maxiter, thresh, icgreset
 902     format(/,18x,35('-'),/,
     $          20x,'Maximum iterations:',7x,i5,/,
     $          20x,'Converge threshold:',e12.2,/,
     $          20x,'Conjugacy resets:',9x,i5,/,
     $          18x,35('-'))
       endif
c
c  Create occupation and orbital energy vectors
c
       if (.not.ma_push_get(MT_DBL, nbf, 'eps', l_evals, k_evals))
     $     call errquit('rohf: cannot allocate',0)
       if (.not.ma_push_get(MT_DBL, nbf, 'MO occ', l_occ, k_occ))
     $     call errquit('rohf: cannot allocate',0)
c
c  Allocate 1- & 2-PDM
c
       if (.not.ma_push_get(MT_DBL, (nact*nact*nact*nact),
     $                      '2P density', l_dm2, k_dm2))
     $     call errquit('rohf: cannot allocate MO density',0)
       if (.not.ma_push_get(MT_DBL, (nact*nact),
     $                      '1P density', l_dm1, k_dm1))
     $     call errquit('rohf: cannot allocate MO density',0)
c
c  Initialize some stuff
c
       call int_init(rtdb, 1, basis)
       call int_normalize(basis)
       call schwarz_init(geom, basis)
       call scf_get_fock_param( rtdb, tol2e )
c
c  Create/get MO vectors 
c
       if (.not.ga_create(MT_DBL,nbf,nbf,'MO vectors',nbf,1,g_movecs))
     $   call errquit('mcscf: cannot allocate MO vectors',0)
       call mcscf_getmovecs( rtdb, geom, basis, nbf, dbl_mb(k_occ),
     $                       dbl_mb(k_evals), g_movecs )
c
c  Get orbital symmetries
c 
       if (.not.ma_push_get(MT_INT, nbf, 'MO occ', l_sym, k_sym))
     $     call errquit('mcscf: cannot allocate symmetry',0)
       call sym_movecs_adapt( basis, 1.d-8, g_movecs,
     $                        int_mb(k_sym), nmixed )
       if (nmixed .ne. 0) call errquit(
     $   'mcscf: symmetry contamination in starting MOs', nmixed)
       if (util_print('orbital symmetry',print_default)) then
         write(6,887) 
 887     format(/,'Orbital symmetry')
         write(6,888) (int_mb(k_sym+i),i=0,nbf-1)
 888     format(16i3)
       endif
       if (.not.ga_duplicate(g_movecs,g_hcore,'1e hamiltonian'))
     $   call errquit('mcscf: cannot dup MO vectors',0)
c
c  Create conjugacy vectors
c
       if (.not.ga_create(MT_DBL,orlen,1,'Conj',1,1,g_conj))
     $      call errquit('mcscf: cannot allocate',0)
       if (.not.ga_create(MT_DBL,orlen,1,'Old Gradient',1,1,g_gradold))
     $      call errquit('mcscf: cannot allocate',0)
       call ga_zero(g_conj)
       call ga_zero(g_gradold)
       pfac = 1.d0
       oconj = .false.
c
c  Allocate operator matrices
c  Memory test required here!
c
       mo_lo = nclosed + 1
       mo_hi = nclosed + nact
       noper = (nact*(nact+1))/2
       g_coul = ga_create_JKblocked(noper,nbf,nbf,'Coulomb Oper')
       g_exch = ga_create_JKblocked(noper,nbf,nbf,'X Oper')
c
c
c  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c                    MCSCF Macroiteractions
c
c  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
       do iter=1,maxiter
         if ((ga_nodeid().eq.0).and.(oprinfo)) write(6,996) iter
 996     format(/,10x,'%%%%% Iteration ',i3,' %%%%%',/)
c
c  4-Index Tranform 
c
         if (ga_nodeid().eq.0) PRINT*,'Calling 4index'
         call moints_build_6x( basis, osym,
     $                         mo_lo, mo_lo, mo_hi, 1, nbf,
     $                         g_movecs, g_coul, .true.,
     $                         g_exch, .true., 16 )
         call mcscf_fcore( basis, nbf, nclosed, nact, g_movecs,
     $                     g_coul, g_exch, g_hcore )
c
c  CI solver
c
         if (iter.lt.1000) then
           if (ga_nodeid().eq.0) print*,'Calling CIsolve'
           call mcscf_cisolve( rtdb, geom, basis, nbf, nclosed, nact,
     $                         dbl_mb(k_evals), g_hcore, g_coul,
     $                         dbl_mb(k_dm1), dbl_mb(k_dm2), .true.,
     $                         .false. )
         endif
c
c  Orbital solver
c
         oconj = .true.
         if ((oconj).and.(mod(iter,icgreset).eq.0)) then
C           PRINT*,'=== Conjugacy Reset ==='
           call ga_zero(g_conj)
           call ga_zero(g_gradold)
           pfac = 1.d0
         endif
         if (oorbsolve) then
           if (ga_nodeid().eq.0) PRINT*,'Calling orbsolve'
           call mcscf_orbsolve( rtdb, geom, basis, nbf, nclosed, nact,
     $                          osym, int_mb(k_sym), dbl_mb(k_dm1),
     $                          dbl_mb(k_dm2), g_coul, g_exch, oconj,
     $                          pfac, g_conj, g_gradold, g_hcore,
     $                          dbl_mb(k_evals), eone, etwo, gnorm,
     $                          g_movecs )
         else
           if (ga_nodeid().eq.0) PRINT*,'Orbital solve skipped'
         endif
         energy = eone + etwo + enrep
         write(6,927) iter,energy,gnorm
 927     format(//,10x,'####',i5,5x,f20.12,5x,e12.5)
         if (gnorm.lt.thresh) goto 27
       enddo
c
c  Convert to natural orbital representation
c
 27    continue
       write(6,928)
 928   format(////,5x,30('='),/,12x,'MCSCF Converged',/,5x,30('='))
       call mcscf_naturalorb( nbf, nclosed, nact, dbl_mb(k_dm1),
     $                        dbl_mb(k_occ), g_movecs )
       if (ga_nodeid().eq.0) then
         write(6,880)
 880     format(///,10x,'Natural Orbitals and Occupation')
         if (.not.ma_push_get(MT_DBL, (nbf*(nclosed+nact)), 'tmp',
     $                        l_tmp, k_tmp))
     $     call errquit('mcscf: cannot allocate local MO',0)
         call ga_get(g_movecs, 1, nbf, 1, (nclosed+nact),
     $               dbl_mb(k_tmp), nbf)
         chi = 0
 33      clo = chi + 1
         chi = min((clo + 7),(nact+nclosed))
         write(6,*)
         write(6,881) (dbl_mb(k_occ+i-1),i=clo,chi)
 881     format(8f12.6)
         write(6,*)
         do i=1,nbf
           write(6,881) (dbl_mb(k_tmp+(j-1)*nbf+i-1),j=clo,chi)
         enddo
         if (chi.ne.(nact+nclosed)) goto 33
         if (.not.ma_pop_stack(l_tmp))
     $     call errquit('mcscf: cannot pop local MO',0)
       endif
c
c  Resolve secular eqn with natural orbitals
c
       call moints_build_6x( basis, osym,
     $                       mo_lo, mo_lo, mo_hi, 1, nbf,
     $                       g_movecs, g_coul, .true.,
     $                       g_exch, .true., 16 )
       call mcscf_fcore( basis, nbf, nclosed, nact, g_movecs,
     $                   g_coul, g_exch, g_hcore )
       if (ga_nodeid().eq.0) write(6,882) 
 882   format(//,10x,' Active Space Ground State CI vector')
       call mcscf_cisolve( rtdb, geom, basis, nbf, nclosed, nact,
     $                     dbl_mb(k_evals), g_hcore, g_coul,
     $                     dbl_mb(k_dm1), dbl_mb(k_dm2), .false.,
     $                     .true. )
c
c  Cleanup
c
       if (.not.ga_destroy(g_hcore))
     $      call errquit('mcscf: cannot destroy 1e hamiltonian',0)
       if (.not.ga_destroy(g_conj))
     $      call errquit('mcscf: cannot destroy conjugacy',0)
       if (.not.ga_destroy(g_gradold))
     $      call errquit('mcscf: cannot destroy old gradient',0)
c
c
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c                         MCSCF -> ROHF Test Section
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c  ROHF occupation
c
c$$$       call mcscf_occ2int( nbf, dbl_mb(k_occ), mclosed, mopen )
c$$$       write(6,903) mclosed, mopen
c$$$ 903   format(' ROHF Occupation',6x,'(closed):',i3,5x,'(open):',i3)
c
c  Trivial ROHF density matrix 
c
c$$$       call mcscf_rohf_modens(mopen,nact,dbl_mb(k_dm1),dbl_mb(k_dm2))
c
c
c
c

c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c         
c              MCSCF Hessian vector product Test Section
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c  Allocate MO matrices
c
       if (.not.ga_create(MT_DBL,nbf,nbf,'Act Fock',nbf,1,g_afock))
     $      call errquit('mcscf: cannot allocate active Fock',0)
       if (.not.ga_create(MT_DBL,nbf,nbf,'In Fock',nbf,1,g_ifock))
     $      call errquit('mcscf: cannot allocate inactive Fock',0)
       if (.not.ga_create(MT_DBL,nbf,nbf,'Gen Fock',nbf,1,g_gfock))
     $      call errquit('mcscf: cannot allocate general Fock',0)
c
c  Allocate MCSCF Gradient vector
c
       if (.not.ga_create(MT_DBL,orlen,1,'Gradient',1,1,g_grad))
     $      call errquit('rohf_head: cannot allocate',0)
       if (.not.ga_create(MT_DBL,orlen,1,'Product',1,1,g_prod))
     $      call errquit('rohf_head: cannot allocate',0)
       if (.not.ga_create(MT_DBL,orlen,1,'Arg vec',1,1,g_x))
     $      call errquit('rohf_head: cannot allocate',0)
       call ga_zero(g_grad)
       call ga_zero(g_x)
       call ga_zero(g_prod)
       call moints_build_6x( basis, osym,
     $                       mo_lo, mo_lo, mo_hi, 1, nbf,
     $                       g_movecs, g_coul, .true.,
     $                       g_exch, .true., 16 )
c
c
       call mcscf_fcore( basis, nbf, nclosed, nact, g_movecs,
     $                   g_coul, g_exch, g_afock )
c
c  Print 1-PDM
c 
       if (util_print('density matrix',print_debug)) then
         if (ga_nodeid().eq.0) then
           write(6,671)
 671       format(/,'<<<<<<< 1pdm density matrix >>>>>>>>>')
           call moints_matprint( nact, nact, dbl_mb(k_dm1) )
           write(6,672)
 672       format(/,'<<<<<<< 2pdm density matrix >>>>>>>>>')
           call mcscf_twopdm_print(nact,dbl_mb(k_dm2))
           write(6,673)
 673       format(/,'<<<<<<< symm. 2pdm density matrix >>>>>>>>>')
           call mcscf_symmetrize_2pdm( nact, dbl_mb(k_dm2), crap )
           call mcscf_twopdm_print(nact,crap)
         endif
       endif
c
c  Energy trace routine
c     
       call mcscf_etrace( geom, basis, nbf, nclosed, nact,
     $                    oskel, tol2e, dbl_mb(k_dm1), dbl_mb(k_dm2),
     $                    g_movecs, g_coul, eone, etwo )
       energy = eone + etwo + enrep
       if (util_print('trace energy',print_default)) then
         if (ga_nodeid().eq.0) write(6,674) energy
 674     format(/,20x,'<<<<<<<< trace energy: ',f20.10,' >>>>>>>',/)
       endif
c
c  Finite difference gradient
c
c$$$       call mcscf_fdiff_grad( geom, basis, nbf, nclosed, nact,
c$$$     $                        oskel, tol2e, dbl_mb(k_dm1),
c$$$     $                        dbl_mb(k_dm2), g_movecs, g_coul, g_grad)
c$$$       call ga_rowprint( ' ==== Finite Diff Gradient ====', g_grad )
c$$$       gnorm = sqrt(ga_ddot(g_grad,g_grad))
c$$$       write(6,773) gnorm
c$$$ 773   format(/,10x,'Finite diff gradient norm:',e12.4,/)
c
c  All MCSCF Fock matrices
c
       call mcscf_fock( geom, basis, nbf, nclosed, nact,
     $                  oskel, tol2e, dbl_mb(k_dm1), dbl_mb(k_dm2),
     $                  g_movecs, g_coul, eone, etwo,
     $                  g_ifock, g_afock, g_gfock )
       call mcscf_gfock2grad( nbf, nclosed, nact, g_gfock, g_grad )
       call mcscf_vsymm_chk( nbf, nclosed, nact, int_mb(k_sym),
     $                       g_grad, .true., .true., xx)

       call ga_rowprint( ' ==== Gradient ====', g_grad )
       gnorm = sqrt(ga_ddot(g_grad,g_grad))
       energy = eone + etwo + enrep
       if (ga_nodeid().eq.0) then
         write(6,772) gnorm
 772     format(/,10x,'Final gradient norm:',e12.4,/)
       endif
       if (util_print('fock energy',print_default)) then
         if (ga_nodeid().eq.0) write(6,675) energy
 675     format(/,20x,'<<<<<<<< fock energy: ',f20.10,' >>>>>>>',/)
       endif

c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c         
c        Explicit Hessian contruction and Eigenvalues Section
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
       call makeJK( basis, nbf, nclosed, nact, g_movecs )
       call dfill((orlen*orlen),0.d0,hh,1)
       call hmat( nbf, nclosed, nact, orlen, dbl_mb(k_dm1),
     $            dbl_mb(k_dm2), g_ifock, g_afock, g_gfock,
     $            g_coul, g_exch, hh)

       write(6,951) 
 951   format('Hessian diagonal')
       write(6,222) (hh((i-1)*orlen+i),i=1,orlen)
 222   format(10f12.6)
       call dsyev( 'V', 'L', orlen, hh, orlen, ev, scr, ssiz, info )
       write(6,953)
 953   format('Hessian eigenvalues')
       write(6,222) (ev(i),i=1,orlen)
c
c
c  Make explicit Hessian matrix from vector products
c
c$$$       call dfill((orlen*orlen),0.d0,hh,1)
c$$$       call mcscf_hessmake( geom, basis, nbf, nclosed, nact,
c$$$     $                      oskel, orlen, g_movecs, dbl_mb(k_dm1),
c$$$     $                      dbl_mb(k_dm2), g_ifock, g_afock,
c$$$     $                      g_gfock, g_coul, g_exch, g_x, g_prod,
c$$$     $                      hh )
c
c
 666   continue
c
c  Finite difference Hessian
c

c$$$       call mcscf_fdiff_hess( geom, basis, nbf, nclosed, nact,
c$$$     $                        oskel, tol2e, dbl_mb(k_dm1),
c$$$     $                        dbl_mb(k_dm2), g_movecs, g_coul,
c$$$     $                        g_grad )
c
c  Single call to Hessian vector product
c
c$$$       call mcscf_hessv( geom, basis, nbf, nclosed, nact,
c$$$     $                   oskel, dbl_mb(k_dm1), dbl_mb(k_dm2),
c$$$     $                   g_movecs, g_ifock, g_afock, g_gfock,
c$$$     $                   g_coul, g_exch, g_grad, g_prod )
c$$$       call ga_rowprint( '==== Product ====', g_prod )
c
c
c Free Coulomb and exchange operators
c
       if (.not.ga_destroy(g_exch))
     $      call errquit('rohf_head: cannot destroy exchange',0)
       if (.not.ga_destroy(g_coul))
     $      call errquit('rohf_head: cannot destroy Coulomb',0)
c
c Free MCSCF Gradient vectors and Fock matrices
c
       if (.not.ga_destroy(g_grad))
     $      call errquit('rohf_head: cannot destroy gradient',0)
       if (.not.ga_destroy(g_prod))
     $      call errquit('rohf_head: cannot destroy product',0)
       if (.not.ga_destroy(g_x))
     $      call errquit('rohf_head: cannot destroy product',0)
       if (.not.ga_destroy(g_afock))
     $      call errquit('rohf_head: cannot destroy MO matrix',0)
       if (.not.ga_destroy(g_ifock))
     $      call errquit('rohf_head: cannot destroy MO matrix',0)
       if (.not.ga_destroy(g_gfock))
     $      call errquit('rohf_head: cannot destroy MO matrix',0)
c
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c         
c              ROHF Test Section
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c  Allocate MCSCF Gradient vectors
c
       if (.not.ga_create(MT_DBL,orlen,1,'Gradient',1,1,g_grad))
     $      call errquit('rohf_head: cannot allocate',0)
       if (.not.ga_create(MT_DBL,orlen,1,'Product',1,1,g_prod))
     $      call errquit('rohf_head: cannot allocate',0)
       if (.not.ga_create(MT_DBL,orlen,1,'Arg vec',1,1,g_x))
     $      call errquit('rohf_head: cannot allocate',0)
c
c  Fock build and 1e-Hessian vector product
c
       if (.not.ga_create(MT_DBL,nbf,nbf,'CV Fock',nbf,1,g_fcv))
     $      call errquit('mcscf: cannot allocate active Fock',0)
       if (.not.ga_create(MT_DBL,nbf,nbf,'PV Fock',nbf,1,g_fpv))
     $      call errquit('mcscf: cannot allocate active Fock',0)
       if (.not.ga_create(MT_DBL,nbf,nbf,'CP Fock',nbf,1,g_fcp))
     $      call errquit('mcscf: cannot allocate active Fock',0)
c
c  ROHF Gradient
c
c$$$       call rohf_fock( geom, basis, nclosed, nopen, tol2e, g_movecs,
c$$$     $                 eone, etwo, g_fcv, g_fpv, g_fcp, oskel )
c$$$       call rohf_fock2grad( nbf, nclosed, nact,
c$$$     $                      g_fcv, g_fpv, g_fcp, g_grad )
c$$$       call ga_rowprint( '==== ROHF Gradient ====', g_grad )
c
c  ROHF Hessian vector product
c
c$$$       call ga_copy(g_grad, g_x)
c$$$       call ga_zero(g_prod)
c$$$       pflg = 2
c$$$       lshift = 0.d0
c$$$       call rohf_hessv_xx( basis, geom, nbf, nclosed, nact, pflg,
c$$$     $                     g_movecs, oskel, g_fcv, g_fpv, g_fcp,
c$$$     $                     tol2e, lshift, g_x, g_prod )
c$$$
c$$$       call ga_rowprint( '==== ROHF Product', g_prod )
c
c  Make ROHF Hessian
c
c$$$       call rohf_hessmake( basis, geom, nbf, nclosed, nact, 
c$$$     $                 g_movecs, oskel, g_fcv, g_fpv, g_fcp,
c$$$     $                 tol2e, g_x, g_prod )
c$$$
c
c  Create ROHF Hessian diagonal to compare with
c
c$$$       call rohf_hxxx( nbf, nclosed, nact, 0.d0, g_fcv, g_fpv,
c$$$     $                 g_fcp, g_prod )
c
c  Deallocate stuff
c
       if (.not.ga_destroy(g_fcv))
     $      call errquit('mcscf: cannot destroy MO vectors',0)
       if (.not.ga_destroy(g_fpv))
     $      call errquit('mcscf: cannot destroy MO vectors',0)
       if (.not.ga_destroy(g_fcp))
     $      call errquit('mcscf: cannot destroy MO vectors',0)
       if (.not.ga_destroy(g_grad))
     $      call errquit('mcscf: cannot destroy gradient',0)
       if (.not.ga_destroy(g_prod))
     $      call errquit('mcscf: cannot destroy product',0)
       if (.not.ga_destroy(g_x))
     $      call errquit('mcscf: cannot destroy product',0)
c
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c         
c              Cleanup
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
       call scf_tidy()
       call schwarz_tidy()
       call int_terminate()
       if (.not.bas_destroy(basis))
     $   call errquit('cannot free basis handle',0)
c
       if (.not.ma_pop_stack(l_sym))
     $      call errquit('mcscf: cannot pop stack?',0)
       if (.not.ma_pop_stack(l_dm1))
     $      call errquit('mcscf: cannot pop stack?',0)
       if (.not.ma_pop_stack(l_dm2))
     $      call errquit('mcscf: cannot pop stack?',0)
       if (.not.ma_pop_stack(l_occ))
     $      call errquit('mcscf: cannot pop stack?',0)
       if (.not.ma_pop_stack(l_evals))
     $      call errquit('mcscf: cannot pop stack?',0)
       if (.not.ga_destroy(g_movecs))
     $      call errquit('mcscf: cannot destroy MO vectors',0)
c
c
c
       call util_print_pop
       if (.not. context_pop('mcscf'))
     $   call errquit('context_pop failed',0)
       return
       end







       subroutine mcscf_getinfo( rtdb, geom, basis, nbf, nclosed, nact )
       implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
c
       integer rtdb
       integer geom
       integer basis
       integer nbf
       integer nclosed
       integer nact
c
       integer nactel
       integer nelec
       double precision charge
       double precision nuclear_charge
c
c  Create geom and basis handles
c
       if (.not. geom_create(geom, 'geometry'))
     $     call errquit('get_db_handles: geom_create?', 0)
       if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('get_db_handles: no geometry ', 0)
       if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('get_db_handles: bas_create?', 0)
       if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     $     call errquit('get_db_handles: no ao basis ', 0)
c
c  Basis info
c
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('scf_init: basis info',0)
c
c  Find number of electrons
c
      if (.not. rtdb_get(rtdb, 'charge', MT_DBL, 1, charge))
     $     charge = 0.0d0
      if (.not. geom_nuc_charge(geom, nuclear_charge))
     $     call errquit('scf: geom_nuc_charge failed', 0)
      nelec = nint(nuclear_charge - charge)
      if (nelec .le. 0) call errquit
     $     ('scf: negative no. of electrons ?', nelec)
      if (abs(nuclear_charge - charge - dble(nelec)) .gt. 1.0d-8)
     $     call errquit('mcscf: charge neglected?', 0)
       write(6,882) nelec
 882   format('rohf: got electrons info',i5)
       call util_flush(6)

c
c  Determine no. of closed and active shells
c  **** Have to correctly determine number of active from MCSCF input ****
c
       if (.not. rtdb_get(rtdb, 'mcscf:nact', MT_INT, 1, nact))
     $    call errquit('mcscf:no active orbitals',0)
       if (nact .gt. nelec) call errquit('mcscf: nact>nelec ', nact)
       if (.not. rtdb_get(rtdb, 'mcscf:nactelec', MT_INT, 1, nactel))
     $    call errquit('mcscf:no active electrons',0)
       nclosed = (nelec-nactel) / 2

       return
       end









       subroutine mcscf_hessmake( geom, basis, nbf, nclosed, nact,
     $                            oskel, orlen, g_movecs, dm1, dm2,
     $                            g_ifock, g_afock, g_gfock,
     $                            g_coul, g_exch,
     $                            g_x, g_prod, hh )
       implicit none
#include "mafdecls.fh"
       integer geom, basis
       integer nbf, nclosed, nact
       logical oskel
       integer orlen
       double precision dm1(*), dm2(*)
       integer g_movecs, g_ifock, g_afock, g_gfock
       integer g_coul, g_exch
       integer g_x, g_prod
       double precision hh(orlen,orlen)
c
c
C       INTEGER L_HD, K_HD
       integer l_hv, k_hv, l_hy, k_hy
       integer nvir, incr, xoff, i, j, ii, pflg
       double precision xx, lshift
       data pflg/2/
       data lshift/0.d0/
c
c
c
       if (.not.ma_push_get(MT_DBL, orlen, 'H', l_hv, k_hv))
     $     call errquit('rohf: cannot allocate',0)
       if (.not.ma_push_get(MT_DBL, orlen, 'Hy', l_hy, k_hy))
     $     call errquit('rohf: cannot allocate',0)
       nvir = nbf - nclosed - nact
C       goto 10000
       open(unit=88,file='hess',form='unformatted',
     $      status='unknown')
       write(88) nclosed,nact,nvir
       write(88) (nclosed*nvir),((nclosed+nact)*nvir),orlen
       xx = 1.d0
       incr = 0
c
c
c
       do j=1,nclosed+nact
         do i=1,nvir
           ii = (j-1)*nvir + i
           call ga_zero(g_x)
           call ga_put(g_x,ii,ii,1,1,xx,1)
           call mcscf_hessv( geom, basis, nbf, nclosed, nact,
     $                       .false., pflg, lshift, dm1, dm2,
     $                       g_movecs, g_ifock, g_afock, g_gfock,
     $                       g_coul, g_exch, g_x, g_prod )
           call ga_get(g_prod,1,orlen,1,1,dbl_mb(k_hv),orlen)
           call hv_writev(88,orlen,dbl_mb(k_hv))
           call dcopy(orlen,dbl_mb(k_hv),1,hh(1,ii),1)
         enddo
       enddo
c
       xoff = nvir*(nclosed+nact)
       do j=1,nclosed
         do i=1,nact
           ii = xoff + (j-1)*nact + i
           call ga_zero(g_x)
           call ga_put(g_x,ii,ii,1,1,xx,1)
           call mcscf_hessv( geom, basis, nbf, nclosed, nact,
     $                       .false., pflg, lshift, dm1, dm2,
     $                       g_movecs, g_ifock, g_afock, g_gfock,
     $                       g_coul, g_exch, g_x, g_prod )
           call ga_get(g_prod,1,orlen,1,1,dbl_mb(k_hv),orlen)
           call hv_writev(88,orlen,dbl_mb(k_hv))
           call dcopy(orlen,dbl_mb(k_hv),1,hh(1,ii),1)
         enddo
       enddo
c
c
c
c$$$       if (.not.ma_push_get(MT_DBL, orlen, 'Hy', l_hd, k_hd))
c$$$     $     call errquit('rohf: cannot allocate',0)
c$$$c
c$$$c
c$$$c
c$$$       open(unit=11,file='hessian.ascii',form='formatted',
c$$$     $      status='unknown')
c$$$           call hv_write(11,incr,orlen,dbl_mb(k_hv))
c$$$           call hv_writev(12,orlen,dbl_mb(k_hv))
c$$$       open(unit=12,file='hessian',form='unformatted',
c$$$     $      status='unknown')
c$$$       close(11)
c$$$       close(12)

c$$$10000  continue
c$$$       xx = 1.d0
c$$$       do i=1,orlen
c$$$         call ga_zero(g_x)
c$$$         call ga_put(g_x,i,i,1,1,xx,1)
c$$$         call mcscf_hessv( geom, basis, nbf, nclosed, nact,
c$$$     $                     oskel, dm1, dm2, g_movecs,
c$$$     $                     g_ifock, g_afock, g_gfock,
c$$$     $                     g_coul, g_exch, g_x, g_prod )
c$$$         call ga_get(g_prod,i,i,1,1,dbl_mb(k_hd+i-1),1)
c$$$       enddo
c$$$       write(6,900) 
c$$$ 900   format('Exact Hessian diagonal')
c$$$       write(6,901) (dbl_mb(k_hd+i-1),i=1,orlen)
c$$$ 901   format(10f12.6)

c$$$       if (.not.ma_pop_stack(l_hd))
c$$$     $      call errquit('rohf: damn',0)
c
c  Clean up
c
       if (.not.ma_pop_stack(l_hy))
     $      call errquit('rohf: damn',0)
       if (.not.ma_pop_stack(l_hv))
     $      call errquit('rohf: damn',0)

       return
       end



      
c
c  Make J and K integrals for closed + active operators
c  and dump to disk
c
       subroutine makeJK( basis, nbf, nclosed, nact, g_movecs )
       implicit none
#include "global.fh"
       integer basis, nbf, nclosed, nact, g_movecs
       integer mo_lo, mo_hi, noper, nn, i
       integer g_coul, g_exch
       double precision tmp(1000)
       logical osym
       data osym/.false./
       integer ga_create_JKblocked
       external ga_create_JKblocked

       mo_lo = 1
       mo_hi = (nclosed+nact)
       noper = (mo_hi*(mo_hi+1))/2
       nn = nbf*nbf
       g_coul = ga_create_JKblocked(noper,nbf,nbf,'Coulomb Oper')
       g_exch = ga_create_JKblocked(noper,nbf,nbf,'X Oper')
       call moints_build_6x( basis, osym,
     $                       mo_lo, mo_lo, mo_hi, 1, nbf,
     $                       g_movecs, g_coul, .true.,
     $                       g_exch, .true., 16 )
       open(unit=88,file='JKints',status='unknown',form='unformatted')
       do i=1,noper
         call ga_get(g_coul,1,nn,i,i,tmp,1)
         call hv_writev(88,nn,tmp)
       enddo
       do i=1,noper
         call ga_get(g_exch,1,nn,i,i,tmp,1)
         call hv_writev(88,nn,tmp)
       enddo
       close(88)
       if (.not. ga_destroy(g_coul)) call errquit('mcscf: ga?',0)
       if (.not. ga_destroy(g_exch)) call errquit('mcscf: ga?',0)
       return
       end









