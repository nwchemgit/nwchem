      subroutine grad_force ( rtdb, basis, geom )
c$Id: grad_force.F,v 1.13 1996-11-02 00:32:58 vg038 Exp $
C calculate energy gradients with respect to nuclear coordinates
C------------------------------------------------------------------------------
C         ___                 ___                         ___
C dE      \         dh(i,j)   \             d(mn|op)      \          dS(i,j)
C -- = 2   > D(i,j) ------- +  > P(m,n,o,p) -------- - 2   > Dw(i,j) -------
C dA      /           dA      /                dA         /            dA
C         ---                 ---                         ---
C         i,j                ijkl                         i,j
C
C
C        dV(nuc-nuc)
C      + ----------
C            dA
C------------------------------------------------------------------------------
C RHF: (H.Schlegel & S. Wolfe, JCP 63, p3632)
C                                1
C P(i,j,k,l) = [2 D(i,j)D(k,l) - - ( D(i,k)D(j,l) + D(i,l)D(j,k) ]
C                                2
C------------------------------------------------------------------------------
C UHF: (M. Dupuis & H.F. King, JCP 68, p3998)
C                                  1 
C P(i,j,k,l) = [2 D+(i,j)D+(k,l) - - ( D+(i,k)D+(j,l) + D-(i,k)D-(j,l) 
C                                  2
C                                    + D+(i,l)D+(j,k) + D-(i,l)D-(j,k) ) ]
C    D == D+
C------------------------------------------------------------------------------
C ROHF: (Y. Yamaguchi, Y. Osamura, J.D. Goddard, H.F. Schaefer:
C    'A New Dimension to Quantum Chemistry', Oxford University Press '94, p74)
C
C            -         I    J               I    J       I    J
C P(ijkl) =  > [ a(IJ)D(ij)D(kl) + b(IJ) ( D(ik)D(jl) + D(il)D(jk) ) ]
C            -
C            IJ
C     I,J = open, closed
C------------------------------------------------------------------------------
C UMP2: (M.J. Frisch, M. Head-Gordon, J.A. Pople, CPL 166, Nr. 3, p275,
C explicit spin formalism (alpha/beta) from Robert
C             1
C P(ijkl) = [ - ( (Da(ij)+Db(ij))*(Pa(kl)+Pb(kl))  
C             4
C               + (Pa(ij)+Pb(ij))*(Da(kl)+Db(kl)) )
C             1
C           - - ( Da(il)*Pa(jk) + Db(il)*Pb(jk) + Pa(il)*Da(jk) + Pb(il)*Db(jk)
C             8
C               + Da(jl)*Pa(ik) + Db(jl)*Pb(ik) + Pa(jl)*Da(ik) + Pb(jl)*Db(ik)
C
C               ) ]
C------------------------------------------------------------------------------
      implicit none

#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
#include "util.fh"
#include "sym.fh"
#include "stdio.fh"
#include "msgids.fh"

      integer rtdb, basis, geom

#ifdef OLDTIMERS
#ifdef SUN
      real t0, tim1(3), tim2(3), dtime
#endif
#ifdef KSR
      real tim1, tim2
#endif
#ifdef __PARAGON__
      double precision tim1, tim2
#endif
#endif
      double precision cpu_tim(2), wall_tim(2)

      integer ga_create_atom_blocked
      external ga_create_atom_blocked

      logical movecs_read_header, movecs_read
      external movecs_read_header, movecs_read

      integer g_dens,     ! density matrices (up to 6)
     $     g_wdens,       ! energy weighted density
     $     g_eigen_diag,  ! eigenvalue diagonal matrix or lagrangian
     $     g_force        ! forces on atoms (cartesian)

      character*255 movecs_in

      integer nproc, nat, nbf, max_sh_bf, max_at_bf, nopen, nclosed, 
     $     ndens, nbf2, nocc, max1e, max2e, mscratch_1e, mscratch_2e,
     $     lbuf, lscratch, lsqa, i, j

      integer lforce, l_force, k_force, l_evals, k_evals, l_occ, k_occ,
     $     l_sch, k_sch, l_act, k_act, nactive, l_actlist, k_actlist,
     $     l_buf, k_buf, l_scr, k_scr, l_dens, k_dens, l_wdens, k_wdens,
     $     k_frc_nuc, k_frc_kin, k_frc_2el, k_frc_wgh, 
     $     l_frc_nuc, l_frc_kin, l_frc_2el, l_frc_wgh,
     $     l_frc_mp2, k_frc_mp2

      integer ld_ij, ld_kl, lh_ik, ld_ik, lh_jl, ld_jl, lh_il, ld_il,
     $     lh_jk, ld_jk, 
     $     lh_ij2, ld_ij2, lh_kl2, ld_kl2, lh_ik2, ld_ik2,
     $     lh_jl2, ld_jl2, lh_il2, ld_il2, lh_jk2, ld_jk2,
     $     lh_ij3, ld_ij3, lh_kl3, ld_kl3, lh_ik3, ld_ik3,
     $     lh_jl3, ld_jl3, lh_il3, ld_il3, lh_jk3, ld_jk3,
     $     lh_ij4, ld_ij4, lh_kl4, ld_kl4, lh_ik4, ld_ik4,
     $     lh_jl4, ld_jl4, lh_il4, ld_il4, lh_jk4, ld_jk4

      double precision crd,  ! atomic coordinates
     $                 charge, charge_nuc, tol2e, q, gnorm, gmax

      double precision grad_norm, grad_max
      external grad_norm, grad_max

      dimension crd(3), g_dens(6) ! max # of density matrices

      character*16 tag
      character*20 task
      character*4 scftype
      character*32 rtdb_string
      logical oskel, omp2, status

      nproc = ga_nnodes()

C     get information about basis set
      status = geom_ncent(geom,nat)
      if ( .not.status ) call errquit('grad: could not get natoms',0)
      status = bas_numbf(basis,nbf)
      if ( .not.status ) call errquit('grad: could not get nbf',0)
      status = bas_nbf_cn_max(basis,max_sh_bf)
      if ( .not.status ) call errquit('grad: could not get max_sh_bf',0)
      status = bas_nbf_ce_max(basis,max_at_bf)
      if ( .not.status ) call errquit('grad: could not get max_at_bf',0)

      status = rtdb_get(rtdb, 'charge', MT_DBL,1, charge)
      if ( .not.status ) charge = 0.D0
      status = geom_nuc_charge(geom, charge_nuc)
      if ( .not.status ) call errquit('grad: no nuclear charge',0)

C     get SCF MO vectors for density
      if (.not. rtdb_cget(rtdb, 'scf:input vectors', 1, movecs_in))
     $     call errquit('gradients: SCF MO vectors not defined',0)

C     get information about type of calculation

C     scftype: RHF, ROHF or UHF
      if (.not. rtdb_get(rtdb, 'scf:nopen', MT_INT, 1, nopen))
     $     nopen = 0
      if (.not. rtdb_get(rtdb, 'scf:nclosed', MT_INT, 1, nclosed))
     $     call errquit ('gradients: nopen > nelec ', nopen)

      if (.not. rtdb_cget(rtdb, 'scf:scftype', 1, scftype))
     $     call errquit('gradients: scftype not defined',0)
      call inp_ucase(scftype)

C     is it MP2?
      status = rtdb_cget(rtdb, 'task:theory', 1, task)
      if (      (task .eq. 'mp2') 
     $     .or. (task .eq. 'semi_dir_mp2')
     $     .or. (task .eq. 'direct_mp2') ) then
        omp2 = .true.
      else
        omp2 = .false.
      endif

      if ( .not. omp2 ) then
        if ( ga_nodeid() .eq. 0 ) then
          if (util_print('information', print_medium)) then
            write(luout,1) charge, scftype
 1          format(/
     $           '  charge          = ', f6.2/
     $           '  wavefunction    =   ', a/)
            call util_flush(luout)
          endif
        endif
      endif

C     # of eigenvalues and density matrices
      if ( omp2 ) then
        if ( scftype .eq. 'UHF' ) then
          nbf2 = 2 * nbf
          ndens = 5
        else if ( scftype .eq. 'RHF' ) then
          nbf2 = nbf
          ndens = 3
        else    ! ROHF
          nbf2 = nbf
          ndens = 3
          call errquit ('no ROMP2 gradients yet', 110)
        endif
      else   ! SCF
        if ( scftype .eq. 'UHF' ) then
          nbf2 = 2 * nbf
          ndens = 2
        else if ( scftype .eq. 'RHF' ) then
          nbf2 = nbf
          ndens = 1
        else    ! ROHF
          nbf2 = nbf
          ndens = 3
        endif
      endif

C  allocate and initialize global and local memory

C     energy weighted density (NxN)
      g_wdens = ga_create_atom_blocked(geom, basis, 'weighted density')

C     forces on atoms (3xnat)
      status = ga_create ( MT_DBL, 3, nat, 'forces', 3, 1, g_force )
      call ga_zero ( g_force )

C     local replication (separate for the different pieces)
      lforce = nat * 3
      status = MA_alloc_get(MT_DBL, lforce, 'forces', l_force, k_force)
      if (.not.status) call errquit('could not allocate l_force',1)
      status = MA_alloc_get(MT_DBL, lforce, 'forces', l_frc_nuc, 
     $     k_frc_nuc)
      if (.not.status) call errquit('could not allocate l_force',1)
      status = MA_alloc_get(MT_DBL, lforce, 'forces', l_frc_kin, 
     $     k_frc_kin)
      if (.not.status) call errquit('could not allocate l_force',1)
      status = MA_alloc_get(MT_DBL, lforce, 'forces', l_frc_wgh, 
     $     k_frc_wgh)
      if (.not.status) call errquit('could not allocate l_force',1)
      status = MA_alloc_get(MT_DBL, lforce, 'forces', l_frc_2el, 
     $     k_frc_2el)
      if (.not.status) call errquit('could not allocate l_force',1)
      if ( omp2 ) then
        status = MA_alloc_get(MT_DBL, lforce, 'forces', l_frc_mp2, 
     $       k_frc_mp2)
        if (.not.status) call errquit('could not allocate l_force',1)
      endif

C     eigenvalue diagonal matrix
      nocc = nopen + nclosed
      status = ga_create ( MT_DBL, nocc, nocc,'eigen_diag', 1, nocc, 
     $                       g_eigen_diag )
      call ga_zero ( g_eigen_diag )

C     eigenvalues
      status = MA_alloc_get(MT_DBL, nbf2,'MO evals', l_evals, k_evals)
      if (.not.status) 
     $     call errquit('grad: could not allocate l_evals',1)
      
C     occupation numbers
      status = MA_alloc_get(MT_DBL, nbf2,'occ. numbers', l_occ, k_occ)
      if (.not.status) 
     $     call errquit('grad: could not allocate l_occ',1)

C     global density
      do i=1, ndens
        g_dens(i) = ga_create_atom_blocked ( geom, basis, 
     $              'density matrix' )
      enddo

C     scratch array for Schwarz screening
      status = MA_alloc_get(MT_DBL, nat, 'schwarz limits', l_sch, k_sch)
      if (.not.status) 
     $        call errquit('grad: could not allocate l_sch',1)

C     lookup table and list of active atoms
      status = MA_alloc_get(MT_LOG, nat, 'active atoms', l_act, k_act)
      if (.not.status) 
     $        call errquit('grad: could not allocate l_act',1)

C     symmetry or not
      if (.not. rtdb_get(rtdb, 'gradients:use symmetry', mt_log, 1,
     $     oskel)) then
         if (.not. rtdb_get(rtdb, 'scf:skeleton', mt_log, 1,
     $        oskel)) then
            oskel = sym_number_ops(geom) .gt. 0
         endif
      endif
c
      if (oskel .and. ga_nodeid().eq.0 .and. 
     $     util_print('information', print_default)) then
         write(luout,*) ' Using symmetry'
         call util_flush(luout)
      endif
c
      status = rtdb_get ( rtdb, 'geometry:nactive', MT_INT, 1, nactive )
      if ( status ) then ! there is a list of active atoms
        if ( nactive.le.0 .or. nactive.gt.nat )
     $       call errquit('invalid number of active atoms',110)
        status = MA_alloc_get(MT_INT, nactive, 'active list', l_actlist, 
     $       k_actlist)
        if (.not.status) 
     $       call errquit('grad: could not allocate l_act',1)
        status = rtdb_get ( rtdb, 'geometry:actlist', MT_INT, nactive, 
     $       int_mb(k_actlist) )
        if (.not.status) 
     $       call errquit('active atom list missing',110)
C       fill lookup table
        do 10, i=0, nat-1
          log_mb(k_act+i) = .false.
 10     continue
        do 12, i=0, nactive-1
          log_mb ( k_act + int_mb(k_actlist + i) - 1 ) = .true.
 12     continue
        status = MA_free_heap ( l_actlist )
      else      ! all atoms are active
        do 15, i=0, nat-1
          log_mb(k_act+i) = .true.
 15     continue
      endif

      call grad_dens (geom, basis, g_dens, g_wdens, g_eigen_diag, 
     $                dbl_mb(k_evals), dbl_mb(k_occ), ndens, nbf, nbf2, 
     $                nopen, nclosed, scftype, movecs_in, omp2)

C     free temporary arrays
      status = ga_destroy ( g_eigen_diag )
      if (.not.status) 
     $     call errquit('grad: could not destroy g_eigen_diag',1)
      status = MA_free_heap ( l_occ )
      status = MA_free_heap ( l_evals )

c
C  initialize for integral gradients
      call int_init ( rtdb, 1, basis )
      call schwarz_init ( geom, basis )
      call int_terminate()
      call intd_init ( rtdb, 1, basis )
      call dfill ( lforce, 0.0D0, dbl_mb(k_force), 1 )
      call dfill ( lforce, 0.0D0, dbl_mb(k_frc_nuc), 1 )
      call dfill ( lforce, 0.0D0, dbl_mb(k_frc_kin), 1 )
      call dfill ( lforce, 0.0D0, dbl_mb(k_frc_wgh), 1 )
      call dfill ( lforce, 0.0D0, dbl_mb(k_frc_2el), 1 )
      call int_mem ( max1e, max2e, mscratch_1e, mscratch_2e )
      lbuf = max ( max1e, max2e )
      lscratch = max ( mscratch_1e, mscratch_2e )

C  one-electron contribution
C     buffers for one electron integral derivatives
      status = MA_alloc_get(MT_DBL, lbuf, 'deriv buffer', l_buf, k_buf)
      if (.not.status) call errquit('could not allocate buffer',1)
      status = MA_alloc_get(MT_DBL, lscratch, 'deriv scratch', 
     $         l_scr, k_scr)
      if (.not.status) 
     $        call errquit('grad: could not allocate scratch',1)

C     local density matrix block
      lsqa = max_at_bf * max_at_bf
      status = MA_alloc_get(MT_DBL,lsqa,'local_density',l_dens,k_dens)
      if (.not.status) call errquit('could not allocate l_dens',1)
      status = MA_alloc_get(MT_DBL,lsqa,'local_w_density',l_wdens, 
     $                      k_wdens)
      if (.not.status) call errquit('could not allocate l_wdens',1)

#ifdef OLDTIMERS
#ifdef SUN
      tim1(3) = dtime( tim1 )
      t0 = real(tcgtime())
#endif
#ifdef KSR
      tim1 = all_seconds()
#endif
#ifdef __PARAGON__
      tim1 = dclock()
#endif
#endif
C      goto 112
      cpu_tim(1)  = util_cpusec()
      wall_tim(1) = util_wallsec()
      call grad1( dbl_mb(k_buf), lbuf, dbl_mb(k_scr), lscratch,
     $     dbl_mb(k_dens), dbl_mb(k_wdens), dbl_mb(k_frc_nuc),
     $     dbl_mb(k_frc_kin), dbl_mb(k_frc_wgh), g_force,
     $     g_dens, g_wdens, basis, geom, nproc, nat, max_at_bf,
     $     rtdb, oskel)
      cpu_tim(1)  = util_cpusec() - cpu_tim(1)
      wall_tim(1) = util_wallsec() - wall_tim(1)
 112  continue
#ifdef OLDTIMERS
#ifdef SUN
      tim1(3) = dtime( tim1 )
      tim1(3) = real(tcgtime()) - t0
#endif
#ifdef KSR
      tim1 = all_seconds() - tim1
#endif
#ifdef __PARAGON__
      tim1 = dclock() - tim1
#endif
#endif

C  two-electron contribution
C  allocate arrays for two-electron integral stuff (or rename old ones)
C  blocks of density matrix:
      ld_ij = k_dens
      ld_kl = k_wdens
      status = MA_alloc_get(MT_DBL, lsqa, 'ld_ik', lh_ik, ld_ik)
      if (.not.status) call errquit('could not allocate ld_ik',1)
      status = MA_alloc_get(MT_DBL, lsqa, 'ld_jl', lh_jl, ld_jl)
      if (.not.status) call errquit('could not allocate ld_jl',1)
      status = MA_alloc_get(MT_DBL, lsqa, 'ld_il', lh_il, ld_il)
      if (.not.status) call errquit('could not allocate ld_il',1)
      status = MA_alloc_get(MT_DBL, lsqa, 'ld_jk', lh_jk, ld_jk)
      if (.not.status) call errquit('could not allocate ld_jk',1)

      if ( scftype .ne. 'RHF' .or. omp2 ) then ! UHF or ROHF (or MP2)
        if ( scftype .eq. 'ROHF' .or. omp2 ) then
          status = MA_alloc_get(MT_DBL, lsqa, 'ld_ij2', lh_ij2, ld_ij2)
          if (.not.status) call errquit('could not allocate ld_ij2',1)
          status = MA_alloc_get(MT_DBL, lsqa, 'ld_kl2', lh_kl2, ld_kl2)
          if (.not.status) call errquit('could not allocate ld_kl2',1)
        endif          
        status = MA_alloc_get(MT_DBL, lsqa, 'ld_ik2', lh_ik2, ld_ik2)
        if (.not.status) call errquit('could not allocate ld_ik2',1)
        status = MA_alloc_get(MT_DBL, lsqa, 'ld_jl2', lh_jl2, ld_jl2)
        if (.not.status) call errquit('could not allocate ld_jl2',1)
        status = MA_alloc_get(MT_DBL, lsqa, 'ld_il2', lh_il2, ld_il2)
        if (.not.status) call errquit('could not allocate ld_il2',1)
        status = MA_alloc_get(MT_DBL, lsqa, 'ld_jk2', lh_jk2, ld_jk2)
        if (.not.status) call errquit('could not allocate ld_jk2',1)
      endif

      if ( scftype .eq. 'UHF' .and. omp2 ) then
        status = MA_alloc_get(MT_DBL, lsqa, 'ld_ij3', lh_ij3, ld_ij3)
        if (.not.status) call errquit('could not allocate ld_ij3',1)
        status = MA_alloc_get(MT_DBL, lsqa, 'ld_kl3', lh_kl3, ld_kl3)
        if (.not.status) call errquit('could not allocate ld_jl3',1)
        status = MA_alloc_get(MT_DBL, lsqa, 'ld_ik3', lh_ik3, ld_ik3)
        if (.not.status) call errquit('could not allocate ld_ik3',1)
        status = MA_alloc_get(MT_DBL, lsqa, 'ld_jl3', lh_jl3, ld_jl3)
        if (.not.status) call errquit('could not allocate ld_jl3',1)
        status = MA_alloc_get(MT_DBL, lsqa, 'ld_il3', lh_il3, ld_il3)
        if (.not.status) call errquit('could not allocate ld_il3',1)
        status = MA_alloc_get(MT_DBL, lsqa, 'ld_jk3', lh_jk3, ld_jk3)
        if (.not.status) call errquit('could not allocate ld_jk3',1)
        status = MA_alloc_get(MT_DBL, lsqa, 'ld_ij4', lh_ij4, ld_ij4)
        if (.not.status) call errquit('could not allocate ld_ij4',1)
        status = MA_alloc_get(MT_DBL, lsqa, 'ld_kl4', lh_kl4, ld_kl4)
        if (.not.status) call errquit('could not allocate ld_kl4',1)
        status = MA_alloc_get(MT_DBL, lsqa, 'ld_ik4', lh_ik4, ld_ik4)
        if (.not.status) call errquit('could not allocate ld_ik4',1)
        status = MA_alloc_get(MT_DBL, lsqa, 'ld_jl4', lh_jl4, ld_jl4)
        if (.not.status) call errquit('could not allocate ld_jl4',1)
        status = MA_alloc_get(MT_DBL, lsqa, 'ld_il4', lh_il4, ld_il4)
        if (.not.status) call errquit('could not allocate ld_il4',1)
        status = MA_alloc_get(MT_DBL, lsqa, 'ld_jk4', lh_jk4, ld_jk4)
        if (.not.status) call errquit('could not allocate ld_jk4',1)
      endif

C  define threshold for Schwarz screening (same as in SCF)
      status = rtdb_get(rtdb, 'scf:tol2e', MT_DBL, 1, tol2e)
      if (.not.status) tol2e = 1.0d-9

#ifdef OLDTIMERS
#ifdef SUN
      tim2(3) = dtime( tim2 )
      t0 = real(tcgtime())
#endif
#ifdef KSR
      tim2 = all_seconds()
#endif
#ifdef __PARAGON__
      tim2 = dclock()
#endif
#endif

C      goto 111
      cpu_tim(2)  = util_cpusec()
      wall_tim(2) = util_wallsec()
      call grad2(dbl_mb(ld_ij),dbl_mb(ld_kl),dbl_mb(ld_ik),
     $           dbl_mb(ld_jl),dbl_mb(ld_il),dbl_mb(ld_jk), 
     $           dbl_mb(ld_ij2),dbl_mb(ld_kl2),  !    only
     $           dbl_mb(ld_ik2),dbl_mb(ld_jl2),  !    for
     $           dbl_mb(ld_il2),dbl_mb(ld_jk2),  !  UHF,ROHF
     $           dbl_mb(ld_ij3),dbl_mb(ld_kl3),    !
     $           dbl_mb(ld_ik3),dbl_mb(ld_jl3),    !
     $           dbl_mb(ld_il3),dbl_mb(ld_jk3),    ! only
     $           dbl_mb(ld_ij4),dbl_mb(ld_kl4),    ! UMP2
     $           dbl_mb(ld_ik4),dbl_mb(ld_jl4),    !
     $           dbl_mb(ld_il4),dbl_mb(ld_jk4),    !
     $           dbl_mb(k_frc_2el), g_dens, g_force, max_at_bf,
     $           rtdb, geom, basis, nproc, nat,
     $           dbl_mb(k_scr), lscratch, dbl_mb(k_buf), lbuf, tol2e,
     $           dbl_mb(k_sch), log_mb(k_act), oskel, scftype, omp2 )
      cpu_tim(2)  = util_cpusec()  - cpu_tim(2)
      wall_tim(2) = util_wallsec() - wall_tim(2)
 111  continue
#ifdef OLDTIMERS
#ifdef SUN
      tim2(3) = dtime( tim2 )
      tim2(3) = real(tcgtime()) - t0
#endif
#ifdef KSR
      tim2 = all_seconds() - tim2
#endif
#ifdef __PARAGON__
      tim2 = dclock() - tim2
#endif
#endif

c
c
c     terminate integral scope
c
      call schwarz_tidy ()
      call intd_terminate()

C     print forces
      call ga_dgop ( msg_grad_nuc, dbl_mb(k_frc_nuc), 3*nat, '+' )
      call ga_dgop ( msg_grad_wgh, dbl_mb(k_frc_wgh), 3*nat, '+' )
      call ga_dgop ( msg_grad_kin, dbl_mb(k_frc_kin), 3*nat, '+' )
      call ga_sync()
      status = rtdb_parallel (.false.)
      if ( ga_nodeid() .eq. 0 ) then
        call ga_get ( g_force, 1, 3, 1, nat, dbl_mb(k_frc_2el), 3 )

        do i=0, 3*nat-1
          dbl_mb(k_force+i) = dbl_mb(k_frc_2el+i) + dbl_mb(k_frc_nuc+i) 
     $                      + dbl_mb(k_frc_wgh+i) + dbl_mb(k_frc_kin+i)
        enddo

        if ( omp2 ) then
          status = rtdb_get ( rtdb, 'mp2:nonseparable gradient', mt_dbl,
     $         3*nat, dbl_mb(k_frc_mp2) )
          if (.not. status ) 
     $         call errquit ('grad_force: no nonseparable gradient',110)
          do i=0, 3*nat-1
            dbl_mb(k_force+i) = dbl_mb(k_force+i) + dbl_mb(k_frc_mp2+i) 
          enddo
        endif

C       zero force contributions on inactive atoms
        call zero_forces ( dbl_mb(k_force), log_mb(k_act), nat )
        
C       symmetrize
        if (oskel) then
          call sym_grad_symmetrize(geom, dbl_mb(k_force))
        endif

        if ( util_print('forces', print_high) ) then
          write (luout,2200) 'nuclear repulsion gradient',
     $         ((dbl_mb(k_frc_nuc+i-1+3*(j-1)),i=1,3),j=1,nat)
          write (luout,2200) 'weighted density gradient',
     $         ((dbl_mb(k_frc_wgh+i-1+3*(j-1)),i=1,3),j=1,nat)
          write (luout,2200) 'kinetic energy gradient',
     $         ((dbl_mb(k_frc_kin+i-1+3*(j-1)),i=1,3),j=1,nat)
          write (luout,2200) '2-electron gradient',
     $         ((dbl_mb(k_frc_2el+i-1+3*(j-1)),i=1,3),j=1,nat)
          if ( omp2 ) then
            write (luout,2200) 'non-separable gradient',
     $         ((dbl_mb(k_frc_mp2+i-1+3*(j-1)),i=1,3),j=1,nat)
            write (luout,2200) 
     $           'total MP2 gradient',
     $           ((dbl_mb(k_force+i-1+3*(j-1)),i=1,3),j=1,nat)
          else
            write (luout,2200) 'total SCF gradient',
     $           ((dbl_mb(k_force+i-1+3*(j-1)),i=1,3),j=1,nat)
          endif
 2200     format(A,/,1000(3(1x,F12.6),/))
          call util_flush(luout)
        endif

C       store in rtdb
        if ( omp2 ) then
          rtdb_string = 'mp2:gradient'
        else  ! scf
          rtdb_string = 'scf:gradient'
        endif
        status = rtdb_put ( rtdb, rtdb_string, MT_DBL, 3*nat, 
     $       dbl_mb(k_force))
        if ( .not. status ) 
     $       call errquit('could not store gradients',1)
      endif
      status = rtdb_parallel (.true.)

C     default: print the total forces 
      if ( ga_nodeid() .eq. 0 
     $     .and. util_print('forces', print_low) ) then
        write (luout,1000)'x','y','z','x','y','z'
        do 30, i=1, nat
          status = geom_cent_get ( geom, i, tag, crd, q )
          write (luout,2000) i, tag, (crd(j),j=1,3),
     $         (dbl_mb(k_force+3*(i-1)+j),j=0,2)
 30     continue
        write(luout,*)
 1000   format (/,/,30X,'ENERGY GRADIENTS',/,/,4X,'atom',15X,
     $       'coordinates',
     $       24X,'gradient',/,6X,2(1X,(3(10X,A1))))
 2000   format (1X,I3,1X,A4,2(1X,3(1X,F10.6)))
        call util_flush(luout)
      endif

C     forces for ECCE
      if ( ga_nodeid() .eq. 0 ) then
        gnorm = grad_norm ( dbl_mb(k_force), nat )
        gmax  = grad_max  ( dbl_mb(k_force), nat )
        call ecce_print2 ('gradients', mt_dbl, dbl_mb(k_force), 3,3,nat)
        call ecce_print1 ('gradient norm', mt_dbl, 1, gnorm )
        call ecce_print1 ('gradient max', mt_dbl, 1, gmax )
      endif

C     print timing information
      if (ga_nodeid().eq.0 .and.
     $     util_print('timing', print_default)) then
        write(luout,03000)cpu_tim,wall_tim
03000   format(17x,40('-'),/,
     &      17x,'|  Time  |  1-e (secs)  |  2-e (secs)  |',/,
     &      17x,40('-'),/,
     &      17x,'|  CPU   |',f11.2,3x,'|',f11.2,3x,'|',/,
     &      17x,40('-'),/,
     &      17x,'|  WALL  |',f11.2,3x,'|',f11.2,3x,'|',/,
     &      17x,40('-'))
#ifdef OLDTIMERS
#ifdef SUN
        write (luout,3000) '1-el',(tim1(j),j=1,3),'2-el',
     &   (tim2(k),k=1,3)
 3000   format (/,20X,'TIMING INFORMATION',/,/,23X,'CPU',15X,
     $        'SYS',15X,'elapsed',/,2(A10,3(3X,F15.2),/))
#endif
#if defined(KSR) || defined(__PARAGON__)
        write (luout,3001) '1-el', tim1, '2-el', tim2
 3001   format (/,/,7X,'ELAPSED TIME (SECS)',/,/,2(A10,F15.2,/))
#endif
#endif
        call util_flush(luout)
      endif

C  free memory

      status = MA_free_heap ( l_act )
      if ( scftype .ne. 'RHF' .or. omp2 ) then
        if ( scftype .eq. 'ROHF' .or. omp2 ) then
          status = MA_free_heap ( lh_kl2 )
          status = MA_free_heap ( lh_ij2 )
        endif
        status = MA_free_heap ( lh_jk2 )
        status = MA_free_heap ( lh_il2 )
        status = MA_free_heap ( lh_jl2 )
        status = MA_free_heap ( lh_ik2 )
        if ( scftype .eq. 'ROHF' ) then
          status = MA_free_heap ( lh_ij2 )
          status = MA_free_heap ( lh_kl2 )
        endif
      endif
      if ( scftype .eq. 'UHF' .and. omp2 ) then
        status = MA_free_heap ( lh_ij3 )
        status = MA_free_heap ( lh_kl3 )
        status = MA_free_heap ( lh_ik3 )
        status = MA_free_heap ( lh_jl3 )
        status = MA_free_heap ( lh_il3 )
        status = MA_free_heap ( lh_jk3 )
        status = MA_free_heap ( lh_ij4 )
        status = MA_free_heap ( lh_kl4 )
        status = MA_free_heap ( lh_ik4 )
        status = MA_free_heap ( lh_jl4 )
        status = MA_free_heap ( lh_il4 )
        status = MA_free_heap ( lh_jk4 )
      endif
      status = MA_free_heap ( lh_jk )
      status = MA_free_heap ( lh_il )
      status = MA_free_heap ( lh_jl )
      status = MA_free_heap ( lh_ik )
      status = MA_free_heap ( l_wdens )
      status = MA_free_heap ( l_dens )
      status = MA_free_heap ( l_scr )
      status = MA_free_heap ( l_buf )
      status = MA_free_heap ( l_sch )
      status = MA_free_heap ( l_frc_2el )
      status = MA_free_heap ( l_frc_wgh )
      status = MA_free_heap ( l_frc_kin )
      status = MA_free_heap ( l_frc_nuc )
      status = MA_free_heap ( l_force )
      if ( omp2 ) then
        status = MA_free_heap ( l_frc_mp2 )
      endif
      call ga_sync()
      status = ga_destroy ( g_force )
      status = ga_destroy ( g_wdens )
      do i=1, ndens
        status = ga_destroy ( g_dens(i) )
        if ( .not. status ) call errquit ('error destroying density', 1)
      enddo

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function grad_norm ( frc, nat )
      implicit none

      integer nat
      double precision frc(3*nat)

      integer i
      double precision g2

      g2 = 0
      do i=1, 3*nat
        g2 = g2 + frc(i)*frc(i)
      enddo

      grad_norm = sqrt(g2)

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function grad_max (frc, nat )
      implicit none

      integer nat
      double precision frc(3,nat)

      integer i
      double precision norm, mxnorm

      mxnorm = 0
      do i=1, nat
        norm = frc(1,i)*frc(1,i) + frc(2,i)*frc(2,i) + frc(3,i)*frc(3,i) 
        mxnorm = max(norm,mxnorm)
      enddo

      grad_max = sqrt(mxnorm)

      return
      end
