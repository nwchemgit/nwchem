      subroutine grad_force ( rtdb, basis, geom, movecs )
c$Id: grad_force.F,v 1.6 1996-08-29 01:36:40 vg038 Exp $
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
      implicit none

#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
#include "util.fh"
#include "sym.fh"
#include "stdio.fh"

      integer rtdb, basis, geom
      character*(*) movecs

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

      integer g_dens,     ! density matrices (up to 3)
     $     g_wdens,       ! energy weighted density
     $     g_eigen_diag,  ! eigenvalue diagonal matrix or lagrangian
     $     g_force        ! forces on atoms (cartesian)

      integer nproc, nat, nbf, max_sh_bf, max_at_bf, nopen, nclosed, 
     $     ndens, nbf2, nocc, max1e, max2e, mscratch_1e, mscratch_2e,
     $     lbuf, lscratch, lsqa, i, j

      integer lforce, l_force, k_force, l_evals, k_evals, l_occ, k_occ,
     $     l_sch, k_sch, l_act, k_act, nactive, l_actlist, k_actlist,
     $     l_buf, k_buf, l_scr, k_scr, l_dens, k_dens, l_wdens, k_wdens

      integer ld_ij, ld_kl, lh_ik, ld_ik, lh_jl, ld_jl, lh_il, ld_il,
     $     lh_jk, ld_jk, lh_ij2, ld_ij2, lh_kl2, ld_kl2, lh_ik2, ld_ik2,
     $     lh_jl2, ld_jl2, lh_il2, ld_il2, lh_jk2, ld_jk2

      double precision crd,  ! atomic coordinates
     $                 charge, charge_nuc, tol2e, q

      dimension crd(3), g_dens(3) ! max # of density matrices

      character*16 tag
      character*4 scftype
      logical oskel, status

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

C     get information about type of calculation
      if (.not. rtdb_get(rtdb, 'scf:nopen', MT_INT, 1, nopen))
     $     nopen = 0
      if (.not. rtdb_get(rtdb, 'scf:nclosed', MT_INT, 1, nclosed))
     $     call errquit ('scf: nopen > nelec ', nopen)

      if (.not. rtdb_cget(rtdb, 'scf:scftype', 1, scftype)) then
         if (nopen .eq. 0) then
            scftype = 'RHF'
         else
            scftype = 'ROHF'
         endif
      endif
      call inp_ucase(scftype)

      if ( ga_nodeid() .eq. 0 ) then
        if (util_print('information',print_medium)) then
          write(6,1) charge, scftype
 1        format(/
     $         '  charge          = ', f6.2/
     $         '  wavefunction    =   ', a/)
          call util_flush(6)
        endif
      endif
C     # of eigenvalues and density matrices
      if ( scftype .eq. 'UHF' ) then
        nbf2 = 2 * nbf
        ndens = 2
      else if ( scftype .eq. 'RHF' ) then
        nbf2 = nbf
        ndens = 2
      else    ! ROHF
        nbf2 = nbf
        ndens = 3
      endif

C  allocate and initialize global and local memory

C     energy weighted density (NxN)
      g_wdens = ga_create_atom_blocked(geom, basis, 'weighted density')

C     forces on atoms (3xnat)
      status = ga_create ( MT_DBL, 3, nat, 'forces', 3, 1, g_force )
      call ga_zero ( g_force )

C     local replication
      lforce = nat * 3
      status = MA_alloc_get(MT_DBL, lforce, 'forces', l_force, k_force)
      if (.not.status) call errquit('could not allocate l_force',1)

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
     $     util_print('information',print_default)) then
         write(6,*) ' Using symmetry'
         call util_flush(6)
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
     $                nopen, nclosed, scftype, movecs)

C     free temporary arrays
      status = ga_destroy ( g_eigen_diag )
      if (.not.status) 
     $     call errquit('grad: could not destroy g_eigen_diag',1)
      status = MA_free_heap ( l_occ )
      status = MA_free_heap ( l_evals )

C  initialize for integral gradients
      call dfill ( lforce, 0.0D0, dbl_mb(k_force), 1 )
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
      cpu_tim(1)  = util_cpusec()
      wall_tim(1) = util_wallsec()
      call grad1( dbl_mb(k_buf), lbuf, dbl_mb(k_scr), lscratch,
     $     dbl_mb(k_dens), dbl_mb(k_wdens), dbl_mb(k_force), 
     $     g_dens, g_wdens, basis, geom, nproc, nat, max_at_bf,
     $     rtdb, oskel)
      cpu_tim(1)  = util_cpusec() - cpu_tim(1)
      wall_tim(1) = util_wallsec() - wall_tim(1)
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

      if ( scftype .ne. 'RHF' ) then ! UHF or ROHF
        if ( scftype .eq. 'ROHF' ) then
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

      cpu_tim(2)  = util_cpusec()
      wall_tim(2) = util_wallsec()
      call grad2(dbl_mb(ld_ij),dbl_mb(ld_kl),dbl_mb(ld_ik),
     $           dbl_mb(ld_jl),dbl_mb(ld_il),dbl_mb(ld_jk), 
     $           dbl_mb(ld_ij2),dbl_mb(ld_kl2),  !    only
     $           dbl_mb(ld_ik2),dbl_mb(ld_jl2),  !    for
     $           dbl_mb(ld_il2),dbl_mb(ld_jk2),  !  UHF,ROHF
     $           dbl_mb(k_force), g_dens, g_force, max_at_bf,
     $           rtdb, geom, basis, nproc, nat,
     $           dbl_mb(k_scr), lscratch, dbl_mb(k_buf), lbuf, tol2e,
     $           dbl_mb(k_sch), log_mb(k_act), oskel, scftype )
      cpu_tim(2)  = util_cpusec()  - cpu_tim(2)
      wall_tim(2) = util_wallsec() - wall_tim(2)
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

C     print forces
      if (ga_nodeid().eq.0 .and.
     $     util_print('forces',print_low)) then
         write (luout,1000)'x','y','z','x','y','z'
         do 30, i=1, nat
            status = geom_cent_get ( geom, i, tag, crd, q )
            write (luout,2000) i, tag, (crd(j),j=1,3),
     $           (dbl_mb(k_force+3*(i-1)+j),j=0,2)
 30      continue
         write(luout,*)
 1000    format (/,/,30X,'ENERGY GRADIENTS',/,/,4X,'atom',15X,
     $        'coordinates',
     $        24X,'gradient',/,6X,2(1X,(3(10X,A1))))
 2000    format (1X,I3,1X,A4,2(1X,3(1X,F10.6)))
      endif

C     print timing information
      if (ga_nodeid().eq.0 .and.
     $     util_print('timing',print_default)) then
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
      endif

C  free memory

      status = MA_free_heap ( l_act )
      if ( scftype .ne. 'RHF' ) then
        status = MA_free_heap ( lh_jk2 )
        status = MA_free_heap ( lh_il2 )
        status = MA_free_heap ( lh_jl2 )
        status = MA_free_heap ( lh_ik2 )
        if ( scftype .eq. 'ROHF' ) then
          status = MA_free_heap ( lh_ij2 )
          status = MA_free_heap ( lh_kl2 )
        endif
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
      status = MA_free_heap ( l_force )
      status = ga_destroy ( g_force )
      status = ga_destroy ( g_wdens )
      do i=1, ndens
        status = ga_destroy ( g_dens(i) )
      enddo

      return
      end
