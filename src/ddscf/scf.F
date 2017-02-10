      logical function scf(rtdb)
C$Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "cscfps.fh"
#include "sym.fh"
#include "util.fh"
#include "cscf.fh"
#include "case.fh"
#include "frozemb.fh"
c
c     ROHF module.
c
c     Currently performs closed-shell and high-spin open-shell
c     restricted Hartree-Fock calculations
c
      integer rtdb              ! [input] database handle
      logical converged
      logical status, olocal
      double precision energy
      logical oprint
      character*255 name, ri_trans
c
      integer nbf_per_ir(20),i
      integer natoms, Nff, molo, mohi
      integer g_dens(2) ! dummy global array handles for initial guess
      integer maxireps
      parameter (maxireps=20)
      integer nop,nir,class_dim(maxireps),ipnt(maxireps)
      character*8 zir(maxireps),zclass(maxireps)
      character*8 zname
      double precision chars(maxireps*maxireps)
      integer ioff(maxireps),icount(maxireps)
      integer irs,reord
      integer k_reord,l_reord,k_reord2,l_reord2
c     integer reord(10),reord2(10),irs(10)
      integer g_hcore,noper,g_coul,g_exch
      logical rohf, uhf, riscf_init, movecs_converged
      external rohf, uhf, riscf_init, movecs_converged
      logical int_normalize
      external int_normalize
      logical cphf_poliz, do_poliz
      external cphf_poliz
      integer necp,itype
      integer blen
c
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      integer ga_create_JKblocked
      external ga_create_JKblocked
c
      data blen/25/
c
c     Make sure database is parallel, 
c     and initialize performance statistics
c
      status = rtdb_parallel(.true.)
      call util_print_push
      call util_print_rtdb_load(rtdb, 'scf')
      call ecce_print_module_entry('scf')
      oprint = util_print('information', print_low)
c
c     Extract high level info from the data-base setting defaults
c     return parameters here and also shoving stuff into commons
c
c !!! BGJ
c !!! Can't place the following in scf_init since scf_init is called
c !!! by cphf for either HF or DFT
c     Store SCF hamiltonian type as HF for use in BGJ routines
      if (.not. rtdb_put(rtdb, 'bgj:scf_type', MT_INT, 1, 1))
     $     call errquit('scf: put of bgj:scf_type failed',0, RTDB_ERR)
c !!! BGJ
      call scf_init(rtdb)
c
c     Initialize performance statistics
c
      call scf_pstat_init(rtdb)
c
c     Print info
c 
      if (ga_nodeid().eq.0 .and. oprint) then
         call util_print_centered(6, 'NWChem SCF Module', 40, .true.)
         write(6,*)
         write(6,*)
         if (title .ne. ' ') then
            call util_print_centered(6, title, 40, .false.)
            write(6,*)
            write(6,*)
         endif
c
         if(.not. geom_systype_get(geom,itype))
     &        itype = 0
         if(itype.ne.0) call errquit(
     A     'scf: code working only for molecules',itype, INPUT_ERR)

         if (.not. geom_ncent(geom, natoms)) call errquit
     $        ('scf: geom_ncent failed', 0, GEOM_ERR)
c
         write(6,1) trans(1:inp_strlen(trans)), nbf
 1       format(/
     $        '  ao basis        = "',a,'"'/
     $        '  functions       = ', i5)
         if (nriscf .ne. 0) then
            status = bas_numbf (riscf_basis, Nff)
            if (.not. bas_name(riscf_basis, name, ri_trans))
     $           call errquit('scf: bas_name?', 0,
     &       BASIS_ERR)
            write(6,11) ri_trans(1:inp_strlen(ri_trans)), nff
 11         format(
     $           '  ri basis        = "',a,'"'/
     $           '  ri functions    = ', i5)
         endif
         write(6,12) natoms
 12      format('  atoms           = ', i5)
         if (scftype .eq. 'UHF') then
            write(6,121) nalpha, nbeta
 121        format('  alpha electrons = ', i5/
     $           '  beta  electrons = ', i5)
         else
            write(6,122) nclosed, nopen
 122        format('  closed shells   = ', i5/
     $           '  open shells     = ', i5)
         endif
         write(6,123)
     $        charge, scftype,
     $        movecs_in(1:inp_strlen(movecs_in)),
     $        movecs_out(1:inp_strlen(movecs_out)), oskel, oadapt
 123     format(
     $        '  charge          = ', f6.2/
     $        '  wavefunction    = ', a/
     $        '  input vectors   = ', a/
     $        '  output vectors  = ', a/
     $        '  use symmetry    = ', l1/
     $        '  symmetry adapt  = ', l1)
         if (olock) write(6,124) olock
 124     format('  lock orbitals   = ', l1)
         if (util_print('basis summary', print_default)) then
            if (.not.bas_summary_print(basis))
     &           call errquit('scf: basis summary print failed',911,
     &       BASIS_ERR)
         else
            write(6,*)
         endif
         if (oadapt .and. util_print('char table', print_high)) then
            call sym_print_char_table(geom)
         endif
         if (oadapt .and. util_print('basis sym', print_default)) then
            call sym_bas_irreps(basis, .true., nbf_per_ir)
         endif
         call util_flush(6)
      endif
c
c     Print out other info as requested by user
c
      if (ga_nodeid() .eq. 0) then
         if (util_print('geometry',print_high)) then
            if (.not. geom_print(geom))
     $        call errquit('scf: geom_print', 0, GEOM_ERR)
            write(6,*)
         endif
         if (util_print('symmetry', print_debug)) then
            call sym_print_all(geom, .true., .true., .true., 
     $           .true., .true.)
            write(6,*)
         endif
         if (util_print('basis', print_high)) then
            if (.not. bas_print(basis)) 
     $           call errquit('scf: bas_print', 0, BASIS_ERR)
            write(6,*)
         endif
         if (util_print('basis labels',print_high) .and.
     $        ga_nodeid().eq.0) then
            call bas_print_labels(basis)
         endif
         if (util_print('geombas', print_debug)) then
            if (.not. gbs_map_print(basis)) 
     $           call errquit('scf:gmp', 0, BASIS_ERR)
            write(6,*)
         endif
         call util_flush(6)
      endif
c
c     Get and check the Coulomb Attenuation Method (CAM) parameters
      if (.not.rtdb_get(rtdb, 'scf:cam_exch', mt_log, 1,
     &   cam_exch)) cam_exch=.false.
      if (.not.rtdb_get(rtdb, 'scf:cam_omega', mt_dbl, 1,
     &   cam_omega)) cam_omega=0.d0
      if (.not.rtdb_get(rtdb, 'scf:cam_alpha', mt_dbl, 1,
     &   cam_alpha)) cam_alpha=0.d0
      if (.not.rtdb_get(rtdb, 'scf:cam_beta', mt_dbl, 1,
     &   cam_beta)) cam_beta=0.d0
c
      if (.not.int_normalize(rtdb,basis))
     $           call errquit('scf:int_normalize failed', 0, INT_ERR)
c
c     If input and output movecs are same then might be converged MOs
c
      if (movecs_in .eq. movecs_out) then
         converged = movecs_converged(rtdb, basis, 'scf', movecs_in)
         if (converged) then
            if (.not. rtdb_get(rtdb, 'scf:energy', mt_dbl, 1, energy))
     $           call errquit('scf: failed getting converged energy',0,
     &       RTDB_ERR)
            if (ga_nodeid().eq.0 .and. oprint) then
               write(6,1101) energy
 1101          format(/'  The SCF is already converged '//,
     $              '         Total SCF energy =', f20.12/)
               call util_flush(6)
            endif
            goto 3131
         endif
      endif
c
c     initialize for RI and calculate 3c integrals if necessary
c
      if (nriscf.gt.0) then
        status = riscf_init(rtdb, geom, basis, riscf_basis, tol2e)
        if (.not. status) call errquit('scf: riscf_init failed', 0,
     &       UNKNOWN_ERR)
        if (ga_nodeid().eq.0 .and. oprint) then
          if (nriscf.eq.1) then
            call util_print_centered(6, 'RI hessian', 40, .true.)
          else if (nriscf.eq.2) then
            call util_print_centered(6, 'RISCF', 40, .true.)
          else if (nriscf.eq.3) then
            call util_print_centered(6, 'preconverge with RISCF', 40, 
     $           .true.)
          endif
          write(6,*)
        endif
      endif
c
c     Initialize integral library, schwarz screening, cached 1-e integs
c
      call int_init(rtdb,1, basis)
c
      call print_integrals(basis,oskel) ! print integrals, if requested
c 
      call schwarz_init(geom, basis)
      if ((nbf*nbf*3*8)/ga_nnodes() .lt. 2000000) ! Cache 1e if < 2MB/node
     $     call int_1e_cache_ga(basis, oskel)
c
c     Form intial guess vectors
c
      if (ga_nodeid().eq.0 .and. oprint) write(6,711) util_wallsec()
 711  format(/' Forming initial guess at ',f9.1,'s'/)
c      call scf_vectors_guess(rtdb)
c
c     prev call replaced with new call
c
      call scf_vectors_guess(rtdb, tol2e, geom, basis, trans, 
     &                       movecs_in, movecs_out, movecs_guess,
     &                       scftype, nclosed, nopen, nbf, nmo,
     &                       nalpha, nbeta,  k_eval, k_occ, k_irs,
     &                       g_movecs, g_dens, 'vectors', 
     &                       'scf', title, oskel, oadapt, 
     &                       oinitialized) 
c
c     Here branch to the appropriate wavefunction type ... returns
c     true if happily converged, false otherwise (e.g., if restart
c     is needed or problems encountered).
c
      if (ga_nodeid().eq.0 .and. oprint) write(6,712) util_wallsec()
 712  format(/' Starting SCF solution at ',f9.1,'s'/)
c
c     frozen embedding 
      frozemb = .false.
      if (.not.rtdb_get(rtdb, 'scf:frozemb', mt_log, 1, frozemb))
     &   frozemb=.false.
      if (frozemb) then
       g_frozemb = ga_create_atom_blocked(geom, basis, 'rohf: frozemb')
       call ga_zero(g_frozemb)
       call scf_frozemb(rtdb, geom, basis, nbf, g_frozemb)
      end if   ! frozemb
c
      if (scftype .eq. 'UHF') then
         converged = uhf(rtdb, energy)
      else
         converged = rohf(rtdb, energy)
      endif
c
c     reconverge with exact energy/gradient if requested
c
      if (nriscf.eq.3) then
        if (ga_nodeid().eq.0 .and. oprint) then
          call util_print_centered(6,
     $         'switching to exact energy/gradient', 40, .true.)
          write(6,*)
        endif
        nriscf = 1
        if (scftype .eq. 'UHF') then
          converged = uhf(rtdb, energy)
        else
          converged = rohf(rtdb, energy)
        endif
      endif
c
c     Store energy and convergence status ... not sure if better to
c     do this here or have each wavefunction type do it
c
      if (.not. rtdb_put(rtdb, 'scf:energy', MT_DBL, 1, energy))
     $     call errquit
     $     ('scf: failed to store energy in rtdb', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'scf:converged', MT_LOG, 1, converged)) 
     $     call errquit
     $     ('scf: failed to store converged in rtdb', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'scf:alpha irreps', MT_INT, nbf,
     &   int_mb(k_irs)))
     &   call errquit('scf: failed to MO irreps in rtdb', 0, RTDB_ERR)
      if (scftype.eq.'UHF') then
      if (.not. rtdb_put(rtdb, 'scf:beta irreps', MT_INT, nbf,
     &   int_mb(k_irs+nbf)))
     &   call errquit('scf: failed to MO irreps in rtdb', 0, RTDB_ERR)
      endif
c
c     Must store the final MOs AFTER writing the convergence flag
c     since the time it is put into the RTDB is used to verify
c     that the two correspond.
c
      if (.not. rtdb_get(rtdb, 'scf:localize', mt_log, 1, olocal))
     $     olocal = .false.
      if (olocal) CALL TESTLOCALIZE(rtdb)
c
      call scf_movecs_write(rtdb)
c
c     With some better data compression in the ECCE output
c     we can now output all occupied plus the top 20 virtuals.
c     Also, ECCE may currently be confused if molo != 1
c
      molo = 1
      mohi = nmo ! min(nmo,nalpha+20)
c
      if (scftype .eq. 'UHF') then
         call movecs_ecce(nbf, nmo, molo, mohi, dbl_mb(k_eval),
     $        dbl_mb(k_occ), int_mb(k_irs), 
     $        g_movecs(1), scftype, 'alpha')
         call movecs_ecce(nbf, nmo, molo, mohi, dbl_mb(k_eval+nbf),
     $        dbl_mb(k_occ+nbf), int_mb(k_irs+nmo),
     $        g_movecs(2), scftype, 'beta')
      else
         call movecs_ecce(nbf, nmo, molo, mohi, dbl_mb(k_eval),
     $        dbl_mb(k_occ), int_mb(k_irs),
     $        g_movecs(1), scftype, ' ')
      endif
c
c     terminate integrals 
c
      call int_1e_uncache_ga()
      call schwarz_tidy()
      call int_terminate
c 
c     clean up frozen embedding
      if (frozemb) then
        if (.not. ga_destroy(g_frozemb)) call errquit
     $     ('rohf_fock: ga_destroy failed',0, GA_ERR)
      end if
c
c     SCF is done ... tidy up memory etc.
c
      if (nriscf.ne.0) then
        call riscf_tidy ( geom, basis, riscf_basis )
      endif
 3131 call scf_tidy(rtdb)
c
      call scf_pstat_print
c
c     Nothing else for now ... return convergence status
c
      call util_print_pop
c
      if (converged) then
         call ecce_print_module_exit('scf', 'ok')
      else
         call ecce_print_module_exit('scf', 'failed')
      endif
c
      scf = converged
c
c !!! BGJ
      if (.not. rtdb_get(rtdb, 'bgj:poliz', mt_log,
     &     1, do_poliz)) then
         do_poliz = .false.
      endif
      if (do_poliz) then
         write(*,*)'*** scf: calling cphf_poliz'
         if (.not. cphf_poliz(rtdb)) ! Never executed.
     $        call errquit(' cphf_poliz: failed from scf !',0, RTDB_ERR)
      endif
c !!! BGJ
c
      end
