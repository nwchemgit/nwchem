      logical function scf(rtdb)
C$Id: scf.F,v 1.32 1995-10-10 06:11:13 d3g681 Exp $
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "context.fh"
#include "inp.fh"
#include "cscfps.fh"
#include "sym.fh"
#include "util.fh"
#include "cscf.fh"
c
c     ROHF module.
c
c     Context is '...:scf'
c
c     Currently performs closed-shell and high-spin open-shell
c     restricted Hatree-Fock calculations
c
      integer rtdb              ! [input] database handle
      logical converged
      logical status
      double precision energy
      logical oprint
      character*255 name, ri_trans
c
      integer nbf_per_ir(20)
      integer natoms, Nff
      logical rohf, uhf, riscf_init
      external rohf, uhf, riscf_init
c
c     Make sure database is parallel, push context down to SCF,
c     and initialize performance statistics
c
      status = rtdb_parallel(.true.)
      if (.not. context_push('scf'))
     $     call errquit('scf: context_push failed', 0)
      call util_print_push
      call util_print_rtdb_load(rtdb, 'scf')
      oprint = util_print('information', print_low)
c
c     Extract high level info from the data-base setting defaults
c     return parameters here and also shoving stuff into commons
c
      call scf_init(rtdb)
c
c     Initialize performance statistics
c
      call scf_pstat_init(rtdb)
c
c     Print info
c 
      if (ga_nodeid().eq.0 .and. oprint) then
         call util_print_centered(6, 'NWCHEM SCF Module', 40, .true.)
         write(6,*)
         write(6,*)
         if (title .ne. ' ') then
            call util_print_centered(6, title, 40, .false.)
            write(6,*)
            write(6,*)
         endif
c
         if (.not. geom_ncent(geom, natoms)) call errquit
     $        ('scf: geom_ncent failed', 0)
c
         write(6,1) trans(1:inp_strlen(trans)), nbf
 1       format(/
     $        '  ao basis        = "',a,'"'/
     $        '  functions       = ', i5)
         if (nriscf .ne. 0) then
            status = bas_numbf (riscf_basis, Nff)
            if (.not. bas_name(riscf_basis, name, ri_trans))
     $           call errquit('scf: bas_name?', 0)
            write(6,11) ri_trans(1:inp_strlen(ri_trans)), nff
 11         format(
     $           '  ri basis        = "',a,'"'/
     $           '  ri functions    = ', i5)
         endif
         write(6,12) natoms,
     $        nclosed, nopen, nelec, nalpha, nbeta, charge,
     $        scftype,
     $        movecs_in(1:inp_strlen(movecs_in)),
     $        movecs_out(1:inp_strlen(movecs_out)), oskel, oadapt
 12      format(
     $        '  atoms           = ', i5/
     $        '  closed shells   = ', i5/
     $        '  open shells     = ', i5/
     $        '  total electrons = ', i5/
     $        '  alpha electrons = ', i5/
     $        '  beta  electrons = ', i5/
     $        '  charge          = ', f6.2/
     $        '  wavefunction    = ', a/
     $        '  input vectors   = ', a/
     $        '  output vectors  = ', a/
     $        '  use symmetry    = ', l1/
     $        '  symmetry adapt  = ', l1/)
         if (oadapt) then
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
     $        call errquit('scf: geom_print', 0)
            write(6,*)
         endif
         if (util_print('symmetry', print_debug)) then
            call sym_print_all(geom, .true., .true., .true., 
     $           .true., .true.)
            write(6,*)
         endif
         if (util_print('basis', print_high)) then
            if (.not. bas_print(basis)) 
     $           call errquit('scf: bas_print', 0)
            write(6,*)
         endif
         if (util_print('geombas', print_debug)) then
            if (.not. gbs_map_print(basis)) 
     $           call errquit('scf:gmp', 0)
            write(6,*)
         endif
         call util_flush(6)
      endif
c
      call int_normalize(basis)
c
c     initialize for RI and calculate 3c integrals if necessary
c
      if (nriscf.gt.0) then
        status = riscf_init(rtdb, geom, basis, riscf_basis, tol2e)
        if (.not. status) call errquit('scf: riscf_init failed', 0)
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
c     Initialize integral library, schwarz screening
c
      call int_init(rtdb,1, basis)
      call schwarz_init(geom, basis)
c
c     Form intial guess vectors
c
      call scf_vectors_guess(rtdb)
c
c     Here branch to the appropriate wavefunction type ... returns
c     true if happily converged, false otherwise (e.g., if restart
c     is needed or problems encountered).
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
      if (.not. context_prefix('energy', name)) call errquit
     $     ('scf: context_prefix failed', 0)
      if (.not. rtdb_put(rtdb, name, MT_DBL, 1, energy)) call errquit
     $     ('scf: failed to store energy in rtdb', 0)
      if (.not. context_prefix('converged', name)) call errquit
     $     ('scf: context_prefix failed', 0)
      if (.not. rtdb_put(rtdb, name, MT_LOG, 1, converged)) call errquit
     $     ('scf: failed to store converged in rtdb', 0)
c
c     terminate integrals 
c
      call schwarz_tidy()
      call int_terminate
c
c     SCF is done ... tidy up memory etc.
c
      call scf_tidy
      if (nriscf.ne.0) then
        call riscf_tidy
      endif
c
c     Nothing else for now ... return convergence status
c
      call util_print_pop
      if (.not. context_pop('scf')) call errquit('context_pop failed',0)
      scf = converged
c
      call scf_pstat_print
c
      end
