      subroutine dft_cg_init(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "sym.fh"
#include "util.fh"
#include "cscf.fh"
#include "cosmo.fh"
#include "stdio.fh"
c !!! BGJ
#include "bgj.fh"
c !!! BGJ
c     
      integer rtdb              ! database handle
c     
      double precision nuclear_charge
      character*255 name
      integer len_occ
      external cscfdata ! For T3D linker
c
      logical osome
c
      logical hf_job
      character*30 tag
      character*255 theory
      integer mult
      integer cd_basis
c
c     hf_job = .not. bgj_have_xc()
c     if (bgj_print() .gt. 0)
c    &     write(*,*)'*** scf_get_info: hf_job ',hf_job
c
c     if (.not. rtdb_cget(rtdb, 'title', 1, title))
c    $     title = ' '
c
c     load geometry and symmetry info
c
c     if (.not. geom_create(geom, 'geometry'))
c    $     call errquit('dft_cg_init: geom_create?', 0, GEOM_ERR)
c     if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
c    $     call errquit('dft_cg_init: no geometry ', 0, RTDB_ERR)
      if (rtdb_get(rtdb, 'dft:level shift info', MT_DBL, 6,shifts)) then
        ouser_changed_conv = .true.
      else
        call dfill(6, -1.0d0, shifts, 1)
      endif
c
      if (rtdb_get(rtdb, 'dft:full hessian switch', MT_DBL, 1,
     $     nr_gswitch)) then
        ouser_changed_conv = .true.
      else
        nr_gswitch = 0.1d0
      endif
c
      if (.not.rtdb_get(rtdb, 'dft:skeleton',MT_LOG, 1, oskel)) then
         oskel = sym_number_ops(geom) .gt. 0
      endif
c
      if (.not.rtdb_get(rtdb, 'dft:adapt',MT_LOG, 1, oadapt)) then
         oadapt = sym_number_ops(geom) .gt. 0
      endif
c
      if (.not.rtdb_get(rtdb, 'dft:lock',MT_LOG, 1, olock)) then
         olock = .false.
      endif
c
      if (.not. rtdb_get(rtdb, 'dft:tol2e', MT_DBL, 1, tol2e))
     &     tol2e = min(1.0d-7,gnorm_tol*1d-2)
c
      if (rtdb_get(rtdb, 'dft:ri', MT_INT, 1, nriscf)) then
        if (.not. bas_create(riscf_basis, 'riscf basis'))
     &       call errquit('dft_cg_init: bas_create?', 0, BASIS_ERR)
        if (.not. bas_rtdb_load(rtdb, geom, riscf_basis, 'riscf basis'))
     &       call errquit('dft_cg_init: no riscf basis set', 0,
     &                    RTDB_ERR)
      else
        nriscf = 0
        riscf_basis = 0
      endif
c
c     check for CD fitting
c
      if (.not. bas_create(cd_basis, 'cd basis'))
     &     call errquit('dft_cg_init: bas_create?', 0, BASIS_ERR)
       if(bas_rtdb_load(rtdb, geom, cd_basis, 'cd basis')) then
          if(ga_nodeid().eq.0) write(6,*) 
     C         'WARNING: CD fitting not available for cgmin yet'
        if (.not. bas_destroy(cd_basis)) call errquit
     $         ('dftcginfo: basis destroy failed',0, BASIS_ERR)
      else
        if (.not. bas_destroy(cd_basis)) call errquit
     $         ('dftcginfo: basis destroy failed',33, BASIS_ERR)

       endif
c     
c     load the basis set and get info about it
c
c     if (.not. bas_create(basis, 'ao basis'))
c    $     call errquit('dft_cg_init: bas_create?', 0, BASIS_ERR)
c     if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis')) 
c    $        call errquit('dft_cg_init: no ao basis set', 0, RTDB_ERR)
c
c     For debug ... call int_init and do the 2-e
c
      if (util_print('texas init debug',print_never)) then
        call int_init(rtdb, 1, basis)
        write(6,*) ' DONE INIT'
*        call schwarz_init(geom, basis)
*        call schwarz_tidy()
        call int_terminate
      endif
c
      if (.not. bas_name(basis, name, trans))
     $     call errquit('dft_cg_init: bas_name?', 0, BASIS_ERR)
c
c     if (.not. bas_numbf(basis, nbf)) call errquit
c    $     ('dft_cg_init: basis info',0, BASIS_ERR)
c
c     Is RI approximation to be used? If so get fitting basis set.
c
      if (rtdb_get(rtdb, 'dft:ri', MT_INT, 1, nriscf)) then
        if (.not. bas_create(riscf_basis, 'riscf basis'))
     $     call errquit('dft_cg_init: bas_create?', 0, BASIS_ERR)
        if (.not. bas_rtdb_load(rtdb, geom, riscf_basis, 'riscf basis'))
     $     call errquit('dft_cg_init: no riscf basis set', 0, RTDB_ERR)
      else
        nriscf = 0
        riscf_basis = 0
      endif
c
c     Figure input/output MO vectors
c     
c     if (hf_job) then
c        tag = 'scf:input vectors'
c     else
c        tag = 'dft:input vectors'
c     endif
c     if (.not. rtdb_cget(rtdb, tag, 1, movecs_in))
c    $     movecs_in = 'atomic'
c     if (hf_job) then
c        tag = 'scf:output vectors'
c     else
c        tag = 'dft:output vectors'
c     endif
c     if (.not. rtdb_cget(rtdb, tag, 1, movecs_out))
c    $     movecs_out = ' '
c     if (movecs_out.eq.' ') then
c        if (movecs_in.eq.'atomic' .or. movecs_in.eq.'hcore' .or. 
c    $        movecs_in.eq.'project' .or. movecs_in.eq.'fragment'
c    $          .or.movecs_in.eq.'rotate') then
c           call util_file_name('movecs', .false.,.false.,movecs_out)
c        else
c           movecs_out = movecs_in
c        endif
c     endif
c
c     Resolve names of MO files to full paths defaulting to the
c     permanent directory
c
c     if (movecs_in.eq.'atomic' .or. movecs_in.eq.'hcore' .or. 
c    $        movecs_in.eq.'project' .or. movecs_in.eq.'fragment'
c    $          .or.movecs_in.eq.'rotate') then
c        continue
c     else
c        call util_file_name_resolve(movecs_in, .false.)
c     endif
c     call util_file_name_resolve(movecs_out, .false.)
c     
c     Figure out the number of electrons from the required total
c     charge and the sum of nuclear charges
c     
c     if (.not. rtdb_get(rtdb, 'charge', MT_DBL, 1, charge))
c    $     charge = 0.0d0
c     if (.not. geom_nuc_charge(geom, nuclear_charge))
c    $     call errquit('scf: geom_nuc_charge failed', 0, GEOM_ERR)
c     nelec = nint(nuclear_charge - charge)
c     if (nelec .le. 0) call errquit
c    $     ('scf: negative no. of electrons ?', nelec, INPUT_ERR)
c     if (abs(nuclear_charge - charge - dble(nelec)) .gt. 1d-8)
c    $     call errquit('scf: non-integral no. of electrons ?', 0,
c    &       INPUT_ERR)
c     
c     Determine no. of open and closed shells ... default is to run closed
c     shell unless told otherwise
c     
c     if (bgj_print() .gt. 0)
c    &   write(*,*)'*** !!! more stuff to check in scf_get_info.F !!!'
c     if(.not.rtdb_cget(rtdb,'task:theory',1,theory))
c    +     call errquit('task: no task input for theory?',0, INPUT_ERR)
c     if (.not. rtdb_get(rtdb, 'scf:nopen', MT_INT, 1, nopen))
c     $     nopen = 0
c     if(theory .eq. 'dft')then
c        if (.not. rtdb_get(rtdb, 'dft:mult', MT_INT, 1,mult))
c    *        mult = 1
c        nopen = mult - 1
c     endif
c     if (.not. rtdb_get(rtdb, 'scf:nopen', MT_INT, 1, nopen))
c    $     nopen = 0
      if (nopen .gt. nelec) call errquit
     $     ('scf: nopen > nelec ', nopen, INPUT_ERR)
      if (mod(nelec-nopen,2) .ne. 0) call errquit
     $     ('dft: no. of closed-shell electrons is not even!',nopen,
     &       INPUT_ERR)
c     nclosed = (nelec-nopen) / 2
c     
c     if (.not. rtdb_cget(rtdb, 'scf:scftype', 1, scftype)) then
c        if (nopen .eq. 0) then
c           scftype = 'RHF'
c        else
c           scftype = 'ROHF'
c        endif
c     endif
c
      call inp_ucase(scftype)
c
c     Take care of holes in the input routines
c
      if (scftype.eq.'RHF' .and. nopen.gt.0) then
         scftype = 'ROHF'
      else if (scftype.eq.'ROHF' .and. nopen.eq.0) then
         scftype = 'RHF'
      endif
c
      if ( scftype.ne.'ROHF' .and. scftype.ne.'RHF' .and.
     $     scftype.ne.'UHF' ) call errquit
     $     ('scf: only ROHF, RHF, and UHF currently supported', 0,
     &       INPUT_ERR)
c
c     Dump lagrangian?  Yes by default now since if the SCF 
c     has converged for an energy it will not be rerun for the gradient.
c
      if (.not.rtdb_get(rtdb, 'dft:lagrangian',MT_LOG, 1, olagr))
     $   olagr = scftype.eq.'ROHF'
c
      nalpha = nclosed + nopen
      nbeta  = nclosed
c
c   DIIS toggle
c
c     if (.not.rtdb_get(rtdb, 'scf:diis',MT_LOG, 1, odiis))
c    $   odiis = .false.
c
      call ga_sync()
c
c     For now set NMO = NBF, however this may change later when the 
c     linear dependency analysis is done just before the starting guess
c
      nmo = nbf
c
c     Store the derived info in the database for other wavefunction
c     modules and/or restart to access
c
      if (.not. rtdb_cput(rtdb, 'dft:scftype', 1, scftype)) 
     $   call errquit('dft_cg_init: put of scftyp failed',0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'dft:nopen', MT_INT, 1, nopen))
     $   call errquit('dft_cg_init: put of nopen failed',0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'dft:nclosed', MT_INT, 1, nclosed))
     $   call errquit('dft_cg_init: put of nclosed failed',0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'dft:nelec', MT_INT, 1, nelec))
     $   call errquit('dft_cg_init: put of nelec failed',0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'dft:nmo', MT_INT, 1, nmo))
     $   call errquit('dft_cg_init: put of nmo failed',0, RTDB_ERR)
      if (scftype .eq. 'UHF') then
        if (.not. rtdb_put(rtdb, 'dft:nalpha', MT_INT, 1, nalpha))
     $     call errquit('dft_cg_init: put of nalpha failed',0, RTDB_ERR)
        if (.not. rtdb_put(rtdb, 'dft:nbeta', MT_INT, 1, nbeta))
     $     call errquit('dft_cg_init: put of nbeta failed',0, RTDB_ERR)
      endif
c     
c     Allocate persistent local and global arrays ... these may
c     be reallocated later when the dependency analysis is done
c
c     if (scftype .eq. 'UHF') then
c        if (.not. ga_create(MT_DBL, nbf, nmo, 'dft_cg_init: alpha MOs',
c    $        32, 32, g_movecs)) call errquit('dft_cg_init: alpha MOs', 0,
c    &       GA_ERR)
c        if (.not. ga_create(MT_DBL, nbf, nmo, 'dft_cg_init: beta MOs',
c    $        32, 32, g_movecs(2))) call errquit('dft_cg_init: beta MOs',0,
c    &       GA_ERR)
c     else
c        if (.not. ga_create(MT_DBL, nbf, nmo, 'dft_cg_init: MOs',
c    $        32, 32, g_movecs)) call errquit('dft_cg_init: MOs', 0,
c    &       GA_ERR)
c     endif
c     
      len_occ = nmo
      if (scftype .eq. 'UHF') len_occ = nbf * 2
      if (.not. ma_push_get(mt_dbl, len_occ, 'dft_cg_init: mo evals',
     $     l_eval, k_eval)) call errquit
     $     ('dft_cg_init: insufficient memory?', len_occ, MA_ERR)
c     
      if (.not. ma_push_get(mt_dbl, len_occ, 'dft_cg_init: mo occ',
     $     l_occ, k_occ)) call errquit
     $     ('dft_cg_init: insufficient memory?', len_occ, MA_ERR)
c
      if (.not. ma_push_get(mt_int, len_occ, 'dft_cg_init: mo irs',
     $     l_irs, k_irs)) call errquit
     $     ('dft_cg_init: insufficient memory?', len_occ, MA_ERR)
c
      call ifill(len_occ, 1, int_mb(k_irs), 1) ! In case not adapting
c
c     Fill in the SCF convergence info
c     
      call dft_get_conv_info(rtdb)
      call scf_get_fock_param(rtdb, tol2e)
c
c     ----- cosmo initialization ----
c
c      cosmo_last = .false.
c      if ( rtdb_get(rtdb,'slv:cosmo',mt_log,1,cosmo_on)) then
c         if(cosmo_on) then
c            osome = util_print('information', print_low)
c
c            call cosmo_initialize(rtdb,geom,basis,osome)
c
c           Turn cosmo on, we want to run the calculation
c           Start with gas_phase run unless told otherwise
c
c            cosmo_last = .true.
c            cosmo_on = .true.
c            if(.not.rtdb_get(rtdb,'cosmo_phase',mt_int,1,cosmo_phase))
c     >         cosmo_phase = 1
c         endif
c      endif
      oinitialized = .true.
c     
      end
      subroutine dft_cg_tidy(rtdb)
      implicit none
#include "errquit.fh"
#include "cscf.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "bas.fh"
#include "cfock.fh"
#include "rtdb.fh"
#include "cosmo.fh"
#include "cgridfile.fh"
      integer rtdb
c
      logical status
c
      if (oinitialized) then
c        if (.not. geom_destroy(geom)) call errquit
c    $        ('scf_tidy: geom destroy failed', 0, GEOM_ERR)
c        if (.not. bas_destroy(basis)) call errquit
c    $        ('scf_tidy: basis destroy failed',0, BASIS_ERR)
         status = ma_pop_stack(l_irs)
         status = ma_pop_stack(l_occ) .and. status
         status = ma_pop_stack(l_eval) .and. status
         if (.not. status) call errquit
     $        ('scf_tidy: failed to free irs/occupation/evals',0,
     &       MA_ERR)
c        if (.not. ga_destroy(g_movecs)) call errquit
c    $        ('scf_tidy: failed to free movecs',0, GA_ERR)
c        if (scftype .eq. 'UHF') then
c           if (.not. ga_destroy(g_movecs(2))) call errquit
c    $           ('scf_tidy: failed to free beta movecs',0, GA_ERR)
c        endif
         oinitialized = .false.
      endif
c
      call fock_2e_tidy(rtdb)
c
c     The quasi Newton-Raphson solver uses FOCK_XC to calculate the DFT
c     contributions to the Kohn-Sham matrix. FOCK_XC will open the grid
c     file at the beginning of every Fock build and close it again when
c     done. Therefore at this point the grid file is always closed and
c     we need to set GRID_WRITTEN to .FALSE. so that the code will
c     always reopen the grid file before trying anything else.
c
      grid_written = .false.
c
c     ----- cosmo cleanup and reset -----
c
c      if ( rtdb_get(rtdb,'slv:cosmo',mt_log,1,cosmo_on)) then
c         if(cosmo_on) then
c           call cosmo_tidy(rtdb)
c           cosmo_on = .false.
c           cosmo_phase = 1
c         endif
c      endif
c
c     --- ri-scf cleanup ---
c
      if (nriscf.gt.0) then
        if (.not. bas_destroy(riscf_basis))
     &     call errquit('dft_cg_tidy: bas_destroy?',0, BASIS_ERR)
        nriscf = 0
        riscf_basis = 0
      endif
c
      end
c $Id$
