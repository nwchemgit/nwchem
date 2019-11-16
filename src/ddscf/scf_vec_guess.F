      subroutine scf_vectors_guess(rtdb, tol2e, geom, basis, trans, 
     &                             movecs_in, movecs_out, movecs_guess,
     &                             scftype, nclosed, nopen, nbf, nmo,
     &                             nalpha, nbeta,  k_eval, k_occ, k_irs,
     &                             g_movecs, g_dens, vecs_or_dens, 
     &                             calling_module, title, oskel, oadapt, 
     &                             oinitialized) 
c
c     $Id$
c
      implicit none
#include "errquit.fh"
c
c     arguments:
c
c     tol2e         = 2-e integral selection threshold for energy/gradient
c     geom          = Geometry handle
c     basis         = Basis handle
c     trans         = Translated basis set name
c     movecs_in     = Name of input movecs file
c     movecs_out    = Name of output movecs file
c     movecs_guess  = Name of source of guess orbtials (=movecs_in initially)
c     scftype       = type of SCF (UHF/ROHF)
c     nclosed       = no. of closed shells
c     nopen         = no. of open shells
c     nbf           = no. of basis functions
c     nmo           = no. of linearly indep. MOs (= nbf for now)
c     nalpha        = no. of alpha spin electrons for UHF only
c     nbeta         = no. of beta spin electrons for UHF only
c     k_eval        = scf eigen values (ma index)
c     k_occ         = scf occupation (ma index)
c     k_irs         = MO symmetries (2*nmo for UHF) (ma index)
c     g_movecs      = eigen-vectors (global ga handle) (2x array for UHF)
c     g_dens        = AO densisty (global ga handle) (2x array for UHF)
c     vecs_or_dens  = "vectors" or "density" requested
c     calling_module= calculation type; SCF and DFT currently recognized
c     title         = Job Title
c     oskel         = if true use symmetry in fock building
c     oadapt        = if true force symmetry on the MOs
c     oinitialized  = true if initialized
c
      integer geom, basis, nopen, nclosed, nbf, nmo
      integer nalpha, nbeta
      double precision tol2e,ssmear
      logical oskel, oadapt, oinitialized
      integer g_movecs(2), g_dens(2)
      character*255 movecs_in, movecs_out, title, trans
      character*8 movecs_guess
      character*7 vecs_or_dens
      character*4 scftype
      character*3 calling_module
      integer k_eval, k_occ, k_irs
c
#include "global.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "bas.fh"
#include "rtdb.fh"
c#include "cscf.fh"
#include "cfock.fh"
c
c     Form guess vectors on disk for SCF wavefunctions
c
      integer rtdb
c
c     local variables
c
      logical oprint
      integer g_tmp, g_tmp2
      integer i, nprint
      character*255 title_vecs, bas_vecs
      integer nbf_vecs, nsets, nmo_vecs(2)
      logical status, fon
      double precision tolguess, zero, one, two
      character*20 file_scftype
      character*255 string
      integer type, nelem, swap, handle
c
c     externals
c
      integer ga_create_atom_blocked
      logical movecs_read_header, movecs_write, movecs_read
      parameter (zero = 0.d0, one = 1.d0, two = 2.d0)
      external movecs_read_header, movecs_write, movecs_read
      external ga_create_atom_blocked
      logical quickguess
c
      if (.not. oinitialized) call errquit
     &     ('scf_v_g: not initialized?', 0, UNKNOWN_ERR)
c
      if(.not.rtdb_get(rtdb,'tolguess',mt_dbl,1,tolguess)) then
        if (owritefile) then
           tolguess = tol2e
        else
          tolguess = 1d-7        ! 1e-5 broke too many large basis SCFs
        endif
      endif
c
      oprint = util_print('mo guess', print_default) .and.
     &     ga_nodeid().eq.0
*      if (oprint)write(6,*)' Calling module = ',calling_module
c
c     check if smearing is on
c
      if (.not. rtdb_get(rtdb, 'dft:smear_sigma', mt_dbl, 1,
     &   ssmear)) then
         ssmear=0.d0
      endif
c
c     check if fon is on
c
      if (.not. rtdb_get(rtdb, 'dft:fon', mt_log, 1, fon)) then
         fon=.false.
      endif
c
c     Analyze linear dependence ... put eigenvectors of S into g_movecs
c     and if necessary reallocate the mo vectors to reflect the number
c     of independent MO vectors.
c
      call scf_lindep(rtdb, geom, basis, g_movecs, nmo)
      if (nmo .ne. nbf) then
         if (.not. ga_create(MT_DBL, nbf, nmo, 'scf_v_g: tmp MOs',
     &        32, 32, g_tmp)) call errquit('scf_v_g: tmp MOs', 0,
     &       GA_ERR)
         call ga_copy_patch('n', 
     &        g_movecs, 1, nbf, 1, nmo, g_tmp, 1, nbf, 1, nmo)
         if (.not. ga_destroy(g_movecs)) 
     &        call errquit('scf_v_g: ga_destroy?',0, GA_ERR)
         if (scftype .eq. 'UHF') then
            if (.not. ga_destroy(g_movecs(2))) 
     &           call errquit('scf_v_g: ga_destroy(2)?',0, GA_ERR)
            if (.not. ga_create(MT_DBL, nbf, nmo, 'scf_v_g: alpha MOs',
     &           32, 32, g_movecs)) call errquit('scf_v_g: alpha MOs',0,
     &       GA_ERR)
            if (.not. ga_create(MT_DBL, nbf, nmo, 'scf_v_g beta MOs',
     &           32, 32, g_movecs(2)))call errquit('scf_v_g:beta MOs',0,
     &       GA_ERR)
         else
            if (.not. ga_create(MT_DBL, nbf, nmo, 'scf_g_g: MOs',
     &           32,32, g_movecs)) call errquit('scf_v_g: MOs', 0,
     &       GA_ERR)
         endif
         call ga_copy(g_tmp, g_movecs)
         if (.not. ga_destroy(g_tmp)) 
     &        call errquit('scf_v_g: ga_destroy?',0, GA_ERR)
      endif
c
c     Need to force higher accuracy if linear dependence is
c     detected.  Also turn off density screening.  Have to do 
c     this also in the datebase so that subsequent invocations
c     of scf_get_fock_param() do the right thing no matter where
c     they are called from (e.g., scf, mp2, cphf, ...).
c
      if (nmo .ne. nbf) then
         tol2e = min(tol2e, 1d-10)
         dentolmax = 1d-9
         odensityscreen = .false.
         if (.not. rtdb_put(rtdb,'scf:tol2e',mt_dbl,1,tol2e))
     $        call errquit('scf_vec_guess: rtdb failed?',0, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'fock:dentolmax',mt_dbl,1,dentolmax))
     $        call errquit('scf_vec_guess: rtdb failed?',0, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'fock:densityscreen',mt_log,1,
     $        odensityscreen))
     $        call errquit('scf_vec_guess: rtdb failed?',0, RTDB_ERR)
      endif
c
 3212 continue
      if (movecs_in .eq. 'hcore' .or. movecs_in .eq. 'atomic') then
c
c     hcore or atomic guess
c
         if (.not.inp_compare(.false., 'dft', calling_module))then
            g_dens(1) = 
     &         ga_create_atom_blocked(geom, basis, 'Guess density')
         endif
         call ga_zero(g_dens(1))
         if (movecs_in .eq. 'atomic') then
            call guess_dens(rtdb, geom, basis, g_dens(1))
         endif
         if(.not.rtdb_get(rtdb,'quickguess',mt_log,1,quickguess)) 
     .        quickguess=.false.
         if (quickguess) return
c
c     Since the caching of integrals may be turned on automagically
c     we always need to use the full integral precision in the starting 
c     guess
c
c     Modified by RJH so that unless writing of integrals has been
c     forced on, it will do a quick and dirty initial guess and
c     not store any integrals at all.
c
         if (.not. owritefile) then
            call fock_force_direct(rtdb)
         endif
*         write(6,*) ' tolguess', tolguess
         call rhf_dens_to_mo(rtdb, geom, basis, 2*nclosed+nopen,
     &        nbf, nmo, tolguess,
     &        g_dens(1), g_movecs, oskel, dbl_mb(k_eval))
c
         if (.not. oreadfile) call scf_get_fock_param(rtdb,tol2e)
c
         if (scftype .eq. 'UHF') then
            call ga_copy(g_movecs(1),g_movecs(2))
            call dcopy(nbf,dbl_mb(k_eval),1,dbl_mb(k_eval+nbf),1)
         endif
c
         if (.not.inp_compare(.false., 'dft', calling_module))then
            if (.not. ga_destroy(g_dens(1))) call errquit
     &        ('scf_vectors_guess: destroy of density failed',0,
     &       GA_ERR)
         endif
c
      else if (movecs_in .eq. 'project') then
c         if (nmo .ne. nbf) call errquit
c     &        ('projection not yet working with linear dependence ',0)
c
c        Project from a smaller basis set ... must orthog before
c        the symmetry adaption. Also, sym adapt thoroughly
c        since the strong linear dependence in the 
c        random vectors may cause loss of precision.
c
         call movecs_project(rtdb, nmo, geom, basis, g_movecs,1,
     A        calling_module,k_eval)
         call ga_orthog_mos(basis, g_movecs)
         if (oadapt) then
            call scf_movecs_sym_adapt(basis, g_movecs, .false.,
     &           0, '- projected ', .true.,int_mb(k_irs))
            call scf_movecs_sym_adapt(basis, g_movecs, .false.,
     &           0, '- projected ', .true.,int_mb(k_irs))
         endif
         if (scftype .eq. 'UHF') then
         call movecs_project(rtdb, nmo, geom, basis, g_movecs,2,
     A           calling_module,k_eval+nbf)
         call ga_orthog_mos(basis, g_movecs(2))
         if (oadapt) then
            call scf_movecs_sym_adapt(basis, g_movecs(2), .false.,
     &           0, '- projected ', .true.,int_mb(k_irs+nmo))
            call scf_movecs_sym_adapt(basis, g_movecs(2), .false.,
     &           0, '- projected ', .true.,int_mb(k_irs+nmo))
          endif
          endif
c
      else if (movecs_in .eq. 'rotate') then
         nsets = 1
         if (scftype .eq. 'UHF') nsets = 2
         call movecs_rotate(rtdb, geom, basis,nsets, g_movecs,
     .        dbl_mb(k_eval),dbl_mb(k_occ))
c
      else if (movecs_in .eq. 'fragment') then
         if (scftype .ne. 'UHF') then
            call movecs_fragment(rtdb, geom, basis, 
     &           calling_module(1:inp_strlen(calling_module)), 1,
     &           nclosed, nopen, g_movecs(1), dbl_mb(k_eval))
         else
            call movecs_fragment(rtdb, geom, basis, 
     &           calling_module(1:inp_strlen(calling_module)), 1,
     &           0, nalpha, g_movecs(1), dbl_mb(k_eval))
            call movecs_fragment(rtdb, geom, basis, 
     &           calling_module(1:inp_strlen(calling_module)), 2,
     &           0, nbeta, g_movecs(2), dbl_mb(k_eval+nbf))
         endif
      else if (movecs_in .eq. 'molden') then
         nsets = 1
         if (scftype .eq. 'UHF') nsets = 2
         call movecs_molden(rtdb, geom, basis,nsets, nmo, g_movecs,
     .        dbl_mb(k_eval),dbl_mb(k_occ))
         call scf_movecs_sym_adapt(basis, g_movecs, oprint,
     &        0, '- initial ', .true.,int_mb(k_irs))
      else
c
c        Attempt to read vectors from file movecs_in and see if 
c        they make sense
c
         status = movecs_read_header(movecs_in, title_vecs, bas_vecs,
     &        file_scftype, nbf_vecs, nsets, nmo_vecs, 2)
         if (status) then
c     
c        Check that:
c        1) nbf from file = current nbf
c        2) name of basis from file = current basis name
c        3) nmo from file = current nmo (in case of change of linear dep)
c     
            if (nbf.eq.nbf_vecs .and. bas_vecs.eq.trans .and.
     &         nmo.eq.nmo_vecs(1)) then
c
c              If there is linear dependency, save the eigen vectors of
c              S, read in the MOs (perhaps from another geometry), and project
c              them into the linearly independent MO space of this geometry.
c
               if (nbf .ne. nmo) then
                  if (.not. ga_create(MT_DBL, nbf, nmo, 'scf_v_g:tmp',
     &                 32,32, g_tmp)) call errquit('scf_v_g: tmp', 0,
     &       GA_ERR)
                  call ga_copy(g_movecs, g_tmp)
               endif
c
c              Try reading movecs
c
               if (.not. movecs_read(movecs_in, 1,
     &            dbl_mb(k_occ), dbl_mb(k_eval), g_movecs))
     &            call errquit('scf_movecs_read failed ',0, DISK_ERR)
c
               if (scftype .eq. 'UHF') then
                  if (.not. movecs_read(movecs_in, 2,
     &               dbl_mb(k_occ+nbf), dbl_mb(k_eval+nbf),
     &               g_movecs(2))) then
                     if (ga_nodeid() .eq. 0 .and.
     &                  util_print('information', print_default)) then
                        write(6,*) 
     &                  '  Duplicating RHF/ROHF vectors for UHF'
                     endif
                     call ga_copy(g_movecs(1), g_movecs(2))
                     call dcopy(nbf,dbl_mb(k_occ),1,dbl_mb(k_occ+nbf),1)
                     call dcopy(nbf,dbl_mb(k_eval),1,dbl_mb(k_eval+nbf),
     &                          1)
                     if (inp_compare(.false., 'dft', 
     &                    calling_module)) 
     &                    call dscal(2*nbf,0.5d0,dbl_mb(k_occ),1)
                     
                  endif
               endif
c               call scf_movecs_read
*
               if (nbf .ne. nmo) then
                  if (.not. ga_create(MT_DBL, nmo, nmo, 'scf_v_g:tmp2',
     &                 32,32, g_tmp2)) call errquit('scf_v_g: tmp2', 0,
     &       GA_ERR)
                  call ga_dgemm('t', 'n', nmo, nmo, nbf, 1.0d0, 
     &                 g_tmp, g_movecs, 0.0d0, g_tmp2)
                  call ga_dgemm('n', 'n', nbf, nmo, nmo, 1.0d0,
     &                 g_tmp, g_tmp2, 0.0d0, g_movecs)
                  if (scftype .eq. 'UHF') then
                     call ga_dgemm('t', 'n', nmo, nmo, nbf, 1.0d0,
     &                    g_tmp, g_movecs(2), 0.0d0, g_tmp2)
                     call ga_dgemm('n', 'n', nbf, nmo, nmo, 1.0d0, 
     &                    g_tmp, g_tmp2, 0.0d0, g_movecs(2))
                  endif
                  if (.not. ga_destroy(g_tmp))
     &                 call errquit('scf_v_g: ga_destroy?',0, GA_ERR)
                  if (.not. ga_destroy(g_tmp2))
     &                 call errquit('scf_v_g: ga_destroy?',0, GA_ERR)
               endif
            else
               if (ga_nodeid() .eq. 0 .and. 
     $              util_print('warning',print_low)) then
                  write(6,111) nbf, nbf_vecs, nmo, nmo_vecs(1),
     &              trans(1:inp_strlen(trans)),
     &              bas_vecs(1:inp_strlen(bas_vecs))
 111              format(/' !! nbf/nmo/basis-name mismatch '/
     &                 10x, ' nbf=', i4, ' nbf_file=',i4/
     &                 10x, ' nmo=', i4, ' nmo_file=',i4/
     &                 10x, ' basis="',a,'"', ' basis_file="',a,'"'/)
                  write(6,*) ' Either an incorrect movecs file was ',
     &                 'specified, or linear dependence has changed,'
                  write(6,*) ' or the basis name was changed. '
                  call util_flush(6)
               endif
               status = .false.
            endif
            if (oprint) then
               write(6,77) title_vecs(1:inp_strlen(title_vecs))
 77            format(/' Loading old vectors from job with title :'//a/)
               call util_flush(6)
            endif
         endif
         if (.not. status) then
*     Rather then failing simply force atomic guess
*     call errquit('rhf: failed to load vectors',0)
            if (ga_nodeid() .eq. 0 .and. 
     $           util_print('warning',print_low)) then
               write(6,*) 
               write(6,*) ' Load of old vectors failed. ',
     &              '  Forcing atomic density guess'
               write(6,*) 
            endif
            movecs_in = 'atomic'
            goto 3212
         endif
      endif
c
c     Apply any swaps of orbitals
c
      call movecs_swap(rtdb, 
     &     calling_module(1:inp_strlen(calling_module)), scftype, 
     &     g_movecs, dbl_mb(k_occ), dbl_mb(k_eval))
c
c
c    Apply any swaps of nuclei in the mapping of mo vectors
c
       if ( rtdb_ma_get(rtdb, 'reorder', type, nelem, handle)) then
         call movecs_swap_nuclei(rtdb,basis,geom, scftype,g_movecs)
       endif
c
c
c     Optionally symmetry adapt the initial molecular orbitals
c
      if (oadapt) then
         if (scftype .eq. 'UHF') then
c     
c     This may not be a good idea for UHF, but we will do it
c     if we are asked to!
c     
            nprint = min(nalpha+10,nmo)
            if (util_print('all vector symmetries', print_high))
     &           nprint = nmo
            call scf_movecs_sym_adapt(basis, g_movecs, oprint,
     &           nprint, '- initial alpha',.true.,int_mb(k_irs))
            call scf_movecs_sym_adapt(basis, g_movecs(2), oprint,
     &           nprint, '- initial beta', .true.,int_mb(k_irs+nmo))
         else
            nprint = min(nclosed+nopen+10,nmo)
            if (util_print('all vector symmetries', print_high))
     &           nprint = nmo
            call scf_movecs_sym_adapt(basis, g_movecs, oprint,
     &           nprint, '- initial ', .true.,int_mb(k_irs))
         endif
      endif
c
c     Orthogonalize all input/guess vectors ... do not mix the occupied
c     with the virtuals.  For the hcore/atomic guess we will have done a
c     diagonalization so there is no need to orthog. again.
c     (also if the symmetry projection found mixing an orthog is done)
c
c     Also record if this is a restart or a startup.
c
      if ((movecs_in.eq.'atomic' .or. movecs_in.eq.'hcore')) then
         movecs_guess = 'start'
      else
         movecs_guess = 'restart'
      endif
c
      if (movecs_guess.eq.'restart') then
         call ga_orthog_mos(basis, g_movecs)
         if (scftype .eq. 'UHF') call ga_orthog_mos(basis, g_movecs(2))
      endif
      if(ssmear.lt.1.e-9 .and. (.not. fon)) then
c
c     Generate correct occupation numbers
c
      if (scftype .eq. 'UHF') then
         call dfill(2*nbf, 0.0d0, dbl_mb(k_occ), 1)
         do i = 1, nalpha
            dbl_mb(i-1+k_occ) = 1.0d0
         enddo
         do i = 1, nbeta
            dbl_mb(i-1+nbf+k_occ) = 1.0d0
         enddo
      else
         call dfill(nbf, 0.0d0, dbl_mb(k_occ), 1) ! Current occupations
         do i = 1, nclosed
            dbl_mb(i-1+k_occ) = 2.0d0
         enddo
         do i = nclosed+1,nclosed+nopen
            dbl_mb(i-1+k_occ) = 1.0d0
         enddo
      endif
      endif
c
c     Apply phase conventions
c
      call movecs_fix_phase(g_movecs)
      if (scftype .eq. 'UHF') call movecs_fix_phase(g_movecs(2))
c
      if (util_print('initial vectors', print_debug)) then
         if (ga_nodeid() .eq. 0) then
            write(6,*)
            call util_print_centered(6,'Initial MO vectors',40,.true.)
            write(6,*)
            call util_flush(6)
         endif
         call ga_print(g_movecs)
         if (scftype .eq. 'UHF') call ga_print(g_movecs(2))
      endif
c
c     Save vectors to disk
c
      nmo_vecs(1) = nmo
      if (scftype .eq. 'UHF') then
         nsets = 2
         nmo_vecs(2) = nmo
      else
         nsets = 1
      endif
c
      if (.not. movecs_write(rtdb, basis, movecs_out, 
     &     calling_module(1:inp_strlen(calling_module)), title, 
     &     nbf, nsets, nmo_vecs,
     &     dbl_mb(k_occ), nbf, dbl_mb(k_eval),
     &     nbf, g_movecs)) call errquit
     &     ('scf_v_g: movecs_write failed', 0, DISK_ERR)
c
c      call scf_movecs_write(rtdb)
c     
c     Overwrite the input movecs option with where the final movecs
c     will end up so that a restart automatically picks them up
c
      string = ' '
      write(string,'(a,a)') 
     &   calling_module(1:inp_strlen(calling_module)),':input vectors'
      if (.not. rtdb_cput(rtdb, string(1:inp_strlen(string)),
     &     1, movecs_out)) 
     &     call errquit('scf_v_g:rtdb_cput failed scf:input vectors',0,
     &       RTDB_ERR)
      string = ' '
      write(string,'(a,a)') 
     &   calling_module(1:inp_strlen(calling_module)),':output vectors'
      if (.not. rtdb_cput(rtdb, string(1:inp_strlen(string)),
     &     1, movecs_out)) 
     &     call errquit('scf_v_g:rtdb_cput failed scf:output vectors',0,
     &       RTDB_ERR)
c
      movecs_in = movecs_out
c
c     Print out an analysis if explicitly requested
c
      if (util_print('initial vector analysis', print_debug)) then
         nprint = min(nclosed+nopen+30,nmo)
         if (scftype.eq.'RHF' .or. scftype.eq.'ROHF') then
            call movecs_print_anal(basis, 1,
     &           nprint, 0.15d0, g_movecs, 
     &           'ROHF Initial Molecular Orbital Analysis', 
     &           .true., dbl_mb(k_eval), oadapt, int_mb(k_irs),
     &           .true., dbl_mb(k_occ))
         else
            nprint = min(nalpha+20,nmo)
            call movecs_print_anal(basis, max(1,nbeta-20),
     &           nprint, 0.15d0, g_movecs, 
     &           'UHF Initial Alpha Molecular Orbital Analysis', 
     &           .true., dbl_mb(k_eval), oadapt, int_mb(k_irs),
     &           .true., dbl_mb(k_occ))
            call movecs_print_anal(basis, max(1,nbeta-20),
     &           nprint, 0.15d0, g_movecs(2), 
     &           'UHF Initial Beta Molecular Orbital Analysis', 
     &           .true., dbl_mb(k_eval+nbf), oadapt, int_mb(k_irs+nmo),
     &           .true., dbl_mb(k_occ+nbf))
         endif
      endif
c
c     If 'density' has been requested form it here before returning.
c
      if (inp_compare(.false., 'density', vecs_or_dens))then
         if (scftype .eq. 'UHF') then
            call ga_dgemm('n', 't', nbf, nmo, nalpha, one, g_movecs(1),
     &      g_movecs(1), zero, g_dens(1))
            call ga_dgemm('n', 't', nbf, nmo, nbeta, one, g_movecs(2),
     &      g_movecs(2), zero, g_dens(2))
         else
            call ga_dgemm('n', 't', nbf, nmo, nclosed, two, g_movecs(1),
     &      g_movecs(1), zero, g_dens(1))
         endif
c
c        Symmetrize the guess density
c
         if (oadapt)then
            call sym_symmetrize(geom, basis, .true., g_dens(1))
            if (scftype .eq. 'UHF')
     &         call sym_symmetrize(geom, basis, .true., g_dens(2))
         endif
c
c        Print guess density if debugging
c
         if (util_print('initial density', print_debug)) then
            if (ga_nodeid() .eq. 0) then
               write(6,*)
               call util_print_centered(6,'Initial Density',40,.true.)
               write(6,*)
               call util_flush(6)
            endif
            call ga_print(g_dens)
            if (scftype .eq. 'UHF') call ga_print(g_dens(2))
         endif
c
      endif
      return
      end
