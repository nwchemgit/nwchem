      subroutine scf_vectors_guess(rtdb)
C     $Id: scf_vec_guess.F,v 1.1 1995-12-13 01:26:14 d3g681 Exp $
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "cscf.fh"
c
c     Form guess vectors on disk for SCF wavefunctions
c
      integer rtdb
c
c     local variables
c
      logical oprint
      integer g_aodens          ! Atomic density guess
      integer i, nprint
      character*255 title_vecs, bas_vecs
      integer nbf_vecs, nsets, nmo_vecs(2)
      logical status
c
c     externals
c
      integer ga_create_atom_blocked
      logical movecs_read_header
      external movecs_read_header
      external ga_create_atom_blocked
c
      if (.not. oinitialized) call errquit
     $     ('scf_v_g: not initialized?', 0)
c
      oprint = util_print('mo guess', print_default) .and.
     $     ga_nodeid().eq.0
c
      call ga_zero(g_movecs)
c
      if (movecs_in .eq. 'hcore' .or. movecs_in .eq. 'atomic') then
c
c     hcore or atomic guess
c
         g_aodens = ga_create_atom_blocked(geom, basis, 'Guess density')
         call ga_zero(g_aodens)
         if (movecs_in .eq. 'atomic') then
            call guess_dens(geom, basis, g_aodens)
         endif
         call rhf_dens_to_mo(rtdb, geom, basis, 2*nclosed+nopen,
     $        nbf, tol2e,
     $        g_aodens, g_movecs, oskel)
c
         if (scftype .eq. 'UHF') then
            call ga_copy(g_movecs(1),g_movecs(2))
            call dcopy(nbf,dbl_mb(k_eval),1,dbl_mb(k_eval+nbf),1)
         endif
c
         if (.not. ga_destroy(g_aodens)) call errquit
     $        ('scf_vectors_guess: destroy of density failed',0)
c
      else if (movecs_in .eq. 'project') then
c
c     Project from a smaller basis set ... must orthog before
c     the symmetry adaption
c
         call movecs_project(rtdb, geom, basis, g_movecs)
         call ga_orthog_mos(basis, g_movecs)
         if (scftype .eq. 'UHF') then
            call ga_copy(g_movecs(1),g_movecs(2))
            call dcopy(nbf,dbl_mb(k_eval),1,dbl_mb(k_eval+nbf),1)
         endif
      else
c
c     Attempt to read vectors from file movecs_in and see if 
c     they make sense
c
         status = movecs_read_header(movecs_in, title_vecs, bas_vecs,
     $        nbf_vecs, nsets, nmo_vecs, 2)
         if (status) then
c     
c     Check that:
c     1) nbf from file = current nbf
c     2) name of basis from file = current basis name
c     
            if (nbf .eq. nbf_vecs .and. bas_vecs.eq.trans) then
               call scf_movecs_read
            else
               If (ga_nodeid() .eq. 0 )
     $            write(6,*) ' nbf/basis name mismatch ',
     $              ' nbf=', nbf, ' nbf_file=',nbf_vecs, 
     $              ' trans="',trans(1:inp_strlen(trans)),'"',
     $              ' trans_file="',
     $              bas_vecs(1:inp_strlen(bas_vecs)),'"'
               status = .false.
            endif
            if (oprint) then
               write(6,77) title_vecs(1:inp_strlen(title_vecs))
 77            format(/' Loading old vectors from job with title :'//a/)
               call util_flush(6)
            endif
         endif
         if (.not. status) call errquit('rhf: failed to load vectors',0)
      endif
c
c     Apply any swaps of orbitals
c
      call movecs_swap(rtdb, scftype, g_movecs,
     $     dbl_mb(k_occ), dbl_mb(k_eval))
c
c     Optionally symmetry adapt the initial molecular orbitals
c
      if (oadapt) then
         if (scftype .eq. 'UHF') then
c     
c     This may not be a good idea for UHF, but we will do it
c     if we are asked to!
c     
            nprint = min(nalpha+10,nbf)
            if (util_print('all vector symmetries', print_high))
     $           nprint = nbf
            call scf_movecs_sym_adapt(basis, g_movecs, oprint,
     $           nprint, '- initial alpha')
            call scf_movecs_sym_adapt(basis, g_movecs(2), oprint,
     $           nprint, '- initial beta')
         else
            nprint = min(nclosed+nopen+10,nbf)
            if (util_print('all vector symmetries', print_high))
     $           nprint = nbf
            call scf_movecs_sym_adapt(basis, g_movecs, oprint,
     $           nprint, '- initial ')
         endif
      endif
c
c     Orthogonalize all input/guess vectors ... do not mix the occupied
c     with the virtuals.  For the atomic guess we will have done a
c     diagonalization so there is no need to orthg. again.
c
      if (movecs_in .ne. 'atomic') then
         call ga_orthog_mos(basis, g_movecs)
         if (scftype .eq. 'UHF') call ga_orthog_mos(basis, g_movecs(2))
      endif
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
      call scf_movecs_write
c     
c     Overwrite the input movecs option with where the final movecs
c     will end up so that a restart automatically picks them up
c
      if (.not. rtdb_cput(rtdb, 'scf:input vectors', 1, movecs_out))
     $     call errquit('scf:rtdb_cput failed scf:input vectors',0)
c
      movecs_in = movecs_out
c
      end

