       subroutine cphf_poliz( rtdb, geom, basis, movecs )
       implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "context.fh"
c
c
       integer rtdb
       integer geom
       integer basis
       character*255 movecs
       integer gx, gy, gz
       logical oskel
c
       character*255 title_vecs, basis_vecs
       integer nbf, nsets, nmo
       logical status
c
       write(6,900) movecs
 900   format(/,1x,'Inside CPHF_POLIZ',
     $        /,'Movecs file:',a80)
c
c     read in the mo vectors
c
       status = movecs_read_header(movecs, title_vecs, basis_vecs,
     $      nbf, nsets, nmo, 1)
       if (.not. status) call errquit('cphf: no vectors?',0)
c
      g_vecs = ga_create_atom_blocked(geom, basis,'Molecular orbital')
      if (.not. ma_push_get(mt_dbl, nbf, 'mo evals',l_evals, k_evals))
     $     call errquit('cphf: no space for eigenvalues?', nbf)
      if (.not. ma_push_get(mt_dbl, nbf, 'mo occ',l_occ, k_occ))
     $     call errquit('cphf: no space for occupation?', nbf)
c
      if (.not. status) call errquit('rhf: insufficient memory ', nbf)
       status = movecs_read(movecs, 1, dbl_mb(k_occ),
     $      dbl_mb(k_evals), g_vecs)
c
c     initialize integrals and make sure the AO basis is normalized
c
      call int_init(rtdb,1, basis)
      call int_normalize(basis)
c
c...  collect the dipole moment integrals from int_mpole:
c
c
       g_x = ga_create_atom_blocked(geom, basis, 'cphf: dip x')
       g_y = ga_create_atom_blocked(geom, basis, 'cphf: dip y')
       g_z = ga_create_atom_blocked(geom, basis, 'cphf: dip z')
       call int_dip_ga(basis, basis, g_x, g_y, g_z)
c
c     transform dipole moments into the MO basis
c

c
c     Intialize stuff needed for the fock build
c
      call rhf_fock_param(rtdb)
      call schwarz_init(geom, basis)
c


c

      call schwarz_tidy()
      call int_terminate
c
c

       return
       end
