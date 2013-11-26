      subroutine dft_cg_answer(rtdb,energy,converged)
      implicit none
c     
c     Store energy and convergence status ... must store before
c     write movecs since date of insertion is used.
c     
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "cscf.fh"
c
      integer rtdb
      logical converged       ! .true. if the calculation converged
      double precision energy ! the final total energy
c
      logical  movecs_write
      external movecs_write
c
      integer nsets
      integer nmo_t(2)
c
      nmo_t(1) = nmo
      if (scftype.eq.'UHF') then
        nsets = 2
        nmo_t(2) = nmo
      else
        nsets = 1
      endif
c
c     == needed for dftmp2 ==
      if (.not. rtdb_cput(rtdb,'scf:scftype', 1, scftype))
     &   call errquit('dft_scf: rtdb_cput scftype failed', 0, RTDB_ERR)
      if (.not. rtdb_put(rtdb,'scf:nclosed', mt_int, 1, nclosed))
     &   call errquit('dft_scf: rtdb_put nclosed failed', 0, RTDB_ERR)

c     if (.not. rtdb_put(rtdb,'uhf:coulomb', mt_dbl, 1, Ecoul)) call
c    $     errquit('uhf: writing ecoul failed', 0, RTDB_ERR)
c
      if (.not. rtdb_put(rtdb, 'dft:energy', MT_DBL, 1, energy))
     &   call errquit('dft_scf: failed to store energy in rtdb', 0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'dft:converged', MT_LOG, 1, converged))
     &   call errquit('dft_scf: failed to store converged in rtdb', 0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'dft:alpha irreps', MT_INT, nbf, 
     &   int_mb(k_irs)))
     &   call errquit('dft_scf: failed to MO irreps in rtdb', 0,
     &       RTDB_ERR)
      if (nsets.eq.2) then
        if (.not. rtdb_put(rtdb, 'dft:beta irreps', MT_INT, nbf, 
     &     int_mb(k_irs+nbf)))
     &     call errquit('dft_scf: failed to MO irreps in rtdb', 0,
     &         RTDB_ERR)
      endif
c
c     output energies and eigenvectors to disk
c     
      if (.not.movecs_write(rtdb, basis, movecs_out, 'dft', title,
     &                      nbf, nsets, nmo_t, dbl_mb(k_occ), nbf, 
     &                      dbl_mb(k_eval), nbf, g_movecs))
     &    call errquit('dft_scf: movec_write failed', 0, DISK_ERR)
c
      end
c $Id$
