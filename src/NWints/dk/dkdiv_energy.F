      logical function dkdiv_energy(rtdb)
c $Id$
      implicit none
#include "global.fh"
#include "rtdb.fh"
#include "context.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "geom.fh"
#include "util.fh"
#include "stdio.fh"
#include "inp.fh"
#include "errquit.fh"
c
c     The routine computes the energy for a given basis, geometry
c     and density (read from the rtdb) for the integral difference
c     (Douglas-Kroll - non-relativistic Hamiltonian)
c
      integer rtdb ! [input] handle to database
c
      logical int_normalize, file_read_ga
      external int_normalize, file_read_ga
c
      character*255 dkdiv_name
      integer DKg_a, g_dens
      integer geom, basis
      integer nbf
      double precision energy
c
c     Get basis, geometry, density or densities and additional data
c
      if (.not. geom_create(geom, 'DK geometry'))
     $     call errquit('dkdiv_energy: geom_create failed', 0,
     &       GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('dkdiv_energy: geom_rtdb_load failed', 0,
     &       RTDB_ERR)
      if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('dkdiv_energy: bas_create failed', 0,
     &       BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     $     call errquit('dkdiv_energy: bas_rtdb_load failed', 0,
     &       RTDB_ERR)
c
      call int_init(rtdb,1,basis)
c
      if (.not. int_normalize(rtdb,basis))
     $     call errquit('dkdiv_energy: int_normalize failed',0,
     &       INT_ERR)
      if (.not. bas_numbf(basis,nbf))
     $   call errquit('dkdiv_energy: bas_numbf failed',0, BASIS_ERR)
c
c     Get kinetic and potential energy integral matrix. Multiply
c     by -1 to setup for DK - NR
c
      if (.not. ga_create(MT_DBL,nbf,nbf,'GA dkdiv',1,1,DKg_a))
     $     call errquit('dkdiv_energy:ga_create failed GA dkdiv',0,
     &       GA_ERR)
      call ga_zero(DKg_a)
      call int_1e_oldga(basis,basis,DKg_a,'kinetic',.false.)
      call int_1e_oldga(basis,basis,DKg_a,'potential',.false.)
      call ga_dscal(Dkg_a, -1.0d0)
c
c     Add Douglas-Kroll integrals to GA with - kinetic - potential
c     We are left with the required difference integrals
c
      call int_1edk_ga(basis,DKg_a,3,.false.)
c
c     Get density from file back into a GA
c
      if (.not. ga_create(MT_DBL,nbf,nbf,'GA dens',1,1,g_dens))
     $     call errquit('dkdiv_energy:ga_create failed GA dens',0,
     &       GA_ERR)
      call ga_zero(g_dens)
      call util_file_name('dkdiv_density',.true.,.true.,dkdiv_name)
      if (.not. file_read_ga(dkdiv_name,g_dens))
     $    call errquit('dkdiv_energy:file_read_ga density failed',0,
     &       GA_ERR)
c
c     Compute the energy E = Density * Integrals
c
      energy = 0.0d0
      energy = ga_ddot(DKg_a,g_dens)
c
c     Store energy in database under specific name: grad_dk:dkdiv_energy
c
      if (.not. rtdb_put(rtdb,'grad_dk:dkdiv_energy',MT_DBL,1,energy))
     $   call errquit('dkdiv_energy:rtdb_put energy failed',0,
     &       RTDB_ERR)
      dkdiv_energy = .true.
c
      if (.not.((ga_destroy(g_dens)).and.(ga_destroy(DKg_a))))
     $   call errquit('dkdiv_energy: destroy GAs failed',0, GA_ERR)
      if (.not.((bas_destroy(basis)).and.(geom_destroy(geom))))
     $    call errquit('dkdiv_energy: destroy bas or geom failed',0,
     &       MEM_ERR)
c
      call int_terminate()
c
      return
c
      end
