      subroutine mp2_lai_fock_uhf(geom, basis,
     $     g_pab_a, g_pab_b, g_pij_a, g_pij_b,
     $     g_vecs_a, g_vecs_b,
     $     no_a, no_b, nv_a, nv_b,
     $     g_lai_a, g_lai_b)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer geom, basis        ! [input] 
      integer g_pab_a, g_pab_b, g_pij_a, g_pij_b ! [input] GA densities
      integer g_vecs_a, g_vecs_b ! [input] MO vector handles
      integer no_a, no_b, nv_a, nv_b ! [input] occ/virt dimensions
      integer g_lai_a, g_lai_b  ! [output] Accumulate results into
c
c     g_lai += Pbc [(ia|cb) - (ib|ca)] + Pjk [(kj|ia) - (ka|ij)]
c
c            = CubPbcCvc [(ia|uv) - (iu|av)] + CujPjkCvk [(ia|uv) - (iu|av)]
c
c            = J - K with density from symmetrized total AO density
c
      integer g_dens(2), g_fock(2), nfock, nmo, nbf, gtype
      integer g_tmp1, g_tmp2
      double precision jfac(2), kfac(2), tol2e
      logical oskel
c
      data nfock /2/
      data jfac, kfac /1.0d0, 1.0d0, -1.0d0, -1.0d0/
      data tol2e /1.0d-9/
      data oskel /.false./
c
c     Allocate space for AO fock and density matrices
c
      call ga_create_atom_blocked(geom,basis,'mp2_grad: da',g_dens(1))
      call ga_create_atom_blocked(geom,basis,'mp2_grad: db',g_dens(2))
      call ga_create_atom_blocked(geom,basis,'mp2_grad: fa',g_fock(1))
      call ga_create_atom_blocked(geom,basis,'mp2_grad: fb',g_fock(2))
      call ga_zero(g_dens(1))
      call ga_zero(g_dens(2))
      call ga_zero(g_fock(1))
      call ga_zero(g_fock(2))
c
c     Workspace
c
      call ga_inquire(g_vecs_a, gtype, nbf, nmo)
      if (no_a+nv_a.ne.nmo .or. no_b+nv_b.ne.nmo)
     $     call errquit('mp2_grad: weird nmo?',0)
c      
      if (.not. ga_create(mt_dbl, nmo, nmo, 'mp2_grad: tmp1',
     $     1, 1, g_tmp1)) call errquit('mp2_lai_fock_uhf: GA', nmo*nmo)
      if (.not. ga_create(mt_dbl, nbf, nmo, 'mp2_grad: tmp2',
     $     1, 1, g_tmp2)) call errquit('mp2_lai_fock_uhf: GA', nmo*nmo)
c
c     Copy OO and VV patches into N*N array, transform and symmetrize
c
      call ga_zero(g_tmp1)
      call ga_copy_patch('n',
     $     g_pij_a, 1, no_a, 1, no_a, g_tmp1, 1, no_a, 1, no_a)
      call ga_copy_patch('n',
     $     g_pab_a, 1, nv_a, 1, nv_a, g_tmp1, no_a+1, nmo, no_a+1, nmo)
      call ga_dgemm('n', 'n', nbf, nmo, nmo, 1.0d0, g_vecs_a,
     $     g_tmp1, 0.0d0, g_tmp2)
      call ga_dgemm('n', 't', nbf, nbf, nmo, 1.0d0, g_tmp2,
     $     g_vecs_a, 0.0d0, g_dens(1))
      call ga_symmetrize(g_dens(1))
c
      call ga_zero(g_tmp1)
      call ga_copy_patch('n',
     $     g_pij_b, 1, no_b, 1, no_b, g_tmp1, 1, no_b, 1, no_b)
      call ga_copy_patch('n',
     $     g_pab_b, 1, nv_b, 1, nv_b, g_tmp1, no_b+1, nmo, no_b+1, nmo)
      call ga_dgemm('n', 'n', nbf, nmo, nmo, 1.0d0, g_vecs_b,
     $     g_tmp1, 0.0d0, g_tmp2)
      call ga_dgemm('n', 't', nbf, nbf, nmo, 1.0d0, g_tmp2,
     $     g_vecs_b, 0.0d0, g_dens(2))
      call ga_symmetrize(g_dens(2))
c
      if (.not. ga_destroy(g_tmp1)) call errquit('mp2_lfu: ga?',0)
      if (.not. ga_destroy(g_tmp2)) call errquit('mp2_lfu: ga?',0)
c
c     Make the fock matrices
c
      call fock_2e(geom, basis, nfock, jfac, kfac,
     $     tol2e, oskel, g_dens, g_fock)
c
c     Transform back
c
      g_tmp1 = g_dens(1)
c
      call ga_zero(g_tmp1)
      call ga_dgemm('n', 'n', nbf, nbf, nv_a, 1.0d0, g_fock(1), 
     $     g_vecs_a, 0.0d0, g_tmp1)
      call ga_dgemm('t', 'n', nv_a, no_a, nbf, 1.0d0, g_vecs_a,
     $     g_tmp1, 1.0d0, g_lai_a)
c
      call ga_zero(g_tmp1)
      call ga_dgemm('n', 'n', nbf, nbf, nv_b, 1.0d0, g_fock(1), 
     $     g_vecs_b, 0.0d0, g_tmp1)
      call ga_dgemm('t', 'n', nv_b, no_b, nbf, 1.0d0, g_vecs_b,
     $     g_tmp1, 1.0d0, g_lai_b)
c
      if (.not. ga_destroy(g_dens(1))) call errquit('mp2_lfu: ga?',0)
      if (.not. ga_destroy(g_dens(2))) call errquit('mp2_lfu: ga?',0)
      if (.not. ga_destroy(g_fock(1))) call errquit('mp2_lfu: ga?',0)
      if (.not. ga_destroy(g_fock(2))) call errquit('mp2_lfu: ga?',0)
*
      end

