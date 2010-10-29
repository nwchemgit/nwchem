      subroutine mp2_lai_fock_uhf_prepar(g_pab_a, g_pab_b,
     $     g_pij_a, g_pij_b, no_a, no_b, nv_a, nv_b, 
     $     g_p_a, g_p_b, nmo)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer g_pab_a, g_pab_b,
     $     g_pij_a, g_pij_b, no_a, no_b, nv_a, nv_b, 
     $     g_p_a, g_p_b, nmo 
c
c     Put the OO and VV blocks of P into the appropriate
c     diagonal blocks of a density matrix which is also 
c     allocated here.  This is in preparation for building
c     L3+l4. It used to be part of lai_fock_uhf but has
c     been pulled out so that that routine can be used by the RIMP2.
c
*ga:1:0
      if (.not. ga_create(mt_dbl, nmo, nmo, 'mp2_grad: p_a',
     $     0, 0, g_p_a)) call errquit('mp2_lai_fock_uhf: GA', nmo*nmo,
     &       GA_ERR)
*ga:1:0
      if (.not. ga_create(mt_dbl, nmo, nmo, 'mp2_grad: p_b',
     $     0, 0, g_p_b)) call errquit('mp2_lai_fock_uhf: GA', nmo*nmo,
     &       GA_ERR)
c
      call ga_zero(g_p_a)
      call ga_copy_patch('n',
     $     g_pij_a, 1, no_a, 1, no_a, g_p_a, 1, no_a, 1, no_a)
      call ga_copy_patch('n',
     $     g_pab_a, 1, nv_a, 1, nv_a, g_p_a, no_a+1, nmo, no_a+1, nmo)
c
      call ga_zero(g_p_b)
      call ga_copy_patch('n',
     $     g_pij_b, 1, no_b, 1, no_b, g_p_b, 1, no_b, 1, no_b)
      call ga_copy_patch('n',
     $     g_pab_b, 1, nv_b, 1, nv_b, g_p_b, no_b+1, nmo, no_b+1, nmo)
c
      end
      subroutine mp2_lai_fock_uhf_tidy(g_p_a,g_p_b)
      implicit none
#include "errquit.fh"
#include "global.fh"
      integer g_p_a,g_p_b
      if (.not. ga_destroy(g_p_a)) callerrquit('mp2_l_f_u_t: GA?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_p_b)) callerrquit('mp2_l_f_u_t: GA?',1,
     &       GA_ERR)
      end
      subroutine mp2_lai_fock_uhf(geom, basis,
     $     g_p_a, g_p_b, g_vecs_a, g_vecs_b,
     $     no_a, no_b, nv_a, nv_b,
     $     g_lai_a, g_lai_b, rtdb, tol2e)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "bas.fh"
      integer geom, basis       ! [input] 
      integer rtdb
      double precision tol2e
      integer g_p_a, g_p_b ! [input] GA densities
      integer g_vecs_a, g_vecs_b ! [input] MO vector handles
      integer no_a, no_b, nv_a, nv_b ! [input] occ/virt dimensions
      integer g_lai_a, g_lai_b  ! [output] Accumulate results into
      integer ga_create_atom_blocked
      external ga_create_atom_blocked

c     g_lai += J[Paa] - K[Paa] + J[Pbb]
c     
      integer g_dens(4), g_fock(4), nfock, nmo, nbf, gtype
      integer g_tmp, i
      double precision jfac(4), kfac(4), one, zero, mone
      logical oskel
c     
      data nfock /4/
      data jfac / 1.0d0, 1.0d0, 1.0d0, 1.0d0/
      data kfac /-1.0d0, 0.0d0,-1.0d0, 0.0d0/
      data oskel /.false./
c     
c     
      one = 1.0d0
      zero= 0.0d0
      mone = -1.0d0
c     
c     Allocate space for AO fock and density matrices
c     
      g_dens(1)=ga_create_atom_blocked(geom,basis,'mp2_grad:da')
      g_dens(2)=g_dens(1)
      g_dens(3)=ga_create_atom_blocked(geom,basis,'mp2_grad:db')
      g_dens(4)=g_dens(3)
      g_fock(1)=ga_create_atom_blocked(geom,basis,'mp2_grad:fa')
      g_fock(2)=ga_create_atom_blocked(geom,basis,'mp2_grad:fb')
      g_fock(3)=ga_create_atom_blocked(geom,basis,'mp2_grad:fb')
      g_fock(4)=ga_create_atom_blocked(geom,basis,'mp2_grad:fb')
      do i = 1, 4
         call ga_zero(g_dens(i))
         call ga_zero(g_fock(i))
      end do
      call scf_get_fock_param(rtdb,tol2e)
      call fock_force_direct(rtdb)    ! Force Fock build to be direct
c     
c     Workspace
c     
      call ga_inquire(g_vecs_a, gtype, nbf, nmo)
      if (no_a+nv_a.ne.nmo .or. no_b+nv_b.ne.nmo)
     $     call errquit('mp2_grad: weird nmo?',0, INPUT_ERR)
c     
*ga:1:0
      if (.not. ga_create(mt_dbl, nbf, nmo, 'mp2_grad: tmp',
     $     0, 0, g_tmp)) call errquit('mp2_lai_fock_uhf: GA', nmo*nmo,
     &       GA_ERR)
c     
      call ga_dgemm('n', 'n', nbf, nmo, nmo, one, g_vecs_a,
     $     g_p_a, zero, g_tmp)
      call ga_dgemm('n', 't', nbf, nbf, nmo, one, g_tmp,
     $     g_vecs_a, zero, g_dens(1))
      call ga_symmetrize(g_dens(1))
      call ga_dscal(g_dens(1),mone)
c     
      call ga_dgemm('n', 'n', nbf, nmo, nmo, one, g_vecs_b,
     $     g_p_b, zero, g_tmp)
      call ga_dgemm('n', 't', nbf, nbf, nmo, one, g_tmp,
     $     g_vecs_b, zero, g_dens(3))
      call ga_symmetrize(g_dens(3))
      call ga_dscal(g_dens(3),mone)
c     
      if (.not. ga_destroy(g_tmp)) call errquit('mp2_lfu: ga?',0,
     &       GA_ERR)
c     
c     Make the fock matrices
c     
      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)

      call fock_2e(geom, basis, nfock, jfac, kfac,
     $     tol2e, oskel, g_dens, g_fock,.false. )
c     
      call int_terminate()
      call schwarz_tidy()
      call fock_2e_tidy(rtdb)

      call ga_dadd(one, g_fock(1), one, g_fock(4), g_fock(1))
      call ga_dadd(one, g_fock(3), one, g_fock(2), g_fock(3))
      call ga_symmetrize(g_fock(1))
      call ga_symmetrize(g_fock(3))
c     
c     Transform back
c     
      g_tmp = g_dens(1)
c     
      call ga_zero(g_tmp)
      call ga_matmul_patch('n','n',one, zero, 
     $     g_fock(1), 1, nbf, 1, nbf,
     $     g_vecs_a, 1, nbf, no_a+1, nmo,
     $     g_tmp, 1, nbf, 1, nv_a)
      call ga_matmul_patch('t','n', one, one,
     $     g_vecs_a, 1, no_a, 1, nbf,
     $     g_tmp, 1, nbf, 1, nv_a,
     $     g_lai_a, 1, no_a, 1, nv_a)
c     
      call ga_zero(g_tmp)
      call ga_matmul_patch('n','n',one, zero, 
     $     g_fock(3), 1, nbf, 1, nbf,
     $     g_vecs_b, 1, nbf, no_b+1, nmo,
     $     g_tmp, 1, nbf, 1, nv_b)
      call ga_matmul_patch('t','n', one, one,
     $     g_vecs_b, 1, no_b, 1, nbf,
     $     g_tmp, 1, nbf, 1, nv_b,
     $     g_lai_b, 1, no_b, 1, nv_b)
c     
      if (.not. ga_destroy(g_dens(1))) call errquit('mp2_lfu: ga?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_dens(3))) call errquit('mp2_lfu: ga?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_fock(1))) call errquit('mp2_lfu: ga?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_fock(2))) call errquit('mp2_lfu: ga?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_fock(3))) call errquit('mp2_lfu: ga?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_fock(4))) call errquit('mp2_lfu: ga?',0,
     &       GA_ERR)
c     
      end
      subroutine mp2_wij_fock_uhf(rtdb, geom, basis, tol2e,
     $     g_ppq_a, g_ppq_b, 
     $     no_a, no_b,
     $     g_vecs_a, g_vecs_b,
     $     g_wij_a, g_wij_b)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer rtdb, geom, basis ! [input] 
      double precision tol2e    ! [input]
      integer no_a, no_b        ! [input] No. of occupied
      integer g_ppq_a, g_ppq_b  ! [input] GA densities
      integer g_vecs_a, g_vecs_b ! [input] MO vector handles
      integer g_wij_a, g_wij_b  ! [output] Accumulate results into
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c     
c     w_ij += J[Paa] - K[Paa] + J[Pbb]
c     
      integer g_dens(4), g_fock(4), nfock, nmo, nbf, gtype
      integer g_tmp1, g_tmp2, i
      double precision jfac(4), kfac(4), one, zero, mone
      logical oskel
c     
      data nfock /4/
      data jfac / 1.0d0, 1.0d0, 1.0d0, 1.0d0/
      data kfac /-1.0d0, 0.0d0,-1.0d0, 0.0d0/
      data oskel /.false./
c     
      one = 1.0d0
      zero= 0.0d0
      mone = -1.0d0
c     
c     Allocate space for AO fock and density matrices
c     
      g_dens(1)=ga_create_atom_blocked(geom,basis,'mp2_grad:da')
      g_dens(2)=g_dens(1)
      g_dens(3)=ga_create_atom_blocked(geom,basis,'mp2_grad:db')
      g_dens(4)=g_dens(3)
      g_fock(1)=ga_create_atom_blocked(geom,basis,'mp2_grad:fa')
      g_fock(2)=ga_create_atom_blocked(geom,basis,'mp2_grad:fb')
      g_fock(3)=ga_create_atom_blocked(geom,basis,'mp2_grad:fb')
      g_fock(4)=ga_create_atom_blocked(geom,basis,'mp2_grad:fb')
      do i = 1, 4
         call ga_zero(g_dens(i))
         call ga_zero(g_fock(i))
      end do
      call scf_get_fock_param(rtdb,tol2e)
      call fock_force_direct(rtdb)    ! Force Fock build to be direct
c     
c     Workspace
c     
      call ga_inquire(g_vecs_a, gtype, nbf, nmo)
c     
*ga:1:0
      if (.not. ga_create(mt_dbl, nbf, nmo, 'mp2_grad: tmp2',
     $     0, 0, g_tmp2)) call errquit('mp2_wij_fock_uhf: GA', nmo*nmo,
     &       GA_ERR)
c     
c     Transform the densities and symmetrize
c     
      call ga_dgemm('n', 'n', nbf, nmo, nmo, one, g_vecs_a,
     $     g_ppq_a, zero, g_tmp2)
      call ga_dgemm('n', 't', nbf, nbf, nmo, mone, g_tmp2,
     $     g_vecs_a, zero, g_dens(1))
      call ga_symmetrize(g_dens(1))
c     
      call ga_dgemm('n', 'n', nbf, nmo, nmo, one, g_vecs_b,
     $     g_ppq_b, zero, g_tmp2)
      call ga_dgemm('n', 't', nbf, nbf, nmo, mone, g_tmp2,
     $     g_vecs_b, zero, g_dens(3))
      call ga_symmetrize(g_dens(3))
c     
      if (.not. ga_destroy(g_tmp2)) call errquit('mp2_lfu: ga?',0,
     &       GA_ERR)
c     
c     Make the fock matrices
c     
      call int_init(rtdb,1,basis)
      call schwarz_init(geom,basis)

      call fock_2e(geom, basis, nfock, jfac, kfac,
     $     tol2e, oskel, g_dens, g_fock, .false.)
c     
      call int_terminate()
      call schwarz_tidy()
      call fock_2e_tidy(rtdb)

      call ga_dadd(one, g_fock(1), one, g_fock(4), g_fock(1))
      call ga_dadd(one, g_fock(3), one, g_fock(2), g_fock(3))
      call ga_symmetrize(g_fock(1))
      call ga_symmetrize(g_fock(3))
c     
c     Transform back ... just want wij
c     
      g_tmp1 = g_dens(1)
c     
      call ga_matmul_patch('n','n',one, zero, 
     $     g_fock(1), 1, nbf, 1, nbf,
     $     g_vecs_a, 1, nbf, 1, no_a,
     $     g_tmp1, 1, nbf, 1, no_a)
      call ga_matmul_patch('t','n', one, one,
     $     g_vecs_a, 1, no_a, 1, nbf,
     $     g_tmp1, 1, nbf, 1, no_a,
     $     g_wij_a, 1, no_a, 1, no_a)
c     
      call ga_matmul_patch('n','n',one, zero, 
     $     g_fock(3), 1, nbf, 1, nbf,
     $     g_vecs_b, 1, nbf, 1, no_b,
     $     g_tmp1, 1, nbf, 1, no_b)
      call ga_matmul_patch('t','n', one, one,
     $     g_vecs_b, 1, no_b, 1, nbf,
     $     g_tmp1, 1, nbf, 1, no_b,
     $     g_wij_b, 1, no_b, 1, no_b)
c     
      if (.not. ga_destroy(g_dens(1))) call errquit('mp2_lfu: ga?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_dens(3))) call errquit('mp2_lfu: ga?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_fock(1))) call errquit('mp2_lfu: ga?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_fock(2))) call errquit('mp2_lfu: ga?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_fock(3))) call errquit('mp2_lfu: ga?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_fock(4))) call errquit('mp2_lfu: ga?',0,
     &       GA_ERR)
c     
      end

