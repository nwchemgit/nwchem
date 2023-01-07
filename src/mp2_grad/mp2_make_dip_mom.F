      subroutine mp2_make_dip_mom(rtdb, basis, nbf, nmo, noa, nob,
     $     g_p_tot_a, g_p_tot_b, g_vecs_a, g_vecs_b)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "util.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
      integer basis,rtdb,nbf,nmo
      integer g_p_tot_a,g_p_tot_b
      integer g_vecs_a,g_vecs_b
      integer noa, nob
c
      integer g_x, g_y, g_z
      integer g_scratch, g_dens
      double precision mp2x, scfx
      double precision mp2y, scfy
      double precision mp2z, scfz
      double precision ndip(3)
      logical ga_create_atom_blocked
      external ga_create_atom_blocked
      integer geom
c
      if (.not. bas_cando_mpoles(basis)) return
c
      if (.not. bas_geom(basis, geom)) 
     $     call errquit('mp2_dip: bad basis?',0, BASIS_ERR)
c
c     Generate dipole integrals
c
      g_x = ga_create_atom_blocked(geom, basis,'ao integs x')
      g_y = ga_create_atom_blocked(geom, basis,'ao integs y')
      g_z = ga_create_atom_blocked(geom, basis,'ao integs z')
c
      call int_init(rtdb,1,basis)
      call int_dip_ga(basis, basis, g_x, g_y, g_z)
      call int_terminate()
c
c     Generate SCF AO density
c
*ga:1:0
      g_dens = ga_create_atom_blocked(geom, basis,'ao density')
      call ga_dgemm('n', 't', nbf, nbf, noa, 1.0d0, g_vecs_a,
     $     g_vecs_a, 0.0d0, g_dens)
      call ga_dgemm('n', 't', nbf, nbf, nob, 1.0d0, g_vecs_b,
     $     g_vecs_b, 1.0d0, g_dens)
c
      scfx = -ga_ddot(g_dens, g_x) ! Electrons are negative
      scfy = -ga_ddot(g_dens, g_y)
      scfz = -ga_ddot(g_dens, g_z)
c
c     Form the MP2 AO density
c
*ga:1:0
      g_scratch = ga_create_atom_blocked(geom, basis,'scratch')
      call ga_dgemm('n','t',nmo,nbf,nmo,1.0d0,g_p_tot_a,g_vecs_a,
     $     0.0d0,g_scratch)
      call ga_dgemm('n','n',nbf,nbf,nmo,1.0d0,g_vecs_a,g_scratch,
     $     0.0d0,g_dens)
      call ga_dgemm('n','t',nmo,nbf,nmo,1.0d0,g_p_tot_b,g_vecs_b,
     $     0.0d0,g_scratch)
      call ga_dgemm('n','n',nbf,nbf,nmo,1.0d0,g_vecs_b,g_scratch,
     $     1.0d0,g_dens)
c
      mp2x = -ga_ddot(g_dens, g_x) ! Electrons are negative
      mp2y = -ga_ddot(g_dens, g_y)
      mp2z = -ga_ddot(g_dens, g_z)
c
      if(.not.ga_destroy(g_x))call errquit('mp2_grad: ga destroy?',0,
     &       GA_ERR)
      if(.not.ga_destroy(g_y))call errquit('mp2_grad: ga destroy?',0,
     &       GA_ERR)
      if(.not.ga_destroy(g_z))call errquit('mp2_grad: ga destroy?',0,
     &       GA_ERR)
      if(.not.ga_destroy(g_dens))call errquit('mp2_grad: ga destroy?',0,
     &       GA_ERR)
      if(.not.ga_destroy(g_scratch))
     $     call errquit('mp2_grad: ga destroy?',0, GA_ERR)
c
      if (.not. geom_nuc_dipole(geom,ndip)) 
     $     call errquit('mp2_dip: failed to get nuclear dipole',0,
     &       GEOM_ERR)
c
      if (ga_nodeid().eq.0 .and. util_print('dipole',print_low)) then
         write(6,1) 
         write(6,2) 'X', ndip(1), scfx+ndip(1),mp2x+scfx+ndip(1)
         write(6,2) 'Y', ndip(2), scfy+ndip(2),mp2y+scfy+ndip(2)
         write(6,2) 'Z', ndip(3), scfz+ndip(3),mp2z+scfz+ndip(3)
         write(6,*)
 1       format(
     $        6x,'          Dipole Moments (a.u.)'/
     $        6x,'          --------------'/
     $        6x,'   Nuclear       SCF       SCF+MP2'/
     $        6x,' ----------- ----------- -----------')
 2       format(1x,a1,3x,3f12.6)
         call util_flush(6)
      end if
c
c     Store the MP2 dipole in the database
c
      ndip(1) = mp2x+scfx+ndip(1)
      ndip(2) = mp2y+scfy+ndip(2)
      ndip(3) = mp2z+scfz+ndip(3)
c
      if (.not. rtdb_put(rtdb,'mp2:dipole',mt_dbl,3,ndip))
     $     call errquit('mp2_dipole: database is corrupt',0,
     &       RTDB_ERR)
      call ecce_print1('total dipole', mt_dbl, ndip, 3)
c
      end


