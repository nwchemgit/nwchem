      subroutine mp2_make_tot_dens(g_p_ij_a,g_p_ij_b,g_p_ab_a,
     $     g_p_ab_b,g_rhs,g_p_tot_a,g_p_tot_b,noa,nva,nob,nvb,nbf,
     $     g_vecs_a,g_vecs_b,geom,basis,rtdb)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer nbf,noa,nva,nob,nvb
      integer g_p_ij_a,g_p_ij_b
      integer g_p_ab_a,g_p_ab_b
      integer g_rhs
      integer g_p_tot_a,g_p_tot_b
      integer g_scratch,g_vecs_a,g_vecs_b,g_a_dens,g_b_dens
      integer ga_create_atom_blocked
      integer i,a,indx
      double precision temp
      integer geom,basis,rtdb
      integer g_over,g_vecs_inv
      double precision one,zero
      one=1.0d0
      zero=0.0d0
      call int_init(rtdb,1,basis)
      if(.not.ga_create(mt_dbl,nbf,nbf,'scratch',1,1,g_scratch))
     $     call errquit('scratch',1)
      call ga_dscal(g_rhs,2.0d0)
      call ga_zero(g_p_tot_a)
      call ga_zero(g_p_tot_b)
      call ga_copy_patch('n',g_rhs,1,nva*noa,1,1,g_p_tot_a,
     $     noa+1,nbf,1,noa)
c$$$      call ga_copy_patch('t',g_rhs,1,nva*noa,1,1,g_p_tot_a,
c$$$     $     1,noa,noa+1,nbf)
      do i=1,noa
         do a=1,nva
            indx=(i-1)*nva+a
            call ga_get(g_rhs,indx,indx,1,1,temp,1)
            call ga_put(g_p_tot_a,i,i,a+noa,a+noa,temp,1)
         enddo
      enddo
      call ga_copy_patch('n',g_p_ij_a,1,noa,1,noa,g_p_tot_a,
     $     1,noa,1,noa)
      call ga_copy_patch('n',g_p_ab_a,1,nva,1,nva,g_p_tot_a,
     $     noa+1,nbf,noa+1,nbf)
      call ga_copy_patch('n',g_rhs,(nva*noa)+1,(nva*noa)+(nvb*nob)
     $     ,1,1,g_p_tot_b,nob+1,nbf,1,nob)
      do i=1,nob
         do a=1,nvb
            indx=(i-1)*nvb+a+(nva*noa)
            call ga_get(g_rhs,indx,indx,1,1,temp,1)
            call ga_put(g_p_tot_b,i,i,a+nob,a+nob,temp,1)
         enddo
      enddo
c$$$      call ga_copy_patch('t',g_rhs,(nva*noa)+1,(nva*noa)+(nvb*nob)
c$$$     $     ,1,1,g_p_tot_b,1,nob,nob+1,nbf)
      call ga_copy_patch('n',g_p_ij_b,1,nob,1,nob,g_p_tot_b,
     $     1,nob,1,nob)
      call ga_copy_patch('n',g_p_ab_b,1,nvb,1,nvb,g_p_tot_b,
     $     nob+1,nbf,nob+1,nbf)
      if(.not.ga_destroy(g_scratch))call errquit(
     $     'could not destroy g_scratch',1)
c...make the overlap matrix
c$$$      g_scf_dens=ga_create_atom_blocked(geom,basis,'uhf overlap')
c$$$      call ga_zero(g_scf_dens)
c$$$      call ga_dgemm('n', 't', nbf, nbf, noa, one, g_vecs_a,
c$$$     $     g_vecs_a, zero, g_scf_dens)
c... add in the scf density if desired:
c$$$      do i=1,noa
c$$$         call ga_acc(g_p_tot_a,i,i,i,i,1.0d0,1,1.0d0)
c$$$      enddo
c$$$      do i=1,nob
c$$$         call ga_acc(g_p_tot_b,i,i,i,i,1.0d0,1,1.0d0)
c$$$      enddo
c$$$c... transform the total density to the AO basis
c$$$      g_over = ga_create_atom_blocked(geom, basis, 'over')
c$$$      call two_index_transf(g_p_tot_a,g_vecs_a,g_vecs_a,g_over,
c$$$     $     g_p_tot_a)
c$$$      write(6,*) 'total density in AO basis'
c$$$      call ga_print(g_p_tot_a)
c$$$      call int_1e_ga(basis,basis,g_over,'overlap',.false.)
c$$$      write(6,*) 'overlap in ao basis'
c$$$      call ga_print(g_over)
c$$$      if(.not.ga_create(mt_dbl,nbf,nbf,'inv_mo',1,1,g_vecs_inv))
c$$$     $     call errquit('inv_mo',1)
c$$$      call ga_dgemm('n','n',nbf,nbf,nbf,1.0d0,g_over,g_vecs_a,
c$$$     $     0.0d0,g_vecs_inv)
c$$$      write(6,*) 'g_vecs_inv'
c$$$      call ga_print(g_vecs_inv)
c$$$      call int_terminate()
      end


