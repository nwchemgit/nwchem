      subroutine mp2_make_dip_mom(g_p_tot_a,g_p_tot_b,
     $     basis,rtdb,nbf,g_vecs_a,g_vecs_b)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer basis,rtdb,nbf
      integer g_p_tot_a,g_p_tot_b
      integer g_x_a,g_y_a,g_z_a
      integer g_x_b,g_y_b,g_z_b
      integer g_vecs_a,g_vecs_b
      integer g_scratch
      double precision dipx,dipx_a,dipx_b
      double precision dipy,dipy_a,dipy_b
      double precision dipz,dipz_a,dipz_b
c
      if(.not.ga_create(mt_dbl,nbf,nbf,' ao integs ',
     $     1,1,g_x_a)) call errquit('ga',0)
      if(.not.ga_create(mt_dbl,nbf,nbf,' ao integs ',
     $     1,1,g_x_b)) call errquit('ga',0)
      if(.not.ga_create(mt_dbl,nbf,nbf,' ao integs ',
     $     1,1,g_y_a)) call errquit('ga',0)
      if(.not.ga_create(mt_dbl,nbf,nbf,' ao integs ',
     $     1,1,g_y_b)) call errquit('ga',0)
      if(.not.ga_create(mt_dbl,nbf,nbf,' ao integs ',
     $     1,1,g_z_a)) call errquit('ga',0)
      if(.not.ga_create(mt_dbl,nbf,nbf,' ao integs ',
     $     1,1,g_z_b)) call errquit('ga',0)
      if(.not.ga_create(mt_dbl,nbf,nbf,'uhf_s_e: alpha MOs',
     $     1,1,g_scratch))call errquit('uhf_s_e: MOs',0)
      call int_init(rtdb,1,basis)
      call int_dip_ga(basis,basis,g_x_a,g_y_a,g_z_a)
      call int_dip_ga(basis,basis,g_x_b,g_y_b,g_z_b)
      call two_index_transf(g_x_a,g_vecs_a,
     $     g_vecs_a,g_scratch,g_x_a)
      call two_index_transf(g_x_b,g_vecs_b,
     $     g_vecs_b,g_scratch,g_x_b)
      call two_index_transf(g_y_a,g_vecs_a,
     $     g_vecs_a,g_scratch,g_y_a)
      call two_index_transf(g_y_b,g_vecs_b,
     $     g_vecs_b,g_scratch,g_y_b)
      call two_index_transf(g_z_a,g_vecs_a,
     $     g_vecs_a,g_scratch,g_z_a)
      call two_index_transf(g_z_b,g_vecs_b,
     $     g_vecs_b,g_scratch,g_z_b)
      dipx_a=ga_ddot(g_x_a,g_p_tot_a)
      dipx_b=ga_ddot(g_x_b,g_p_tot_b)
      dipx=dipx_a+dipx_b
      dipy_a=ga_ddot(g_y_a,g_p_tot_a)
      dipy_b=ga_ddot(g_y_b,g_p_tot_b)
      dipy=dipy_a+dipy_b
      dipz_a=ga_ddot(g_z_a,g_p_tot_a)
      dipz_b=ga_ddot(g_z_b,g_p_tot_b)
      dipz=dipz_a+dipz_b
      write(6,*) 'dipx_a = ',dipx_a
      write(6,*) 'dipx_b = ',dipx_b
      write(6,*) 'dipx = ',dipx
      write(6,*) 'dipy_a = ',dipy_a
      write(6,*) 'dipy_b = ',dipy_b
      write(6,*) 'dipy = ',dipy
      write(6,*) 'dipz_a = ',dipz_a
      write(6,*) 'dipz_b = ',dipz_b
      write(6,*) 'dipz = ',dipz
      if(.not.ga_destroy(g_x_a))call errquit(
     $     'could not destroy g_x_a',1)
      if(.not.ga_destroy(g_x_b))call errquit(
     $     'could not destroy g_x_b',1)
      if(.not.ga_destroy(g_y_a))call errquit(
     $     'could not destroy g_y_a',1)
      if(.not.ga_destroy(g_y_b))call errquit(
     $     'could not destroy g_y_b',1)
      if(.not.ga_destroy(g_z_a))call errquit(
     $     'could not destroy g_z_a',1)
      if(.not.ga_destroy(g_z_b))call errquit(
     $     'could not destroy g_z_b',1)
      if(.not.ga_destroy(g_scratch))call errquit(
     $     'could not destroy g_scratch',1)
c      call int_terminate()
      end


