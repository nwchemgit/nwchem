      subroutine tce_ccsd_1prdm(d_t1,d_t2,d_lambda1,d_lambda2,
     1                    d_hh,d_pp,
     1                    k_t1_offset,k_t2_offset,
     1                    k_l1_offset,k_l2_offset,
     1                    k_hh_offset,k_pp_offset)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "rtdb.fh"
c
      integer d_t1,k_t1_offset
      integer d_t2,k_t2_offset
      integer d_lambda1,k_l1_offset
      integer d_lambda2,k_l2_offset
      integer d_hh,k_hh_offset
      integer d_pp,k_pp_offset
c
c     integer k_force
c
      character*255 filename
c
c     ====
c     1PDM
c     ====
      integer d_ph,l_ph_offset,k_ph_offset,size_ph
      integer d_hp,l_hp_offset,k_hp_offset,size_hp
      integer d_1pdm,l_1pdm_offset,k_1pdm_offset,size_1pdm
c     ===========================
c     atomic orbital tiles
c     ===========================
      integer atpart2
      integer nalength2(200)
      integer a2length2(200)
c
      double precision cpu,wall,r1,residual
      integer dummy, iter
c     ===========================
c     GA for 1PDM_AO and 1PWDM_AO
c     ===========================
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      integer  g_1pdm
c
      integer dima,l_a,k_a
c
      logical nodezero
c
      integer l_test,k_test,l_r2,k_r2,size,l_r1,k_r1
      integer g1,g2,g3,g4,g1b,g2b,g3b,g4b
      integer ig1,ig2,ig3,ig4,igg1,igg2,igg3,igg4
      integer k,l,ind1,ind2
      integer dim_2pdm, k_2pdm, l_2pdm
      integer dim_1pdm, k_1pdm, l_1pdm
      integer dim_1pwdm,k_1pwdm,l_1pwdm
      integer k_2eint,l_2eint,k_1eint,l_1eint
c
      double precision res
      integer l_zvec_x,k_zvec_x
      integer m,e,a,b,q,r,s,n,f
c
c     debug
      logical zvec_debug, wdm_debug
c
      logical ao_1prdm_write
      external ao_1prdm_write
c     parallel
      integer next
      integer count
      integer nxtask
      external nxtask
      integer nprocs
c
      nodezero=(ga_nodeid().eq.0)
c
c     =================================================================
c     !!! Attention: frozen orbital option was not considered !!!
c     =================================================================
c
c     ========================
c     1PDM Hole-Particle Block
c     ========================
      call tce_filename('hp',filename)
      call tce_dens_hp_offset(l_hp_offset,k_hp_offset,size_hp)
      call createfile(filename,d_hp,size_hp)
      call ccsd_1pdm_hp_mo(d_hp,d_lambda1,k_hp_offset,k_l1_offset)
      call reconcilefile(d_hp,size_hp)
c     ========================
c     1PDM Particle-Hole Block
c     ========================
      call tce_filename('ph',filename)
      call tce_dens_ph_offset(l_ph_offset,k_ph_offset,size_ph)
      call createfile(filename,d_ph,size_ph)
      call ccsd_1pdm_ph_mo(d_ph,d_t1,d_t2,d_lambda1,d_lambda2,
     1     k_ph_offset,k_t1_offset,k_t2_offset,k_l1_offset,k_l2_offset)
      call reconcilefile(d_ph,size_ph)
c     =====================================================================================================
c     put 1PDM Hole-Hole, Hole-Particle, Particle-Hole, and Particle-Particle Block four pieces in one file
c     =====================================================================================================
      call tce_filename('1pdm',filename)
      call tce_1pdm_offset(l_1pdm_offset,k_1pdm_offset,size_1pdm)
      call createfile(filename,d_1pdm,size_1pdm)
c
      call put_1pdm_hh(d_1pdm,k_1pdm_offset,d_hh,k_hh_offset)
      call put_1pdm_hp(d_1pdm,k_1pdm_offset,d_hp,k_hp_offset)
      call put_1pdm_ph_2(d_1pdm,k_1pdm_offset,d_ph,k_ph_offset)
      call put_1pdm_pp(d_1pdm,k_1pdm_offset,d_pp,k_pp_offset)
      call reconcilefile(d_1pdm,size_1pdm)
c    =================================
c    HF REFERENCE CONTRIBUTION TO 1PDM
c    =================================
      nprocs=ga_nnodes( )
      count=0
      next=nxtask(nprocs,1)
c
      do g1b=1,noab
         if(count.eq.next) then
c
         if ((.not.restricted).or.(int_mb(k_spin+g1b-1).ne.2)) then
            dima=int_mb(k_range+g1b-1)*int_mb(k_range+g1b-1)
            if (.not.ma_push_get(mt_dbl,dima,'noname',l_a,k_a))
     1      call errquit('ccsd_gradients: ma problem',2,ma_err)
            call get_hash_block(d_1pdm,dbl_mb(k_a),dima,
     1           int_mb(k_1pdm_offset),(g1b-1)+(g1b-1)*(noab+nvab))
            do g1=1,int_mb(k_range+g1b-1)
               dbl_mb(k_a+(g1-1)+(g1-1)*int_mb(k_range+g1b-1))=
     1         dbl_mb(k_a+(g1-1)+(g1-1)*int_mb(k_range+g1b-1))+1.d0
            enddo
            call put_hash_block(d_1pdm,dbl_mb(k_a),dima,
     1           int_mb(k_1pdm_offset),(g1b-1)+(g1b-1)*(noab+nvab))
            if (.not.ma_pop_stack(l_a))
     1      call errquit('ccsd_gradients: ma problem',3,ma_err)
         endif
         next=nxtask(nprocs,1)
         endif
         count=count+1
      enddo
      next=nxtask(-nprocs,1)
      call ga_sync( )
c     ====================================
c     Do the back transformation
c     ===================================
      g_1pdm = ga_create_atom_blocked(geom, ao_bas_han,
     $         'density matrix')
      call ao_tiles(atpart2,nalength2,a2length2,30)
      call btrans1(d_1pdm,k_1pdm_offset,g_1pdm,atpart2,nalength2)
      call ga_symmetrize(g_1pdm)
c     ====================================
c     DUMP it to a file
c     ====================================
c      call ga_print(g_1pdm)
       if(.not.ao_1prdm_write(nbf,g_1pdm))
     1   call errquit('tce_ccsd_1prdm: ao_1prdm_write failed',0,0)
c      call ga_zero(g_1pdm)
c      call ga_print(g_1pdm)
c      call ao_1prdm_read(nbf,g_1pdm)
c      call ga_print(g_1pdm)
c     =============================
c     clean up the files and arrays
c     =============================
      if (.not. ga_destroy(g_1pdm)) call errquit
     1   ('tce_ccsd_1prdm: error destroying density', 1, GA_ERR)
c
      call deletefile(d_1pdm)
      if(.not.ma_pop_stack(l_1pdm_offset))
     1  call errquit('tce_ccsd_1prdm: ma problem',17,ma_err)
c
      call deletefile(d_ph)
      if(.not.ma_pop_stack(l_ph_offset))
     1  call errquit('tce_ccsd_1prdm: ma problem',17,ma_err)
c
      call deletefile(d_hp)
      if(.not.ma_pop_stack(l_hp_offset))
     1  call errquit('tce_ccsd_1prdm: ma problem',17,ma_err)
c
      end

c $Id$
