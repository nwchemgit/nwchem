      subroutine tce_eomccsd_1prdm(d_hh,d_pp,d_t1,d_t2,d_x0,d_x1,d_x2,
     1                         d_y1,d_y2,k_hh_offset,k_pp_offset,
     1                         k_t1_offset,k_t2_offset,
     1                         k_x0_offset,k_x1_offset,k_x2_offset,
     3                         k_y1_offset,k_y2_offset)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "rtdb.fh"
      integer d_hh,k_hh_offset
      integer d_pp,k_pp_offset
      integer d_t1,k_t1_offset
      integer d_t2,k_t2_offset
      integer d_x0,k_x0_offset
      integer d_x1,k_x1_offset
      integer d_x2,k_x2_offset
      integer d_y1,k_y1_offset
      integer d_y2,k_y2_offset
c
      character*255 filename
c     ====
c     1PDM
c     ====
      integer d_hp,l_hp_offset,k_hp_offset,size_hp
      integer d_ph,l_ph_offset,k_ph_offset,size_ph
      integer d_1pdm,l_1pdm_offset,k_1pdm_offset,size_1pdm
c     ====================
c     atomic orbital tiles
c     ====================
c      integer atpart
c      integer nalength(200)
c      integer a2length(200)
c
      double precision cpu,wall,r1,residual
      integer dummy, iter
c     ==============
c     GA for 1PDM_AO
c     ==============
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
      integer  g_1pdm
c
      logical  nodezero
c
      integer g1b, g1
      integer dima, l_a, k_a
c
c     parallel
      integer next
      integer count
      integer nxtask
      external nxtask
      integer nprocs
c
      nodezero=(ga_nodeid().eq.0)
c
c     ============================================================
c     !!! Attention: frozen orbital option was not considered !!!
c     ============================================================
c
c     ==========================
c     1PDM Hole-Particle Block
c     ==========================
      call tce_filename('hp',filename)
      call tce_dens_hp_offset(l_hp_offset,k_hp_offset,size_hp)
      call createfile(filename,d_hp,size_hp)
      call eomccsd_1pdm_hp_mo(d_hp,d_x0,d_x1,d_y1,d_y2,
     1      k_hp_offset,k_x0_offset,k_x1_offset,
     2      k_y1_offset,k_y2_offset)
      call reconcilefile(d_hp,size_hp)
c     ===========================
c     1PDM Particle-Hole Block
c     ==========================
      call tce_filename('ph',filename)
      call tce_dens_ph_offset(l_ph_offset,k_ph_offset,size_ph)
      call createfile(filename,d_ph,size_ph)
      call eomccsd_1pdm_ph_mo(d_ph,d_t1,d_t2,
     1              d_x0,d_x1,d_x2,d_y1,d_y2,
     1     k_ph_offset,k_t1_offset,k_t2_offset,
     1     k_x0_offset,k_x1_offset,k_x2_offset,
     1     k_y1_offset,k_y2_offset)
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
     1      call errquit('tce_eomccsd_1prdm: ma problem',2,ma_err)
            call get_hash_block(d_1pdm,dbl_mb(k_a),dima,
     1           int_mb(k_1pdm_offset),(g1b-1)+(g1b-1)*(noab+nvab))
            do g1=1,int_mb(k_range+g1b-1)
               dbl_mb(k_a+(g1-1)+(g1-1)*int_mb(k_range+g1b-1))=
     1         dbl_mb(k_a+(g1-1)+(g1-1)*int_mb(k_range+g1b-1))+1.d0
            enddo
            call put_hash_block(d_1pdm,dbl_mb(k_a),dima,
     1           int_mb(k_1pdm_offset),(g1b-1)+(g1b-1)*(noab+nvab))
            if (.not.ma_pop_stack(l_a))
     1      call errquit('tce_eomccsd_1prdm: ma problem',3,ma_err)
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
c      call ga_zero(g_1pdm)
      call ao_tiles(atpart,nalength,a2length,30)
      call btrans1(d_1pdm,k_1pdm_offset,g_1pdm,atpart,nalength)
c      call ga_symmetrize(g_1pdm)
c     ====================================
c     DUMP it to a file
c     ====================================
c      call ga_zero(g_1pdm)
c      call ga_print(g_1pdm)
      call ao_1prdm_write(nbf,g_1pdm)
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
