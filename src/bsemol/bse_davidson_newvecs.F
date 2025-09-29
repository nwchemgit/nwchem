      subroutine bse_davidson_newvecs(ntrials,ndim,npoles,g_x,g_y,
     $               g_ym,g_s,g_ks,g_mks,g_w,g_kw,g_mkw,lvkv,vkv,map,
     $               nblocks)

      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"

      integer ntrials,ndim,npoles,nblocks
      integer g_x,g_y,g_ym,g_s,g_ks,g_mks
      integer g_w,g_kw,g_mkw,lvkv,vkv
      integer map(nblocks),grid(2)

      logical stat

      grid = (/ga_nnodes(),1/)

      ! Save working vectors
      call ga_copy(g_s,g_x)
      call ga_copy(g_ks,g_y)
      call ga_copy(g_mks,g_ym)

      ! Destroy working vectors
      stat = ga_destroy(g_s).and.ga_destroy(g_ks).and.ga_destroy(g_mks)
      if(.not.stat)
     $    call errquit('bse_davidson_newvecs: destroy failed',1,GA_ERR)

      ! Create new working vectors
      g_s = ga_create_handle()
      call ga_set_data(g_s,2,(/npoles,ntrials+ndim/),mt_dbl)
      call ga_set_array_name(g_s,'s')
      call ga_set_tiled_irreg_proc_grid(g_s,map,(/nblocks,1/),grid)
      stat = ga_allocate(g_s)

      stat = stat.and.ga_duplicate(g_s,g_ks,'ks').and.
     $                ga_duplicate(g_s,g_mks,'mks')
      if(.not.stat)
     $  call errquit('bse_davidson_newvecs: create failed',1,GA_ERR)


      ! Copy back old working vectors plus new vectors
      call ga_copy_patch('n',g_x,1,npoles,1,ntrials,
     $                       g_s,1,npoles,1,ntrials)
      call ga_copy_patch('n',g_w,1,npoles,1,ndim,
     $                       g_s,1,npoles,ntrials+1,ntrials+ndim)

      call ga_copy_patch('n',g_y,1,npoles,1,ntrials,
     $                       g_ks,1,npoles,1,ntrials)
      call ga_copy_patch('n',g_kw,1,npoles,1,ndim,
     $                       g_ks,1,npoles,ntrials+1,ntrials+ndim)

      call ga_copy_patch('n',g_ym,1,npoles,1,ntrials,
     $                       g_mks,1,npoles,1,ntrials)
      call ga_copy_patch('n',g_mkw,1,npoles,1,ndim,
     $                       g_mks,1,npoles,ntrials+1,ntrials+ndim)

      ! Update ntrials
      ntrials = ntrials + ndim

      ! Destroy old solutions
      stat = ga_destroy(g_x).and.ga_destroy(g_y).and.ga_destroy(g_ym)
      if(.not.stat)
     $    call errquit('bse_davidson_newvecs: destroy failed',2,GA_ERR)

      ! Create new solutions
      stat = ga_duplicate(g_s,g_x,'x').and.
     $       ga_duplicate(g_s,g_y,'y').and.
     $       ga_duplicate(g_s,g_ym,'ym')
      if(.not.stat)
     $  call errquit('bse_davidson_newvecs: create failed',2,GA_ERR)

      ! Reallocate VKV
      if(.not.ma_chop_stack(lvkv))
     $  call errquit('bse_davidson_restart: alloc failed',3,MA_ERR)
      if(.not.ma_push_get(mt_dbl,ntrials*ntrials,'vkv',lvkv,vkv))
     $  call errquit('bse_davidson_restart: dealloc failed',3,MA_ERR)

      end subroutine
