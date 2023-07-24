      subroutine bse_davidson_restart(ntrials,nroots,npoles,g_x,g_y,
     $                                g_ym,g_s,g_ks,g_mks,lvkv,vkv)

      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"

      integer ntrials,nroots,npoles
      integer g_x,g_y,g_ym,g_s,g_ks,g_mks,lvkv,vkv

      logical stat

      ! Reset the number of trials
      ntrials = nroots

      ! Copy current solutions
      call ga_copy_patch('n',g_x,1,npoles,1,ntrials,
     $                       g_s,1,npoles,1,ntrials)
      call ga_copy_patch('n',g_y,1,npoles,1,ntrials,
     $                       g_ks,1,npoles,1,ntrials)
      call ga_copy_patch('n',g_ym,1,npoles,1,ntrials,
     $                       g_mks,1,npoles,1,ntrials)

      ! Destroy current solutions
      stat = ga_destroy(g_x).and.ga_destroy(g_y).and.ga_destroy(g_ym)
      if (.not.stat)
     $   call errquit('bse_davidson_restar: failed destroy',1,GA_ERR)

      ! Re-create new solutions
      stat = ga_create(mt_dbl,npoles,ntrials,'x',0,ntrials,g_x) .and.
     $       ga_create(mt_dbl,npoles,ntrials,'y',0,ntrials,g_y) .and. 
     $       ga_create(mt_dbl,npoles,ntrials,'ym',0,ntrials,g_ym)
      if (.not.stat)
     $   call errquit('bse_davidson_restar: failed create',1,GA_ERR)

      ! Copy back solutions
      call ga_copy_patch('n',g_s,1,npoles,1,ntrials,
     $                       g_x,1,npoles,1,ntrials)
      call ga_copy_patch('n',g_ks,1,npoles,1,ntrials,
     $                       g_y,1,npoles,1,ntrials)
      call ga_copy_patch('n',g_mks,1,npoles,1,ntrials,
     $                       g_ym,1,npoles,1,ntrials)

      ! Destroy working vectors
      stat = ga_destroy(g_s).and.ga_destroy(g_ks).and.ga_destroy(g_mks)
      if (.not.stat)
     $   call errquit('bse_davidson_restar: failed destroy',2,GA_ERR)

      ! Re-create new solutions
      stat = ga_create(mt_dbl,npoles,ntrials,'s',0,ntrials,g_s) .and.
     $       ga_create(mt_dbl,npoles,ntrials,'ks',0,ntrials,g_ks) .and. 
     $       ga_create(mt_dbl,npoles,ntrials,'mks',0,ntrials,g_mks)
      if (.not.stat)
     $   call errquit('bse_davidson_restar: failed create',2,GA_ERR)

      ! Copy current solutions
      call ga_copy(g_x,g_s)
      call ga_copy(g_y,g_ks)
      call ga_copy(g_ym,g_mks)

      ! Reallocate VKV
      if(.not.ma_chop_stack(lvkv))
     $  call errquit('bse_davidson_restart: alloc failed',3,MA_ERR)
      if(.not.ma_push_get(mt_dbl,ntrials*ntrials,'vkv',lvkv,vkv))
     $  call errquit('bse_davidson_restart: dealloc failed',3,MA_ERR)

      end subroutine
