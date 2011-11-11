      subroutine mp2_make_nos(rtdb, basis, nbf, nmo, noa, nob,
     $     g_p_tot_a, g_p_tot_b, g_vecs_a, g_vecs_b, scftype, nfrozen,
     $     title)
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
#include "inp.fh"
      integer basis,rtdb,nbf,nmo
      integer g_p_tot_a,g_p_tot_b
      integer g_vecs_a,g_vecs_b
      integer noa, nob
      character*(*) scftype, title
      integer nfrozen
      character*255 filename
      character*8 string
c
c     Make the MP2 natural orbitals
c
      integer l_eval, k_eval, g_dens, g_nos(2), 
     $     g_work, g_vecs(2), g_p(2), no(2)
      integer iset, i, nsets, ioff, nmos(2)
      logical oprint
      logical movecs_write
      external movecs_write
c
      call util_file_name('mp2nos',.false.,.false.,filename)
c
      oprint = util_print('no occupation',print_default)
c
      g_vecs(1) = g_vecs_a
      g_vecs(2) = g_vecs_b
      g_p(1)    = g_p_tot_a
      g_p(2)    = g_p_tot_b
      no(1)     = noa
      no(2)     = nob
      nmos(1)   = nmo
      nmos(2) = nmo
c
      if(.not.ma_push_get(mt_dbl, 2*nbf,'mp2 nos', l_eval,k_eval))      
     $     call errquit('mp2_nos: allocating eval?', nbf, MA_ERR)
      call dfill(2*nbf, 0.0d0, dbl_mb(k_eval), 1)
c
*ga:1:0
      if(.not.ga_create(mt_dbl, nmo, nmo, 'mp2 nos a', nmo, 0,
     $     g_dens)) call errquit('mp2_nos: ga_create a failed',0,
     &       GA_ERR)
*ga:1:0
      if(.not.ga_create(mt_dbl, nmo, nmo, 'mp2 nos b', nmo, 0,
     $     g_work)) call errquit('mp2_nos: ga_create b failed',0,
     &       GA_ERR)
*ga:1:0
      if(.not.ga_create(mt_dbl, nbf, nmo, 'mp2 nos c', nbf, 0,
     $     g_nos(1))) call errquit('mp2_nos: ga_create c failed',0,
     &       GA_ERR)
      if (scftype .eq. 'UHF') then
*ga:1:0
         if(.not.ga_create(mt_dbl, nbf, nmo, 'mp2 nos d', nbf, 0,
     $        g_nos(2))) call errquit('mp2_nos: ga_create d failed',0,
     &       GA_ERR)
      endif
c
      if (scftype .eq. 'RHF') then
         nsets = 1
      else
         nsets = 2
      endif
c
      do iset = 1, nsets
         ioff = (iset-1)*nbf
c
c     Make the total density
c
         call ga_copy(g_p(iset), g_dens)
         if (ga_nodeid() .eq. 0) then
c
c     Put distinct -ve values into the frozen core so that they
c     are not mixed and the order is preserved ... undo this later.
c     Also add in the SCF occupation.
c
            do i = 1, nfrozen
               call ga_put(g_dens, i,i,i,i, dble(2*nfrozen-i), 1)
            enddo
            do i = nfrozen+1,no(iset)
               call ga_acc(g_dens, i,i,i,i, 1.0d0, 1, 1.0d0)
            enddo
         endif
         call ga_sync
c
c     This is dumb.  In order to get the NOs in order most strongly
c     to most weakly occupied, negate the density.
c
         call ga_scale(g_dens,-1.0d0)
c
#ifdef PARALLEL_DIAG
#ifdef SCALAPACK
         call ga_pdsyev  (g_dens, g_work, dbl_mb(k_eval+ioff), 0)
#else
         call ga_diag_std(g_dens, g_work, dbl_mb(k_eval+ioff))
#endif
#else
         call ga_diag_std_seq(g_dens, g_work, dbl_mb(k_eval+ioff))
#endif
c
c     Undo sign munging
c
         do i = 1, nmo
            dbl_mb(k_eval+i-1+ioff) = -dbl_mb(k_eval+i-1+ioff)
         enddo
c
c     Undo munging to stop core from mixing
c
         do i = 1, nfrozen      
            dbl_mb(k_eval+i-1+ioff) = 1.0d0
         enddo
c         
         if (scftype .eq. 'RHF') then
            do i = 1, nmo
               dbl_mb(k_eval+i-1+ioff) = dbl_mb(k_eval+i-1+ioff)*2.0d0
            enddo
         endif
c
         call ga_dgemm('n', 'n', nbf, nmo, nmo, 1.0d0,
     $        g_vecs(iset), g_work, 0.0d0, g_nos(iset))
c
         if (oprint) then
            if (ga_nodeid() .eq. 0) then
               if (scftype .eq. 'RHF') then
                  string = 'RHF'
               else
                  if (iset .eq. 1) then
                     string = 'alpha'
                  else
                     string = 'beta'
                  endif
               endif
               write(6,1) string
 1             format(/' MP2 ',a,' natural orbital occupation',
     $              ' numbers'/)
               call output(dbl_mb(k_eval+ioff),1,nmo,1,1,nmo,1,1)
               call util_flush(6)
            endif
            if (util_print('natural orbitals',print_debug)) 
     $           call ga_print(g_nos(iset))
         endif
c
      enddo
c
      call ga_sync
c
      call util_file_name('mp2nos',.false.,.false.,filename)
c
      if (.not. movecs_write(rtdb, basis, filename, 
     $     'mp2', title,
     $     nbf, nsets, nmos, dbl_mb(k_eval), nbf, 
     $     dbl_mb(k_eval), nbf, g_nos))
     $     call errquit('mp2_nos: failed writing nos',0, DISK_ERR)
c
      if (oprint .and. ga_nodeid().eq.0) then
         write(6,*)
         write(6,*) ' MP2 natural orbitals output to ',
     $        filename(1:inp_strlen(filename))
         write(6,*)
         call util_flush(6)
      endif
c
      if (.not. ma_pop_stack(l_eval)) call errquit('mp2_nos: ma?',0,
     &       MA_ERR)
      if (.not. ga_destroy(g_dens)) call errquit('mp2_nos: ga?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_work)) call errquit('mp2_nos: ga?',0,
     &       GA_ERR)
      if (.not. ga_destroy(g_nos(1))) call errquit('mp2_nos: ga?',0,
     &       GA_ERR)
      if (scftype .eq. 'UHF') then
         if (.not. ga_destroy(g_nos(2))) call errquit('mp2_nos: ga?',0,
     &       GA_ERR)
      endif
c
      end

