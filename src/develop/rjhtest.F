      logical function rjhtest(rtdb)
*     
*     $Id: rjhtest.F,v 1.15 1997-12-22 20:10:58 d3g681 Exp $
*
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
      integer rtdb
      integer basis, geom, nbf, g_q, g_v, nsh, i, 
     $ k_i, k_j, k_k, k_l, l_i, l_j, l_k, l_l,
     $ maxg2, maxs2, k_buf, l_buf, k_scr, l_scr, l_eval, k_eval
c
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('scf_init: geom_create?', 0)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('scf_init: no geometry ', 0)
c     
c     load the basis set and get info about it
c
      if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('scf_init: bas_create?', 0)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis')) 
     $        call errquit('scf_init: no ao basis set', 0)
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('scf_init: basis info',0)
      if ( .not. bas_numcont(basis, nsh) ) call errquit(
     $     'ao_fock_2e: problem with call to bas_numcont', basis)
      if (.not. ga_create(mt_dbl,nbf*nbf,nbf*nbf,'Q',1,1,g_q))
     $   call errquit('ga_create failed q',0)
      if (.not. ga_create(mt_dbl,nbf*nbf,nbf*nbf,'Q',1,1,g_v))
     $   call errquit('ga_create failed v',0)
c
      call int_init(rtdb, 1, basis)
      call int_mem_2e4c(maxg2,maxs2)
      if (.not. ma_push_get(mt_dbl,nbf*nbf,'eval',l_eval, k_eval))
     $   call errquit('ma eval',nbf*nbf)
      if (.not. ma_push_get(mt_dbl,maxs2,'scr',l_scr, k_scr))
     $   call errquit('ma scr',maxg2)
      if (.not. ma_push_get(mt_dbl,maxg2,'buf',l_buf, k_buf))
     $   call errquit('ma buf',maxg2)
      if (.not. ma_push_get(mt_int,maxg2,'i',l_i, k_i))
     $   call errquit('ma ibuf',maxg2)
      if (.not. ma_push_get(mt_int,maxg2,'j',l_j, k_j))
     $   call errquit('ma jbuf',maxg2)
      if (.not. ma_push_get(mt_int,maxg2,'k',l_k, k_k))
     $   call errquit('ma kbuf',maxg2)
      if (.not. ma_push_get(mt_int,maxg2,'l',l_l, k_l))
     $   call errquit('ma lbuf',maxg2)
c
      call rjh_fill_buf(basis, dbl_mb(k_buf), dbl_mb(k_scr),
     $   int_mb(k_i), int_mb(k_j), int_mb(k_k), int_mb(k_l),
     $   maxg2, maxs2, nbf, nsh, .false., g_q)
      call ga_diag_std(g_q, g_v, dbl_mb(k_eval))
c
      if (ga_nodeid().eq.0) then
        write(6,*) ' Coulomb '
        write(6,1) (i,dbl_mb(k_eval+i-1),i=1,nbf*nbf)
  1     format(i6,2x,1p,d16.8)
        call util_flush(6)
      endif
      call ga_sync
c
      call rjh_fill_buf(basis, dbl_mb(k_buf), dbl_mb(k_scr),
     $   int_mb(k_i), int_mb(k_j), int_mb(k_k), int_mb(k_l),
     $   maxg2, maxs2, nbf, nsh, .true., g_q)
      call ga_diag_std(g_q, g_v, dbl_mb(k_eval))
c
      if (ga_nodeid().eq.0) then
        write(6,*) ' Exchange '
        write(6,1) (i,dbl_mb(k_eval+i-1),i=1,nbf*nbf)
        call util_flush(6)
      endif
      call ga_sync
c
      rjhtest = .true.
c
      end
      subroutine rjh_fill_buf(basis, buf, scr, ilab, jlab, klab,
     $  llab, maxg2, maxs2, nbf, nsh, oexch, g_q)
      implicit none
#include "global.fh"
      integer basis, nbf, nsh, maxg2, maxs2, g_q
      logical oexch
      double precision buf(maxg2), scr(maxs2)
      integer ilab(maxg2), jlab(maxs2), klab(maxs2), llab(maxs2)
      integer i, j, k ,l, ish, jsh, ksh, lsh
      integer ij, ik, kl, jl
      integer count, nint, ijkl, me, nproc
      double precision zerotol
c
      call ga_zero(g_q)
      zerotol = 1d-12
      me = ga_nodeid()
      nproc = ga_nnodes()
c  
      count = 0
      do ish = 1, nsh
         do jsh = 1, nsh
            do ksh = 1, nsh
               count = count + 1
               if (mod(count,nproc) .eq. me) then
                  do lsh = 1, nsh
                     call int_l2e4c(basis, ish, jsh, basis, ksh, lsh,
     &                    zerotol, .false., maxg2, buf, nint, 
     $                    ilab, jlab, klab, llab, maxs2, scr)
                     do ijkl = 1, nint
                        i = ilab(ijkl)
                        j = jlab(ijkl)
                        k = klab(ijkl)
                        l = llab(ijkl)
                        if (.not. oexch) then
                           ij = (i-1)*nbf + j
                           kl = (k-1)*nbf + l
                           call ga_put(g_q, ij,ij,kl,kl,buf(ijkl),1)
                        else
                           ik = (i-1)*nbf + k
                           jl = (j-1)*nbf + l
                           call ga_put(g_q, ik,ik,jl,jl,buf(ijkl),1)
                        endif
                     enddo
                  enddo
               endif
            enddo
         enddo
      enddo
      call ga_sync
      call ga_symmetrize(g_q)
c
      end
