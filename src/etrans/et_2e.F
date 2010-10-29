c     $Id$
      subroutine et_2e(rtdb,nbf,basis,geom,max2e,mem2,tol2e,
     &                 g_pa,g_pb,omega2)
      implicit none
#include "errquit.fh"
c
c     calculates the two-electron contribution to Vba, omega2
c
c     nbf, basis, geom                                          [input]
c     max2e, mem2 (for int_2e4c)                                [input]
c     g_pa, g_pb (handles for alpha and beta densities)         [input]
c     omega2 (two electron contrib to Vab)                      [output]
c
#include "rtdb.fh"
#include "schwarz.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "apiP.fh"
#include "bas.fh"
#include "geom.fh"
#include "msgids.fh"   
#include "stdio.fh"   
c
      integer lpa_ij, kpa_ij, lpa_kl, kpa_kl, lpa_ji, kpa_ji
      integer lpa_lk, kpa_lk, lpa_il, kpa_il, lpa_kj, kpa_kj
      integer lpa_jl, kpa_jl, lpa_ki, kpa_ki, lpa_ik, kpa_ik
      integer lpa_lj, kpa_lj, lpa_jk, kpa_jk, lpa_li, kpa_li
      integer lpb_ij, kpb_ij, lpb_kl, kpb_kl, lpb_ji, kpb_ji
      integer lpb_lk, kpb_lk, lpb_il, kpb_il, lpb_kj, kpb_kj
      integer lpb_jl, kpb_jl, lpb_ki, kpb_ki, lpb_ik, kpb_ik
      integer lpb_lj, kpb_lj, lpb_jk, kpb_jk, lpb_li, kpb_li
      integer rtdb
      integer g_pa,g_pb
      integer nbf, geom, basis,max2e, mem2
      integer l_gJ, k_gJ, l_scr, k_scr, l_den, k_den
      integer ish, jsh, ksh, lsh, nsh
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
      integer nproc, blen, max_sh_bf
      integer maxlsh, shellnum, ng, numi, numj
      double precision omega2, omega2_val
      integer nxtask, next
      external nxtask
      double precision smax, sij, skl
      double precision norm, G(10000), dij
      double precision tol2e
      double precision dabsmax, adenmax,bdenmax
      external  dabsmax
c     -------------------------------------------     
c     find the maximum value of the alpha density
c     -------------------------------------------     
c
      call ga_distribution(g_pa,ga_nodeid(),ilo,ihi,jlo,jhi)   
      numi=ihi-ilo+1
      numj=jhi-jlo+1
c
        if (.not.ma_push_get(MT_DBL,numi*numj,'temp dens',l_den,k_den))
     $  call errquit('et_2e:ma_push_get failed for temp dens',k_den,
     &       MA_ERR)
        call dfill(numi*numj, 0.0d0, dbl_mb(k_den), 1)
c
        call ga_get(g_pa, ilo, ihi, jlo, jhi, dbl_mb(k_den),1)
        adenmax = dabsmax(numi*numj,dbl_mb(k_den))
c
        call ga_sync()
c
        call ga_dgop(msg_et_adens,adenmax,1,'max')
c
c     -------------------------------------------     
c     find the maximum value of the beta density
c     -------------------------------------------     
c
      call ga_distribution(g_pb,ga_nodeid(),ilo,ihi,jlo,jhi)   
      numi=ihi-ilo+1
      numj=jhi-jlo+1
c
        call dfill(numi*numj, 0.0d0, dbl_mb(k_den), 1)
c
        call ga_get(g_pb, ilo, ihi, jlo, jhi, dbl_mb(k_den),1)
        bdenmax = dabsmax(numi*numj,dbl_mb(k_den))
c
        if (.not. ma_pop_stack(l_den)) call errquit('et_2e: pop',1,
     &       MA_ERR)
c
        call ga_sync()
c
        call ga_dgop(msg_et_bdens,bdenmax,1,'max')
c
c     -------------------------------------------------
c     initialize omega2 and establish 2e and scr arrays
c     -------------------------------------------------
c     
      omega2 = 0.0d0
c     
      if (.not.ma_push_get(MT_DBL, max2e, '2e buffer', l_gJ, k_gJ))
     $     call errquit('et_2e:ma_push_get failed for 2e J buffer',k_gJ,
     &       MA_ERR)
c     
      if (.not.ma_push_get(MT_DBL, mem2, '2e scr', l_scr, k_scr))
     $     call errquit('et_2e:ma_push_get failed for 2e scratch',k_scr,
     &       MA_ERR)
c     
c     ------------------
c     get basis set info
c     ------------------
c     
      if (.not. bas_numcont(basis, nsh)) call errquit
     $     ('et_2e: bas_numcont failed', basis, BASIS_ERR)
c     
      if (.not. bas_nbf_cn_max(basis,max_sh_bf))
     $     call errquit('et_2e: could not get max_sh_bf',0, BASIS_ERR)

      blen = min(nbf,36,6*max_sh_bf) ! d(6)*6=36, 6**4=1296 quartets
c
c     ----------------------------------
c     set and print screening parameters
c     ----------------------------------
c
      call schwarz_init(geom, basis)
      smax = schwarz_max() 
      dij = max(adenmax*adenmax,bdenmax*bdenmax)
c
 1212 format (1x,a30,1x,a1,1pE14.4)
      if (ga_nodeid().eq.0) then
      write(LuOut,*)
      write(LuOut,*)'Start two-electron contribution,  H2(RP)...'
      write(LuOut,1212)'Density screening tolerance   ',':',dij
      write(LuOut,1212)'Largest two-electron integral ',':',smax
      write(LuOut,1212)'Integral screening threshold  ',':',tol2e
      write(LuOut,*)
      endif
c
c     ---------------------------------
c     allocate temporary density arrays
c     ---------------------------------
c
      if (.not. ma_push_get(mt_dbl, blen*blen, 'pa_ij', lpa_ij, kpa_ij))
     $     call errquit('et_2e:could not allocate lpa_ij',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, blen*blen, 'pa_kl', lpa_kl, kpa_kl))
     $     call errquit('et_2e:could not allocate lpa_kl',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, blen*blen, 'pa_ji', lpa_ji, kpa_ji))
     $     call errquit('et_2e:could not allocate lpa_ji',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, blen*blen, 'pa_lk', lpa_lk, kpa_lk))
     $     call errquit('et_2e:could not allocate lpa_lk',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, blen*blen, 'pa_il', lpa_il, kpa_il))
     $     call errquit('et_2e:could not allocate lpa_il',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, blen*blen, 'pa_kj', lpa_kj, kpa_kj))
     $     call errquit('et_2e:could not allocate lpa_kj',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, blen*blen, 'pa_jl', lpa_jl, kpa_jl))
     $     call errquit('et_2e:could not allocate lpa_jl',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, blen*blen, 'pa_ki', lpa_ki, kpa_ki))
     $     call errquit('et_2e:could not allocate lpa_ki',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, blen*blen, 'pa_ik', lpa_ik, kpa_ik))
     $     call errquit('et_2e:could not allocate lpa_ik',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, blen*blen, 'pa_lj', lpa_lj, kpa_lj))
     $     call errquit('et_2e:could not allocate lpa_lj',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, blen*blen, 'pa_jk', lpa_jk, kpa_jk))
     $     call errquit('et_2e:could not allocate lpa_jk',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, blen*blen, 'pa_li', lpa_li, kpa_li))
     $     call errquit('et_2e:could not allocate lpa_li',0, MA_ERR)
c
      if (.not. ma_push_get(mt_dbl, blen*blen, 'pb_ij', lpb_ij, kpb_ij))
     $     call errquit('et_2e:could not allocate lpb_ij',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, blen*blen, 'pb_kl', lpb_kl, kpb_kl))
     $     call errquit('et_2e:could not allocate lpb_kl',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, blen*blen, 'pb_ji', lpb_ji, kpb_ji))
     $     call errquit('et_2e:could not allocate lpb_ji',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, blen*blen, 'pb_lk', lpb_lk, kpb_lk))
     $     call errquit('et_2e:could not allocate lpb_lk',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, blen*blen, 'pb_il', lpb_il, kpb_il))
     $     call errquit('et_2e:could not allocate lpb_il',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, blen*blen, 'pb_kj', lpb_kj, kpb_kj))
     $     call errquit('et_2e:could not allocate lpb_kj',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, blen*blen, 'pb_jl', lpb_jl, kpb_jl))
     $     call errquit('et_2e:could not allocate lpb_jl',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, blen*blen, 'pb_ki', lpb_ki, kpb_ki))
     $     call errquit('et_2e:could not allocate lpb_ki',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, blen*blen, 'pb_ik', lpb_ik, kpb_ik))
     $     call errquit('et_2e:could not allocate lpb_ik',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, blen*blen, 'pb_lj', lpb_lj, kpb_lj))
     $     call errquit('et_2e:could not allocate lpb_lj',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, blen*blen, 'pb_jk', lpb_jk, kpb_jk))
     $     call errquit('et_2e:could not allocate lpb_jk',0, MA_ERR)
      if (.not. ma_push_get(mt_dbl, blen*blen, 'pb_li', lpb_li, kpb_li))
     $     call errquit('et_2e:could not allocate lpb_li',0, MA_ERR)
c
c     ----------------------     
c     begin loop over shells
c     ----------------------     
c     
      call ga_sync()
      shellnum = 0
      nproc = ga_nnodes()
      next = nxtask(nproc,1)
c
      do ish = nsh, 1, -1
         if (.not. bas_cn2bfr(basis, ish, ilo, ihi))
     $        call errquit('et_2e: bas_cn2bfr', ish, BASIS_ERR)
c     
         do jsh = ish, 1, -1
            if (.not. bas_cn2bfr(basis, jsh, jlo, jhi))
     $           call errquit('et_2e: bas_cn2bfr', jsh, BASIS_ERR)
c
              sij = schwarz_shell(ish,jsh)
                if (sij*smax*dij .gt. tol2e) then
c     
               do ksh = ish, 1, -1
                  if (.not. bas_cn2bfr(basis, ksh, klo, khi))
     $                 call errquit('et_2e: bas_cn2bfr', ksh, BASIS_ERR)
c     
                  maxlsh = ksh
                  if (ksh.EQ.ish) maxlsh = jsh
                  do lsh = maxlsh, 1, -1
                     if (.not. bas_cn2bfr(basis, lsh, llo, lhi))
     $                 call errquit('et_2e: bas_cn2bfr', lsh, BASIS_ERR)
c
                   if (shellnum.eq.next) then
c
                     skl = schwarz_shell(ksh,lsh)      
c
                   call ga_get(g_pa,ilo,ihi,jlo,jhi,dbl_mb(kpa_ij),blen)
                   call ga_get(g_pa,klo,khi,llo,lhi,dbl_mb(kpa_kl),blen)
                   call ga_get(g_pa,jlo,jhi,ilo,ihi,dbl_mb(kpa_ji),blen)
                   call ga_get(g_pa,llo,lhi,klo,khi,dbl_mb(kpa_lk),blen)
                   call ga_get(g_pa,ilo,ihi,llo,lhi,dbl_mb(kpa_il),blen)
                   call ga_get(g_pa,klo,khi,jlo,jhi,dbl_mb(kpa_kj),blen)
                   call ga_get(g_pa,jlo,jhi,llo,lhi,dbl_mb(kpa_jl),blen)
                   call ga_get(g_pa,klo,khi,ilo,ihi,dbl_mb(kpa_ki),blen)
                   call ga_get(g_pa,ilo,ihi,klo,khi,dbl_mb(kpa_ik),blen)
                   call ga_get(g_pa,llo,lhi,jlo,jhi,dbl_mb(kpa_lj),blen)
                   call ga_get(g_pa,jlo,jhi,klo,khi,dbl_mb(kpa_jk),blen)
                   call ga_get(g_pa,llo,lhi,ilo,ihi,dbl_mb(kpa_li),blen)
                   call ga_get(g_pb,ilo,ihi,jlo,jhi,dbl_mb(kpb_ij),blen)
                   call ga_get(g_pb,klo,khi,llo,lhi,dbl_mb(kpb_kl),blen)
                   call ga_get(g_pb,jlo,jhi,ilo,ihi,dbl_mb(kpb_ji),blen)
                   call ga_get(g_pb,llo,lhi,klo,khi,dbl_mb(kpb_lk),blen)
                   call ga_get(g_pb,ilo,ihi,llo,lhi,dbl_mb(kpb_il),blen)
                   call ga_get(g_pb,klo,khi,jlo,jhi,dbl_mb(kpb_kj),blen)
                   call ga_get(g_pb,jlo,jhi,llo,lhi,dbl_mb(kpb_jl),blen)
                   call ga_get(g_pb,klo,khi,ilo,ihi,dbl_mb(kpb_ki),blen)
                   call ga_get(g_pb,ilo,ihi,klo,khi,dbl_mb(kpb_ik),blen)
                   call ga_get(g_pb,llo,lhi,jlo,jhi,dbl_mb(kpb_lj),blen)
                   call ga_get(g_pb,jlo,jhi,klo,khi,dbl_mb(kpb_jk),blen)
                   call ga_get(g_pb,llo,lhi,ilo,ihi,dbl_mb(kpb_li),blen)
c
                        call dens_prod(nbf,
     $                    ilo,ihi,jlo,jhi,klo,khi,llo,lhi, 
     $                    ish,jsh,ksh,lsh, norm,G,ng,blen,
     $                    dbl_mb(kpa_ij),dbl_mb(kpa_kl),dbl_mb(kpa_ji),
     $                    dbl_mb(kpa_lk),dbl_mb(kpa_il),dbl_mb(kpa_kj),
     $                    dbl_mb(kpa_jl),dbl_mb(kpa_ki),dbl_mb(kpa_ik),
     $                    dbl_mb(kpa_lj),dbl_mb(kpa_jk),dbl_mb(kpa_li),
     $                    dbl_mb(kpb_ij),dbl_mb(kpb_kl),dbl_mb(kpb_ji),
     $                    dbl_mb(kpb_lk),dbl_mb(kpb_il),dbl_mb(kpb_kj),
     $                    dbl_mb(kpb_jl),dbl_mb(kpb_ki),dbl_mb(kpb_ik),
     $                    dbl_mb(kpb_lj),dbl_mb(kpb_jk),dbl_mb(kpb_li))
c
                        if (sij*skl*smax*norm.gt.tol2e) then
c
                           call int_2e4c(basis, ish, jsh, basis, 
     $                          ksh, lsh, mem2, dbl_mb(k_scr), max2e, 
     $                          dbl_mb(k_gJ))
c     
                           call omega2_mult(ilo,ihi,jlo,jhi,
     &                          klo,khi,llo,lhi,
     &                          ish,jsh,ksh,lsh,
     &                          G,ng,dbl_mb(k_gJ),omega2_val)
c
                           omega2 = omega2 + omega2_val
c     
                       endif  !norm screening
                   next = nxtask(nproc,1)
                   endif   !shellnum
                   shellnum = shellnum + 1
                  enddo
               enddo
            endif  !sij x smax
         enddo
      enddo
c     
      next = nxtask(-nproc,1)
      call ga_sync()
c     
      call ga_dgop(msg_et_2e,omega2,1,'+')
c     
c     destroy arrays
c     ==============
      call schwarz_tidy()
c     
      if (.not. ma_pop_stack(lpb_li)) call errquit('et_2e: pop',1,
     &       MA_ERR)
      if (.not. ma_pop_stack(lpb_jk)) call errquit('et_2e: pop',2,
     &       MA_ERR)
      if (.not. ma_pop_stack(lpb_lj)) call errquit('et_2e: pop',3,
     &       MA_ERR)
      if (.not. ma_pop_stack(lpb_ik)) call errquit('et_2e: pop',4,
     &       MA_ERR)
      if (.not. ma_pop_stack(lpb_ki)) call errquit('et_2e: pop',5,
     &       MA_ERR)
      if (.not. ma_pop_stack(lpb_jl)) call errquit('et_2e: pop',6,
     &       MA_ERR)
      if (.not. ma_pop_stack(lpb_kj)) call errquit('et_2e: pop',7,
     &       MA_ERR)
      if (.not. ma_pop_stack(lpb_il)) call errquit('et_2e: pop',8,
     &       MA_ERR)
      if (.not. ma_pop_stack(lpb_lk)) call errquit('et_2e: pop',9,
     &       MA_ERR)
      if (.not. ma_pop_stack(lpb_ji)) call errquit('et_2e: pop',10,
     &       MA_ERR)
      if (.not. ma_pop_stack(lpb_kl)) call errquit('et_2e: pop',11,
     &       MA_ERR)
      if (.not. ma_pop_stack(lpb_ij)) call errquit('et_2e: pop',12,
     &       MA_ERR)
      if (.not. ma_pop_stack(lpa_li)) call errquit('et_2e: pop',13,
     &       MA_ERR)
      if (.not. ma_pop_stack(lpa_jk)) call errquit('et_2e: pop',14,
     &       MA_ERR)
      if (.not. ma_pop_stack(lpa_lj)) call errquit('et_2e: pop',15,
     &       MA_ERR)
      if (.not. ma_pop_stack(lpa_ik)) call errquit('et_2e: pop',16,
     &       MA_ERR)
      if (.not. ma_pop_stack(lpa_ki)) call errquit('et_2e: pop',17,
     &       MA_ERR)
      if (.not. ma_pop_stack(lpa_jl)) call errquit('et_2e: pop',18,
     &       MA_ERR)
      if (.not. ma_pop_stack(lpa_kj)) call errquit('et_2e: pop',19,
     &       MA_ERR)
      if (.not. ma_pop_stack(lpa_il)) call errquit('et_2e: pop',20,
     &       MA_ERR)
      if (.not. ma_pop_stack(lpa_lk)) call errquit('et_2e: pop',21,
     &       MA_ERR)
      if (.not. ma_pop_stack(lpa_ji)) call errquit('et_2e: pop',22,
     &       MA_ERR)
      if (.not. ma_pop_stack(lpa_kl)) call errquit('et_2e: pop',23,
     &       MA_ERR)
      if (.not. ma_pop_stack(lpa_ij)) call errquit('et_2e: pop',24,
     &       MA_ERR)
      if (.not. ma_pop_stack(l_scr)) call errquit('et_2e: pop',25,
     &       MA_ERR)
      if (.not. ma_pop_stack(l_gJ))  call errquit('et_2e: pop',26,
     &       MA_ERR)
c
c
      return
      end



