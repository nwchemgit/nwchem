c     
c
c  $Id: moints_6x.F,v 1.8 1995-10-27 23:25:18 rg240 Exp $
c
c
c  This routine returns the Coulomb and exchange integral
c  operator matrices for the range of MO-indices as mo_indx_hi, mo_indx_lo
c  The g_coul, g_exch global arrays are ordered as
c
c
c
c               ij
c  (ij|ab) = ( J  )  = g_coul[ ij : (a-1)*N2 + b ] = g_coul [ ij : (b-1)*N2 + a ]
c                  ab
c
c               ij
c  (ia|jb) = ( K  )  = g_exch[ ij : (a-1)*N2 + b ]
c                  ab
c
c -----------------
c Note: This routine differs from the standard in that
c it computes the integrals 6 times as opposed to 2 times
c in the standard. However, it has lower communication scaling
c for K, O(N^4), compared to O(N^5) in the standard. Therefore
c this routine is more appropriate for very high parallelism
c -----------------
c
c
c
       subroutine moints_build_6x(basis, osym,
     $                            occ_start, mo1_lo, mo1_hi,
     $                            mo2_lo, mo2_hi,
     $                            g_movecs,
     $                            g_coul, ocoul,
     $                            g_exch, oexch,
     $                            blen )
       implicit none
#include "tcgmsg.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "util.fh"
#include "schwarz.fh"
#include "msgids.fh"
c
c Arguments
c
       integer basis                          ! Basis handle
       logical osym                           ! Symmetry toggle
       integer occ_start                      ! Offset from frozen core
       integer mo1_lo, mo1_hi                 ! 1st Pair Index range
       integer mo2_lo, mo2_hi                 ! 2nd Pair Index range
       integer g_movecs                       ! MO coefficients
       integer g_coul                         ! Coulomb operator
       integer g_exch                         ! Exchange operator
       logical ocoul,oexch                    ! Type selection
       integer blen                           ! Blocking length
c
c Local variables
c
       integer nmo1, nmo2, nbf, nsh, maxbfsh
       integer nblk, bsize, blocksize
       integer l_gmap, k_gmap
       integer ish0, jsh0, ish, jsh, ilen, jlen
       integer ibflo,ibfhi,jbflo,jbfhi,kbflo,kbfhi,lbflo,lbfhi
       integer kshlo, kshhi, lshlo, lshhi
       integer kblen, lblen
       integer kb, lb, i, nmixed
       integer l_ssbb, k_ssbb, l_ssbbt, k_ssbbt
       integer l_hlp, k_hlp, l_ssni,k_ssni
       integer l_eri, k_eri, l_iscr,k_iscr
       integer l_mo, k_mo
       integer l_shmap, k_shmap
       integer l_sym, k_sym
       integer n_ssbb, n_ssbb1, n_ssni, n_hlp, n_ijni
       integer mem2, max2e
       integer num_nodes, ploop, next
       double precision tol2e, scale, scale0, schw_ij
       double precision tz, thalf, tint, t1jidx, t2jidx, t1kidx, t2kidx
       double precision t34kidx, t34jidx, thalf0, ttotal, tsynch
       double precision ttask, ttaskmaxJ, ttaskminJ, ttaskaggJ
       double precision ttaskmaxK, ttaskminK, ttaskaggK
       double precision flop1, flop2k, flop2j
       integer tottaskJ, tottaskK
       logical status
       logical oprint
       logical osym
#include "moints_stats.fh"
c
c
c
       integer moints_numgr, gr_len
       external moints_numgr, gr_len
c
c
c
       data tol2e/1.d-12/
       data oprint/.true./
c
c
c  General basis info
c
       ttotal = tcgtime()
       num_nodes = ga_nnodes()
       status = bas_numbf(basis,nbf)
       status = status.and.bas_numcont(basis,nsh)
       status = status.and.bas_nbf_cn_max(basis,maxbfsh)
       if (.not.status) call errquit('moints: cannot get basis info',0)
       nmo1 = mo1_hi - occ_start + 1
       nmo2 = mo2_hi - mo2_lo + 1
c     
c     Symmetry adapt the MOs and renumber irreps to start at zero
c     
       if (osym) then
          if (.not. ma_push_get(MT_INT, nbf, 'movec syms',l_sym, k_sym))
     $         call errquit('moints_6x: no memory for syms?',0)
          call sym_movecs_adapt(basis, 1d-8, g_movecs, int_mb(k_sym),
     $                          nmixed)
          if (nmixed .ne. 0)
     $         call errquit('moints_6x: symmetry contamination', nmixed)
          do i =0, nbf-1
             int_mb(k_sym+i) = int_mb(k_sym+i) - 1
          enddo
          if (util_print('orbital symmetry',print_debug)) then
            write(6,887) 
 887        format('Symmetry of MOs')
            write(6,888) (int_mb(k_sym+i),i=0,nbf-1)
 888        format(16i3)
          endif
       endif
c
c  Local MO coefficients
c
       status = status.and.ma_push_get(MT_DBL,(nbf*nbf),
     $                                 'movecs cols',l_mo,k_mo)
       call ga_get(g_movecs,1,nbf,1,nbf,dbl_mb(k_mo),nbf)
c
c  Integrals allocation
c
       call int_mem_2e4c(max2e, mem2)
       status = ma_push_get(MT_DBL, max2e,'moints: buf', l_eri, k_eri)
       status = ma_push_get(MT_DBL, mem2, 'moints: scr', l_iscr, k_iscr)
c
c  Shell group mapping
c
       blocksize = blen
       nblk = moints_numgr( basis, blocksize )
       status = ma_push_get(MT_INT,(nblk*4),'shell group map',
     $                      l_gmap, k_gmap)
       call moints_grmap( basis, blocksize, nblk, int_mb(k_gmap))
c
c  Reorder shells by descending shell-length
c
       status = ma_push_get(MT_INT,nsh,'shell order map',
     $                      l_shmap, k_shmap)
       call moints_shorder( basis, nsh, maxbfsh, int_mb(k_shmap) )
c
c  Temporary partially-transformed arrays
c     
       bsize = max(blocksize,maxbfsh)
       n_ssbb = maxbfsh*maxbfsh*bsize*bsize
       n_ssbb1 = max((nmo1*nmo1),n_ssbb)
       n_hlp = max((bsize*maxbfsh*maxbfsh*nmo1),(maxbfsh*nbf))
       n_ssni = maxbfsh*maxbfsh*nbf*nmo1
       status = ma_push_get(MT_DBL,n_ssbb1,'ssbb block',l_ssbb,k_ssbb)
       status = ma_push_get(MT_DBL,n_ssbb,'ssbbt block',l_ssbbt,k_ssbbt)
       status = ma_push_get(MT_DBL,n_hlp,'hlp block',l_hlp,k_hlp)
       status = ma_push_get(MT_DBL,n_ssni,'ssni block',l_ssni,k_ssni)
       if (.not.(status)) call errquit('cannot allocate local memory',0)
c
c Initialize
c
       flop1 = 0.d0
       flop2k = 0.d0
       flop2j = 0.d0
       tint = 0.d0
       t1jidx = 0.d0
       t2jidx = 0.d0
       t1kidx = 0.d0
       t2kidx = 0.d0
       t34jidx = 0.d0
       t34kidx = 0.d0
       tottaskJ = 0
       tottaskK = 0
       ttaskmaxJ = 0.d0
       ttaskminJ = 1.d24
       ttaskaggJ = 0.d0
       ttaskmaxK = 0.d0
       ttaskminK = 1.d24
       ttaskaggK = 0.d0
       thalf0 = tcgtime()
       ploop = 0
       next = -1
       if (ocoul) call ga_zero(g_coul)
       if (oexch) call ga_zero(g_exch)
c
c  Exchange part
c  Four fold shell loop
c
c
       if (oexch) then
         if (next.lt.0) next = nxtval(num_nodes)
         do ish0=1,nsh
           ish = int_mb(k_shmap+ish0-1)
           status = bas_cn2bfr(basis,ish,ibflo,ibfhi)
           ilen = ibfhi - ibflo + 1
           do jsh0=1,ish0
             jsh = int_mb(k_shmap+jsh0-1)
             status = bas_cn2bfr(basis,jsh,jbflo,jbfhi)
             jlen = jbfhi - jbflo + 1
             scale0 = 1.d0
             if (ish.eq.jsh) scale0 = 0.5d0
             if (next.eq.ploop) then
               ttask = tcgtime()
               tottaskK = tottaskK + 1
               call dfill(n_ssni,0.d0,dbl_mb(k_ssni),1)
               do kb=1,nblk
                 call gr_sh_idx(nblk,int_mb(k_gmap),kb,kshlo,kshhi)
                 kblen = gr_len(nblk,int_mb(k_gmap),kb,kbflo,kbfhi)
                 do lb=1,kb
                   call gr_sh_idx(nblk,int_mb(k_gmap),lb,lshlo,lshhi)
                   lblen = gr_len(nblk,int_mb(k_gmap),lb,lbflo,lbfhi)
                   tz = tcgtime()
                   scale = 1.d0
                   if (kb.eq.lb) scale = 0.5d0
                   call moints6x_gblkk( basis, ish, jsh, kshlo, kshhi,
     $                                  lshlo, lshhi, schw_ij, tol2e,
     $                                  max2e, dbl_mb(k_eri), mem2,
     $                                  dbl_mb(k_iscr), ilen, jlen,
     $                                  kbflo, kbfhi, lbflo, lbfhi,
     $                                  dbl_mb(k_ssbb),
     $                                  dbl_mb(k_ssbbt) )
                   tint = tint + tcgtime() - tz
                   flop1 = flop1 + 4*lblen*ilen*jlen*kblen*nmo1
                   tz = tcgtime()
                   call moints6x_trf1K( nbf, occ_start, mo1_hi, ilen,
     $                                  jlen, kbflo, kbfhi, lbflo,
     $                                  lbfhi, scale, dbl_mb(k_ssbb),
     $                                  dbl_mb(k_ssbbt), dbl_mb(k_mo),
     $                                  dbl_mb(k_ssni), dbl_mb(k_hlp))
                   t1kidx = t1kidx + tcgtime() - tz
                 enddo
               enddo
               tz = tcgtime()
               call moints6x_trf2K( nbf, occ_start, mo1_lo, mo1_hi,
     $                              ibflo, ibfhi, jbflo, jbfhi,
     $                              dbl_mb(k_ssni),
     $                              dbl_mb(k_hlp), dbl_mb(k_ssbb),
     $                              dbl_mb(k_hlp), dbl_mb(k_ssbb),
     $                              dbl_mb(k_mo), g_exch )
               t2kidx = t2kidx + tcgtime() - tz
               next = nxtval(num_nodes)
               ttask = tcgtime() - ttask
               ttaskmaxK = max(ttaskmaxK,ttask)
               ttaskminK = min(ttaskminK,ttask)
               ttaskaggK = ttaskaggK + ttask
             endif
             ploop = ploop + 1
           enddo
         enddo
       endif
c
c
c
c$$$       if (oexch) then
c$$$         if (util_print('exchange half integral',print_debug)) then
c$$$           call moints_op_print(occ_start,mo1_lo,mo1_hi,nbf,g_exch)
c$$$         endif
c$$$       endif
c
c  Coulomb part
c  4-fold shell loop
c
       if (ocoul) then
         if (next.lt.0) next = nxtval(num_nodes)
         do ish0=1,nsh
           do jsh0=1,ish0
             ish = max(int_mb(k_shmap+ish0-1),int_mb(k_shmap+jsh0-1))
             jsh = min(int_mb(k_shmap+ish0-1),int_mb(k_shmap+jsh0-1))
             status = bas_cn2bfr(basis,ish,ibflo,ibfhi)
             status = bas_cn2bfr(basis,jsh,jbflo,jbfhi)
             ilen = ibfhi - ibflo + 1
             jlen = jbfhi - jbflo + 1
             schw_ij = schwarz_shell(ish,jsh)
             scale = 1.d0
             if (ish.eq.jsh) scale = scale*0.5d0
             if (schw_ij*schwarz_max().ge.tol2e) then
               if (next.eq.ploop) then
                 ttask = tcgtime()
                 tottaskJ = tottaskJ + 1
                 n_ijni = ilen*jlen*nbf*nmo1
                 call dfill(n_ijni,0.d0,dbl_mb(k_ssni),1)
                 do kb=1,nblk
                   call gr_sh_idx(nblk,int_mb(k_gmap),kb,kshlo,kshhi)
                   kblen = gr_len(nblk,int_mb(k_gmap),kb,kbflo,kbfhi)
                   do lb=1,kb
                     call gr_sh_idx(nblk,int_mb(k_gmap),lb,lshlo,lshhi)
                     lblen = gr_len(nblk,int_mb(k_gmap),lb,lbflo,lbfhi)
                     tz = tcgtime()
                     call moints_gblk( basis, ish, jsh, kshlo, kshhi,
     $                                 lshlo, lshhi, schw_ij, tol2e,
     $                                 max2e, dbl_mb(k_eri), mem2,
     $                                 dbl_mb(k_iscr), ilen, jlen,
     $                                 kbflo, kbfhi, lbflo, lbfhi,
     $                                 dbl_mb(k_ssbb), osym )
                     tint = tint + tcgtime() - tz
                     flop1 = flop1 + 4*lblen*ilen*jlen*kblen*nmo1
                     if (lb.ne.kb) then
                       tz = tcgtime()
                       call moints_blktr( ilen, jlen, kblen, lblen,
     $                                    dbl_mb(k_ssbb),
     $                                    dbl_mb(k_ssbbt))
                       call moints_trf1( nbf, mo1_lo, mo1_hi, ilen,
     $                                   jlen, kbflo, kbfhi, lbflo,
     $                                   lbfhi, 1.d0, dbl_mb(k_ssbb),
     $                                   dbl_mb(k_ssbbt),
     $                                   dbl_mb(k_mo), dbl_mb(k_ssni),
     $                                   dbl_mb(k_hlp) )
                     else
                       tz = tcgtime()
                       call moints_trf1( nbf, mo1_lo, mo1_hi, ilen,
     $                                   jlen, kbflo, kbfhi, lbflo,
     $                                   lbfhi, 0.5d0, dbl_mb(k_ssbb),
     $                                   dbl_mb(k_ssbb), dbl_mb(k_mo),
     $                                   dbl_mb(k_ssni),
     $                                   dbl_mb(k_hlp) )
                     endif
                     t1jidx = t1jidx + tcgtime() - tz
                   enddo
                 enddo
                 flop2j = flop2j+2*nbf*ilen*jlen*((nmo1*(nmo1+1))/2)
                 tz = tcgtime()
                 call moints_trf2J( nbf, occ_start, mo1_lo, mo1_hi,
     $                              ibflo, ibfhi, jbflo, jbfhi,
     $                              dbl_mb(k_ssni), dbl_mb(k_hlp),
     $                              dbl_mb(k_ssbb), dbl_mb(k_mo),
     $                              g_coul )
                 t2jidx = t2jidx + tcgtime() - tz
                 next = nxtval(num_nodes)
                 ttask = tcgtime() - ttask
                 ttaskmaxJ = max(ttaskmaxJ,ttask)
                 ttaskminJ = min(ttaskminJ,ttask)
                 ttaskaggJ = ttaskaggJ + ttask
               endif
               ploop = ploop + 1
             endif
           enddo
         enddo
       endif
c
c
c$$$       if (ocoul) then
c$$$         if (util_print('coulomb half integral',print_debug)) then
c$$$           call moints_op_print(occ_start,mo1_lo,mo1_hi,nbf,g_coul)
c$$$         endif
c$$$       endif
c
c Synchronize
c
       tz = tcgtime()
       next = nxtval(-num_nodes)
       tsynch = tcgtime() - tz
       thalf = tcgtime() - thalf0
c
c Clean-up
c
       if (.not. ma_pop_stack(l_ssni))
     $     call errquit('moints: failed to pop', l_ssni)
       if (.not. ma_pop_stack(l_hlp))
     $     call errquit('moints: failed to pop', l_hlp)
       if (.not. ma_pop_stack(l_ssbbt))
     $     call errquit('moints: failed to pop', l_ssbbt)
       if (.not. ma_pop_stack(l_ssbb))
     $     call errquit('moints: failed to pop', l_ssbb)
       if (.not. ma_pop_stack(l_shmap))
     $     call errquit('moints: failed to pop', l_shmap)
       if (.not. ma_pop_stack(l_gmap))
     $     call errquit('moints: failed to pop', l_gmap)
       if (.not. ma_pop_stack(l_iscr))
     $     call errquit('moints: failed to pop', l_iscr)
       if (.not. ma_pop_stack(l_eri))
     $     call errquit('moints: failed to pop', l_eri)
c
c 
c
       status = ma_push_get(MT_DBL,(nbf*nbf),'hlp',l_hlp,k_hlp)
c
c
c
       if (oexch) then
         tz = tcgtime()
         call moints_Ktrf34(g_exch, .false., occ_start, mo1_lo, mo1_hi,
     $                       mo2_lo, mo2_hi, nbf, dbl_mb(k_mo),
     $                       dbl_mb(k_hlp), osym, int_mb(k_sym) )
         t34kidx = tcgtime() - tz
         if (util_print('exchange integral',print_debug)) then
           call moints_op_print(occ_start,mo1_lo,mo1_hi,nbf,g_exch)
         endif
       endif
c
c
c
       if (ocoul) then
         tz = tcgtime()
         call moints_Jtrf34( g_coul, occ_start, mo1_lo, mo1_hi,
     $                       mo2_lo, mo2_hi, nbf, dbl_mb(k_mo),
     $                       dbl_mb(k_hlp), osym, int_mb(k_sym) )
         if (util_print('coulomb integral',print_debug)) then
           call moints_op_print(occ_start,mo1_lo,mo1_hi,nbf,g_coul)
         endif
         t34jidx = tcgtime() - tz
       endif
       call ga_sync()
c
c  Clean-up
c
       if (.not. ma_pop_stack(l_hlp))
     $     call errquit('moints: failed to pop', l_hlp)
       if (.not. ma_pop_stack(l_mo))
     $     call errquit('moints: failed to pop', l_mo)
       if (osym) then
          if (.not. ma_pop_stack(l_sym))
     $         call errquit('moints_2x: memory corrupt',0)
       endif
c
c
c
#ifdef BLOCK_TRANSF
       if (ga_nodeid().eq.0) write(6,965)
 965   format(/,10x,'**** BLOCK TRANSFER ENABLED ****')
#endif
#ifdef NOCOMMS
       if (ga_nodeid().eq.0) write(6,334)
 334   format(/,10x,'**** COMMUNICATION DISABLED ****')
#endif


c
c Timings and Statistics bookkeeping
c
       flop1 = flop1*1.d-6
       flop2k = flop2k*1.d-6
       flop2j = flop2j*1.d-6
       ttotal = tcgtime() - ttotal
       mi_npass = mi_npass + 1.d0
       mi_ttotal = mi_ttotal + ttotal
       mi_thalf = mi_thalf + thalf
       mi_tint = mi_tint + tint
       mi_t1j = mi_t1j + t1jidx
       mi_t2j = mi_t2j + t2jidx
       mi_t1k = mi_t1k + t1kidx
       mi_t2k = mi_t2k + t2kidx
       mi_t34k = mi_t34k + t34kidx
       mi_t34j = mi_t34j + t34jidx
       mi_flop1 = mi_flop1 + flop1
       mi_synch = mi_synch + tsynch
       mi_maxsynch = max(mi_maxsynch,tsynch)
       mi_minsynch = min(mi_minsynch,tsynch)
       mi_aggsynch = mi_aggsynch + tsynch
       mi_nsynchs = mi_nsynchs + 1
       if (ocoul) then
         mi_maxtaskJ = max(mi_maxtaskJ,ttaskmaxJ)
         mi_mintaskJ = min(mi_mintaskJ,ttaskminJ)
         mi_aggtaskJ = mi_aggtaskJ + ttaskaggJ
         mi_ntasksJ = mi_ntasksJ + tottaskJ
       else
         mi_maxtaskJ = 0.d0
         mi_mintaskJ = 0.d0
         mi_aggtaskJ = 0.d0
         mi_ntasksJ = 1.d0
       endif
       if (oexch) then
         mi_maxtaskK = max(mi_maxtaskK,ttaskmaxK)
         mi_mintaskK = min(mi_mintaskK,ttaskminK)
         mi_aggtaskK = mi_aggtaskK + ttaskaggK
         mi_ntasksK = mi_ntasksK + tottaskK
       else
         mi_maxtaskK = 0.d0
         mi_mintaskK = 0.d0
         mi_aggtaskK = 0.d0
         mi_ntasksK = 1.d0
       endif

       return
       end










      subroutine moints6x_gblkk( basis, ish, jsh, kshlo, kshhi,
     $                         lshlo, lshhi, schw_ij, tol2e,
     $                         max2e, eri, mem2, iscr, ilen,
     $                         jlen, kblo, kbhi, lblo, lbhi,
     $                         sbsb, sbsbt )
      implicit none
#include "bas.fh"
#include "schwarz.fh"
      integer basis, ish, jsh, kshlo, kshhi, lshlo, lshhi
      integer max2e, mem2
      double precision schw_ij, tol2e, eri(*), iscr(*)
      integer ilen, jlen, kblo, kbhi, lblo, lbhi
      double precision sbsb(lblo:lbhi,kblo:kbhi,jlen,ilen)
      double precision sbsbt(kblo:kbhi,lblo:lbhi,jlen,ilen)
c
      integer ksh, lsh, kbflo, kbfhi, lbflo, lbfhi, ltop
      integer klen, llen, kblen, lblen, bsize, k1, k2, j, i
      logical status

 
      kblen = kbhi - kblo + 1
      lblen = lbhi - lblo + 1
      bsize = kblen*lblen*ilen*jlen
      call dfill(bsize,0.d0,sbsb,1)
      call dfill(bsize,0.d0,sbsbt,1)
      do ksh=kshlo,kshhi
        status = bas_cn2bfr(basis,ksh,kbflo,kbfhi)
        klen = kbfhi - kbflo + 1
        ltop = lshhi
        if (kshlo.eq.lshlo) ltop = ksh
        do lsh=lshlo,ltop
C        do lsh=lshlo,lshhi
          status = bas_cn2bfr(basis,lsh,lbflo,lbfhi)
          llen = lbfhi - lbflo + 1
          if (schwarz_shell(ish,ksh)*
     $        schwarz_shell(jsh,lsh).ge.tol2e) then
            call int_2e4c(basis, ish, ksh, basis, jsh, lsh,
     $                    mem2, iscr, max2e, eri )
            call moints6x_eri2blkk(ilen, jlen, klen, llen,
     $                             eri, sbsb(lbflo,kbflo,1,1),
     $                             lblen, kblen )
          endif
          if (schwarz_shell(ish,lsh)*
     $        schwarz_shell(jsh,ksh).ge.tol2e) then
            call int_2e4c(basis, ish, lsh, basis, jsh, ksh,
     $                    mem2, iscr, max2e, eri )
            call moints6x_eri2blkk(ilen, jlen, llen, klen,
     $                             eri, sbsbt(kbflo,lbflo,1,1),
     $                             kblen, lblen )
          endif
        enddo
      enddo
      if (kshlo.eq.lshlo) then
        do i=1,ilen
          do j=1,jlen
            do k1=kblo,kbhi
              do k2=kblo,k1-1
                sbsb(k1,k2,j,i) = sbsbt(k1,k2,j,i)
                sbsbt(k2,k1,j,i) = sbsb(k2,k1,j,i)
              enddo
            enddo
          enddo
        enddo
      endif

      return
      end











      subroutine moints6x_eri2blkk( ilen, jlen, klen, llen,
     $                              eri, blk, lblen, kblen )
      implicit none
      integer ilen, jlen, klen, llen, lblen, kblen
      double precision blk(lblen, kblen, jlen, ilen )
      double precision eri(llen, jlen, klen, ilen )
      integer k,l,i,j

      do i=1,ilen
        do j=1,jlen
          do k=1,klen
            do l=1,llen
              blk(l,k,j,i) = eri(l,j,k,i)
            enddo
          enddo
        enddo
      enddo
      
      return
      end










      subroutine moints6x_pushK( nbf, ilen, jlen, klo, khi,
     $                           nmo, x, ssni )
      implicit none
      integer nbf, ilen, jlen, klo, khi, nmo
      double precision x(klo:khi,jlen,ilen,nmo)
      double precision ssni(nbf,jlen,ilen,nmo)
      integer a,i,j,k
      
      do a=1,nmo
        do i=1,ilen
          do j=1,jlen
            do k=klo,khi
              ssni(k,j,i,a) = ssni(k,j,i,a) + x(k,j,i,a)
            enddo
          enddo
        enddo
      enddo
      return
      end







      subroutine moints6x_trf1K( nbf, molo, mohi, ilen, jlen,
     $                           klo, khi, llo, lhi, scale,
     $                           ssbb, ssbbt, c, ssni, hlp )
      implicit none
      integer nbf, molo, mohi, ilen, jlen, klo, khi, llo, lhi
      double precision scale
      double precision ssbb(llo:lhi,klo:khi,jlen,ilen)
      double precision ssbbt(klo:khi,llo:lhi,jlen,ilen)
      double precision c(nbf,nbf)
      double precision ssni(nbf,jlen,ilen,molo:mohi)
      double precision hlp(*)
c
c
      integer nmo, llen, klen, kjilen, ljilen
c
c      
      nmo = mohi - molo + 1
      llen = lhi - llo + 1
      klen = khi - klo + 1
      kjilen = klen*jlen*ilen
      ljilen = llen*jlen*ilen
      call dgemm( 't', 'n', kjilen, nmo, llen, scale, ssbb, llen,
     $            c(llo,molo), nbf, 0.d0, hlp, kjilen )
      call moints6x_pushK( nbf, ilen, jlen, klo, khi,
     $                     nmo, hlp, ssni )
      call dgemm( 't', 'n', ljilen, nmo, klen, scale, ssbbt, klen,
     $            c(klo,molo), nbf, 0.d0, hlp, ljilen )
      call moints6x_pushK( nbf, ilen, jlen, llo, lhi,
     $                     nmo, hlp, ssni )

      return
      end









      subroutine moints6x_trf2K( nbf, ostart, olo, ohi, ilo, ihi,
     $                           jlo, jhi, ssni, h1, h2, g1, g2,
     $                           c, g_exch )
      implicit none
      integer nbf, ostart, olo, ohi, ilo, ihi, jlo, jhi
      double precision ssni(nbf,jlo:jhi,ilo:ihi,ostart:ohi)
#ifdef SCATTER_TRANSF
      double precision h1(*)
      integer iv(100),jv(100)
#endif
      double precision h1(nbf,ilo:ihi)
      double precision h2(jlo:jhi,ilo:ihi)
      double precision g1(nbf,jlo:jhi)
      double precision g2(ilo:ihi,jlo:jhi)
      double precision c(nbf,nbf)
      integer g_exch
c
      integer ilen, jlen, ijlen
      integer nni, nnj, ijlo, ijhi, jilo, jihi
      integer aoff, ofroz, aa, bb, a, b, ab, i, j
c
      ofroz = ostart - 1
      aoff = ((olo-ofroz)*(olo-ofroz-1))/2
      ilen = ihi - ilo + 1
      jlen = jhi - jlo + 1
      ijlen = ilen*jlen
      nni = ilen*nbf
      nnj = jlen*nbf

#ifdef BLOCK_TRANSF
      ijlo = (ilo-1)*nbf + 1
      ijhi = ihi*nbf
#endif
      do a=olo,ohi
        do b=ostart,a
          call dgemm('t','n',ijlen,1,nbf,1.d0,ssni(1,jlo,ilo,b),
     $               nbf,c(1,a),nbf,0.d0,h2,ijlen)
#ifndef NOCOMMS
          aa = a - ofroz
          bb = b - ofroz
          ab = (aa*(aa-1))/2 + bb - aoff
#ifdef BLOCK_TRANSF
          call dfill(nni,0.d0,h1,1)
          do i=ilo,ihi
            do j=jlo,jhi
              h1(j,i) = h2(j,i)
            enddo
          enddo
          call ga_acc(g_exch,ijlo,ijhi,ab,ab,h1,nni,1.d0)
#else
#if SCATTER_TRANSF
          ij = 0
          do i=ilo,ihi
            ii = (i-1)*nbf
            do j=jlo,jhi
              ij = ij + 1
              h1(ij) = h2(j,i)
              iv(ij) = ii + j
              jv(ij) = ab
            enddo
          enddo
          call ga_dscatter(g_exch,h1,iv,jv,ijlen)
#else
          do i=ilo,ihi
            ijlo = (i-1)*nbf + jlo
            ijhi = (i-1)*nbf + jhi
            call ga_acc(g_exch,ijlo,ijhi,ab,ab,h2(jlo,i),ilen,1.d0)
          enddo
#endif
#endif
#endif
        enddo
      enddo



      if (ilo.eq.jlo) goto 120
#ifdef BLOCK_TRANSF
      jilo = (jlo-1)*nbf + 1
      jihi = jhi*nbf
#endif
      do a=olo,ohi
        do b=ostart,a
          call dgemm('t','n',ijlen,1,nbf,1.d0,ssni(1,jlo,ilo,a),
     $               nbf,c(1,b),nbf,0.d0,h2,ijlen)
#ifndef NOCOMMS
          aa = a - ofroz
          bb = b - ofroz
          ab = (aa*(aa-1))/2 + bb - aoff
#ifdef BLOCK_TRANSF
          call dfill(nnj,0.d0,g1,1)
          do j=jlo,jhi
            do i=ilo,ihi
              g1(i,j) = h2(j,i)
            enddo
          enddo
          call ga_acc(g_exch,jilo,jihi,ab,ab,g1,nnj,1.d0)
#else
#if SCATTER_TRANSF
          ij = 0
          do j=jlo,jhi
            jj = (j-1)*nbf
            do i=ilo,ihi
              ij = ij + 1
              h1(ij) = h2(j,i)
              iv(ij) = jj + i
              jv(ij) = ab
            enddo
          enddo
          call ga_dscatter(g_exch,h1,iv,jv,ijlen)
#else
          do j=jlo,jhi
            do i=ilo,ihi
              g1(i,j) = h2(j,i)
            enddo
          enddo
          do j=jlo,jhi
            jilo = (j-1)*nbf + ilo
            jihi = (j-1)*nbf + ihi
            call ga_acc(g_exch,jilo,jihi,ab,ab,g1(ilo,j),ilen,1.d0)
          enddo
#endif
#endif
#endif
        enddo
      enddo

 120  continue
      return
      end

