       subroutine moints_build_2xz(basis, ohalf,
     $                             occ_start, mo1_lo, mo1_hi,
     $                             mo2_lo, mo2_hi,
     $                             g_movecs,
     $                             g_coul, ocoul,
     $                             g_exch, oexch,
     $                             g_int1t,
     $                             blen )
       implicit none
#include "tcgmsg.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "schwarz.fh"
#include "msgids.fh"
c
c Arguments
c
       integer basis                          ! Basis handle
       logical ohalf                          ! Toggle alternative algorithm for half-ints intermediates
       integer occ_start                      ! Offset from frozen core
       integer mo1_lo, mo1_hi                 ! 1st Pair Index range
       integer mo2_lo, mo2_hi                 ! 2nd Pair Index range
       integer g_movecs                       ! MO coefficients
       integer g_coul                         ! Coulomb operator
       integer g_exch                         ! Exchange operator
       integer g_int1t                        ! Exchange operator with ao index
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
       integer kb, lb
       integer l_ssbb, k_ssbb, l_ssbbt, k_ssbbt
       integer l_hlp, k_hlp, l_ssni,k_ssni
       integer l_hlp2, k_hlp2
       integer l_eri, k_eri, l_iscr,k_iscr
       integer l_mo,k_mo
       integer l_shmap, k_shmap
       integer n_ssbb, n_ssbb1, n_ssni, n_ssai, n_hlp, n_hlp2
       integer mem2, max2e
       integer n_ijni
       integer num_nodes, ploop, next
       double precision tol2e, scale, schw_ij
       double precision tz, thalf, tint, t1idx, t2kidx, t3kidx, t4kidx
       double precision t34jidx, ttotal, t2jidx, tsynch
       double precision ttask, ttaskmax, ttaskmin, ttaskagg
       double precision flop1
       integer tottask
       logical status
       logical oprint
#ifdef MOINTS_STATS
#include "moints_stats.fh"
#endif
c
c
c
       integer moints_numgr, gr_len
       external moints_numgr, gr_len
c
c
c These are options that should be set
c by input deck
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
       n_hlp2 = maxbfsh*maxbfsh*nmo2
       n_ssni = maxbfsh*maxbfsh*nbf*nmo1
       n_ssai = maxbfsh*maxbfsh*nmo2*nmo1
       status = ma_push_get(MT_DBL,n_ssbb1,'ssbb block',l_ssbb,k_ssbb)
       status = ma_push_get(MT_DBL,n_ssbb,'ssbbt block',l_ssbbt,k_ssbbt)
       status = ma_push_get(MT_DBL,n_hlp,'hlp block',l_hlp,k_hlp)
       status = ma_push_get(MT_DBL,n_hlp2,'hlp block',l_hlp2,k_hlp2)
       status = ma_push_get(MT_DBL,n_ssni,'ssni block',l_ssni,k_ssni)
       if (.not.(status)) call errquit('cannot allocate local memory',0)
c
c Initialize
c
       flop1 = 0.d0
       tint = 0.d0
       t1idx = 0.d0
       t2jidx = 0.d0
       t34jidx = 0.d0
       t2kidx = 0.d0
       t3kidx = 0.d0
       t4kidx = 0.d0
       tottask = 0
       ttaskmax = 0.d0
       ttaskmin = 1.d24
       ttaskagg = 0.d0
       if (oexch) call ga_zero(g_exch)
       if (oexch) call ga_zero(g_int1t)
       if (ocoul) call ga_zero(g_coul)
       ploop = 0
       next = nxtval(num_nodes)
       thalf = tcgtime()
c
c  4-fold shell loop
c
       do ish0=1,nsh
         do jsh0=1,ish0
           ish = max(int_mb(k_shmap+ish0-1),int_mb(k_shmap+jsh0-1))
           jsh = min(int_mb(k_shmap+ish0-1),int_mb(k_shmap+jsh0-1))
           status = bas_cn2bfr(basis,ish,ibflo,ibfhi)
           ilen = ibfhi - ibflo + 1
           status = bas_cn2bfr(basis,jsh,jbflo,jbfhi)
           jlen = jbfhi - jbflo + 1
           schw_ij = schwarz_shell(ish,jsh)
           scale = 1.d0
           if (ish.eq.jsh) scale = scale*0.5d0
           if (schw_ij*schwarz_max().ge.tol2e) then
             if (next.eq.ploop) then
               tottask = tottask + 1
               ttask = tcgtime()
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
     $                               lshlo, lshhi, schw_ij, tol2e,
     $                               max2e, dbl_mb(k_eri), mem2,
     $                               dbl_mb(k_iscr), ilen, jlen,
     $                               kbflo, kbfhi, lbflo, lbfhi,
     $                               dbl_mb(k_ssbb) )
                   tint = tint + tcgtime() - tz
                   flop1 = flop1 + 4*lblen*ilen*jlen*kblen*nmo1
                   if (lb.ne.kb) then
                     tz = tcgtime()
                     call moints_blktr( ilen, jlen, kblen, lblen,
     $                                  dbl_mb(k_ssbb),
     $                                  dbl_mb(k_ssbbt))
                     call moints_trf1( nbf, mo1_lo, mo1_hi,
     $                                 ilen, jlen, kbflo, kbfhi,
     $                                 lbflo, lbfhi, 1.d0,
     $                                 dbl_mb(k_ssbb),
     $                                 dbl_mb(k_ssbbt), dbl_mb(k_mo),
     $                                 dbl_mb(k_ssni), dbl_mb(k_hlp) )

                   else
                     tz = tcgtime()
                     call moints_trf1( nbf, mo1_lo, mo1_hi,
     $                                 ilen, jlen, kbflo, kbfhi,
     $                                 lbflo, lbfhi, 0.5d0,
     $                                 dbl_mb(k_ssbb),
     $                                 dbl_mb(k_ssbb), dbl_mb(k_mo),
     $                                 dbl_mb(k_ssni), dbl_mb(k_hlp) )
                   endif
                   t1idx = t1idx + tcgtime() - tz
                 enddo
               enddo
               if (oexch) then
                 tz = tcgtime()
                 if (ohalf) then
                   call moints2x_trf2K( nbf, occ_start, mo1_lo, mo1_hi,
     $                                  ibflo, ibfhi, jbflo, jbfhi,
     $                                  scale, dbl_mb(k_ssni),
     $                                  dbl_mb(k_hlp), dbl_mb(k_hlp),
     $                                  dbl_mb(k_mo), g_exch )
                   t2kidx = t2kidx + tcgtime() - tz
                 else
                   call moints2xv_trf23K( nbf, occ_start,
     $                                  mo1_lo, mo1_hi, nmo2,
     $                                  mo2_lo, mo2_hi, ibflo, ibfhi,
     $                                  jbflo, jbfhi, scale,
     $                                  dbl_mb(k_ssni),dbl_mb(k_hlp2),
     $                                  dbl_mb(k_hlp), dbl_mb(k_hlp),
     $                                  dbl_mb(k_mo), g_exch )
                   t3kidx = t3kidx + tcgtime() - tz
                 endif
               endif
               if (ocoul) then
                 tz = tcgtime()
                 call moints_trf2J( nbf, occ_start, mo1_lo, mo1_hi,
     $                              ibflo, ibfhi, jbflo, jbfhi,
     $                              dbl_mb(k_ssni), dbl_mb(k_hlp),
     $                              dbl_mb(k_ssbb), dbl_mb(k_mo),
     $                              g_coul )
                 t2jidx = t2jidx + tcgtime() - tz
               endif
               next = nxtval(num_nodes)
               ttask = tcgtime() - ttask
               ttaskmax = max(ttaskmax,ttask)
               ttaskmin = min(ttaskmin,ttask)
               ttaskagg = ttaskagg + ttask
             endif
             ploop = ploop + 1
           endif
         enddo
       enddo
       tsynch = tcgtime()
       next = nxtval(-num_nodes)
       tsynch = tcgtime() - tsynch
       thalf = tcgtime() - thalf
c
c Clean-up
c
       if (.not. ma_pop_stack(l_ssni))
     $     call errquit('moints: failed to pop', l_ssni)
       if (.not. ma_pop_stack(l_hlp2))
     $     call errquit('moints: failed to pop', l_hlp2)
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
         if (ohalf) then
c
c Half transformed ints could be optionally dumped here
c or go ahead and transform the last 2 indices
c
c
c$$$          call moints_print_oper_z(occ_start, mo1_lo, mo1_hi,
c$$$     $                             nmo2, g_exch)
           call moints_Ktrf34z( g_exch, g_int1t, .true., occ_start,
     $                         mo1_lo, mo1_hi, mo2_lo, mo2_hi,
     $                         nbf, dbl_mb(k_mo), dbl_mb(k_hlp) )
         else
           call moints2xv_Ktrf4( g_exch, occ_start, mo1_lo, mo1_hi,
     $                           mo2_lo, mo2_hi, nbf, dbl_mb(k_mo),
     $                           dbl_mb(k_hlp) )
           t4kidx = t4kidx + tcgtime() - tz
c$$$          call moints_print_oper_z(occ_start, mo1_lo, mo1_hi,
c$$$     $                             nmo2, g_exch)
         endif
       endif
c
c
c
       if (ocoul) then
         tz = tcgtime()
         call moints_Jtrf34( g_coul, occ_start, mo1_lo, mo1_hi,
     $                       mo2_lo, mo2_hi, nbf, dbl_mb(k_mo),
     $                       dbl_mb(k_hlp) )
c$$$         call moints_print_oper_z(occ_start, mo1_lo, mo1_hi,
c$$$     $                            nmo2, g_coul)
         t34jidx = tcgtime() - tz
       endif
       call ga_sync()
c
c
c
       if (.not. ma_pop_stack(l_hlp))
     $     call errquit('moints: failed to pop', l_hlp)
       if (.not. ma_pop_stack(l_mo))
     $     call errquit('moints: failed to pop', l_mo)
c
c
c
#ifdef NOCOMMS
       if (ga_nodeid().eq.0) write(6,334)
 334   format(/,10x,'**** COMMUNICATION DISABLED ****')
#endif

c
c
c
#ifdef MOINTS_STATS
c
c Timings and Statistics bookkeeping
c
       flop1 = flop1*1.d-6
       ttotal = tcgtime() - ttotal
       mi_npass = mi_npass + 1.d0
       mi_ttotal = mi_ttotal + ttotal
       mi_thalf = mi_thalf + thalf
       mi_tint = mi_tint + tint
       mi_t1 = mi_t1 + t1idx
       mi_t2k = mi_t2k + t2kidx
       mi_t3k = mi_t3k + t3kidx
       mi_t4k = mi_t4k + t4kidx
       mi_t2j = mi_t2j + t2jidx
       mi_t34j = mi_t34j + t34jidx
       mi_flop1 = mi_flop1 + flop1
       mi_synch = mi_synch + tsynch
       mi_maxsynch = max(mi_maxsynch,tsynch)
       mi_minsynch = min(mi_minsynch,tsynch)
       mi_aggsynch = mi_aggsynch + tsynch
       mi_nsynchs = mi_nsynchs + 1
       mi_maxtaskJ = max(mi_maxtaskJ,ttaskmax)
       mi_mintaskJ = min(mi_mintaskJ,ttaskmin)
       mi_aggtaskJ = mi_aggtaskJ + ttaskagg
       mi_ntasksJ = mi_ntasksJ + tottask
       mi_ntasksK = 1.d0
       mi_mintaskK = 0.d0
#endif

       return
       end

c
c
c
c
c
c
c
c


      subroutine moints_Ktrf34z( g_exch, g_int1t, otrp, ostart, olo,
     $                           ohi, alo, ahi, nbf, c, tmp )
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer g_exch, g_int1t
      logical otrp
      integer nbf, ostart, olo, ohi, alo, ahi
      double precision c(nbf,nbf)
      double precision tmp(nbf,nbf)

      integer nact, clo, chi, rlo, rhi, my_id
      integer ofroz, a1, a2, aa1, aa2, aa, ab, ba, aoff, k_local, ld

      ofroz = ostart - 1
      nact = ahi - alo + 1
      my_id = ga_nodeid()
      call ga_distribution(g_exch, my_id, rlo, rhi, clo, chi )
      aoff = ((olo-ofroz)*(olo-ofroz-1))/2
      do a1=olo,ohi
        aa1 = a1 - ofroz
        do a2=ostart,a1
          aa2 = a2 - ofroz
          aa = (aa1*(aa1-1))/2 + aa2 - aoff
crk note aoff may not yet work for core for ab and ba
          ab = (aa1-1)*(ohi-olo+1) + aa2
          ba = (aa2-1)*(ohi-olo+1) + aa1
          if ((aa.ge.clo).and.(aa.le.chi)) then
            call ga_access(g_exch,rlo,rhi,aa,aa,k_local,ld)
            call moints_Ktrf34az(g_int1t,nbf,nact,ostart,c,
     &                           dbl_mb(k_local),tmp,ab,ba)
            call ga_release(g_exch,rlo,rhi,aa,aa)
          endif
        enddo
      enddo
      return
      end





c
c
c Note: May require transpose the half-transformed matrix
c       to conform to convention for storage
c
      subroutine moints_Ktrf34az(g_int1t,nbf,nact,ostart,c,x,tmp,ab,ba)
      implicit none
      integer g_int1t
      integer a,b,ab,ba
      integer nbf,nact,ostart
      double precision c(nbf,nbf)
      double precision x(nbf,nbf)
      double precision tmp(nbf,nbf)
c
c     pause here to push out one ao index integrals to global array
c
      call dgemm('t','n',nbf,nbf,nbf,1.d0,c,nbf,x,nbf,
     $           0.d0,tmp,nbf)
      do a=1,nbf
        call ga_put(g_int1t,(ba-1)*nact+1,(ba-1)*nact+nact,a,a,
     &              tmp(ostart,a),nact)
      enddo
      call dgemm('t','t',nbf,nbf,nbf,1.d0,c,nbf,x,nbf,
     $           0.d0,tmp,nbf)
      do a=1,nbf
        call ga_put(g_int1t,(ab-1)*nact+1,(ab-1)*nact+nact,a,a,
     &              tmp(ostart,a),nact)
      enddo
c
      call dgemm('n','n',nbf,nbf,nbf,1.d0,tmp,nbf,c,nbf,
     $           0.d0,x,nbf)

      return
      end







