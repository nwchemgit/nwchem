c
c
c  $Id: moints_trp.F,v 1.3 1995-07-27 17:39:01 rg240 Exp $
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
c
c
c

       subroutine moints_trp(basis, ohalf,
     $                            occ_start, mo1_lo, mo1_hi,
     $                            mo2_lo, mo2_hi,
     $                            g_coul, ocoul,
     $                            g_exch, oexch,
     $                            blen, cmo, ncor, ndel, nocc, nvir)
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
capr       integer g_movecs                       ! MO coefficients
       integer g_coul                         ! Coulomb operator
       integer g_exch                         ! Exchange operator
       logical ocoul,oexch                    ! Type selection
       integer blen                           ! Blocking length
       integer ncor, ndel, nocc, nvir
       double precision cmo(*)
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
       double precision half
       data half/0.5d00/
c
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
capr       status = status.and.ma_push_get(MT_DBL,(nbf*nbf),
capr     $                                 'movecs cols',l_mo,k_mo)
capr       call ga_get(g_movecs,1,nbf,1,nbf,dbl_mb(k_mo),nbf)
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
       if (ocoul) call ga_zero(g_coul)
       ploop = 0
       next = nxtval(num_nodes)
       thalf = tcgtime()
c
c  4-fold shell loop
c
       do ish0=1,nsh
         ish = int_mb(k_shmap+ish0-1)
         status = bas_cn2bfr(basis,ish,ibflo,ibfhi)
         ilen = ibfhi - ibflo + 1
         do jsh0=1,ish0
           jsh = int_mb(k_shmap+jsh0-1)
           status = bas_cn2bfr(basis,jsh,jbflo,jbfhi)
           jlen = jbfhi - jbflo + 1
           schw_ij = schwarz_shell(ish,jsh)
           scale = 1.d0
           if (ish.eq.jsh) scale = scale*half
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
     $                                 dbl_mb(k_ssbbt), cmo,
     $                                 dbl_mb(k_ssni), dbl_mb(k_hlp) )

                   else
                     tz = tcgtime()
                     call moints_trf1( nbf, mo1_lo, mo1_hi,
     $                                 ilen, jlen, kbflo, kbfhi,
     $                                 lbflo, lbfhi, half,
     $                                 dbl_mb(k_ssbb),
     $                                 dbl_mb(k_ssbb), cmo,
     $                                 dbl_mb(k_ssni), dbl_mb(k_hlp) )
                   endif
                   t1idx = t1idx + tcgtime() - tz
                 enddo
               enddo
               if (oexch) then
                 tz = tcgtime()
                   call mo_trp_trf23K( nbf, occ_start,
     $                                  mo1_lo, mo1_hi, nmo2,
     $                                  mo2_lo, mo2_hi, ibflo, ibfhi,
     $                                  jbflo, jbfhi, scale,
     $                                  dbl_mb(k_ssni),dbl_mb(k_hlp2),
     $                                  dbl_mb(k_hlp), dbl_mb(k_hlp),
     $                                  cmo, g_exch, ncor, nocc)
                   t3kidx = t3kidx + tcgtime() - tz
               endif
               if (ocoul) then
                 tz = tcgtime()
                 call mo_trp_trf2J( nbf, occ_start, mo2_lo, mo2_hi,
     $                              ibflo, ibfhi, jbflo, jbfhi,
     $                              dbl_mb(k_ssni), dbl_mb(k_hlp),
     $                              dbl_mb(k_ssbb), cmo,
     $                              g_coul, ncor, nocc)
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
           call mo_trp_Ktrf4( g_exch, occ_start, mo1_lo, mo1_hi,
     $                           mo2_lo, mo2_hi, nbf, cmo,
     $                           dbl_mb(k_hlp), ncor, nocc)
           t4kidx = t4kidx + tcgtime() - tz
       endif
c
c
c
       if (ocoul) then
         tz = tcgtime()
         call mo_trp_Jtrf34( g_coul, occ_start, mo1_lo, mo1_hi,
     $                       mo2_lo, mo2_hi, nbf, cmo,
     $                       dbl_mb(k_hlp), ncor, nocc, nvir)
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
capr       if (.not. ma_pop_stack(l_mo))
capr     $     call errquit('moints: failed to pop', l_mo)
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
       subroutine mo_trp_trf23K( nbf, ostart, olo, ohi, nvir,
     $                              vlo, vhi, ilo, ihi, jlo, jhi,
     $                              scale, ssni, ssv, h1, h2, c,
     $                              g_exch, ncor, nocc )
       implicit none
       integer nbf, ostart, ohi, olo, nvir, vlo, vhi
       integer ilo, ihi, jlo, jhi
       double precision scale
       double precision ssni(nbf,jlo:jhi,ilo:ihi,olo:ohi)
       double precision ssv(nvir,jlo:jhi,ilo:ihi)
       double precision h1(nvir,ilo:ihi), h2(nvir,jlo:jhi)
       double precision c(nbf,nbf)
       integer g_exch,ncor,nocc

c
       integer aoff,nni,nnij,nnj,ab,ivlo,ivhi,jvlo,jvhi,ofroz
       integer ssij
       integer aa,bb,a,b,i,j,v
       double precision cc
c
       ofroz = ostart - 1
       ssij = (ihi - ilo + 1)*(jhi - jlo + 1)
       nni = (ihi - ilo + 1)*nvir
       nnij = (jhi - jlo + 1)*nni
       ivlo = (ilo-1)*nvir + 1
       ivhi = ihi*nvir
       aoff = ((olo-ofroz)*(olo-ofroz-1))/2
       nnj = (jhi - jlo + 1)*nvir
       jvlo = (jlo-1)*nvir + 1
       jvhi = jhi*nvir

       do a=olo,ohi
         call dgemm('t', 'n', nvir, ssij, nbf, 1.d0,
     $               c(1,vlo), nbf, ssni(1,jlo,ilo,a),
     $               nbf, 0.d0, ssv, nvir )
         do b=ostart,ostart+nocc-1
capr         do b=ostart,a
           call dfill(nni,0.d0,h1,1)
           do j=jlo,jhi
             cc = c(j,b)
             do i=ilo,ihi
               do v=1,nvir
                 h1(v,i) = h1(v,i) + scale*cc*ssv(v,j,i)
               enddo
             enddo
           enddo
#ifndef NOCOMMS
capr           aa = a - ofroz
capr           bb = b - ofroz
capr           ab = (aa*(aa-1))/2 + bb - aoff
           ab = (a-olo)*nocc+b-ostart+1
           call ga_acc(g_exch,ivlo,ivhi,ab,ab,h1,nni,1.d0)
#endif
         enddo
         do b=ostart,ostart+nocc-1
capr         do b=ostart,a
           call dfill(nnj,0.d0,h2,1)
           do i=ilo,ihi
            cc = c(i,b)
            do j=jlo,jhi
              do v=1,nvir
                h2(v,j) = h2(v,j) + scale*cc*ssv(v,j,i)
              enddo
            enddo
          enddo
#ifndef NOCOMMS
capr          aa = a - ofroz
capr          bb = b - ofroz
capr          ab = (aa*(aa-1))/2 + bb - aoff
          ab = (a-olo)*nocc+b-ostart+1
          call ga_acc(g_exch,jvlo,jvhi,ab,ab,h2,nnj,1.d0)
#endif
         enddo
       enddo
c
c
       return
       end
      subroutine mo_trp_Ktrf4( g_exch, ostart, olo, ohi,
     $                           vlo, vhi, nbf, c, tmp ,ncor, nocc)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer g_exch
      integer nbf, ostart, olo, ohi, vlo, vhi,ncor,nocc
      double precision c(nbf,nbf)
      double precision tmp(nbf,nbf)

      integer nvir, clo, chi, rlo, rhi, my_id
      integer ofroz, a1, a2, aa1, aa2, aa, aoff, k_local, ld

      ofroz = ostart - 1
      nvir = vhi - vlo + 1
      my_id = ga_nodeid()
      call ga_distribution(g_exch, my_id, rlo, rhi, clo, chi )
      aoff = ((olo-ofroz)*(olo-ofroz-1))/2
      do a1=olo,ohi
        aa1 = a1 - ofroz
capr        do a2=ostart,a1
        do a2=ostart,ostart+nocc-1
capr          aa2 = a2 - ofroz
capr          aa = (aa1*(aa1-1))/2 + aa2 - aoff
          aa = (a1-olo)*nocc+a2-ostart+1
          if ((aa.ge.clo).and.(aa.le.chi)) then
            call ga_access(g_exch,rlo,rhi,aa,aa,k_local,ld)
            call moints2xv_Ktrf4a( vlo, vhi, nbf, c,
     $                            dbl_mb(k_local), tmp )
            call dcopy((nvir*nvir),tmp,1,dbl_mb(k_local),1)
            call ga_release(g_exch,rlo,rhi,aa,aa)
          endif
        enddo
      enddo
      return
      end
      subroutine mo_trp_trf2J( nbf, ostart, olo, ohi, ilo, ihi,
     $                         jlo, jhi, ssni, h1, h2, c, g_coul,
     &                         ncor, nocc )
      implicit none
      integer nbf, ostart, olo, ohi, ilo, ihi, jlo, jhi
      double precision ssni(nbf,jlo:jhi,ilo:ihi,olo:ohi)
      double precision h1(nbf,ilo:ihi)
      double precision h2(jlo:jhi,ilo:ihi)
      double precision c(nbf,nbf)
      integer g_coul,ncor,nocc
c
      integer nni, ijlo, ijhi, ilen, jlen, ijlen,lenij
      integer ab, aa, bb, a, b, i, j, aoff, ofroz, jtop
c
c
      ofroz = ostart - 1
      aoff = ((olo-ofroz)*(olo-ofroz-1))/2
      ilen = ihi - ilo + 1
      jlen = jhi - jlo + 1
      nni = ilen*nbf
      ijlen = ilen*jlen
#ifdef BLOCK_TRANSF
      ijlo = (ilo-1)*nbf + 1
      ijhi = ihi*nbf
#endif
      do a=olo,ohi
capr        do b=ostart,a
        do b=ostart,ostart+nocc-1
          call dgemm('t','n',ijlen,1,nbf,1.d0,ssni(1,jlo,ilo,a),
     $               nbf,c(1,b),nbf,0.d0,h2,ijlen)
#ifndef NOCOMMS
capr          aa = a - ofroz
capr          bb = b - ofroz
capr          ab = (aa*(aa-1))/2 + bb - aoff
          ab = (a-olo)*nocc+b-ostart+1
#ifdef BLOCK_TRANSF
          call dfill(nni,0.d0,h1,1)
          do i=ilo,ihi
            jtop = jhi
            if (jhi.eq.ihi) jtop = i
            do j=jlo,jtop
              h1(j,i) = h2(j,i)
            enddo
          enddo
          call ga_acc(g_coul,ijlo,ijhi,ab,ab,h1,nni,1.d0)
#else
CRK       do i=ilo,ihi
CRK         ijlo = i*(i-1)/2 + jlo
CRK         ijhi = i*(i-1)/2 + jhi
CRK         if (jhi.eq.ihi) ijhi = i*(i-1)/2 + i
CRK         lenij=ijhi-ijlo+1
CRK         write(6,*)'ijlo,ijhi',ijlo,ijhi
CRK         write(6,*)'ab,h2(jlo,i)',ab,jlo,i,(h2(jlo+j,i),j=0,lenij-1)
CRK         call ga_acc(g_coul,ijlo,ijhi,ab,ab,h2(jlo,i),1,1.d0)
CRK       enddo
          if(jlo.le.ilo)then
            do i=ilo,ihi
            do j=jlo,min(i,jhi)
              ijlo = i*(i-1)/2 + j
              ijhi = i*(i-1)/2 + j
              lenij=ijhi-ijlo+1
              call ga_acc(g_coul,ijlo,ijhi,ab,ab,h2(j,i),1,1.d0)
            enddo
            enddo
          else
            do j=jlo,jhi
            do i=ilo,min(j,ihi)
              ijlo = j*(j-1)/2 + i
              ijhi = j*(j-1)/2 + i
              lenij=ijhi-ijlo+1
              call ga_acc(g_coul,ijlo,ijhi,ab,ab,h2(j,i),1,1.d0)
            enddo
            enddo
          endif

CRK2      do i=ilo,ihi
CRK2      do j=jlo,min(i,jhi)
CRK2        ijlo = (i-1)*nbf + jlo
CRK2        ijhi = (i-1)*nbf + jhi
CRK2        if (jhi.eq.ihi) ijhi = (i-1)*nbf + i
CRK2        write(6,*)'ab,h2(j,i)',ab,j,i,h2(j,i)
CRK2        call ga_acc(g_coul,ijlo,ijhi,ab,ab,h2(j,i),1,1.d0)
CRK2      enddo
CRK2      enddo
#endif
#endif
        enddo
      enddo
      return
      end
      subroutine mo_trp_Jtrf34( g_coul, ostart, olo, ohi, vlo, vhi,
     $                          nbf, c, tmp, ncor, nocc, nvir)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer g_coul
      integer ostart, olo, ohi, vlo, vhi, nbf, ncor, nocc, nvir
      double precision c(nbf,nbf)
      double precision tmp(nbf,nbf)

      integer clo, chi, rlo, rhi, my_id
      integer ofroz, aoff, a1, a2, aa1, aa2, aa, k_local, ld

      ofroz = ostart - 1
      aoff = ((olo-ofroz)*(olo-ofroz-1))/2
      my_id = ga_nodeid()
      call ga_distribution(g_coul, my_id, rlo, rhi, clo, chi )
      do a1=olo,ohi
capr        aa1 = a1 - ofroz
        do a2=ostart,ostart+nocc-1
capr        do a2=ostart,a1
capr          aa2 = a2 - ofroz
capr          aa = (aa1*(aa1-1))/2 + aa2 - aoff
          aa = (a1-olo)*nocc+a2-ostart+1
          if ((aa.ge.clo).and.(aa.le.chi)) then
            call ga_access(g_coul,rlo,rhi,aa,aa,k_local,ld)
CRK         do aa1=0,nbf*(nbf+1)/2-1
CRK            write(6,*)'read x',aa,dbl_mb(k_local+aa1)
CRK         enddo
            call mo_trp_Jtrf34a( nbf, c, dbl_mb(k_local), tmp, 
     &                           ncor,nocc,nvir)
CRK         do aa1=0,nvir*nvir-1
CRK            write(6,*)'x',dbl_mb(k_local+aa1)
CRK         enddo
            call ga_release(g_coul,rlo,rhi,aa,aa)
          endif
        enddo
      enddo
      return
      end
      subroutine mo_trp_Jtrf34a( nbf, c, x, tmp, ncor, nocc, nvir)
      implicit none
      integer nbf,ncor,nocc,nvir
      double precision c(nbf,nbf)
      double precision x(nbf*nbf)
      double precision tmp(nbf,nbf)

      integer i, j, icnt

      icnt=0
      do i=1,nbf
        do j=1,i
          icnt=icnt+1
          tmp(j,i) = x(icnt)
          tmp(i,j) = x(icnt)
CRK       write(6,*)'tmp',tmp(i,j)
        enddo
      enddo

      call dgemm('t','n',nvir,nbf,nbf,1.d0,c(1,ncor+nocc+1),nbf,tmp,nbf,
     $           0.d0,x,nvir)
      call dgemm('n','n',nvir,nvir,nbf,1.d0,x,nvir,c(1,ncor+nocc+1),nbf,
     $           0.d0,tmp,nvir)

      call dcopy(nvir*nvir,tmp,1,x,1)

      return
      end

