c
c
c  $Id: moints_2x.F,v 1.12 1995-10-27 23:25:12 rg240 Exp $
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
       subroutine moints_build_2x(basis, ohalf, osym,
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
#include "sym.fh"
c
c  Arguments
c
       integer basis                          ! Basis handle
       logical ohalf                          ! Toggle alternative algorithm for half-ints intermediates
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
c  Local variables
c
       integer nmo1, nmo2, nbf, nsh, maxbfsh
       integer nblk, bsize, blocksize
       integer l_gmap, k_gmap
       integer ish0, jsh0, ish, jsh, ilen, jlen
       integer ibflo,ibfhi,jbflo,jbfhi,kbflo,kbfhi,lbflo,lbfhi
       integer kshlo, kshhi, lshlo, lshhi
       integer kblen, lblen
       integer kb, lb, i
       integer l_ssbb, k_ssbb, l_ssbbt, k_ssbbt
       integer l_hlp, k_hlp, l_ssni,k_ssni
       integer l_hlp2, k_hlp2
       integer l_eri, k_eri, l_iscr,k_iscr
       integer l_mo, k_mo
       integer l_shmap, k_shmap
       integer k_sym, l_sym
       integer n_ssbb, n_ssbb1, n_ssni, n_ssai, n_hlp, n_hlp2, n_ijni
       integer mem2, max2e
       integer num_nodes, ploop, next, itmp
       integer nmixed
       double precision tol2e, scale, schw_ij
       double precision tz, thalf, tint, t1idx, t2kidx, t3kidx, t4kidx
       double precision t34jidx, ttotal, t2jidx, tsynch
       double precision ttask, ttaskmax, ttaskmin, ttaskagg
       double precision flop1, q2
       integer tottask
       integer geom
       logical status, oprint
       logical odoit
c
#include "moints_stats.fh"
c
       integer moints_numgr, gr_len
       external moints_numgr, gr_len
c
c  These are options that should be set
c  by input deck
c
       data tol2e/1.d-12/
       data oprint/.true./
c
c
c  General basis info
c
       ttotal = tcgtime()
       num_nodes = ga_nnodes()
       if (.not. bas_geom(basis, geom)) call errquit
     $      ('moints: cannot get geometry', 0)
       status = bas_numbf(basis,nbf)
       status = status.and.bas_numcont(basis,nsh)
       status = status.and.bas_nbf_cn_max(basis,maxbfsh)
       if (.not.status) call errquit('moints: cannot get basis info',0)
       nmo1 = mo1_hi - occ_start + 1
       nmo2 = mo2_hi - mo2_lo + 1
c     
c  Symmetry adapt the MOs and renumber irreps to start at zero
c     
       if (osym) then
          if (.not. ma_push_get(MT_INT, nbf, 'movec syms',l_sym, k_sym))
     $         call errquit('moints_2x: no memory for syms?',0)
          call sym_movecs_adapt(basis, 1d-8, g_movecs, int_mb(k_sym),
     $                          nmixed)
          if (nmixed .ne. 0)
     $         call errquit('moints_2x: symmetry contamination', nmixed)
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
c  Initialize
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
           odoit = schw_ij*schwarz_max().ge.tol2e
           if (odoit .and. osym) then
              odoit = sym_shell_pair(basis, ish, jsh, q2)
           endif
           if (odoit) then
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
     $                               dbl_mb(k_ssbb), osym )
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
c  Clean-up
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
c  2nd half transformation
c
       status = ma_push_get(MT_DBL,(nbf*nbf),'hlp',l_hlp,k_hlp)
c
c  Exchange term
c
       if (oexch) then
         tz = tcgtime()
         if (ohalf) then
c$$$           if (util_print('exchange half integral',print_debug))
c$$$     $       call moints_op_print(occ_start,mo1_lo,mo1_hi,nbf,g_exch)
           call moints_Ktrf34( g_exch, .true., occ_start,
     $                         mo1_lo, mo1_hi, mo2_lo, mo2_hi,
     $                         nbf, dbl_mb(k_mo), dbl_mb(k_hlp),
     $                         osym, int_mb(k_sym))
           itmp = nbf
         else
           call moints2xv_Ktrf4( g_exch, occ_start, mo1_lo, mo1_hi,
     $                           mo2_lo, mo2_hi, nbf, dbl_mb(k_mo),
     $                           dbl_mb(k_hlp), osym, int_mb(k_sym) )
           t4kidx = t4kidx + tcgtime() - tz
           itmp = nmo2
         endif
         if (util_print('exchange integral',print_debug))
     $      call moints_op_print(occ_start,mo1_lo,mo1_hi,itmp,g_exch)
       endif
c
c  Coulomb term
c
       if (ocoul) then
         tz = tcgtime()
c$$$         if (util_print('coulomb half integral',print_debug))
c$$$     $      call moints_op_print(occ_start,mo1_lo,mo1_hi,nbf,g_coul)
         call moints_Jtrf34( g_coul, occ_start, mo1_lo, mo1_hi,
     $                       mo2_lo, mo2_hi, nbf, dbl_mb(k_mo),
     $                       dbl_mb(k_hlp), osym, int_mb(k_sym) )

         t34jidx = tcgtime() - tz
         itmp = nmo2
         if (ohalf) itmp = nbf
         if (util_print('coulomb integral',print_debug))
     $      call moints_op_print(occ_start,mo1_lo,mo1_hi,itmp,g_coul)
       endif
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
#ifdef NOCOMMS
       if (ga_nodeid().eq.0) write(6,334)
 334   format(/,10x,'**** COMMUNICATION DISABLED ****')
#endif

c
c  Timings and Statistics bookkeeping
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
c
c
c
       return
       end









c
c  Exchange 2nd Index transformation to occupied
c  Used when half-transformed integrals
c  are required.
c
       subroutine moints2x_trf2K( nbf, ostart, olo, ohi, ilo, ihi,
     $                           jlo, jhi, scale, ssni, h1, h2,
     $                           c, g_exch )
       implicit none
       integer nbf, ostart, olo, ohi, ilo, ihi, jlo, jhi
       double precision scale
       double precision ssni(nbf,jlo:jhi,ilo:ihi,olo:ohi)
       double precision h1(nbf,ilo:ihi), h2(nbf,jlo:jhi)
       double precision c(nbf,nbf)
       integer g_exch
c
       integer nni,nnj,ab,ijlo,ijhi,ofroz
       integer aoff,a,b,i,j,iz
c
       ofroz = ostart - 1
       aoff = ((olo-ofroz)*(olo-ofroz-1))/2
       nni = (ihi - ilo + 1)*nbf
       nnj = (jhi - jlo + 1)*nbf
       do a=olo,ohi
         do b=ostart,a
           call dfill(nni,0.d0,h1,1)
           do i=ilo,ihi
             do j=jlo,jhi
               do iz=1,nbf
                 h1(iz,i) = h1(iz,i) + scale*c(j,b)*ssni(iz,j,i,a)
               enddo
             enddo
           enddo
#ifndef NOCOMMS
           ab = ((a-ofroz)*(a-ofroz-1))/2 + b - ofroz - aoff
           ijlo = (ilo-1)*nbf + 1
           ijhi = ihi*nbf
           call ga_acc(g_exch,ijlo,ijhi,ab,ab,h1,nni,1.d0)
#endif
         enddo
       enddo
c
c
       do a=olo,ohi
         do b=ostart,a
           call dfill(nnj,0.d0,h2,1)
           do j=jlo,jhi
             do i=ilo,ihi
               do iz=1,nbf
                 h2(iz,j) = h2(iz,j) + scale*c(i,b)*ssni(iz,j,i,a)
               enddo
             enddo
           enddo
#ifndef NOCOMMS
           ab = ((a-ofroz)*(a-ofroz-1))/2 + b - ofroz - aoff
           ijlo = (jlo-1)*nbf + 1
           ijhi = jhi*nbf
           call ga_acc(g_exch,ijlo,ijhi,ab,ab,h2,nnj,1.d0)
#endif
         enddo
       enddo
       return
       end







c
c  Exchange 2nd index & 3rd index transformations
c  Note that the 2nd index is transformed to virtuals
c  initially.
c
       subroutine moints2xv_trf23K( nbf, ostart, olo, ohi, nvir,
     $                              vlo, vhi, ilo, ihi, jlo, jhi,
     $                              scale, ssni, ssv, h1, h2, c,
     $                              g_exch )
       implicit none
       integer nbf, ostart, ohi, olo, nvir, vlo, vhi
       integer ilo, ihi, jlo, jhi
       double precision scale
       double precision ssni(nbf,jlo:jhi,ilo:ihi,olo:ohi)
       double precision ssv(nvir,jlo:jhi,ilo:ihi)
       double precision h1(nvir,ilo:ihi), h2(nvir,jlo:jhi)
       double precision c(nbf,nbf)
       integer g_exch

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
         do b=ostart,a
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
           aa = a - ofroz
           bb = b - ofroz
           ab = (aa*(aa-1))/2 + bb - aoff
           call ga_acc(g_exch,ivlo,ivhi,ab,ab,h1,nni,1.d0)
#endif
         enddo
         do b=ostart,a
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
          aa = a - ofroz
          bb = b - ofroz
          ab = (aa*(aa-1))/2 + bb - aoff
          call ga_acc(g_exch,jvlo,jvhi,ab,ab,h2,nnj,1.d0)
#endif
         enddo
       enddo
c
c
       return
       end







c
c  Transform last index of K to virtual
c
      subroutine moints2xv_Ktrf4( g_exch, ostart, olo, ohi,
     $                           vlo, vhi, nbf, c, tmp, osym, syms )
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer g_exch
      integer nbf, ostart, olo, ohi, vlo, vhi
      double precision c(nbf,nbf)
      double precision tmp(nbf,nbf)
      logical osym
      integer syms(nbf)

      integer nvir, clo, chi, rlo, rhi, my_id
      integer ofroz, a1, a2, aa1, aa2, aa, aoff, k_local, ld

      ofroz = ostart - 1
      nvir = vhi - vlo + 1
      my_id = ga_nodeid()
      call ga_distribution(g_exch, my_id, rlo, rhi, clo, chi )
      aoff = ((olo-ofroz)*(olo-ofroz-1))/2
      do a1=olo,ohi
        aa1 = a1 - ofroz
        do a2=ostart,a1
          aa2 = a2 - ofroz
          aa = (aa1*(aa1-1))/2 + aa2 - aoff
          if ((aa.ge.clo).and.(aa.le.chi)) then
            call ga_access(g_exch,rlo,rhi,aa,aa,k_local,ld)
            call moints2xv_Ktrf4a( vlo, vhi, nbf, c,
     $                            dbl_mb(k_local), tmp )
            if (osym) call moints_sym_zero( nbf, vlo, vhi, syms,
     $                                      a1, a2, tmp)
            call dcopy((nvir*nvir),tmp,1,dbl_mb(k_local),1)
            call ga_release(g_exch,rlo,rhi,aa,aa)
          endif
        enddo
      enddo
      return
      end




      subroutine moints2xv_Ktrf4a( vlo, vhi, nbf, c, x, tmp )
      implicit none
      integer vlo, vhi, nbf
      double precision c(nbf,nbf)
      double precision x(vlo:vhi,nbf)
      double precision tmp(vlo:vhi,vlo:vhi)
      integer nvir
      
      nvir = vhi - vlo + 1
      call dgemm('t', 't', nvir, nvir, nbf, 1.d0, c(1,vlo), nbf,
     $           x, nvir, 0.d0, tmp, nvir)
      return
      end

















