c
c
c  $Id: moints_semi.F,v 1.1 1995-10-17 16:39:09 rg240 Exp $
c
c  This routine writes or returns
c  half-transformed integrals
c
c
c
       subroutine moints_semi( basis, tol2e, ocoul, oexch,
     $                         olo, ohi, vlo, vhi, g_movecs )
       implicit none
#include "tcgmsg.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "bas.fh"
c
       integer basis                          ! [input]  Basis handle
       double precision tol2e                 ! [input]  Integral tolerance
       logical ocoul                          ! [input]  Toggle Coulomb integrals
       logical oexch                          ! [input]  Toggle exchange integrals
       integer olo, ohi                       ! [input]  Occupied Index range for this segment
       integer vlo, vhi                       ! [input]  Virtual Index range
       integer g_movecs                       ! [input]  MO coefficients
c
       integer junit, kunit
       parameter(junit=10,kunit=11)
c
       integer l_gmap, k_gmap
       integer l_shmap, k_shmap
       integer l_gloc, k_gloc
       integer l_iopass, k_iopass
       integer l_rlen, k_rlen
       integer l_rloc, k_rloc
       integer nbf, nsh, maxbfsh, nocc, nvir
       integer ngrp, blksiz, nnbf
       integer balo, bahi, bilo, bihi
       integer myid, numnodes, rlo, rhi, clo, chi
       integer mxrlen, itmp, ierr, vtmp(100)
       integer lmemreq, gmem, len
       integer g_jbuf
       integer g_kbuf
       logical status
       double precision tsemi
       character*256 fname
       data blksiz/16/                         ! Affects performance only
c
       integer moints_numgr
       integer moints_lmem
       external moints_numgr
       external moints_lmem
c
c  General info
c
       tsemi = tcgtime()
       myid = ga_nodeid()
       numnodes = ga_nnodes()
       nocc = ohi - olo + 1
       nvir = vhi - vlo + 1
       status = bas_numbf(basis,nbf)
       status = status.and.bas_numcont(basis,nsh)
       status = status.and.bas_nbf_cn_max(basis,maxbfsh)
       if (.not.status) call errquit('moints: cannot get basis info',0)
c
c  Memory arithmetic
c
       lmemreq = moints_lmem(basis, nocc, nvir, ocoul, oexch, blksiz)
       if ((.not.ga_uses_ma()).and.(.not.(ga_memory_limited())))
     $   call errquit('cannot determine memory limit',0)
       gmem = ga_memory_avail() - numnodes*lmemreq
       itmp = 0
       if (ocoul) itmp = nocc
       if (oexch) itmp = itmp + nvir
       mxrlen = gmem/(8*nocc*itmp)
       mxrlen = min(mxrlen,(nbf*nbf))
       MXRLEN = 512
c
       if (ga_nodeid().eq.0) WRITE(6,900) LMEMREQ, GMEM, MXRLEN
 900   FORMAT(/,'Local Memory required:  ',i10,' bytes '
     $        /,'Global Memory remaining:',i10,' bytes ',
     $        /,'IO Buffer length:       ',i10)
c
c  Create GA coulomb buffer
c
       if (ocoul) then
         itmp = (nocc/numnodes + min(mod(nocc,numnodes),1))*nocc
         if (.not.(ga_create( MT_DBL, mxrlen, nocc*nocc, 'transp Jbuf',
     $                        mxrlen, nocc, g_jbuf )))
     $           call errquit('moints_semi: cannot allocate J buffer')
         call ga_distribution(g_jbuf, myid, rlo, rhi, clo, chi )
         if ((clo.eq.0).and.(chi.eq.-1)) then
           bilo = 0
           bihi = -1
         else
           bilo = olo - (clo-1)/nocc
           bihi = olo - (chi/nocc) - 1
           if ((mod(chi,nocc).ne.0).or.(mod(clo,nocc).ne.1))
     $       call errquit('moints_semi:wrong distrib. for J transp',0)
         endif
       endif
c
c  Create GA exchange buffer
c
       if (oexch) then
         itmp = (nvir/numnodes + min(mod(nvir,numnodes),1))*nocc
         if (.not.(ga_create( MT_DBL, mxrlen, nocc*nvir, 'transp Kbuf',
     $                        mxrlen, itmp, g_kbuf )))
     $           call errquit('moints_semi: cannot allocate K buffer')
         call ga_distribution(g_kbuf, myid, rlo, rhi, clo, chi )
         if ((clo.eq.0).and.(chi.eq.-1)) then
           balo = 0
           bahi = -1
         else
           balo = vlo + (clo-1)/nocc
           bahi = vlo + (chi/nocc) - 1
           if ((mod(chi,nocc).ne.0).or.(mod(clo,nocc).ne.1))
     $       call errquit('moints_semi:wrong distrib. for K transp',0)
         endif
         WRITE(6,882) GA_NODEID(),CLO,CHI,ITMP,BALO,BAHI
 882     FORMAT('ME:',I5,'  CLO:CHI ',2I5,' YCHUNK:',i5,' VRANGE:',2I5)
         CALL UTIL_FLUSH(6)
         CALL GA_SYNC()
       endif
c
c  Reorder shells by descending shell-length
c  and group shells by blocksize
c
       status = ma_push_get(MT_INT,nsh,'sh map',l_shmap, k_shmap)
       call moints_shorder( basis, nsh, maxbfsh, int_mb(k_shmap) )
       ngrp = moints_numgr( basis, blksiz )
       status = ma_push_get(MT_INT,(ngrp*4),'grp map', l_gmap, k_gmap)
       call moints_grmap( basis, blksiz, ngrp, int_mb(k_gmap))
c
c  Generate locator and reverse map
c
       itmp = (nsh*(nsh+1))/2+numnodes
       status = ma_push_get(MT_INT, (nbf*nbf), 'loc', l_gloc, k_gloc )
       status = ma_push_get(MT_INT, itmp, 'io', l_iopass, k_iopass )
       status = ma_push_get(MT_INT, itmp, 'rec len', l_rlen, k_rlen )
       itmp = (nbf*(nbf+1))/2
       itmp = itmp + mod(itmp,2)
       status = ma_push_get(MT_INT, itmp, 'rev loc', l_rloc, k_rloc )
       call ifill(itmp, 0, int_mb(k_rloc), 1 )
       call moints_locmap( basis, nsh, nbf, tol2e, int_mb(k_shmap),
     $                     g_kbuf, int_mb(k_gloc), nnbf, int_mb(k_rloc),
     $                     int_mb(k_iopass), int_mb(k_rlen))
c
c  Open local file and write initial info
c
       call ifill(10, 0, vtmp, 1 )
       itmp = nnbf+mod(nnbf,2)
       vtmp(1) = ma_sizeof(MT_INT,(10 + itmp),MT_DBL) + 1
       vtmp(2) = nnbf
       if (oexch) then
         call moints_localfname( 'KHalf', fname, len )          
         call wopen(kunit, fname, len, 1, 0, ierr )
         vtmp(3) = balo
         vtmp(4) = bahi
         call putwa( kunit, vtmp, 1,
     $               ma_sizeof(MT_INT, 10, MT_DBL), ierr )
         call putwa( kunit, int_mb(k_rloc),
     $               ma_sizeof(MT_INT, 12, MT_DBL),
     $               ma_sizeof(MT_INT, itmp, MT_DBL), ierr )
       endif
       if (ocoul) then
         call moints_localfname( 'JHalf', fname, len )          
         call wopen(junit, fname, len, 1, 0, ierr )
         vtmp(3) = bilo
         vtmp(4) = bihi
         call putwa( junit, vtmp, 1,
     $               ma_sizeof(MT_INT, 10, MT_DBL), ierr )
         call putwa( junit, int_mb(k_rloc),
     $               ma_sizeof(MT_INT, 12, MT_DBL),
     $               ma_sizeof(MT_INT, itmp, MT_DBL), ierr )
       endif
       if (.not. ma_pop_stack(l_rloc))
     $     call errquit('moints: failed to pop', l_rloc)
c
c  Call the real stuff
c
       call moints_semi_a( basis, nbf, nsh, maxbfsh, tol2e,
     $                     olo, ohi, vlo, vhi, ocoul, oexch,
     $                     int_mb(k_shmap), blksiz, ngrp,
     $                     int_mb(k_gmap), int_mb(k_gloc),
     $                     nnbf, int_mb(k_iopass), int_mb(k_rlen),
     $                     g_movecs, g_jbuf, g_kbuf )
c
c  Clean up
c
       if (.not. ma_pop_stack(l_rlen))
     $     call errquit('moints: failed to pop', l_rlen)
       if (.not. ma_pop_stack(l_iopass))
     $     call errquit('moints: failed to pop', l_iopass)
       if (.not. ma_pop_stack(l_gloc))
     $     call errquit('moints: failed to pop', l_gloc)
       if (.not. ma_pop_stack(l_gmap))
     $     call errquit('moints: failed to pop', l_gmap)
       if (.not. ma_pop_stack(l_shmap))
     $     call errquit('moints: failed to pop', l_shmap)
       if ((ocoul).and.(.not. ga_destroy(g_jbuf)))
     $     call errquit('moints: failed to destroy g_jbuf',g_jbuf)
       if ((oexch).and.(.not. ga_destroy(g_kbuf)))
     $     call errquit('moints: failed to destroy g_kbuf',g_kbuf)
c
       if (oexch) call wclose(kunit, ierr )
       if (ocoul) call wclose(junit, ierr )
c
       call ga_sync()
       tsemi = tcgtime() - tsemi
       if (ga_nodeid().eq.0) write(6,881) tsemi
 881   format(/,'Time for semi-direct half-transformation:',f12.3)
       return
       end
       








c
c
c
c
c
       subroutine moints_semi_a( basis, nbf, nsh, maxbfsh, tol2e,
     $                           olo, ohi, vlo, vhi, ocoul, oexch,
     $                           shmap, blksiz, ngrp, gmap,
     $                           gloc, nnbf, iopass, rlen,
     $                           g_movecs, g_jbuf, g_kbuf )
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
       integer nbf                            ! Basis functions
       integer nsh                            ! Shells
       integer maxbfsh                        ! Largest shell
       double precision tol2e                 ! Integral tolerance
       integer olo, ohi                       ! 1st Pair Index range
       integer vlo, vhi                       ! 2nd Pair Index range
       logical ocoul                          ! Toggle Coulomb
       logical oexch                          ! Toggle exchange
       integer shmap(nsh)                     ! Map for shells
       integer blksiz                         ! Blocksize
       integer ngrp                           ! Number of groups of shells
       integer gmap(*)                        ! Map of shell-group
       integer gloc(nbf*nbf)                  ! bf -> GA memory map
       integer nnbf                           ! Number of screened basf pairs <= (nbf*(nbf+1))/2
       integer iopass(0:*)                    ! IO pass for each task
       integer rlen(*)                        ! Record lengths for IO writes
       integer g_movecs                       ! MO coefficients
       integer g_jbuf                         ! Buffer for half-trans Coulomb
       integer g_kbuf                         ! Buffer for half-trans exchange
c
       integer junit, kunit
       parameter(junit=10,kunit=11)
c
c Local variables
c
       integer nmo1, nmo2
       integer ish0, jsh0, ish, jsh, ilen, jlen
       integer ibflo, ibfhi, jbflo, jbfhi
       integer kbflo, kbfhi, lbflo, lbfhi
       integer kshlo, kshhi, lshlo, lshhi
       integer kblen, lblen, kb, lb
       integer l_ssbb, k_ssbb, l_ssbbt, k_ssbbt
       integer l_hlp, k_hlp, l_ssni,k_ssni
       integer l_hlp2, k_hlp2
       integer l_eri, k_eri, l_iscr,k_iscr
       integer l_mo, k_mo
       integer n_ssbb, n_ssbb1, n_ssni, n_hlp, n_hlp2
       integer bsize
       integer mem2, max2e, jopass, iocnt
       integer num_nodes, ploop, next
       double precision tol2e, scale, schw_ij
       double precision tpass
       logical status
c
       integer gr_len
       external gr_len
c
c
c
       num_nodes = ga_nnodes()
       nmo1 = ohi - olo + 1
       nmo2 = vhi - vlo + 1
c
c  Local MO coefficients
c
       status = ma_push_get(MT_DBL,(nbf*nbf),'movecs cols',l_mo,k_mo)
       call ga_get(g_movecs,1,nbf,1,nbf,dbl_mb(k_mo),nbf)
c
c  Integrals allocation
c
       call int_mem_2e4c(max2e, mem2)
       status = ma_push_get(MT_DBL, max2e,'moints: buf', l_eri, k_eri)
       status = ma_push_get(MT_DBL, mem2, 'moints: scr', l_iscr, k_iscr)
c
c  Temporary partially-transformed arrays
c     
       bsize = max(blksiz,maxbfsh)
       n_ssbb = maxbfsh*maxbfsh*bsize*bsize
       n_ssbb1 = max((nmo1*nmo1),n_ssbb)
       n_hlp = max((bsize*maxbfsh*maxbfsh*nmo1),(maxbfsh*nbf))
       n_hlp2 = maxbfsh*maxbfsh*nmo2
       n_ssni = maxbfsh*maxbfsh*nbf*nmo1
       status = ma_push_get(MT_DBL,n_ssbb1,'ssbb block',l_ssbb,k_ssbb)
       status = ma_push_get(MT_DBL,n_ssbb,'ssbbt block',l_ssbbt,k_ssbbt)
       status = ma_push_get(MT_DBL,n_hlp,'hlp block',l_hlp,k_hlp)
       status = ma_push_get(MT_DBL,n_hlp2,'hlp block',l_hlp2,k_hlp2)
       status = ma_push_get(MT_DBL,n_ssni,'ssni block',l_ssni,k_ssni)
       if (.not.(status)) call errquit('cannot allocate local memory',0)
c
c Initialize
c
       ploop = 0
       jopass = 1
       iocnt = ma_sizeof(MT_INT,(10+nnbf+mod(nnbf,2)),MT_DBL) + 1    ! start counting from initial segment
       next = nxtval(num_nodes)
       tpass = tcgtime()
c
c  4-fold shell loop
c
       do ish0=1,nsh
         do jsh0=1,ish0
           ish = max(shmap(ish0),shmap(jsh0))
           jsh = min(shmap(ish0),shmap(jsh0))
           status = bas_cn2bfr(basis,ish,ibflo,ibfhi)
           status = bas_cn2bfr(basis,jsh,jbflo,jbfhi)
           ilen = ibfhi - ibflo + 1
           jlen = jbfhi - jbflo + 1
           schw_ij = schwarz_shell(ish,jsh)
           scale = 1.d0
           if (ish.eq.jsh) scale = scale*0.5d0
           if (schw_ij*schwarz_max().ge.tol2e) then
             if (next.eq.ploop) then
c
c  -------------
c  Parallel task
c
c
c  Checkpoint for IO
c
               if (iopass(next).gt.jopass) then
                 CALL UTIL_FLUSH(6)
                 call ga_sync()
                 if (ocoul) then
	           call moints_wrbuf( junit, olo, ohi, olo, ohi, nnbf,
     $                                iocnt, rlen(jopass), g_jbuf )
                   call ga_zero(g_jbuf)
                 endif
                 if (oexch) then
                   call moints_wrbuf( kunit, olo, ohi, vlo, vhi, nnbf,
     $                                iocnt, rlen(jopass), g_kbuf )
                   call ga_zero(g_kbuf)
                 endif
                 tpass = tcgtime() - tpass
                 write(6,901) ga_nodeid(), next, iopass(next), jopass,
     $                        rlen(jopass), tpass
 901             format('Me:',i5,5x,'Next:',i5,5x,'IOPass(next):',i5,5x,
     $                  'Pass:',i5,5x,'Length:',i6, 5x,'Time:',f12.3)
	         CALL UTIL_FLUSH(6)
                 iocnt = iocnt + rlen(jopass)
                 jopass = jopass + 1
                 tpass = tcgtime()
               endif
c
c  Half-tranformed Integral generation
c
               call dfill((ilen*jlen*nbf*nmo1),0.d0,dbl_mb(k_ssni),1)
               do kb=1,ngrp
                 call gr_sh_idx(ngrp,gmap,kb,kshlo,kshhi)
                 kblen = gr_len(ngrp,gmap,kb,kbflo,kbfhi)
                 do lb=1,kb
                   call gr_sh_idx(ngrp,gmap,lb,lshlo,lshhi)
                   lblen = gr_len(ngrp,gmap,lb,lbflo,lbfhi)
                   call moints_gblk( basis, ish, jsh, kshlo, kshhi,
     $                               lshlo, lshhi, schw_ij, tol2e,
     $                               max2e, dbl_mb(k_eri), mem2,
     $                               dbl_mb(k_iscr), ilen, jlen,
     $                               kbflo, kbfhi, lbflo, lbfhi,
     $                               dbl_mb(k_ssbb) )
                   if (lb.ne.kb) then
                     call moints_blktr( ilen, jlen, kblen, lblen,
     $                                  dbl_mb(k_ssbb),
     $                                  dbl_mb(k_ssbbt))
                     call moints_trf1( nbf, olo, ohi,
     $                                 ilen, jlen, kbflo, kbfhi,
     $                                 lbflo, lbfhi, 1.d0,
     $                                 dbl_mb(k_ssbb),
     $                                 dbl_mb(k_ssbbt), dbl_mb(k_mo),
     $                                 dbl_mb(k_ssni), dbl_mb(k_hlp) )

                   else
                     call moints_trf1( nbf, olo, ohi,
     $                                 ilen, jlen, kbflo, kbfhi,
     $                                 lbflo, lbfhi, 0.5d0,
     $                                 dbl_mb(k_ssbb),
     $                                 dbl_mb(k_ssbb), dbl_mb(k_mo),
     $                                 dbl_mb(k_ssni), dbl_mb(k_hlp) )
                   endif
                 enddo
               enddo
               if (ocoul)
     $           call moints_trf2Ky( nbf, olo, ohi, olo, ohi,
     $                               ibflo, ibfhi, jbflo, jbfhi,
     $                               gloc, dbl_mb(k_ssni),
     $                               dbl_mb(k_hlp), dbl_mb(k_hlp2),
     $                               dbl_mb(k_mo), g_jbuf )
               if (oexch)
     $           call moints_trf2Ky( nbf, olo, ohi, vlo, vhi,
     $                               ibflo, ibfhi, jbflo, jbfhi,
     $                               gloc, dbl_mb(k_ssni),
     $                               dbl_mb(k_hlp), dbl_mb(k_hlp2),
     $                               dbl_mb(k_mo), g_kbuf )
               next = nxtval(num_nodes)
c
c  End parallel task
c  -----------------
c
             endif
             ploop = ploop + 1
           endif
         enddo
       enddo
       call ga_sync()
       if (ocoul) then
         call moints_wrbuf( junit, olo, ohi, olo, ohi, nnbf,
     $                      iocnt, rlen(jopass), g_jbuf )
       endif
       if (oexch) then
         call moints_wrbuf( kunit, olo, ohi, vlo, vhi, nnbf,
     $                      iocnt, rlen(jopass), g_kbuf )
       endif
       tpass = tcgtime() - tpass
       write(6,901) ga_nodeid(), next, iopass(next), jopass,
     $              rlen(jopass), tpass
       CALL UTIL_FLUSH(6)
       next = nxtval(-num_nodes)
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
       if (.not. ma_pop_stack(l_iscr))
     $     call errquit('moints: failed to pop', l_iscr)
       if (.not. ma_pop_stack(l_eri))
     $     call errquit('moints: failed to pop', l_eri)
c
c 
c  Complete in-core section
c
c$$$       status = ma_push_get(MT_DBL,(nbf*nbf),'hlp',l_hlp,k_hlp)
c$$$       status = ma_push_get(MT_DBL,(nbf*nbf),'hlp',l_hlp2,k_hlp2)
c$$$       call moints_trf34Ky( nbf, ostart, olo, ohi,
c$$$     $                      vlo, vhi, nnbf, int_mb(k_rloc),
c$$$     $                      dbl_mb(k_hlp), dbl_mb(k_hlp2),
c$$$     $                      dbl_mb(k_mo), g_buf )
c$$$       if (.not. ma_pop_stack(l_hlp2))
c$$$     $     call errquit('moints: failed to pop', l_hlp2)
c$$$       if (.not. ma_pop_stack(l_hlp))
c$$$     $     call errquit('moints: failed to pop', l_hlp)
c
c
c
       if (.not. ma_pop_stack(l_mo))
     $     call errquit('moints: failed to pop', l_mo)
c
c
c
       return
       end







c
c 
c
c
c
      subroutine moints_trf2Ky( nbf, olo, ohi, vlo, vhi,
     $                          ilo, ihi, jlo, jhi, gloc,
     $                          ssni, h, h2, c, g_buf )
      implicit none
#include "global.fh"
      integer nbf, olo, ohi, vlo, vhi
      integer ilo, ihi, jlo, jhi
      integer gloc(nbf,nbf)
      double precision ssni(nbf,jlo:jhi,ilo:ihi,olo:ohi)
      double precision h(jlo:jhi,ilo:ihi),h2(*)
      double precision c(nbf,nbf)
      integer g_buf
      integer ilen, jlen, ijlen
      integer glo, ghi, i, j, ij
      integer o, v, vo
c
      ilen = ihi - ilo + 1
      jlen = jhi - jlo + 1
      ijlen = ilen*jlen
      glo = gloc(ilo,jlo)
      ghi = gloc(ihi,jhi)
c$$$      WRITE(6,881) ILO,IHI,JLO,JHI,GLO,GHI,GA_NODEID()
c$$$ 881  FORMAT('&&',2(2I4,2X), 'RANGE:',I4,1X,'-',1X,I4,10x,'[',i3,']')

      do o=olo,ohi
        do v=vlo,vhi
          vo = (v-vlo)*(ohi-olo+1) + o
          call dgemm('t','n',ijlen,1,nbf,1.d0,ssni(1,jlo,ilo,o),
     $                nbf,c(1,v),nbf,0.d0,h,ijlen)
          if (ihi.ne.jhi) then
            call ga_put(g_buf,glo,ghi,vo,vo,h,1)
          else
            ij = 0
            do i=ilo,ihi
              do j=ilo,i
                ij =   ij + 1
                h2(ij) = h(i,j) 
              enddo
            enddo
            call ga_put(g_buf,glo,ghi,vo,vo,h2,1)
          endif
        enddo
      enddo
      return
      end






c
c  *** For completely in-core 4-index ***
c

      subroutine moints_trf34Ky( nbf, ostart, olo, ohi, vlo, vhi,
     $                            nnbf, rloc, h1, h2, c, g_buf )
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer nbf, ostart, ohi, olo, vlo, vhi
      integer nnbf, rloc(nnbf)
      double precision h1(*), h2(*)
      double precision c(nbf,nbf)
      integer g_buf
c     
      integer nocc, nvir, v, o, vo
      integer k_local, ld, my_id, rlo, rhi, clo, chi
      
      nocc = ohi - olo + 1
      nvir = vhi - vlo + 1
      my_id = ga_nodeid()
      call ga_distribution(g_buf, my_id, rlo, rhi, clo, chi )
      do v=vlo,vhi
        do o=olo,ohi
          vo = (v-1)*(ohi-olo+1) + o
          if ((vo.ge.clo).and.(vo.le.chi)) then
            call dfill( (nbf*nbf), 0.d0, h1, 1 )
            call ga_access(g_buf, rlo, rhi, vo, vo, k_local, ld )
            call scatter(nnbf, h1, rloc, dbl_mb(k_local) )
            call ga_release(g_buf, rlo, rhi, vo, vo )
            call upper2square(nbf,h1,h1)
            
            call dgemm( 'n', 'n', nbf, nocc, nbf, 1.d0, h1, nbf,
     $                   c(1,olo), nbf, 0.d0, h2, nbf )
            call dgemm( 't', 'n', nvir, nocc, nbf, 1.d0, c(1,vlo), nbf,
     $                   h2, nbf, 0.d0, h1, nvir )

C$$$             WRITE(6,911) V,O
C$$$ 911         FORMAT(//,5X,'V :',I5,5X,'O :',I5)
C$$$             CALL MOINTS_MATPRINT(NVIR,NOCC,H1)

          endif
        enddo
      enddo
      return
      end





      subroutine moints_wrbuf( munit, olo, ohi, plo, phi, nnbf,
     $                         iocnt, rlen, g_buf )
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer munit
      integer olo, ohi, plo, phi
      integer nnbf, iocnt, rlen
      integer g_buf
c
      integer faddr, ierr
      integer o, p, po, nocc, np
      integer myid, rlo, rhi, clo, chi, k_local, ld
c
      nocc = ohi - olo + 1
      np = phi - plo + 1
      myid = ga_nodeid()
      call ga_distribution(g_buf, myid, rlo, rhi, clo, chi )
c$$$      WRITE(6,922) GA_NODEID()
c$$$ 922  FORMAT(' $$$$ WR_BUF:',I5)
c$$$      CALL UTIL_FLUSH(6)
      do o=olo,ohi
        do p=plo,phi
          po = (p-plo)*nocc + o - olo + 1
          faddr = ((o-olo)*np + p - plo)*nnbf +  iocnt
          if ((po.ge.clo).and.(po.le.chi)) then
            call ga_access(g_buf, rlo, rhi, po, po, k_local, ld )
            call putwa( munit, dbl_mb(k_local), faddr, rlen, ierr )
            call ga_release(g_buf, rlo, rhi, po, po )
c$$$            WRITE(6,901) O,P,FADDR
c$$$ 901        FORMAT(' [O,P]:',2I4,5X,'ADDR:',I6)
          endif
        enddo
      enddo
      return
      end
          





c
c  Driver routine to illustrate how to read back
c  half-transformed integrals from local node file.
c
      subroutine moints_readintJ( nbf, olo, ohi, vlo, vhi, c, orbe )
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer nbf, olo, ohi, vlo, vhi
      double precision c(nbf,nbf)
      double precision orbe(nbf)
c
      integer junit
      parameter(junit=10)
c
      integer nocc, nnbf, ioff, p, o
      integer bihi, bilo, itmp
      integer l_v, k_v, l_t, k_t, l_x, k_x
      logical status

      call wopen(junit, 'Jhalfints', 8, 1, 0, itmp )
      PRINT*,'Start readintJ ',itmp
      nocc = ohi - olo + 1
      call moints_vrange( junit, bilo, bihi, nnbf, ioff )
      status = ma_push_get(MT_INT,(nnbf+mod(nnbf,2)),'scatter',l_v,k_v)
      call moints_getscattv( junit, nnbf, int_mb(k_v) )
c
c  Loop over occupied/virtual pairs
c
      itmp = max(nnbf,(nbf*nocc))
      status = ma_push_get(MT_DBL,itmp,'read buff',l_t,k_t)
      status = ma_push_get(MT_DBL,(nbf*nbf),'U half',l_x,k_x)
      do o=olo,ohi
        do p=bilo,bihi
          call dfill((nbf*nbf), 0.d0, dbl_mb(k_x), 1 )
          call moints_rdhfint( junit, olo, ohi, o, p, nnbf,
     $                         ioff, dbl_mb(k_t))
          call scatter( nnbf, dbl_mb(k_x), int_mb(k_v), dbl_mb(k_t) )
          call upper2square( nbf, dbl_mb(k_x), dbl_mb(k_x) )

c$$$          WRITE(6,900) O,P
c$$$ 900      FORMAT(5X,' J READ [O,P]:',2I5)
c$$$          CALL MOINTS_MATPRINT(NBF,NBF,DBL_MB(K_X))

          call dgemm( 'n', 'n', nbf, nocc, nbf, 1.d0, dbl_mb(k_x), nbf,
     $                c(1,olo), nbf, 0.d0, dbl_mb(k_t), nbf )
          call dgemm( 't', 'n', nocc, nocc, nbf, 1.d0, c(1,olo), nbf,
     $                dbl_mb(k_t), nbf, 0.d0, dbl_mb(k_x), nocc )
        enddo
      enddo
c
c  Clean up
c
      if (.not. ma_pop_stack(l_x))
     $  call errquit('moints_readint: failed to pop', l_x)
      if (.not. ma_pop_stack(l_t))
     $  call errquit('moints_readint: failed to pop', l_t)
      if (.not. ma_pop_stack(l_v))
     $  call errquit('moints_readint: failed to pop', l_v)
c
      call wclose(junit, itmp )
      return
      end











      subroutine moints_readintK( nbf, olo, ohi, vlo, vhi, c, orbe )
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
      integer MSG_SEMIMP2_SUM
      parameter(MSG_SEMIMP2_SUM=10241)
      integer nbf, olo, ohi, vlo, vhi
      double precision c(nbf,nbf)
      double precision orbe(nbf)
c
      integer kunit
      parameter(kunit=11)
c
      integer nocc, nvir, nnbf, ioff, v, o, balo, bahi, itmp
      integer l_v, k_v, l_t, k_t, l_x, k_x
      integer j, jj, vvlo, vvhi, k_local, ld, myid
      integer rlo, rhi, clo, chi, len
      integer g_exch
      double precision tmp2
      INTEGER G_PAIR, IJ, NNOCC
      character*255 fname
      double precision xx, denom, e2
      logical status
      double precision moints_epair_eval
      external moints_epair_eval

      tmp2 = tcgtime()
      call moints_localfname( 'KHalf', fname, len )          
      call wopen(kunit, fname, len, 1, 0, itmp )
c$$$      PRINT*,'Start readintK',itmp
      myid = ga_nodeid()
      e2 = 0.d0
      nocc = ohi - olo + 1
      nvir = vhi - vlo + 1
      call moints_vrange( kunit, balo, bahi, nnbf, ioff )
      status = ma_push_get(MT_INT,(nnbf+mod(nnbf,2)),'scatter',l_v,k_v)
      call moints_getscattv( kunit, nnbf, int_mb(k_v) )
      if (.not.ga_create( MT_DBL, (nvir*nvir), nocc, 'exch',
     $                    (nvir*nvir), 1, g_exch))
     $    call errquit('moints_semimp2: cannot allocate exch',0)
      call ga_distribution(g_exch, myid, rlo, rhi, clo, chi )

      NNOCC = (NOCC*(NOCC+1))/2
      STATUS = GA_CREATE( MT_DBL, NNOCC, 1, 'pair', 1, 1, G_PAIR)
c
c$$$      WRITE(6,934) GA_NODEID(),BALO,BAHI,CLO,CHI
c$$$ 934  FORMAT('ME:',I5,5x,'VRANGE:',2I5,5x,'OCC RANGE FOR K:',2I5)
c$$$      CALL UTIL_FLUSH(6)
c
c  Loop over occupied/virtual pairs
c
      itmp = max(nnbf,(nbf*nocc))
      status = ma_push_get(MT_DBL,itmp,'read buff',l_t,k_t)
      status = ma_push_get(MT_DBL,(nbf*nbf),'U half',l_x,k_x)
      do o=olo,ohi
c$$$        do v=balo,bahi
        call ga_zero(g_exch)
        do v=vlo,vhi
          if ((v.ge.balo).and.(v.le.bahi)) then
            call dfill((nbf*nbf), 0.d0, dbl_mb(k_x), 1 )
            call moints_rdhfint( kunit, vlo, vhi, o, v, nnbf,
     $                           ioff, dbl_mb(k_t))
            call scatter( nnbf, dbl_mb(k_x), int_mb(k_v), dbl_mb(k_t) )
            call upper2square( nbf, dbl_mb(k_x), dbl_mb(k_x) )

c$$$            WRITE(6,900) O,V,BALO,BAHI
c$$$ 900        FORMAT(5X,' READ [O,V]:',2I5,5X,'VRANGE:',2I5)
c$$$            CALL MOINTS_MATPRINT(NBF,NBF,DBL_MB(K_X))
c$$$            CALL UTIL_FLUSH(6)

            call dgemm( 'n', 'n', nbf, nocc, nbf, 1.d0, dbl_mb(k_x),
     $                 nbf, c(1,olo), nbf, 0.d0, dbl_mb(k_t), nbf )
            call dgemm( 't', 'n', nvir, nocc, nbf, 1.d0, c(1,vlo),
     $                 nbf, dbl_mb(k_t), nbf, 0.d0, dbl_mb(k_x), nvir )

            vvlo = (v-vlo)*nvir + 1
            vvhi = (v-vlo+1)*nvir
            do j=1,nocc
               jj = (j-1)*nvir
               call ga_put( g_exch, vvlo, vvhi, j, j,
     $                      dbl_mb(k_x+jj), nvir)
            enddo
          endif
          CALL GA_SYNC()
        enddo
        call ga_sync()
        do j=1,o
          if ((j.ge.clo).and.(j.le.chi)) then
            denom = orbe(o) + orbe(j+olo-1)
            call ga_access(g_exch, rlo, rhi, j, j, k_local, ld )

c$$$            WRITE(6,920) O,J
c$$$ 920        FORMAT(' Operator: [',i4,',',i4,']')
c$$$            CALL MOINTS_MATPRINT(NVIR,NVIR,DBL_MB(K_LOCAL))
c$$$            CALL UTIL_FLUSH(6)

            xx = moints_epair_eval( nvir, 0, nocc, nvir,
     $                              dbl_mb(k_local), orbe, denom )
            call ga_release(g_exch, rlo, rhi, j, j )
            if (o.eq.j) xx = xx*0.5d0
            e2 = e2 + xx
	    
            IJ = (O*(O-1))/2 + J
            CALL GA_PUT(G_PAIR,IJ,IJ,1,1,XX,1)

c$$$            WRITE(6,922) GA_NODEID(),O,J,XX
c$$$ 922        FORMAT(I3,' %%%%%% ',2I5,5X,F16.10)
c$$$	    CALL UTIL_FLUSH(6)

          endif
c$$$          CALL GA_SYNC()
        enddo
      enddo
c
c
c$$$      CALL GA_SYNC()
c$$$      WRITE(6,967) GA_NODEID(), E2
c$$$ 967  FORMAT('ME:',I5,'  MP2 Contribution:',f16.10)
c$$$      CALL UTIL_FLUSH(6)
c$$$      CALL GA_PRINT(G_PAIR)
      call ga_sync()
      call ga_dgop(MSG_SEMIMP2_SUM,e2,1,'+')
      if (ga_nodeid().eq.0) then
        write(6,681) e2
 681    format(/,10x,'MP2 correction:',5x,f16.10)
        call util_flush(6)
      endif
c
c  Clean up
c
      STATUS = GA_DESTROY(G_PAIR)
      if (.not.ga_destroy(g_exch))
     $  call errquit('moints_semimp2: failed to destroy',g_exch)
      if (.not. ma_pop_stack(l_x))
     $  call errquit('moints_readint: failed to pop', l_x)
      if (.not. ma_pop_stack(l_t))
     $  call errquit('moints_readint: failed to pop', l_t)
      if (.not. ma_pop_stack(l_v))
     $  call errquit('moints_readint: failed to pop', l_v)
c
      call wclose(kunit, itmp )
c
      call ga_sync()
      tmp2 = tcgtime() - tmp2
      if (ga_nodeid().eq.0) write(6,881) tmp2
 881  format(/,'Time for semi-direct energy:',f12.3)
      
      return
      end




c
c
c =========================================================================
c =========================================================================
c
c  Driver routine to illustrate how to read back
c  half-transformed integrals from local node file.
c
c =========================================================================
c =========================================================================
c
c
#ifdef UNUSED
      subroutine moints_readintK2( nbf, olo, ohi, vlo, vhi, c, orbe )
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer nbf, olo, ohi, vlo, vhi
      double precision c(nbf,nbf)
      double precision orbe(nbf)
c
      integer kunit
      parameter(kunit=11)
c
      integer nocc, nvir, nnbf, ioff, v, o, balo, bahi, itmp
      integer l_v, k_v, l_t, k_t, l_x, k_x
      integer j, jj, vv, vvlo, vvhi, k_local, ld, myid, nnvo
      integer rlo, rhi, clo, chi
      integer g_exch
      character*256 fname
      double precision xx, denom, e2, sp1fact, sp2fact
      logical status
      double precision moints_epair_eval
      external moints_epair_eval
      data sp1fact,sp2fact/4.d0,-2.d0/
      
      call moints_localfname( 'KHalf', fname, len )          
      call wopen(kunit, fname, len, 1, 0, itmp )
      e2 = 0.d0
      nocc = ohi - olo + 1
      nvir = vhi - vlo + 1
      nnvo = nocc*nvir
      call moints_vrange( kunit, balo, bahi, nnbf, ioff )
      status = ma_push_get(MT_INT,(nnbf+mod(nnbf,2)),'scatter',l_v,k_v)
      call moints_getscattv( kunit, nnbf, int_mb(k_v) )
      if (.not.ga_create( MT_DBL, (nvir*nocc), nvir, 'exch',
     $                    (nvir*nocc), 1, g_exch))
     $    call errquit('moints_semimp2: cannot allocate exch',0)
      call ga_distribution(g_exch, myid, rlo, rhi, clo, chi )
c
c  Loop over occupied/virtual pairs
c
      itmp = max(nnbf,(nbf*nocc))
      status = ma_push_get(MT_DBL,itmp,'read buff',l_t,k_t)
      status = ma_push_get(MT_DBL,(nbf*nbf),'U half',l_x,k_x)
      do o=olo,ohi
        call ga_zero(g_exch)
        do v=balo,bahi
          call dfill((nbf*nbf), 0.d0, dbl_mb(k_x), 1 )
          call moints_rdhfint( kunit, vlo, vhi, o, v, nnbf,
     $                         ioff, dbl_mb(k_t))
          call scatter( nnbf, dbl_mb(k_x), int_mb(k_v), dbl_mb(k_t) )
          call upper2square( nbf, dbl_mb(k_x), dbl_mb(k_x) )

          call dgemm( 'n', 'n', nbf, nocc, nbf, 1.d0, dbl_mb(k_x), nbf,
     $                c(1,olo), nbf, 0.d0, dbl_mb(k_t), nbf )
          call dgemm( 't', 'n', nocc, nvir, nbf, 1.d0, dbl_mb(k_t), nbf,
     $                c(1,vlo), nbf, 0.d0, dbl_mb(k_x), nocc )

          vv = v - vlo + 1
          call ga_acc(g_exch,1,nnvo,vv,vv,dbl_mb(k_x),nnvo,sp1fact)
          rlo = (vv-1)*nocc + 1
          rhi = vv*nocc
          call ga_acc(g_exch,rlo,rhi,1,nvir,dbl_mb(k_x),nocc,sp2fact)

c$$$          WRITE(6,900) O,V
c$$$ 900      FORMAT(5X,' READ [O,V]:',2I5)
c$$$          CALL MOINTS_MATPRINT(NVIR,NOCC,DBL_MB(K_X))

        enddo
        call ga_sync()
        call ga_print(g_exch)
c$$$        do j=1,o
c$$$          if ((j.ge.clo).and.(j.le.chi)) then
c$$$            denom = orbe(o) + orbe(j+olo-1)
c$$$            call ga_access(g_exch, rlo, rhi, j, j, k_local, ld )
c$$$            xx = moints_epair_eval( nvir, 0, nocc, nvir,
c$$$     $                              dbl_mb(k_local), orbe, denom )
c$$$c$$$       PRINT*,'I,J: ',O,J,'   DENOM: ',denom,'    XX: ',XX,'  LD:',ld
c$$$c$$$       CALL MOINTS_MATPRINT(NVIR,NVIR,DBL_MB(K_LOCAL))
c$$$            if (o.eq.j) xx = xx*0.5d0
c$$$            e2 = e2 + xx
c$$$          endif
c$$$        enddo
      enddo
c$$$      PRINT*,' MP2 correction:',e2
c
c  Clean up
c
      if (.not.ga_destroy(g_exch))
     $  call errquit('moints_semimp2: failed to destroy',g_exch)
      if (.not. ma_pop_stack(l_x))
     $  call errquit('moints_readint: failed to pop', l_x)
      if (.not. ma_pop_stack(l_t))
     $  call errquit('moints_readint: failed to pop', l_t)
      if (.not. ma_pop_stack(l_v))
     $  call errquit('moints_readint: failed to pop', l_v)
c
      call wclose(kunit, itmp )
      return
      end
#endif      /* UNUSED */






      subroutine moints_rdhfint( munit, vlo, vhi, o, v, nnbf, ioff, t )
      implicit none
      integer munit
      integer vlo, vhi, o, v, nnbf, ioff
      double precision t(nnbf)
      integer faddr, nvir, ierr

      nvir = vhi - vlo + 1
      faddr = ((o-1)*nvir + v - vlo)*nnbf + ioff
      call getwa( munit, t, faddr, nnbf, ierr )
      
      return
      end






      subroutine moints_vrange( munit, balo, bahi, nnbf, ioff )
      implicit none
#include "mafdecls.fh"
      integer munit
      integer balo, bahi, nnbf, ioff
      integer vtmp(10), i, ierr

      call getwa( munit, vtmp, 1,
     $            ma_sizeof(MT_INT, 10, MT_DBL), ierr )
      ioff = vtmp(1)
      nnbf = vtmp(2)
      balo = vtmp(3)
      bahi = vtmp(4)

c$$$      PRINT*,'READING BACK INITIAL PARAMETERS'
c$$$      WRITE(6,771) (VTMP(I),I=1,10)
c$$$ 771  FORMAT(16I4)

      return
      end


      subroutine moints_getscattv( munit, nnbf, v )
      implicit none
#include "mafdecls.fh"
      integer munit
      integer nnbf, v(*)
      integer itmp, ierr, i
      
      itmp = nnbf+mod(nnbf,2)
      call getwa( munit, v,
     $            ma_sizeof(MT_INT, 12, MT_DBL),
     $            ma_sizeof(MT_INT, itmp, MT_DBL), ierr )

c$$$      call getwa(munit, v, 6, itmp, ierr )
c$$$      PRINT*,'READ BACK SCATTER ARRAY'
c$$$      WRITE(6,771) (V(I),I=1,NNBF)
c$$$ 771  FORMAT(16I4)

      return
      end







c
c  Units in bytes 
c  Local memory requirements are independent of (ocoul, oexch)?
c

       integer function moints_lmem( basis, nocc, nvir,
     $                               ocoul, oexch, blksiz )
 193   implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
       integer basis, nocc, nvir, blksiz
       logical ocoul, oexch
c
       integer nsh, nbf, maxbfsh
       integer bsize, ngrp, imax2e, imem2, memi, memd
       logical status
       integer moints_numgr
       external moints_numgr

       status = bas_numbf(basis,nbf)
       status = status.and.bas_numcont(basis,nsh)
       status = status.and.bas_nbf_cn_max(basis,maxbfsh)
       bsize = max(blksiz,maxbfsh)
       ngrp = moints_numgr( basis, blksiz )
       call int_mem_2e4c(imax2e, imem2)
       memi = nsh + 4*ngrp + nbf*nbf + nsh*(nsh+1) +
     $        2*ga_nnodes() + (nbf*(nbf+1))/2
       memd = nbf*nbf + imax2e + imem2 + 
     $        maxbfsh*maxbfsh*bsize*bsize +
     $        max((nocc*nocc),(maxbfsh*maxbfsh*bsize*bsize)) +
     $        max((bsize*maxbfsh*maxbfsh*nocc),(maxbfsh*nbf)) +
     $        maxbfsh*maxbfsh*nvir +
     $        maxbfsh*maxbfsh*nbf*nocc
       moints_lmem = ma_sizeof(MT_DBL, memd, MT_BYTE) +
     $               ma_sizeof(MT_INT, memi, MT_BYTE)
       return
       end







       subroutine moints_localfname( prefix, fname, len )
       implicit none
#include "global.fh"       
       character*(*) prefix
       character*(*) fname
       integer len
       integer plen, i, ii, jj
       integer inp_strlen
       external inp_strlen

       plen = inp_strlen( prefix )
       do i=1,plen
          fname(i:i) = prefix(i:i)
       enddo
       fname(plen+1:plen+1) = '.'
       jj = ga_nodeid()
       do i=1,3
          ii = jj/(10**(3-i))
          fname(plen+i+1:plen+i+1) = char(ichar('0')+ii)
          jj = jj - ii*(10**(3-i))
       enddo
       len = inp_strlen( fname )
       return
       end
       
