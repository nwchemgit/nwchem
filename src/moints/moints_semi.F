c==========================================================
c  This stupid section is to try an localize the
c  interface differences between CRAY and the rest
c  on word-addressable I/O.
c==========================================================
c  Put at the beginning for possible inlining
c
c
      subroutine waopen( unit, fname, nblk, istat, ierr )
      implicit none
      integer unit, nblk, istat, ierr
      character*(*) fname
#ifndef CRAY_T3D
      integer len
      integer inp_strlen
      external inp_strlen

      len = inp_strlen(fname)
      call wopen( unit, fname, len, nblk, istat, ierr )
#else
      integer i
      double precision fdummy
      character*8 aunit
      equivalence (fdummy,aunit)
c
      double precision dunits
      common/crayunits/dunits(100)
      
      do i=1,8
        aunit(i:i) = fname
      enddo
      call wopen( fdummy, nblk, istat, ierr )
      dunits(unit) = fdummy
#endif
      return
      end


      subroutine waclose( unit, ierr )
      implicit none
      integer unit, ierr
#ifndef CRAY_T3D
      call wclose(unit,ierr)
#else
      double precision dunits
      common/crayunits/dunits(100)
      double precision xu
      call wclose(xu,ierr)
#endif
      return
      end



      subroutine waput( unit, v, addr, len, ierr )
      implicit none
      integer unit, addr, len, ierr
      double precision v(*)

#ifndef CRAY_T3D
      call putwa( unit, v, addr, len, ierr )
#else
      double precision dunits
      common/crayunits/dunits(100)
      double precision xu

      xu = dunits(unit)
      call putwa( xu, v, addr, len, ierr )
#endif
      return
      end




      subroutine waget( unit, v, addr, len, ierr )
      implicit none
      integer unit, addr, len, ierr
      double precision v(*)

#ifndef CRAY_T3D
      call getwa( unit, v, addr, len, ierr )
#else
      double precision dunits
      common/crayunits/dunits(100)
      double precision xu

      xu = dunits(unit)
      call getwa( xu, v, addr, len, ierr )
#endif
      return
      end


c
c  $Id: moints_semi.F,v 1.13 1996-01-29 23:55:28 mg141 Exp $
c
c  Generates and writes out half-transformed 
c  Coulomb and exchange integrals for some multipassed
c  segment.
c
c  Half-transformed integrals are transposed in GA memory 
c  buffers and flushed when full. Size of buffer is determined
c  by available memory.
c
c  Flat file is written using word-addressable I/O. Prefaced
c  with some indexing info. File layout of integrals is:
c
c    nnbf
c  |<---->
c  !       nnbf*nvirl
c  |<------------------->
c  |                        nnbf*nvirl*nseg
c  |<------------------------------------------------------------->
c
c
c  nnbf:  triangle of basis functions (with sparsity) <= (nbf*(nbf+1))/2
c  nvirl: number of virtual indices this processor owns ~nvir/nproc
c  nseg:  length of occupied segment
c    
c
c
c
       subroutine moints_semi( basis, tol2e, ocoul, oexch,
     $                         oseglo, oseghi,
     $                         olo, ohi, vlo, vhi, g_movecs )
       implicit none
#include "tcgmsg.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "sym.fh"
#include "bas.fh"
c
       integer basis                          ! [input]  Basis handle
       double precision tol2e                 ! [input]  Integral tolerance
       logical ocoul                          ! [input]  Toggle Coulomb integrals
       logical oexch                          ! [input]  Toggle exchange integrals
       integer oseglo, oseghi                 ! [input]  Occupied segment range
       integer olo, ohi                       ! [input]  Correlated occupied index range
       integer vlo, vhi                       ! [input]  Virtual index range
       integer g_movecs                       ! [input]  MO coefficients
c
       integer l_gmap, k_gmap
       integer l_shmap, k_shmap
       integer l_gloc, k_gloc
       integer l_iopass, k_iopass
       integer l_rlen, k_rlen
       integer l_rloc, k_rloc
       integer nbf, nsh, maxbfsh, nseg, nocc, nvir
       integer ngrp, blksiz, nnbf
       integer balo, bahi, bilo, bihi
       integer myid, numnodes, rlo, rhi, clo, chi
       integer mxrlen, itmp, ierr, vtmp(100)
       integer lmemreq, gmem, g_dummy
       integer g_jbuf, g_kbuf
       logical status
       logical osym
       character*256 fnameK, fnameJ
       integer junit, kunit
       parameter (junit=10,kunit=11)
c
       integer moints_numgr
       integer moints_lmem
       external moints_numgr
       external moints_lmem
       data blksiz/16/                              ! Affects performance only
       data osym/.false./
c
c  General info
c  Compare definition for nseg      : batch range
c                         nocc      : number of occupieds to highest i
c
       myid = ga_nodeid()
       numnodes = ga_nnodes()
       nocc = ohi - olo + 1
       nseg = oseghi - oseglo + 1
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
       mxrlen = gmem/(8*nseg*itmp)
       mxrlen = min(mxrlen,(nbf*nbf))
       MXRLEN = 512
c
c$$$       IF (GA_NODEID().EQ.0) WRITE(6,900) LMEMREQ, GMEM, MXRLEN
c$$$ 900   FORMAT(/,'Local Memory required:  ',i10,' bytes '
c$$$     $        /,'Global Memory remaining:',i10,' bytes ',
c$$$     $        /,'IO Buffer length:       ',i10)
c
c  Create GA coulomb buffer
c
       if (ocoul) then
         itmp = (nocc/numnodes + min(mod(nocc,numnodes),1))*nseg
         if (.not.(ga_create( MT_DBL, mxrlen, nseg*nocc, 'transp Jbuf',
     $                        mxrlen, itmp, g_jbuf )))
     $           call errquit('moints_semi: cannot allocate J buffer',0)
         call ga_distribution(g_jbuf, myid, rlo, rhi, clo, chi )
         if ((clo.eq.0).and.(chi.eq.-1)) then
           bilo = 0
           bihi = -1
         else
           if ((nseg.ne.1).and.((mod(chi,nseg).ne.0).or.
     $         (mod(clo,nseg).ne.1))) then
             write(6,221) ga_nodeid(),clo,chi
 221         format('ME:',i3,5x,' Range:',i3,'-',i3)
             call errquit('moints_semi:wrong distrib. for J transp',0)
           endif
           bilo = olo + (clo-1)/nseg
           bihi = olo + (chi/nseg) - 1
         endif
c$$$         WRITE(6,882) 'COULOMB',GA_NODEID(),CLO,CHI,ITMP,BILO,BIHI
       endif
c
c  Create GA exchange buffer
c
       if (oexch) then
         itmp = (nvir/numnodes + min(mod(nvir,numnodes),1))*nseg
         if (.not.(ga_create( MT_DBL, mxrlen, nseg*nvir, 'transp Kbuf',
     $                        mxrlen, itmp, g_kbuf )))
     $           call errquit('moints_semi: cannot allocate K buffer',0)
         call ga_distribution(g_kbuf, myid, rlo, rhi, clo, chi )
         if ((clo.eq.0).and.(chi.eq.-1)) then
           balo = 0
           bahi = -1
         else
           if ((nseg.ne.1).and.((mod(chi,nseg).ne.0).or.
     $         (mod(clo,nseg).ne.1))) then
             write(6,221) ga_nodeid(),clo,chi
             call errquit('moints_semi:wrong distrib. for K transp',0)
           endif
           balo = vlo + (clo-1)/nseg
           bahi = vlo + (chi/nseg) - 1
         endif
c$$$         WRITE(6,882) 'EXCHANGE',GA_NODEID(),CLO,CHI,ITMP,BALO,BAHI
c$$$ 882     FORMAT(A8,5X,'ME:',I5,'  CLO:CHI ',2I5,' YCHUNK:',i5,
c$$$     $          ' VRANGE:',2I5)
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
       g_dummy = g_kbuf
       if (.not.(oexch).and.(ocoul)) g_dummy = g_jbuf
       call moints_locmap( basis, nsh, nbf, tol2e, int_mb(k_shmap),
     $                     g_dummy, int_mb(k_gloc), nnbf,
     $                     int_mb(k_rloc), int_mb(k_iopass),
     $                     int_mb(k_rlen))
c
c  Open local file and write initial info
c  Note complications to do with writing ints in units of double
c
       call ifill(10, 0, vtmp, 1 )
       itmp = nnbf+mod(nnbf,2)
       vtmp(1) = ma_sizeof(MT_INT,12,MT_DBL) + 
     $      ma_sizeof(MT_INT,nnbf+mod(nnbf,2),MT_DBL)
*       vtmp(1) = ma_sizeof(MT_INT,(10 + itmp),MT_DBL) + 1
       vtmp(2) = nnbf
       if (oexch) then
         Call Util_PName('KHalf', fnameK)
         call waopen( kunit, fnameK, 1, 0, ierr )
         vtmp(3) = balo
         vtmp(4) = bahi
         call putwa( kunit, vtmp, 1,
     $               ma_sizeof(MT_INT, 10, MT_DBL), ierr )
         call putwa( kunit, int_mb(k_rloc),
     $               ma_sizeof(MT_INT, 12, MT_DBL),
     $               ma_sizeof(MT_INT, itmp, MT_DBL), ierr )
       endif
       if (ocoul) then
         Call Util_PName('JHalf', fnameJ)
         call waopen(junit, fnameJ, 1, 0, ierr )
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
     $                     oseglo, oseghi, olo, ohi, vlo, vhi,
     $                     ocoul, oexch, osym, kunit, junit,
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

       if (ocoul) then
         if (.not. ga_destroy(g_jbuf))
     $     call errquit('moints: failed to destroy g_jbuf',g_jbuf)
       endif

       if (oexch) then
         if (.not. ga_destroy(g_kbuf))
     $     call errquit('moints: failed to destroy g_kbuf',g_kbuf)
       endif
c
       if (oexch) call waclose( kunit, ierr )
       if (ocoul) call waclose( junit, ierr )
c
       call ga_sync()
       return
       end
       








c
c
c
c
c
       subroutine moints_semi_a( basis, nbf, nsh, maxbfsh, tol2e,
     $                           oseglo, oseghi, olo, ohi, vlo, vhi,
     $                           ocoul, oexch, osym, kunit, junit,
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
       integer basis                          ! [input] Basis handle
       integer nbf                            ! [input] Basis functions
       integer nsh                            ! [input] Shells
       integer maxbfsh                        ! [input] Largest shell
       double precision tol2e                 ! [input] Integral tolerance
       integer oseglo, oseghi                 ! [input] Occupied segment range
       integer olo, ohi                       ! [input] Occupied index range
       integer vlo, vhi                       ! [input] Virtual index range
       logical ocoul                          ! [input] Toggle Coulomb
       logical oexch                          ! [input] Toggle exchange
       logical osym                           ! [input] Toggle symmetry
       integer junit, kunit                   ! [input] Unit numbers for IO
       integer shmap(nsh)                     ! [input] Map for shells
       integer blksiz                         ! [input] Blocksize
       integer ngrp                           ! [input] Number of groups of shells
       integer gmap(*)                        ! [input] Map of shell-group
       integer gloc(nbf*nbf)                  ! [input] bf -> GA memory map
       integer nnbf                           ! [input] Number of screened basf pairs <= (nbf*(nbf+1))/2
       integer iopass(0:*)                    ! [input] IO pass for each task
       integer rlen(*)                        ! [input] Record lengths for IO writes
       integer g_movecs                       ! [input] MO coefficients
       integer g_jbuf                         ! [scratch] Buffer for half-trans Coulomb
       integer g_kbuf                         ! [scratch] Buffer for half-trans exchange
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
       double precision scale, schw_ij
       double precision tpass
       logical status
c
       integer gr_len
       external gr_len
c
       num_nodes = ga_nnodes()
       nmo1 = oseghi - oseglo + 1
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
       status = ma_push_get(MT_DBL,n_hlp2,'hlp2 block',l_hlp2,k_hlp2)
       status = ma_push_get(MT_DBL,n_ssni,'ssni block',l_ssni,k_ssni)
       if (.not.(status)) call errquit('cannot allocate local memory',0)
c
c Initialize
c
       ploop = 0
       jopass = 1
*       iocnt = ma_sizeof(MT_INT,(10+nnbf+mod(nnbf,2)),MT_DBL) + 1    ! start counting from initial segment
c       write(6,*) ' old iocnt', iocnt
       iocnt = ma_sizeof(MT_INT,12,MT_DBL) + 
     $      ma_sizeof(MT_INT,nnbf+mod(nnbf,2),MT_DBL)
c       write(6,*) ' new iocnt', iocnt
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
                 call ga_sync()
                 if (ocoul) then
	           call moints_wrbuf( junit, oseglo, oseghi,
     $                                olo, ohi, nnbf, iocnt,
     $                                rlen(jopass), g_jbuf )
                   call ga_zero(g_jbuf)
                 endif
                 if (oexch) then
                   call moints_wrbuf( kunit, oseglo, oseghi,
     $                                vlo, vhi, nnbf, iocnt,
     $                                rlen(jopass), g_kbuf )
                   call ga_zero(g_kbuf)
                 endif
                 tpass = tcgtime() - tpass
c$$$                 WRITE(6,901) GA_NODEID(), NEXT, IOPASS(NEXT), JOPASS,
c$$$     $                        RLEN(JOPASS), TPASS
c$$$ 901             FORMAT('Me:',i5,2x,'Next:',i5,2x,'IOPass:',i5,2x,
c$$$     $                  'Pass:',i5,2x,'Length:',i6,2x,'Time:',f10.3)
c$$$                 CALL UTIL_FLUSH(6)
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
     $                               dbl_mb(k_ssbb), osym )
                   if (lb.ne.kb) then
                     call moints_blktr( ilen, jlen, kblen, lblen,
     $                                  dbl_mb(k_ssbb),
     $                                  dbl_mb(k_ssbbt))
                     call moints_trf1( nbf, oseglo, oseghi,
     $                                 ilen, jlen, kbflo, kbfhi,
     $                                 lbflo, lbfhi, 1.d0,
     $                                 dbl_mb(k_ssbb),
     $                                 dbl_mb(k_ssbbt), dbl_mb(k_mo),
     $                                 dbl_mb(k_ssni), dbl_mb(k_hlp) )

                   else
                     call moints_trf1( nbf, oseglo, oseghi,
     $                                 ilen, jlen, kbflo, kbfhi,
     $                                 lbflo, lbfhi, 0.5d0,
     $                                 dbl_mb(k_ssbb),
     $                                 dbl_mb(k_ssbb), dbl_mb(k_mo),
     $                                 dbl_mb(k_ssni), dbl_mb(k_hlp) )
                   endif
                 enddo
               enddo
               if (ocoul)
     $           call moints_trf2Ky( nbf, oseglo, oseghi, olo, ohi,
     $                               ibflo, ibfhi, jbflo, jbfhi,
     $                               gloc, dbl_mb(k_ssni),
     $                               dbl_mb(k_hlp), dbl_mb(k_hlp2),
     $                               dbl_mb(k_mo), g_jbuf )
               if (oexch)
     $           call moints_trf2Ky( nbf, oseglo, oseghi, vlo, vhi,
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
         call moints_wrbuf( junit, oseglo, oseghi, olo, ohi, nnbf,
     $                      iocnt, rlen(jopass), g_jbuf )
       endif
       if (oexch) then
         call moints_wrbuf( kunit, oseglo, oseghi, vlo, vhi, nnbf,
     $                      iocnt, rlen(jopass), g_kbuf )
       endif
       tpass = tcgtime() - tpass
c$$$       WRITE(6,901) GA_NODEID(), NEXT, IOPASS(NEXT), JOPASS,
c$$$     $              RLEN(JOPASS), TPASS
c$$$       CALL GA_SYNC()
c$$$       CALL UTIL_FLUSH(6)
c$$$       PRINT*,'ME:',GA_NODEID(),' End AO integral loop and final write'
c$$$       CALL UTIL_FLUSH(6)
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
      subroutine moints_trf2Ky( nbf, oseglo, oseghi, vlo, vhi,
     $                          ilo, ihi, jlo, jhi, gloc,
     $                          ssni, h, h2, c, g_buf )
      implicit none
#include "global.fh"
      integer nbf, oseglo, oseghi, vlo, vhi
      integer ilo, ihi, jlo, jhi
      integer gloc(nbf,nbf)
      double precision ssni(nbf,jlo:jhi,ilo:ihi,oseglo:oseghi)
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
      do o=oseglo,oseghi
        do v=vlo,vhi
          vo = (v-vlo)*(oseghi-oseglo+1) + o - oseglo + 1
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





      subroutine moints_wrbuf( munit, oseglo, oseghi, plo, phi, nnbf,
     $                         iocnt, rlen, g_buf )
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer munit
      integer oseglo, oseghi, plo, phi
      integer nnbf, iocnt, rlen
      integer g_buf
c
      integer faddr, ierr
      integer o, p, po, nocc, np
      integer myid, rlo, rhi, clo, chi, k_local, ld
c
      nocc = oseghi - oseglo + 1
      np = phi - plo + 1
      myid = ga_nodeid()
      call ga_distribution(g_buf, myid, rlo, rhi, clo, chi )
      do o=oseglo,oseghi
        do p=plo,phi
          po = (p-plo)*nocc + o - oseglo + 1
          faddr = ((o-oseglo)*np + p - plo)*nnbf +  iocnt
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
      subroutine moints_readintJ( nbf, oseglo, oseghi, olo, ohi,
     $                            vlo, vhi, c, orbe )
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer nbf, oseglo, oseghi
      integer olo, ohi, vlo, vhi
      double precision c(nbf,nbf)
      double precision orbe(nbf)
c
      integer nocc, nnbf, ioff, p, o, oo
      integer bihi, bilo, itmp, len
      integer l_v, k_v, l_t, k_t, l_x, k_x
      logical status
      character*256 fname
      integer junit
      parameter(junit=10)

      Call Util_PName('JHalf', fname)
      call waopen(junit, fname, 1, 0, itmp )
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
      do o=oseglo,oseghi
        oo = o - oseglo + 1
        do p=bilo,bihi
          call dfill((nbf*nbf), 0.d0, dbl_mb(k_x), 1 )
          call moints_rdhfint( junit, olo, ohi, oo, p, nnbf,
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
      call waclose(junit, itmp )
      return
      end











      subroutine moints_readintK( nbf, oseglo, oseghi, olo, ohi,
     $                            vlo, vhi, c, orbe, g_epair )
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
      integer MSG_SEMIMP2_SUM
      parameter(MSG_SEMIMP2_SUM=10241)
      integer nbf, oseglo, oseghi, olo, ohi, vlo, vhi
      double precision c(nbf,nbf)
      double precision orbe(nbf)
      integer g_epair
c
      integer nseg, nocc, nvir, nnbf, ioff, v, o, oo
      integer balo, bahi, itmp
      integer l_v, k_v, l_t, k_t, l_x, k_x
      integer j, oj, jj, vvlo, vvhi, k_local, ld, myid
      integer rlo, rhi, clo, chi, len
      integer g_exch
      double precision tmp2
      character*256 fname
      integer kunit
      parameter (kunit=11)
      double precision xx, denom, e2
      logical status
      double precision moints_epair_eval
      external moints_epair_eval
       integer i

      tmp2 = tcgtime()
      Call Util_PName('KHalf', fname)
      call waopen(kunit, fname, 1, 0, itmp )
      myid = ga_nodeid()
      e2 = 0.d0
      nocc = ohi - olo + 1
      nvir = vhi - vlo + 1
      nseg = oseghi - oseglo + 1
      call moints_vrange( kunit, balo, bahi, nnbf, ioff )
      status = ma_push_get(MT_INT,(nnbf+mod(nnbf,2)),'scatter',l_v,k_v)
      call moints_getscattv( kunit, nnbf, int_mb(k_v) )
c$$$       do i=1,nnbf
c$$$         write(6,*) 'k_v',int_mb(k_v+i-1)
c$$$       enddo
      if (.not.ga_create( MT_DBL, (nvir*nvir), nocc, 'exch',
     $                    (nvir*nvir), 1, g_exch))
     $    call errquit('moints_semimp2: cannot allocate exch',0)
      call ga_distribution(g_exch, myid, rlo, rhi, clo, chi )
c
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
      do o=oseglo,oseghi
        oo = o - oseglo + 1
        call ga_zero(g_exch)
        do v=vlo,vhi
          if ((v.ge.balo).and.(v.le.bahi)) then
            call dfill((nbf*nbf), 0.d0, dbl_mb(k_x), 1 )
            call moints_rdhfint( kunit, vlo, vhi, oo, v, nnbf,
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
            do j=olo,o
              oj = j-olo+1
              jj = (j-olo)*nvir
              call ga_put( g_exch, vvlo, vvhi, oj, oj,
     $                     dbl_mb(k_x+jj), nvir)
            enddo
          endif
        enddo
        call ga_sync()
        do j=olo,o
          jj = j-olo+1
          if ((jj.ge.clo).and.(jj.le.chi)) then
            denom = orbe(o) + orbe(j)
            call ga_access(g_exch, rlo, rhi, jj, jj, k_local, ld )

c$$$            WRITE(6,920) O,J
c$$$ 920        FORMAT(/,' Operator: [',i4,',',i4,']')
c$$$            CALL MOINTS_MATPRINT(NVIR,NVIR,DBL_MB(K_LOCAL))
c$$$            CALL UTIL_FLUSH(6)
c$$$
            xx = moints_epair_eval( nvir, 0, nocc, nvir,
     $                              dbl_mb(k_local), orbe, denom )
            call ga_release(g_exch, rlo, rhi, jj, jj )
            if (o.eq.j) xx = xx*0.5d0
            e2 = e2 + xx
            oj = ((o-olo+1)*(o-olo))/2 + j
            call ga_put(g_epair,1,1,oj,oj,xx,1)
	    
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
c$$$      call ga_sync()
c$$$      call ga_dgop(MSG_SEMIMP2_SUM,e2,1,'+')
c$$$      if (ga_nodeid().eq.0) then
c$$$        write(6,681) e2
c$$$ 681    format(/,10x,'MP2 correction:',5x,f16.10)
c$$$        call util_flush(6)
c$$$      endif
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
      call waclose(kunit, itmp )
c
      call ga_sync()
      
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
      subroutine moints_readintK2( nbf, oseglo, oseghi, vlo, vhi, c, orbe )
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer nbf, oseglo, oseghi, vlo, vhi
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
      
      Call Util_PName('KHalf', fname)
      call waopen(kunit, fname, 1, 0, itmp )
      e2 = 0.d0
      nocc = oseghi - oseglo + 1
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
      do o=oseglo,oseghi
        call ga_zero(g_exch)
        do v=balo,bahi
          call dfill((nbf*nbf), 0.d0, dbl_mb(k_x), 1 )
          call moints_rdhfint( kunit, vlo, vhi, o, v, nnbf,
     $                         ioff, dbl_mb(k_t))
          call scatter( nnbf, dbl_mb(k_x), int_mb(k_v), dbl_mb(k_t) )
          call upper2square( nbf, dbl_mb(k_x), dbl_mb(k_x) )

          call dgemm( 'n', 'n', nbf, nocc, nbf, 1.d0, dbl_mb(k_x), nbf,
     $                c(1,oseglo), nbf, 0.d0, dbl_mb(k_t), nbf )
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
c$$$            denom = orbe(o) + orbe(j+oseglo-1)
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
      call waclose(kunit, itmp )
      return
      end
#endif      






      subroutine moints_rdhfint( munit, plo, phi, o, p, nnbf, ioff, t )
      implicit none
      integer munit
      integer plo, phi, o, p, nnbf, ioff
      double precision t(nnbf)
      integer faddr, np, ierr

      np = phi - plo + 1
      faddr = ((o-1)*np + p - plo)*nnbf + ioff
      call getwa( munit, t, faddr, nnbf, ierr )
      
      return
      end






      subroutine moints_vrange( munit, balo, bahi, nnbf, ioff )
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer munit
      integer balo, bahi, nnbf, ioff
      integer vtmp(10), ierr
C      INTEGER I

      call getwa( munit, vtmp, 1,
     $            ma_sizeof(MT_INT, 10, MT_DBL), ierr )
      ioff = vtmp(1)
      nnbf = vtmp(2)
      balo = vtmp(3)
      bahi = vtmp(4)

c$$$      PRINT*,'ME: ',GA_NODEID(),' ERROR:',IERR
c$$$      CALL UTIL_FLUSH(6)
c$$$      IF (GA_NODEID().EQ.0) THEN
c$$$        PRINT*,'READ BACK INITIAL PARAMETERS'
c$$$        WRITE(6,771) (VTMP(I),I=1,10)
c$$$ 771    FORMAT(16I4)
c$$$        CALL UTIL_FLUSH(6)
c$$$      ENDIF
c$$$      CALL GA_SYNC()

      return
      end


      subroutine moints_getscattv( munit, nnbf, v )
      implicit none
#include "mafdecls.fh"
      integer munit
      integer nnbf, v(*)
      integer itmp, ierr
C      INTEGER I
      
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
