
c
c  $Id: moints_semi.F,v 1.47 1997-08-14 18:06:45 rg240 Exp $
c
c  Generates and writes out half-transformed 
c  exchange integrals for some multipassed
c  segment, (oseglo:oseghi)
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
       subroutine moints_semi( basis, tol2e, oseglo, oseghi, 
     $                         olo, ohi, vlo, vhi, g_movecs, oblk )
       implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "eaf.fh"
#include "util.fh"
#include "sym.fh"
#include "bas.fh"
c
       integer basis                          ! [input]  Basis handle
       double precision tol2e                 ! [input]  Integral tolerance
       integer oseglo, oseghi                 ! [input]  Occupied segment range
       integer olo, ohi                       ! [input]  Correlated occupied index range
       integer vlo, vhi                       ! [input]  Virtual index range
       integer g_movecs                       ! [input]  MO coefficients
       logical oblk                           ! [input]  Toggle AO integral blocking
c
       integer l_gmap, k_gmap
       integer l_shmap, k_shmap, l_bfmap, k_bfmap, l_rbfmap, k_rbfmap
       integer l_glo, k_glo, l_ghi, k_ghi
       integer l_gloc, k_gloc
       integer l_iopass, k_iopass
       integer l_rlen, k_rlen
       integer l_rloc, k_rloc
       integer nbf, nsh, maxbfsh, nseg, nocc, nvir
       integer ngrp, blen, nnbf
       integer pvlo, pvhi
       integer polo, pohi
       integer myid, numnodes, rlo, rhi, clo, chi
       integer itmp, vtmp(10), mxrlen
       integer lmemreq
       double precision gmem
       double precision iotmpptr
       integer g_kbuf
       logical status
       logical osym
       character*256 fnameK
       integer kunit
       logical oprint_mem
c
       integer moints_numgr
       integer moints_lmem
       external moints_numgr
       external moints_lmem
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
       oprint_mem = util_print('memory',print_high) .and.
     $      ga_nodeid().eq.0
c
c  Clear 4-index statistics
c
       call moints_stats_clear()
c
c     Initial block length is very generous for efficiency.  Reduce it
c     in while loop (1212) if things do not fit.
c
       blen = min(nbf,16*maxbfsh) ! Initial generous allocation
       blen = 32
 1212  continue
c
c  Memory arithmetic
c
       lmemreq = moints_lmem(basis, nocc, nvir, blen)
       if ((.not.ga_uses_ma()).and.(.not.(ga_memory_limited())))
     $   call errquit('moints_semi: cannot determine memory limit',0)
c
       gmem = ga_memory_avail() / 8 ! Want amount in doubles
       call ga_dgop(313, gmem, 1, 'min')
       if (oprint_mem) then
          write(6,717) blen, gmem, lmemreq
 717      format(/'   Block length           ', i8/
     $            '   Local GA available     ', f10.1/
     $            '   Local memory required  ', i8)
          call util_flush(6)
       endif
c
       if (ga_uses_ma()) gmem = gmem - lmemreq
c
       if (oprint_mem) then
          write(6,718) gmem
 718      format('   Adjusted GA available  ', f10.1)
          call util_flush(6)
       endif
c
       gmem = gmem*numnodes
c
       gmem = gmem * 0.9d0      ! Protect against distribution inequity
c
       if (oprint_mem) then
          write(6,719) gmem
 719      format('   Aggregate GA available ', 1p, d10.2)
          call util_flush(6)
       endif
c
       mxrlen = int(gmem/(nseg*nvir))
       mxrlen = min(
     $      mxrlen, 
     $      min(nbf*nbf,10*maxbfsh*maxbfsh*numnodes))
c
       if (oprint_mem) then
          write(6,720) mxrlen, nseg, nvir, nbf, maxbfsh, numnodes
 720      format( '   Maximum record length  ', i8/
     $            '   No. of segments        ', i8/
     $            '   Virtuals dimension     ', i8/
     $            '   No. of functions       ', i8/
     $            '   Max. functions/shell   ', i8/
     $            '   Number of nodes        ', i8/)
          call util_flush(6)
       endif
c
       call ga_sync()           ! Just to ensure printing completes
c
       if (mxrlen .lt. min(nbf*nbf,maxbfsh*maxbfsh*numnodes)) then
          if (blen .le. maxbfsh) call errquit
     $         ('moints_semi: insufficient GA available',
     $         mxrlen*nseg*nvir)
          blen = max(maxbfsh,blen/2) ! Reduce blocksize and try again
          goto 1212
       endif
c
       if (oprint_mem) 
     $      write(6,900) lmemreq, nint(gmem/numnodes), mxrlen
900   format(/,'Semi-direct integral transformation',
     $     /,' Local Memory required:   ',i10,' words '
     $     /,' Global Memory remaining: ',i10,' words per node',
     $     /,' IO Buffer length:        ',i10,' words')
c
c  Create GA exchange buffer
c
c
      itmp = (nvir/numnodes + min(mod(nvir,numnodes),1))*nseg
**          itmp = ((nvir-1)/numnodes+1)*nseg ! rjh
      if (.not.(ga_create( MT_DBL, mxrlen, nseg*nvir, 'transp Kbuf',
     $                     mxrlen, itmp, g_kbuf )))
     $   call errquit('moints_semi: cannot allocate K buffer',0)
      call ga_distribution(g_kbuf, myid, rlo, rhi, clo, chi )
      if ((clo.eq.0).and.(chi.eq.-1)) then
        pvlo = 0
        pvhi = -1
      else
        if ((nseg.ne.1).and.((mod(chi,nseg).ne.0).or.
     $       (mod(clo,nseg).ne.1))) then
          write(6,221) ga_nodeid(),clo,chi
 221      format('ME:',i3,5x,' Range:',i3,'-',i3)
          call errquit('moints_semi:wrong distrib. for K transp',0)
        endif
        pvlo = vlo + (clo-1)/nseg
        pvhi = vlo + (chi/nseg) - 1
      endif
c
c  Reorder shells by descending shell-length
c  and group shells by blocksize
c
       status = ma_push_get(MT_INT,nsh,'shell order map',
     $                      l_shmap, k_shmap)
       status = ma_push_get(MT_INT,nsh,'group lo', l_glo, k_glo )
       status = ma_push_get(MT_INT,nsh,'group hi', l_ghi, k_ghi)
       status = ma_push_get(MT_INT,nbf,'basis map',
     $                      l_bfmap, k_bfmap)
       status = ma_push_get(MT_INT,nbf,'rev basis map',
     $                      l_rbfmap, k_rbfmap)
       call moints_shorder( basis, nsh, nbf, blen, ngrp,
     $                      int_mb(k_glo), int_mb(k_ghi),
     $                      int_mb(k_shmap),
     $                      int_mb(k_bfmap), int_mb(k_rbfmap) )
c
c  ------- Old stuff
c
#ifdef OLD_AOBLOCK
       ngrp = moints_numgr( basis, blen )
       status = ma_push_get(MT_INT,(ngrp*4),'grp map', l_gmap, k_gmap)
       call moints_grmap( basis, blen, ngrp, int_mb(k_gmap))
#else
       status = ma_push_get(MT_INT,(ngrp*4),'grp map', l_gmap, k_gmap)
#endif
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
     $                     g_kbuf, int_mb(k_gloc), nnbf,
     $                     int_mb(k_rloc), int_mb(k_iopass),
     $                     int_mb(k_rlen))
c
c  Open local file and write initial info
c
       call ifill(10, 0, vtmp, 1 )
       vtmp(1) = ma_sizeof(MT_INT,10+nnbf,MT_BYTE)
       vtmp(2) = nnbf
       vtmp(3) = pvlo
       vtmp(4) = pvhi
       call util_file_name('kh', .true.,.true., fnamek)
       if (eaf_open( fnamek, EAF_RW, kunit).ne.0)
     $   call errquit('moints_semi: failed to open half int file',0)
       if (eaf_write(kunit, 0.d0, vtmp,
     $               ma_sizeof(MT_INT,10,MT_BYTE)).ne.0)
     $   call errquit('moints_semi: failed to write header info',0)
       iotmpptr = ma_sizeof(MT_INT,10,MT_BYTE)
       if (eaf_write(kunit, iotmpptr, int_mb(k_rloc),
     $               ma_sizeof(MT_INT,nnbf,MT_BYTE)).ne.0)
     $   call errquit('moints_semi: failed to write header info',0)
       if (.not. ma_pop_stack(l_rloc))
     $     call errquit('moints: failed to pop', l_rloc)
c
c  Call the real stuff
c
       call moints_semi_a( basis, nbf, nsh, maxbfsh, tol2e,
     $                   oseglo, oseghi, olo, ohi, vlo, vhi,
     $                   polo, pohi, pvlo, pvhi, 
     $                   osym, oblk, kunit, int_mb(k_shmap),
     $                   ngrp, int_mb(k_glo), int_mb(k_ghi),
     $                   int_mb(k_bfmap), int_mb(k_rbfmap),
     $                   blen, int_mb(k_gmap), int_mb(k_gloc), 
     $                   nnbf, int_mb(k_iopass),
     $                   int_mb(k_rlen), g_movecs, g_kbuf )
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
       if (.not. ma_pop_stack(l_rbfmap))
     $     call errquit('moints: failed to pop', l_rbfmap)
       if (.not. ma_pop_stack(l_bfmap))
     $     call errquit('moints: failed to pop', l_bfmap)
       if (.not. ma_pop_stack(l_ghi))
     $     call errquit('moints: failed to pop', l_ghi)
       if (.not. ma_pop_stack(l_glo))
     $     call errquit('moints: failed to pop', l_glo)
       if (.not. ma_pop_stack(l_shmap))
     $     call errquit('moints: failed to pop', l_shmap)
       if (.not. ga_destroy(g_kbuf))
     $      call errquit('moints: failed to destroy g_kbuf',g_kbuf)
c
       if (eaf_close(kunit).ne.0)
     $   call errquit('moints_semi: failed to close file',kunit)
c
       call ga_sync()
       if (util_print('statistics', print_high)) 
     $   call moints_stats_print( 'semi' )

       return
       end
       








c
c
c
c
c
       subroutine moints_semi_a( basis, nbf, nsh, maxbfsh, tol2e,
     $                           oseglo, oseghi, olo, ohi, vlo, vhi,
     $                           polo, pohi, pvlo, pvhi,
     $                           osym, oblk,
     $                           kunit, shmap, 
     $                           ngrp, glo, ghi, bfmap, rbfmap,
     $                           blen, gmap, gloc, nnbf, iopass,
     $                           rlen, g_movecs, g_kbuf )
       implicit none
#include "tcgmsg.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "util.fh"
#include "schwarz.fh"
#include "msgids.fh"
c
       integer MSG_SEMI_IO_RATE
       parameter(MSG_SEMI_IO_RATE=18421)
c
c   Arguments
c
       integer basis                          ! [input] Basis handle
       integer nbf                            ! [input] Basis functions
       integer nsh                            ! [input] Shells
       integer maxbfsh                        ! [input] Largest shell
       double precision tol2e                 ! [input] Integral tolerance
       integer oseglo, oseghi                 ! [input] Occupied segment range
       integer olo, ohi                       ! [input] Occupied index range
       integer vlo, vhi                       ! [input] Virtual index range
       integer polo, pohi                     ! [input] Occupied index range for this processor
       integer pvlo, pvhi                     ! [input] Virtual index range for this processor
       logical osym                           ! [input] Toggle symmetry
       logical oblk                           ! [input] Toggle AO integral blocking
       integer kunit                          ! [input] Unit numbers for IO
       integer shmap(nsh)                     ! [input] Map for shells
       integer ngrp                           ! [input] Number of groups of shells
       integer glo(*)                         ! [input] Group lower shell bound
       integer ghi(*)                         ! [input] Group upper shell bound
       integer bfmap(nbf)                     ! [input] BF map: orig --> new
       integer rbfmap(nbf)                    ! [input] Reverse bf map: new --> orig
       integer blen                           ! [input] Blocksize
       integer gmap(4,ngrp)                   ! [input] Old AOBLOCK shell map
       integer gloc(nbf*nbf)                  ! [input] bf -> GA memory map
       integer nnbf                           ! [input] Number of screened basf pairs <= (nbf*(nbf+1))/2
       integer iopass(0:*)                    ! [input] IO pass for each task
       integer rlen(*)                        ! [input] Record lengths for IO writes
       integer g_movecs                       ! [input] MO coefficients
       integer g_kbuf                         ! [scratch] Buffer for half-trans exchange
c
c   Local variables
c
       integer nmo, nmo1, nmo2
       integer ish0, jsh0, ish, jsh, ilen, jlen
       integer ibflo, ibfhi, jbflo, jbfhi
       integer kbflo, kbfhi, lbflo, lbfhi
       integer kshlo, kshhi, lshlo, lshhi
       integer kblen, lblen
       integer kgr, lgr
       integer qlo, qhi
       integer l_ssbb, k_ssbb
       integer l_ssbbt, k_ssbbt
       integer l_hlp, k_hlp, l_ssni,k_ssni
       integer l_hlp2, k_hlp2
       integer l_eri, k_eri, l_iscr,k_iscr
       integer l_mo, k_mo, l_xmo, k_xmo, l_xmo_t, k_xmo_t
       integer n_ssbb, n_ssbb1, n_ssni, n_hlp, n_hlp2
       integer iz, jz, kz, bsize
       integer mem2, max2e, jopass
       double precision iocnt
       integer num_nodes, ploop, next
       double precision schw_ij
       double precision tpass, tio, ttmp, nwrbytes, iorate
       double precision t0, t1, ttask
       logical st
c
#include "moints_stats.fh"
c       
       integer gr_len
       external gr_len
c
       call ga_zero(g_kbuf)
       num_nodes = ga_nnodes()
       nmo1 = oseghi - oseglo + 1
       nmo2 = vhi - vlo + 1
       qlo = olo
       qhi = vhi
       nmo = qhi - qlo + 1
c
c$$$       WRITE(6,221) GA_NODEID(),PVLO,PVHI,POLO,POHI
c$$$ 221   FORMAT('ME:',I5,5X,' V-RANGE:',2I3,5X,' O-RANGE:',2I3)
c
c  Local MO coefficients
c  Integrals and temporary arrays sizes
c
       call int_mem_2e4c(max2e, mem2)
       bsize = max(blen,maxbfsh)
       n_ssbb = maxbfsh*maxbfsh*bsize*bsize
       n_ssbb1 = max((nmo1*nmo1),n_ssbb)
#ifdef NEW_SPARSE
       n_hlp = max(nmo1,maxbfsh*maxbfsh)*nbf
       n_hlp2 = maxbfsh*maxbfsh
#else
       n_hlp = max((bsize*maxbfsh*maxbfsh*nmo1),(maxbfsh*nbf))
       n_hlp2 = maxbfsh*maxbfsh*nmo2
#endif
       n_ssni = maxbfsh*maxbfsh*nbf*nmo1
       st=ma_push_get(MT_DBL,(nbf*vhi),'r movecs',l_xmo,k_xmo)
       st=st.and.ma_push_get(MT_DBL,(nbf*vhi),'movecs cols',l_mo,k_mo)
       if (.not.(st)) call errquit
     $      ('moints_semi: cannot allocate local memory',0)

       call ga_get(g_movecs, 1, nbf, qlo, qhi, dbl_mb(k_mo), nbf )
#ifdef OLD_AOBLOCK
       call dcopy((nbf*nmo), dbl_mb(k_mo), 1, dbl_mb(k_xmo), 1)
#else
       call row_exch( nbf, nmo, rbfmap, dbl_mb(k_mo), dbl_mb(k_xmo) )
#endif
       if (.not. ma_pop_stack(l_mo))
     $     call errquit('moints: failed to pop', l_mo)

       st=st.and.ma_push_get(MT_DBL, max2e,'ibuf', l_eri, k_eri)
       st=st.and.ma_push_get(MT_DBL, mem2, 'int scr', l_iscr, k_iscr)
       st=st.and.ma_push_get(MT_DBL,n_ssbb1,'ssbb blk',l_ssbb,k_ssbb)
#ifndef NEW_SPARSE
       st=st.and.ma_push_get(MT_DBL,n_ssbb,'ssbbt',l_ssbbt,k_ssbbt)
#endif
       st=st.and.ma_push_get(MT_DBL,n_hlp,'hlp block',l_hlp,k_hlp)
       st=st.and.ma_push_get(MT_DBL,n_hlp2,'hlp2 blk',l_hlp2,k_hlp2)
       st=st.and.ma_push_get(MT_DBL,n_ssni,'ssni blk',l_ssni,k_ssni)
       st=st.and.ma_push_get(MT_DBL,nbf*vhi,'xmot',l_xmo_t,k_xmo_t)
       if (.not.(st)) call errquit
     $      ('moints_semi: cannot allocate local memory',0)

       call moints_xmot(nbf,qhi-qlo+1,dbl_mb(k_xmo),dbl_mb(k_xmo_t))

c     
c
c Initialize
c
       ploop = 0
       jopass = 1
       iocnt = ma_sizeof(MT_INT,10+nnbf,MT_BYTE)
       next = nxtval(num_nodes)
       tio = 0.d0
       nwrbytes = 0.d0
       tpass = util_cpusec()
c
c  4-fold shell loop
c
       do ish0=1,nsh
         do jsh0=1,ish0
           ish = max(shmap(ish0),shmap(jsh0))
           jsh = min(shmap(ish0),shmap(jsh0))
           st = bas_cn2bfr(basis,ish,ibflo,ibfhi)
           st = bas_cn2bfr(basis,jsh,jbflo,jbfhi)
           ilen = ibfhi - ibflo + 1
           jlen = jbfhi - jbflo + 1
           schw_ij = schwarz_shell(ish,jsh)
           if (schw_ij*schwarz_max().ge.tol2e) then
             if (next.eq.ploop) then
c
c  -------------
c  Parallel task
c
c
c  Checkpoint for I/O
c
               mi_ntasks = mi_ntasks + 1
               ttask = util_cpusec()
               if (iopass(next).gt.jopass) then
                 t0 = util_cpusec()
                 call ga_sync()
                 t1 = util_cpusec() - t0
                 mi_minsynch = min(mi_minsynch,t1)
                 mi_maxsynch = max(mi_maxsynch,t1)
                 mi_aggsynch = mi_aggsynch + t1
                 mi_nsynchs = mi_nsynchs + 1
                 ttmp = util_cpusec()
                 nwrbytes = nwrbytes + nmo1*(pvhi-pvlo+1)*rlen(jopass)
                 call moints_wrbuf( kunit, oseglo, oseghi,
     $                              vlo, pvlo,  pvhi, nnbf, iocnt,
     $                              rlen(jopass), g_kbuf )
                 call ga_zero(g_kbuf)
                 tio = tio + util_cpusec() - ttmp
                 tpass = util_cpusec() - tpass
                 iocnt = iocnt + ma_sizeof(MT_DBL,rlen(jopass),MT_BYTE)
                 jopass = jopass + 1
                 tpass = util_cpusec()
               endif
c
c  Half-tranformed Integral generation
c
               call dfill((ilen*jlen*nbf*nmo1),0.d0,dbl_mb(k_ssni),1)
#ifdef OLD_AOBLOCK
               do kgr=1,ngrp
                 call gr_sh_idx(ngrp,gmap,kgr,kshlo,kshhi)
                 kblen = gr_len(ngrp,gmap,kgr,kbflo,kbfhi)
                 do lgr=1,kgr
                   call gr_sh_idx(ngrp,gmap,lgr,lshlo,lshhi)
                   lblen = gr_len(ngrp,gmap,lgr,lbflo,lbfhi)
                  call moints_gblk_old( basis, ish, jsh, kshlo, kshhi,
     $                               lshlo, lshhi, schw_ij, tol2e,
     $                               max2e, dbl_mb(k_eri),
     $                               mem2, dbl_mb(k_iscr),
     $                               ibflo, ibfhi, jbflo, jbfhi,
     $                               kbflo, kbfhi, lbflo, lbfhi,
     $                               dbl_mb(k_ssbb), osym, oblk )
#else
               do kgr=1,ngrp
                 kshlo = glo(kgr)
                 kshhi = ghi(kgr)
                 st = bas_cn2bfr(basis,shmap(kshlo),iz,kz)
                 st = bas_cn2bfr(basis,shmap(kshhi),kz,jz)
                 kbflo = rbfmap(iz)
                 kbfhi = rbfmap(jz)
                 kblen = kbfhi - kbflo + 1
                 do lgr=1,kgr
                   lshlo = glo(lgr)
                   lshhi = ghi(lgr)
                   st = bas_cn2bfr(basis,shmap(lshlo),iz,kz)
                   st = bas_cn2bfr(basis,shmap(lshhi),kz,jz)
                   lbflo = rbfmap(iz)
                   lbfhi = rbfmap(jz)
                   lblen = lbfhi - lbflo + 1
                   t0 = util_cpusec()
                   call moints_gblk( basis, ish, jsh,
     $                               kshlo, kshhi, lshlo, lshhi,
     $                               shmap, rbfmap,
     $                               schw_ij, tol2e, osym, oblk,
     $                               max2e, dbl_mb(k_eri),
     $                               mem2, dbl_mb(k_iscr),
     $                               ibflo, ibfhi, jbflo, jbfhi,
     $                               kbflo, kbfhi, lbflo, lbfhi,
     $                               dbl_mb(k_ssbb) )
                   mi_tint = mi_tint + util_cpusec() - t0
#endif
                   t0 = util_cpusec()
                   mi_flop1 = mi_flop1 + 
     $                   4.d-6*ilen*jlen*kblen*lblen*(oseghi-oseglo+1)
#ifdef NEW_SPARSE
                   call moints_trf1_new(nbf, qlo, qhi, oseglo, oseghi,
     $                  ilen, jlen, kbflo, kbfhi,
     $                  lbflo, lbfhi,
     $                  dbl_mb(k_ssbb),
     $                  dbl_mb(k_xmo_t),
     $                  dbl_mb(k_ssni), dbl_mb(k_hlp))
#else
                   if (lgr.ne.kgr) then
                     call moints_blktr( ilen, jlen, kblen, lblen,
     $                                  dbl_mb(k_ssbb),
     $                                  dbl_mb(k_ssbbt))
                     call moints_trf1( nbf, qlo, qhi, oseglo, oseghi,
     $                                 ilen, jlen, kbflo, kbfhi,
     $                                 lbflo, lbfhi, 1.d0,
     $                                 dbl_mb(k_ssbb),
     $                                 dbl_mb(k_ssbbt), dbl_mb(k_xmo),
     $                                 dbl_mb(k_ssni), dbl_mb(k_hlp) )
                   else
                     call moints_trf1( nbf, qlo, qhi, oseglo, oseghi,
     $                                 ilen, jlen, kbflo, kbfhi,
     $                                 lbflo, lbfhi, 0.5d0,
     $                                 dbl_mb(k_ssbb),
     $                                 dbl_mb(k_ssbb), dbl_mb(k_xmo),
     $                                 dbl_mb(k_ssni), dbl_mb(k_hlp) )
                   endif
#endif
                   mi_t1 = mi_t1 + util_cpusec() - t0
                  enddo
               enddo
               t0 = util_cpusec()
#ifdef NEW_SPARSE
               call moints_trf2Kynew(nbf, qlo, qhi, oseglo, oseghi, 
     $                               vlo, vhi, ibflo, ibfhi, 
     $                               jbflo, jbfhi, gloc, dbl_mb(k_ssni),
     $                               dbl_mb(k_hlp), dbl_mb(k_hlp2),
     $                               dbl_mb(k_xmo_t), g_kbuf )
#else
               call moints_trf2Ky(nbf, qlo, qhi, oseglo, oseghi, 
     $                            vlo, vhi, ibflo, ibfhi, 
     $                            jbflo, jbfhi, gloc, dbl_mb(k_ssni),
     $                            dbl_mb(k_hlp), dbl_mb(k_hlp2),
     $                            dbl_mb(k_xmo), g_kbuf )
#endif
               mi_t2k = mi_t2k + util_cpusec() - t0

               ttask = util_cpusec() - ttask
               mi_mintask = min(mi_mintask,ttask)
               mi_maxtask = max(mi_maxtask,ttask)
               mi_aggtask = mi_aggtask + ttask
c
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
       t0 = util_cpusec()
       call ga_sync()
       t1 = util_cpusec() - t0
       mi_nsynchs = mi_nsynchs + 1
       mi_aggsynch = mi_aggsynch + t1
       mi_minsynch = min(mi_minsynch,t1)
       mi_maxsynch = max(mi_maxsynch,t1)
       ttmp = util_cpusec()
       nwrbytes = nwrbytes + nmo1*(pvhi-pvlo+1)*rlen(jopass)
       call moints_wrbuf( kunit, oseglo, oseghi, vlo, pvlo, pvhi,
     $                    nnbf, iocnt, rlen(jopass), g_kbuf )
       tio = tio + util_cpusec() - ttmp
       tpass = util_cpusec() - tpass
       next = nxtval(-num_nodes)
c
c I/O stats
c
       nwrbytes = nwrbytes*8.e-6
       iorate = nwrbytes/tio
       call ga_dgop(MSG_SEMI_IO_RATE,iorate,1,'+')
       if ((ga_nodeid().eq.0).and.
     $     (util_print('io_stats',print_default))) then
         write(6,886) nwrbytes, tio, iorate
 886     format('Node 0 wrote ',f8.1,' Mb in ',f8.1,' s',5x,
     $          'Agg I/O rate:', f8.1, ' Mb/s')
         call util_flush(6)
       endif
       mi_tio = mi_tio + tio
c
c Clean-up
c
       if (.not. ma_pop_stack(l_xmo_t))
     $     call errquit('moints: failed to pop', l_xmo_t)
       if (.not. ma_pop_stack(l_ssni))
     $     call errquit('moints: failed to pop', l_ssni)
       if (.not. ma_pop_stack(l_hlp2))
     $     call errquit('moints: failed to pop', l_hlp2)
       if (.not. ma_pop_stack(l_hlp))
     $     call errquit('moints: failed to pop', l_hlp)
#ifndef NEW_SPARSE
       if (.not. ma_pop_stack(l_ssbbt))
     $     call errquit('moints: failed to pop', l_ssbbt)
#endif
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
       if (.not. ma_pop_stack(l_xmo))
     $     call errquit('moints: failed to pop', l_xmo)
c
c
c
       return
       end








c
      subroutine moints_trf2Ky( nbf, qlo, qhi, oseglo, oseghi, 
     $                          vlo, vhi, ilo, ihi, jlo, jhi, gloc,
     $                          ssni, h, h2, c, g_buf )
      implicit none
#include "global.fh"
      integer nbf, qlo, qhi, oseglo, oseghi, vlo, vhi
      integer ilo, ihi, jlo, jhi
      integer gloc(nbf,nbf)
      double precision ssni(nbf,jlo:jhi,ilo:ihi,oseglo:oseghi)
      double precision h(jlo:jhi,ilo:ihi),h2(*)
      double precision c(nbf,qlo:qhi)
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
#ifndef NOCOMMS
            call ga_put(g_buf,glo,ghi,vo,vo,h,1)
#endif
          else
            ij = 0
            do i=ilo,ihi
              do j=ilo,i
                ij =   ij + 1
                h2(ij) = h(i,j) 
              enddo
            enddo
#ifndef NOCOMMS
            call ga_put(g_buf,glo,ghi,vo,vo,h2,1)
#endif
          endif
        enddo
      enddo
      return
      end



      subroutine moints_trf2Kynew
     $     ( nbf, qlo, qhi, oseglo, oseghi, 
     $     vlo, vhi, ilo, ihi, jlo, jhi, gloc,
     $     ssni, h, h2, ct, g_buf )
      implicit none
#include "global.fh"
      integer nbf, qlo, qhi, oseglo, oseghi, vlo, vhi
      integer ilo, ihi, jlo, jhi
      integer gloc(nbf,nbf)
      double precision ssni(nbf,jlo:jhi,ilo:ihi,oseglo:oseghi)
      double precision h(vlo:vhi,jlo:jhi,ilo:ihi),h2(*)
      double precision ct(qlo:qhi,nbf),s, sum
      integer g_buf
      integer ilen, jlen, ijlen
      integer glo, ghi, i, j, ij,k
      integer o, v, vo, jtop
c
      ilen = ihi - ilo + 1
      jlen = jhi - jlo + 1
      ijlen = ilen*jlen
      glo = gloc(ilo,jlo)
      ghi = gloc(ihi,jhi)

      do o = oseglo, oseghi
        call dfill((vhi-vlo+1)*ilen*jlen,0.0d0,h,1)
        do i = ilo, ihi
          do j = jlo, jhi
            do k = 1, nbf
              s = ssni(k,j,i,o)
              if (abs(s) .gt. 1d-12) then
                do v = vlo,vhi
                  h(v,j,i) = h(v,j,i) + s*ct(v,k)
                enddo
              endif
            enddo
          enddo
        enddo
        do v=vlo,vhi
          vo = (v-vlo)*(oseghi-oseglo+1) + o - oseglo + 1
          ij = 0
          sum = 0.0d0
          do i=ilo,ihi
            jtop = jhi
            if (ihi.eq.jhi) jtop = i
            do j=jlo,jtop
              ij =   ij + 1
              h2(ij) = h(v,j,i) 
              sum = sum + h(v,j,i)*h(v,j,i)
            enddo
          enddo
c     Note that sum is the sum of squares hence test against tol**2
          if (sum .gt. 1d-20) call ga_put(g_buf,glo,ghi,vo,vo,h2,1)
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





      subroutine moints_wrbuf( munit, oseglo, oseghi, alo, plo, phi,
     $                         nnbf, iocnt, rlen, g_buf )
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "eaf.fh"
      integer munit                       ! I/O unit
      integer oseglo, oseghi              ! Segment range
      integer alo                         ! Lowest index to compute GA index
      integer plo, phi                    ! Index range for this processor
      integer nnbf, rlen            
      double precision iocnt              ! file header offset
      integer g_buf
c
      double precision faddr
      integer o, p, po, nocc, np
      integer myid, rlo, rhi, clo, chi, k_local, ld
c
      nocc = oseghi - oseglo + 1
      np = phi - plo + 1
      myid = ga_nodeid()
      call ga_distribution(g_buf, myid, rlo, rhi, clo, chi )
      do o=oseglo,oseghi
        do p=plo,phi
          po = (p-alo)*nocc + o - oseglo + 1
          faddr = iocnt + 
     $            ma_sizeof(MT_DBL,((o-oseglo)*np+p-plo)*nnbf,MT_BYTE)
          if ((po.ge.clo).and.(po.le.chi)) then
            call ga_access(g_buf, rlo, rhi, po, po, k_local, ld )
            if (eaf_write(munit, faddr, dbl_mb(k_local), 
     $                    ma_sizeof(MT_DBL,rlen,MT_BYTE)).ne.0)
     $        call errquit('moints_wrbuf: failed on write',0)
            call ga_release(g_buf, rlo, rhi, po, po )
c$$$            WRITE(6,901) O,P,FADDR,FADDR+MA_SIZEOF(MT_DBL,RLEN,MT_BYTE)
c$$$ 901        FORMAT(' [O,P]:',2I4,5X,'ADDR:',F10.0,'  EOR:',F10.0)
          else
            call errquit('moints_semi: wrong distrib for GA buffer',0)
          endif
        enddo
      enddo
      return
      end
          






      subroutine moints_readintK( nbf, oseglo, oseghi, olo, ohi,
     $                            vlo, vhi, c, orbe, g_epair )
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "eaf.fh"
#include "util.fh"
      integer MSG_SEMIMP2_SUM
      parameter(MSG_SEMIMP2_SUM=10241)
      integer nbf, oseglo, oseghi, olo, ohi, vlo, vhi
      double precision c(nbf,nbf)
      double precision orbe(nbf)
      integer g_epair
c
      integer nseg, nocc, nvir, nnbf, ioff, v, o, oo
      integer pvlo, pvhi, itmp
      integer l_v, k_v, l_t, k_t, l_x, k_x
      integer j, oj, jj, vvlo, vvhi, k_local, ld, myid
      integer rlo, rhi, clo, chi
      integer gtype, dim1, oomax
      integer g_exch
      double precision tmp2
      character*256 fname
      integer kunit
      double precision xx, denom, e2
      logical status
      double precision moints_epair_eval
      external moints_epair_eval

      tmp2 = util_cpusec()
      call util_file_name('kh',.true.,.true.,fname)
      if (eaf_open(fname, EAF_RW, kunit).ne.0)
     $  call errquit('moints_readintK: cannot open half int file',0)
      myid = ga_nodeid()
      e2 = 0.d0
      nocc = ohi - olo + 1
      nvir = vhi - vlo + 1
      nseg = oseghi - oseglo + 1
      call moints_vrange( kunit, pvlo, pvhi, nnbf, ioff )
      status = ma_push_get(MT_INT,(nnbf+mod(nnbf,2)),'scatter',l_v,k_v)
      call moints_getscattv( kunit, nnbf, int_mb(k_v) )
c$$$       do i=1,nnbf
c$$$         write(6,*) 'k_v',int_mb(k_v+i-1)
c$$$       enddo
      if (.not.ga_create( MT_DBL, (nvir*nvir), nocc, 'exch',
     $                    (nvir*nvir), 1, g_exch))
     $    call errquit('moints_semimp2: cannot allocate exch',0)
      call ga_distribution(g_exch, myid, rlo, rhi, clo, chi )
      call ga_inquire(g_epair,gtype,dim1,oomax)
c
c
c$$$      WRITE(6,934) GA_NODEID(),PVLO,PVHI,CLO,CHI
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
        do v=pvlo,pvhi
          call dfill((nbf*nbf), 0.d0, dbl_mb(k_x), 1 )
          call moints_rdhfint( kunit, pvlo, pvhi, oo, v, nnbf,
     $                        ioff, dbl_mb(k_t))
          call scatter( nnbf, dbl_mb(k_x), int_mb(k_v), dbl_mb(k_t) )
          call upper2square( nbf, dbl_mb(k_x), dbl_mb(k_x) )
          call dgemm( 'n', 'n', nbf, nocc, nbf, 1.d0, dbl_mb(k_x),
     $               nbf, c(1,olo), nbf, 0.d0, dbl_mb(k_t), nbf )
          call dgemm( 't', 'n', nvir, nocc, nbf, 1.d0, c(1,vlo),
     $               nbf, dbl_mb(k_t), nbf, 0.d0, dbl_mb(k_x), nvir )
          vvlo = (v-vlo)*nvir + 1
          vvhi = (v-vlo+1)*nvir
          do j=olo,o
            oj = j-olo+1
            jj = (j-olo)*nvir
            call ga_put( g_exch, vvlo, vvhi, oj, oj,
     $                   dbl_mb(k_x+jj), nvir)
          enddo
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
            xx = moints_epair_eval( nvir, 0, nvir, dbl_mb(k_local),
     $                              orbe(ohi+1), denom )
            call ga_release(g_exch, rlo, rhi, jj, jj )
            if (o.eq.j) xx = xx*0.5d0
            e2 = e2 + xx
            oj = ((o-olo+1)*(o-olo))/2 + jj
            call ga_put(g_epair,1,1,oj,oj,xx,1)
	    
c$$$            WRITE(6,922) GA_NODEID(),O,J,XX
c$$$ 922        FORMAT(I3,' %%%%%% ',2I5,5X,F16.10)
c$$$	        CALL UTIL_FLUSH(6)

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
      if (eaf_close(kunit).ne.0)
     $  call errquit('moints_readintK: failed to close file',kunit)
c
      call ga_sync()
      
      return
      end





c
c  Read half-integrals from local file
c
      subroutine moints_rdhfint( munit, plo, phi, o, p, nnbf, ioff, t )
      implicit none
#include "mafdecls.fh"
#include "eaf.fh"
      integer munit
      integer plo, phi, o, p, nnbf
      integer ioff                       ! header offset in bytes
      double precision t(nnbf)
      integer np, rlen, recnum
      double precision faddr

      np = phi - plo + 1
      recnum = (o-1)*np + p - plo
      rlen = ma_sizeof(MT_DBL,nnbf,MT_BYTE)
      faddr = ioff + recnum*rlen
      if (eaf_read(munit, faddr, t, rlen).ne.0)
     $  call errquit('moints_rdhfint:cannot read integral record',0)

c$$$      WRITE(6,771) RECNUM, FADDR
c$$$ 771  FORMAT(' READ RECORD#:',I8,' @ ',F10.0)

      return
      end




c
c  Read header on the local file and
c  return index ranges, offsets and 
c  number of AO indices
c

      subroutine moints_vrange( munit, plo, phi, nnbf, ioff )
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "eaf.fh"
      integer munit
      integer plo, phi, nnbf, ioff
      integer vtmp(10)
c$$$      INTEGER I

      if (eaf_read(munit, 0.d0, vtmp, 
     $             ma_sizeof(MT_INT,10,MT_BYTE)).ne.0)
     $  call errquit('moints_vrange: cannot read header info',0)
      ioff = vtmp(1)
      nnbf = vtmp(2)
      plo = vtmp(3)
      phi = vtmp(4)

c$$$      IF (GA_NODEID().EQ.0) THEN
c$$$        PRINT*,'READ BACK INITIAL PARAMETERS'
c$$$        WRITE(6,771) (VTMP(I),I=1,10)
c$$$ 771    FORMAT(16I4)
c$$$        CALL UTIL_FLUSH(6)
c$$$      ENDIF
c$$$      CALL GA_SYNC()

      return
      end


c
c  Recover the scatter array from local file header.
c  This scatter array maps the dense compound
c  AO index (mu nu) to square array (mu,nu). 
c  Includes sparsity by Schwarz screening and petite list.
c
      subroutine moints_getscattv( munit, nnbf, v )
      implicit none
#include "mafdecls.fh"
#include "eaf.fh"
      integer munit
      integer nnbf, v(*)
      integer itmp
      double precision ioptr
c$$$      INTEGER I
      
      itmp = nnbf+mod(nnbf,2)

      ioptr = ma_sizeof(MT_INT,10,MT_BYTE)
      if (eaf_read(munit, ioptr, v, 
     $             ma_sizeof(MT_INT,nnbf,MT_BYTE)).ne.0)
     $  call errquit('moints_getscattv: cannot read scatter v',0)

c$$$      PRINT*,'READ BACK SCATTER ARRAY'
c$$$      WRITE(6,771) (V(I),I=1,NNBF)
c$$$ 771  FORMAT(16I4)

      return
      end







c
c  Units in bytes 
c
       integer function moints_lmem( basis, nocc, nvir, blen )
       implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "util.fh"
       integer basis, nocc, nvir, blen
c
       integer nsh, nbf, maxbfsh
       integer bsize, ngrp, imax2e, imem2, memi, memd
       logical status, oprint_mem
       integer moints_numgr
       external moints_numgr
c
       oprint_mem = util_print('memory',print_high) .and.
     $      ga_nodeid().eq.0
c
       status = bas_numbf(basis,nbf)
       status = status.and.bas_numcont(basis,nsh)
       status = status.and.bas_nbf_cn_max(basis,maxbfsh)
       bsize = max(blen,maxbfsh)
       ngrp = moints_numgr( basis, blen )
       call int_mem_2e4c(imax2e, imem2)
       memi = nsh + 4*ngrp + nbf*nbf + nsh*(nsh+1) +
     $        2*ga_nnodes() + (nbf*(nbf+1))/2
       memd = 2*nbf*nbf + imax2e + imem2 + 
     $        max((nocc*nocc),(maxbfsh*maxbfsh*bsize*bsize)) +
     $        max(nocc,maxbfsh*maxbfsh)*nbf +
     $        maxbfsh*maxbfsh +
     $        maxbfsh*maxbfsh*nbf*nocc
c
c     Add on the ubiquitous 10% 
c
       moints_lmem = (memd + ma_sizeof(MT_INT, memi, MT_DBL))*1.1
c
       if (oprint_mem) then
          write(6,*)
          write(6,1) ' Integer Lmemory ',
     $         nsh , 4*ngrp , nbf*nbf , nsh*(nsh+1) ,
     $        2*ga_nnodes() , (nbf*(nbf+1))/2
 1        format(a,1x,10i8)
          write(6,1) ' Real Lmemory    ',
     $        2*nbf*nbf , imax2e , imem2 , 
     $        max((nocc*nocc),(maxbfsh*maxbfsh*bsize*bsize)) ,
     $        max(nocc,maxbfsh*maxbfsh)*nbf +
     $        maxbfsh*maxbfsh ,
     $        maxbfsh*maxbfsh*nbf*nocc
          write(6,1) ' Total Lmemory   ',
     $         int((memd + ma_sizeof(MT_INT, memi, MT_DBL))*1.1)
          call util_flush(6)
       endif
c
       end

      subroutine moints_xmot(nbf,nq,xmo, xmo_t)
      implicit none
      integer nbf, nq
      double precision xmo(nbf,nq), xmo_t(nq,nbf)
      integer i,j
      do i = 1, nbf
         do j = 1, nq
            xmo_t(j,i) = xmo(i,j)
         enddo
      enddo
      end

            



