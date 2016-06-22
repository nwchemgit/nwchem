#ifdef SOLARIS
#define NEW_SPARSE
#endif
#define GAADD 1
c
c  $Id$
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
       subroutine moints_semi( rtdb, basis, tol2e, oseglo, oseghi, 
     $                         olo, ohi, vlo, vhi, g_movecs, oblk ,
     K     k_file_size)
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "eaf.fh"
#include "util.fh"
#include "sym.fh"
#include "bas.fh"
#include "rtdb.fh"
c
       integer rtdb
       integer basis                          ! [input]  Basis handle
       double precision tol2e                 ! [input]  Integral tolerance
       integer oseglo, oseghi                 ! [input]  Occupied segment range
       integer olo, ohi                       ! [input]  Correlated occupied index range
       integer vlo, vhi                       ! [input]  Virtual index range
       integer g_movecs                       ! [input]  MO coefficients
       logical oblk                           ! [input]  Toggle AO integral blocking
c
       integer l_shmap, k_shmap, l_bfmap, k_bfmap, l_rbfmap, k_rbfmap
       integer l_glo, k_glo, l_ghi, k_ghi
       integer l_gloc, k_gloc
       integer l_ionext, k_ionext
       integer l_rlen, k_rlen
       integer l_rloc, k_rloc
       integer nbf, nsh, maxbfsh, nseg, nocc, nvir
       integer ngrp, blen, nnbf
       integer pvlo, pvhi
       integer niofl
       integer myid, numnodes, rlo, rhi, clo, chi
       integer itmp, vtmp(10), mxrlen,mxrlen_in
       integer lmemreq, nvir_node
       double precision gmem, gmem_agg
       double precision iotmpptr
       integer g_kbuf,g_kbuf_trans
       logical status
       logical osym
       character*256 fnameK
       integer kunit
       logical oprint_mem
c
       integer moints_numgr
       integer moints_lmem
       integer eaftype,eaf_size_in_mb,inntsize
       double precision k_file_size
       integer mp2_eaftype
       external mp2_eaftype
       external moints_numgr
       external moints_lmem
       logical mp2cpybck
       data osym/.false./
c
c  General info
c  Compare definition for nseg      : batch range
c                         nocc      : number of occupieds to highest i
c
       myid = ga_nodeid()
       numnodes = ga_nnodes()
       inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
       nocc = ohi - olo + 1
       nseg = oseghi - oseglo + 1
       nvir = vhi - vlo + 1
       status = bas_numbf(basis,nbf)
       status = status.and.bas_numcont(basis,nsh)
       status = status.and.bas_nbf_cn_max(basis,maxbfsh)
       if (.not.status) call errquit('moints: cannot get basis info',0,
     &       BASIS_ERR)
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
*       blen = 32
 1212  continue
c
c  Memory arithmetic
c
       lmemreq = moints_lmem(basis, nseg, nvir, blen)
       if ((.not.ga_uses_ma()).and.(.not.(ga_memory_limited())))
     $   call errquit('moints_semi: cannot determine memory limit',0,
     &       GA_ERR)
c
       gmem = ga_memory_avail() / 8 ! Want amount in doubles
       call ga_dgop(313, gmem, 1, 'min')
       if (oprint_mem) then
          write(6,717) blen, gmem, dble(lmemreq)
 717      format(/'   Block length           ', i8/
     $            '   Local GA available     ',  1p, d10.2/
     $            '   Local memory required  ', 1p, d10.2)
          call util_flush(6)
       endif
c
       if (ga_uses_ma()) gmem = gmem - lmemreq
c
       if (oprint_mem) then
          write(6,718) gmem
 718      format('   Adjusted GA available  ', 1p, d10.2)
          call util_flush(6)
       endif
c
       gmem_agg = gmem*numnodes
c
       if (oprint_mem) then
          write(6,719) gmem_agg
 719      format('   Aggregate GA available ', 1p, d10.2)
          call util_flush(6)
       endif
c
       nvir_node = int(nvir/numnodes)
       if (mod(nvir,numnodes).ne.0) nvir_node=nvir_node+1
c     we have g_kbuf and g_kbuf trans plus ga_add duplicates .. /3
c     we have g_kbuf and g_kbuf trans 
       gmem=gmem/2
       mxrlen = int(gmem/(nseg*nvir_node))
       mxrlen = min(
     $      mxrlen, 
     $      nbf*nbf)
C     $      min(nbf*nbf,maxbfsh*maxbfsh*numnodes))

       if (rtdb_get(rtdb, 'mp2:mxrlen', mt_int, 1, mxrlen_in)) 
     c      mxrlen=max(mxrlen_in*1024,mxrlen)
 1964  continue
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
     $         mxrlen*nseg*nvir, GA_ERR)
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
      itmp = nvir/numnodes + min(mod(nvir,numnodes),1)
      gmem=0
*ga:1:0
      if (.not.(ga_create( MT_DBL, mxrlen*nseg, nvir, 'transp Kbuf',
     $                     (mxrlen*nseg), 0, g_kbuf ))) gmem=-1
       call ga_dgop(319, gmem, 1, 'min')
       call ga_sync()
       if(gmem.lt.0) then
c     cant alloc ga .. reduce mxrlen
         mxrlen=mxrlen/2
         if(ga_nodeid().eq.0) write(6,1962)
 1962    format(//,5x,' Halved record length',//)
         goto 1964
cold    call errquit('moints_semi: cannot allocate K buffer',0,
cold     &       GA_ERR)
         endif
      call ga_distribution(g_kbuf, myid, rlo, rhi, clo, chi )
      if ((clo.eq.0).and.(chi.eq.-1)) then
        pvlo = 0
        pvhi = -1
      else
        pvlo = vlo + clo - 1
        pvhi = vlo + chi - 1
      endif
#ifdef GAADD
      gmem=0
      if (.not.(ga_create( MT_DBL, mxrlen*nseg, nvir, ' Kbuf',
     $                     0, nvir, g_kbuf_trans ))) gmem=-1
       call ga_dgop(321, gmem, 1, 'min')
       call ga_sync()
       if(gmem.lt.0) then
c     cant alloc ga .. reduce mxrlen
         mxrlen=mxrlen/2
         if (.not. ga_destroy(g_kbuf))
     $        call errquit('moints_semi: failed to destroy g_kbuf',
     G        g_kbuf,   GA_ERR)
         if(ga_nodeid().eq.0) write(6,1962) 
         goto 1964
cold     $   call errquit('moints_semi: cannot allocate K buffer',0,
cold     &       GA_ERR)
      endif
#else
      g_kbuf_trans=g_kbuf
#endif
c
c  Reorder shells by descending shell-length
c  and group shells by blocksize
c
       status = .true.
       status = status .and. ma_push_get(MT_INT,nsh,'shell order map',
     $                                   l_shmap, k_shmap)
       status = status .and. ma_push_get(MT_INT,nsh,'group lo',
     $                                   l_glo, k_glo )
       status = status .and. ma_push_get(MT_INT,nsh,'group hi',
     $                                   l_ghi, k_ghi)
       status = status .and. ma_push_get(MT_INT,nbf,'basis map',
     $                                   l_bfmap, k_bfmap)
       status = status .and. ma_push_get(MT_INT,nbf,'rev basis map',
     $                                   l_rbfmap, k_rbfmap)
       if (.not.status) then
         call errquit("moints_semi: failed to allocate buffers shorder",
     $                0,MA_ERR)
       endif
       call moints_shorder( basis, nsh, nbf, blen, ngrp,
     $                      int_mb(k_glo), int_mb(k_ghi),
     $                      int_mb(k_shmap),
     $                      int_mb(k_bfmap), int_mb(k_rbfmap) )
c
c  Generate locator and reverse map
c
       itmp = (nsh*(nsh+1))/2+numnodes
       status = .true.
       status = status .and. ma_push_get(MT_INT, (nbf*nbf), 'loc',
     $                                   l_gloc, k_gloc )
       status = status .and. ma_push_get(MT_INT, itmp, 'io',
     $                                   l_ionext, k_ionext )
       status = status .and. ma_push_get(MT_INT, itmp, 'rec len',
     $                                   l_rlen, k_rlen )
       itmp = (nbf*(nbf+1))/2
       itmp = itmp + mod(itmp,2)
       status = status .and. ma_push_get(MT_INT, itmp, 'rev loc',
     $                                   l_rloc, k_rloc )
       if (.not.status) then
         call errquit("moints_semi: failed to allocate buffers locmap",
     $                0,MA_ERR)
       endif
       call ifill(itmp, 0, int_mb(k_rloc), 1 )
       call moints_locmap( basis, nsh, nbf, tol2e, int_mb(k_shmap),
     $                     mxrlen, int_mb(k_gloc), nnbf,
     $                     int_mb(k_rloc), niofl, int_mb(k_ionext),
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
#ifdef NOIO
       eaftype=mp2_eaftype(rtdb,k_file_size)
       if(ga_nodeid().eq.0) write(6,*) ' mp2_eaf for mointskh  ',eaftype
#else
       eaftype=EAF_RW
#endif
       if (eaf_open( fnamek, eaftype, kunit).ne.0)
     $   call errquit('moints_semi: failed to open half int file',0,
     &       DISK_ERR)
       if (eaf_write(kunit, 0.d0, vtmp,
     $               ma_sizeof(MT_INT,10,MT_BYTE)).ne.0)
     $   call errquit('moints_semi: failed to write header info1',0,
     &       DISK_ERR)
       iotmpptr = ma_sizeof(MT_INT,10,MT_BYTE)
       if (eaf_write(kunit, iotmpptr, int_mb(k_rloc),
     $               ma_sizeof(MT_INT,nnbf,MT_BYTE)).ne.0)
     $   call errquit('moints_semi: failed to write header info2',0,
     &       DISK_ERR)
       if (.not. ma_pop_stack(l_rloc))
     $     call errquit('moints: failed to pop', 1, MA_ERR)
c
c  Call the real stuff
c
        if (.not. rtdb_get(rtdb, 'mp2:copyback', mt_log, 1, mp2cpybck))
     $    mp2cpybck=.false.

       call moints_semi_a( basis, nbf, nsh, maxbfsh, tol2e, oseglo,
     $                   oseghi, olo, ohi, vlo, vhi, pvlo, pvhi,
     $                   mxrlen, osym, oblk, kunit, int_mb(k_shmap), 
     $                   ngrp, int_mb(k_glo), int_mb(k_ghi),
     $                   int_mb(k_bfmap), int_mb(k_rbfmap),
     $                   blen, int_mb(k_gloc), 
     $                   nnbf, niofl, int_mb(k_ionext),
     $                   int_mb(k_rlen), g_movecs, g_kbuf,g_kbuf_trans,
     M                   mp2cpybck)
c
c  Clean up
c
       if (.not. ma_pop_stack(l_rlen))
     $     call errquit('moints_semi: failed to pop', 2, MA_ERR)
       if (.not. ma_pop_stack(l_ionext))
     $     call errquit('moints_semi: failed to pop', 3, MA_ERR)
       if (.not. ma_pop_stack(l_gloc))
     $     call errquit('moints_semi: failed to pop', 4, MA_ERR)
       if (.not. ma_pop_stack(l_rbfmap))
     $     call errquit('moints_semi: failed to pop', 6, MA_ERR)
       if (.not. ma_pop_stack(l_bfmap))
     $     call errquit('moints_semi: failed to pop', 7, MA_ERR)
       if (.not. ma_pop_stack(l_ghi))
     $     call errquit('moints_semi: failed to pop', 8, MA_ERR)
       if (.not. ma_pop_stack(l_glo))
     $     call errquit('moints_semi: failed to pop', 9, MA_ERR)
       if (.not. ma_pop_stack(l_shmap))
     $     call errquit('moints_semi: failed to pop', 10, MA_ERR)
       if (.not. ga_destroy(g_kbuf))
     $      call errquit('moints_semi: failed to destroy g_kbuf',g_kbuf,
     &       GA_ERR)
#ifdef GAADD
       if (.not. ga_destroy(g_kbuf_trans))
     $      call errquit('moints_semi: failed to destroy g_kbuf',g_kbuf,
     &       GA_ERR)
#endif
c
       if (eaf_close(kunit).ne.0)
     $   call errquit('moints_semi: failed to close file',kunit,
     &       DISK_ERR)
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
     $                           pvlo, pvhi, mxrlen, osym, oblk, 
     $                           kunit, shmap, ngrp, glo, ghi, bfmap, 
     $                           rbfmap, blen, gloc, nnbf, 
     $                           niofl, ionext, rlen, g_movecs, g_kbuf,
     G     g_kbuf_trans, mp2cpybck)
       implicit none
#include "errquit.fh"
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
       integer pvlo, pvhi                     ! [input] Virtual index range for this processor
       integer mxrlen                         ! [input] Max record length
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
       integer gloc(nbf*nbf)                  ! [input] bf -> GA memory map
       integer nnbf                           ! [input] Number of screened basf pairs <= (nbf*(nbf+1))/2
       integer niofl                          ! [input] Total IO flushes required
       integer ionext(niofl)                  ! [input] Next value for each IO flush
       integer rlen(niofl)                    ! [input] Record lengths for IO writes
       integer g_movecs                       ! [input] MO coefficients
       integer g_kbuf                         ! [scratch] Buffer for half-trans exchange
       integer g_kbuf_trans                         ! [scratch] Buffer for half-trans exchange
       logical mp2cpybck
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
c$$$       integer l_ssbbt, k_ssbbt
       integer l_hlp, k_hlp, l_ssni,k_ssni
       integer l_hlp2, k_hlp2
#ifndef NEW_SPARSE
       integer l_hlp3, k_hlp3
#endif
       integer l_eri, k_eri, l_iscr,k_iscr
       integer l_xmo, k_xmo, l_xmo_t, k_xmo_t
       integer n_ssbb, n_ssbb1, n_ssni, n_hlp, n_hlp2, n_hlp3,
     I      n_hlpx
       integer iz, jz, kz, bsize
       integer mem2, max2e, jopass
       double precision iocnt
       integer num_nodes, ploop, next, iiofl
       double precision schw_ij
       double precision tpass, nwrbytes, iorate
       double precision t0, t1, ttask
       logical st
c
#include "moints_stats.fh"
c       
       integer gr_len, nxtask
       external gr_len, nxtask
c
       call ga_zero(g_kbuf)
       call ga_zero(g_kbuf_trans)
       num_nodes = ga_nnodes()
       nmo1 = oseghi - oseglo + 1
       nmo2 = vhi - vlo + 1
       qlo = olo
       qhi = vhi
       nmo = qhi - qlo + 1
c
c$$$       WRITE(6,221) GA_NODEID(),PVLO,PVHI
c$$$ 221   FORMAT('ME:',I5,5X,' V-RANGE:',2I3)
c
c  Local MO coefficients
c  Integrals and temporary arrays sizes
c
*     call int_mem_2e4c(max2e, mem2)  ! For non-blocking code
       call intb_mem_2e4c(max2e, mem2) ! Determine mem2 = scratch space
       max2e = max(max2e,min(50*maxbfsh**4,21**4)) ! Enuf room for 1 cartesian H shell
c
       bsize = max(blen,maxbfsh)
       n_ssbb = maxbfsh*maxbfsh*bsize*bsize
       n_ssbb1 = max((nmo1*nmo1),n_ssbb)

       n_hlp = max(nmo1,maxbfsh*maxbfsh)*nbf
       n_hlp2 = maxbfsh*maxbfsh*nmo2

       n_hlp3 = maxbfsh*maxbfsh*nmo2

       n_ssni = maxbfsh*maxbfsh*nbf*nmo1
c
c     Get the MOS, reorder the basis functions, transpose
c
       if (.not. ma_push_get(MT_DBL,nbf*vhi,'xmot',l_xmo_t,k_xmo_t))
     $     call errquit('moints_semi: failed to alloc mosT', nbf*vhi,
     &       MA_ERR)
c
       if (.not. ma_push_get(MT_DBL,(nbf*vhi),'r movecs',l_xmo,k_xmo))
     $     then
      write(6,*) ' avail dbles',ma_inquire_avail(MT_DBL)
       call util_flush(6)
          call errquit('moints_semi: failed to alloc mos', nbf*vhi,
     &       MA_ERR)
      endif
c
#if 0
       call ga_get(g_movecs, 1, nbf, qlo, qhi, dbl_mb(k_xmo_t), nbf )
#else
       call util_mygabcast2(g_movecs,1, nbf, qlo, qhi, 
     D      dbl_mb(k_xmo_t), nbf )
#endif
       call row_exch( nbf, nmo, rbfmap, dbl_mb(k_xmo_t), dbl_mb(k_xmo) )
       call moints_xmot(nbf,qhi-qlo+1,dbl_mb(k_xmo),dbl_mb(k_xmo_t))

c     Only the transposed MOs are needed in the inner loop
c
       if (.not. ma_pop_stack(l_xmo)) 
     $      call errquit('moints_semi: pop failed', 5551212, MA_ERR)

       if (.not. ma_push_get(MT_DBL,n_ssni,'ssni blk',l_ssni,k_ssni))
     $      call errquit('mp2_semi: failed to alloc ssni ', n_ssni,
     &       MA_ERR)
c
c Initialize
c
       ploop = 0
       iiofl = 1
       jopass = 1
       iocnt = ma_sizeof(MT_INT,10+nnbf,MT_BYTE)
       next = nxtask(num_nodes, 1)
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
               mi_ntasks = mi_ntasks + 1
               ttask = util_cpusec()
               call dfill((ilen*jlen*nbf*nmo1),0.d0,dbl_mb(k_ssni),1)
c
c  Half-tranformed Integral generation
c
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

                   if (.not. ma_push_get(MT_DBL,n_ssbb1,'ssbb blk',
     $                  l_ssbb,k_ssbb)) then
                      write(6,*) ga_nodeid(),
     C                     ' avail dbles',ma_inquire_avail(MT_DBL)
                      call util_flush(6)
                      if(ga_nodeid().eq.0)
     M                     call ma_summarize_allocated_blocks()
                      
                      call ga_sync()
                       call errquit('moints_semi: ssbb', n_ssbb1,
     &       MA_ERR)
                    endif
                   if (.not. ma_push_get(MT_DBL, max2e, 'ibuf', 
     $                  l_eri, k_eri)) 
     $                  call errquit('moints_semi: eri', max2e,
     &       MA_ERR)
                   if (.not. ma_push_get(MT_DBL, mem2, 'int scr', 
     $                  l_iscr, k_iscr))
     $                  call errquit('moints_semi: scr', mem2,
     &       MA_ERR)

                   call moints_gblk( basis, ish, jsh,
     $                               kshlo, kshhi, lshlo, lshhi,
     $                               shmap, rbfmap,
     $                               schw_ij, tol2e, osym, oblk,
     $                               max2e, dbl_mb(k_eri),
     $                               mem2, dbl_mb(k_iscr),
     $                               ibflo, ibfhi, jbflo, jbfhi,
     $                               kbflo, kbfhi, lbflo, lbfhi,
     $                               dbl_mb(k_ssbb) )


                  if (.not. ma_pop_stack(l_iscr)) 
     $                 call errquit('moints_semi: pop failed ',14,
     &       MA_ERR)
                  if (.not. ma_pop_stack(l_eri))
     $                 call errquit('moints_semi: pop failed ',15,
     &       MA_ERR)

                   mi_tint = mi_tint + util_cpusec() - t0
c
                   t0 = util_cpusec()
                   mi_flop1 = mi_flop1 + 
     $                   4.d-6*ilen*jlen*kblen*lblen*(oseghi-oseglo+1)

       n_hlpx = (oseghi-oseglo+1)*max(kbfhi-kbflo+1,lbfhi-lbflo+1)
                   if (.not. ma_push_get(MT_DBL,n_hlpx,'hlp block',
     $                  l_hlp,k_hlp)) 
     $                  call errquit('moints_semi: hlp', n_hlp,
     &       MA_ERR)

                   call moints_trf1_new(nbf, qlo, qhi, oseglo, oseghi,
     $                  ilen, jlen, kbflo, kbfhi,
     $                  lbflo, lbfhi,
     $                  dbl_mb(k_ssbb),
     $                  dbl_mb(k_xmo_t),
     $                  dbl_mb(k_ssni), dbl_mb(k_hlp))

                   if (.not. ma_pop_stack(l_hlp)) 
     $                 call errquit('moints_semi: pop failed ',16,
     &       MA_ERR)

                  if (.not. ma_pop_stack(l_ssbb)) 
     $                 call errquit('moints_semi: pop failed ',17,
     &       MA_ERR)

                   mi_t1 = mi_t1 + util_cpusec() - t0
                  enddo
               enddo
c
c  Transform 2nd index and put into buffer
c
               t0 = util_cpusec()
               
               if (.not. ma_push_get(MT_DBL,n_hlp2,'hlp2 block',
     $              l_hlp2,k_hlp2)) 
     $              call errquit('moints_semi: hlp2', n_hlp2, MA_ERR)

#ifdef NEW_SPARSE
               if (.not. ma_push_get(MT_DBL,n_hlp,'hlp block',
     $              l_hlp,k_hlp)) 
     $              call errquit('moints_semi: hlp', n_hlp, MA_ERR)

               call moints_trf2Kynew(nbf, qlo, qhi, oseglo, oseghi, 
     $                               vlo, vhi, ibflo, ibfhi, 
     $                               jbflo, jbfhi, rlen(iiofl), 
     $                               gloc, dbl_mb(k_ssni),
     $                               dbl_mb(k_hlp), dbl_mb(k_hlp2),
     $                               dbl_mb(k_xmo_t), g_kbuf_trans)

               if (.not. ma_pop_stack(l_hlp)) 
     $                 call errquit('moints_semi: pop failed ',18,
     &       MA_ERR)
#else
               if (.not. ma_push_get(MT_DBL,n_hlp3,'hlp3 block',
     $              l_hlp3,k_hlp3)) 
     $              call errquit('moints_semi: hlp3', n_hlp3,
     &       MA_ERR)
               if (.not. ma_push_get(MT_DBL,(nbf*vhi),'r movecs',
     $              l_xmo,k_xmo)) then
      write(6,*) ' avail dbles',ma_inquire_avail(MT_DBL)
       call util_flush(6)
                     call errquit
     $              ('moints_semi: failed to alloc mos3', nbf*vhi,
     &       MA_ERR)
             endif
c
c     Generate the untranposed mos by transposing the transpose
c
               call moints_xmot(qhi-qlo+1,nbf,dbl_mb(k_xmo_t),
     $              dbl_mb(k_xmo))
c
               call moints_trf2Ky(nbf, qlo, qhi, oseglo, oseghi, 
     $              vlo, vhi, ibflo, ibfhi, 
     $              jbflo, jbfhi, rlen(iiofl), gloc, 
     $              dbl_mb(k_ssni),
     $              dbl_mb(k_hlp3), dbl_mb(k_hlp2),
     $              dbl_mb(k_xmo), g_kbuf_trans)

               if (.not. ma_pop_stack(l_xmo)) 
     $                 call errquit('moints_semi: pop failed ',180,
     &       MA_ERR)
               if (.not. ma_pop_stack(l_hlp3)) 
     $                 call errquit('moints_semi: pop failed ',18,
     &       MA_ERR)
#endif
               if (.not. ma_pop_stack(l_hlp2)) 
     $                 call errquit('moints_semi: pop failed ',19,
     &       MA_ERR)
               
               
               mi_t2k = mi_t2k + util_cpusec() - t0
               
               ttask = util_cpusec() - ttask
               mi_mintask = min(mi_mintask,ttask)
               mi_maxtask = max(mi_maxtask,ttask)
               mi_aggtask = mi_aggtask + ttask
c
               next = nxtask(num_nodes, 1)
c
c  End parallel task
c  -----------------
c
             endif
             ploop = ploop + 1
c
c  Checkpoint for I/O
c
             if (ploop.eq.ionext(iiofl)) then
               t0 = util_cpusec()
               call ga_mask_sync(.true.,.false.)
               call ga_sync()
               t1 = util_cpusec() - t0
               mi_nsynchs = mi_nsynchs + 1
               call summaxmin(t1, mi_aggsynch, mi_maxsynch, mi_minsynch)
               t0 = util_cpusec()
#ifdef GAADD
               if (mp2cpybck) then
                  call mp2_copyback(g_kbuf,g_kbuf_trans)
               else
                  call ga_copy(g_kbuf_trans,g_kbuf)
               endif
C               call ga_mask_sync(.true.,.false.)
               call ga_sync()
c               call ga_add(1d0,g_kbuf_trans,1d0,g_kbuf,g_kbuf)

#endif
               call moints_wrbuf( kunit, oseglo, oseghi,
     $                            vlo, pvlo, pvhi, nnbf, iocnt,
     $                            rlen(iiofl), g_kbuf)
               mi_tio = mi_tio + util_cpusec() - t0
               nwrbytes = nwrbytes + nmo1*(pvhi-pvlo+1)*rlen(iiofl)
               call ga_mask_sync(.false.,.false.)
               call ga_zero(g_kbuf)
               call ga_mask_sync(.false.,.true.)
               call ga_zero(g_kbuf_trans)
               iocnt = iocnt + ma_sizeof(MT_DBL,rlen(iiofl),MT_BYTE)
               iiofl = iiofl + 1
             endif
           endif
         enddo
       enddo
       t0 = util_cpusec()
cedo       call ga_sync
       next = nxtask(-num_nodes, 1)
       t1 = util_cpusec() - t0
       mi_nsynchs = mi_nsynchs + 1
       call summaxmin(t1, mi_aggsynch, mi_maxsynch, mi_minsynch)
       tpass = util_cpusec() - tpass

c
c I/O stats
c
       nwrbytes = nwrbytes*8.e-6
       if(mi_tio.eq.0) then
          iorate=0d0
       else
          iorate = nwrbytes/mi_tio
       endif
       call ga_dgop(MSG_SEMI_IO_RATE,iorate,1,'+')
       if ((ga_nodeid().eq.0).and.
     $     (util_print('io_stats',print_default))) then
         write(6,886) nwrbytes, mi_tio, iorate
 886     format('Node 0 wrote ',f8.1,' Mb in ',f8.1,' s',5x,
     $          'Agg I/O rate:', f8.1, ' Mb/s')
         call util_flush(6)
       endif
c
c Clean-up
c
       if (.not. ma_pop_stack(l_ssni))
     $      call errquit('moints_semi: pop failed ',20, MA_ERR)
c 
c  Complete in-core section
c
       if (.not. ma_pop_stack(l_xmo_t))
     $      call errquit('moints_semi: pop failed ',22, MA_ERR)
c
       end



      subroutine moints_trf2Ky( nbf, qlo, qhi, oseglo, oseghi, 
     $                          vlo, vhi, ilo, ihi, jlo, jhi, rlen,
     $                          gloc, ssni, h, h2, c, g_buf )
      implicit none
#include "global.fh"
      integer nbf, qlo, qhi, oseglo, oseghi, vlo, vhi
      integer ilo, ihi, jlo, jhi
      integer rlen
      integer gloc(nbf,nbf)
      double precision ssni(nbf,jlo:jhi,ilo:ihi,oseglo:oseghi)
      double precision h(jlo:jhi,ilo:ihi,vlo:vhi),h2(*)
      double precision c(nbf,qlo:qhi)
      integer g_buf
      integer ilen, jlen, ijlen, nvir
      integer glo, ghi, i, j, ij, ijolo, ijohi
      integer o, v
      integer nij
c
      ilen = ihi - ilo + 1
      jlen = jhi - jlo + 1
      nij = (ilen*(ilen+1))/2
      ijlen = ilen*jlen
      glo = gloc(ilo,jlo)
      ghi = gloc(ihi,jhi)
      nvir = vhi - vlo + 1
      do o=oseglo,oseghi
        call dgemm('t', 'n', ijlen, nvir, nbf, 1.d0, ssni(1,jlo,ilo,o),
     $              nbf, c(1,vlo), nbf, 0.d0, h(jlo,ilo,vlo), ijlen)
        ijolo = (o - oseglo)*rlen + glo
        ijohi = (o - oseglo)*rlen + ghi
        if (ihi.eq.jhi) then
          ij = 0
          do v=vlo,vhi
            do i=ilo,ihi
              do j=ilo,i
                ij =   ij + 1
                h2(ij) = h(j,i,v)
              enddo
            enddo
          enddo
#ifndef NOCOMMS
          call ga_put( g_buf, ijolo, ijohi, 1, nvir, h2, nij )
#endif
        else
#ifndef NOCOMMS
          call ga_put( g_buf, ijolo, ijohi, 1, nvir, h, ijlen )
#endif
        endif
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






      subroutine moints_wrbuf( munit, oseglo, oseghi, vlo, pvlo, pvhi,
     $                         nnbf, iocnt, rlen, g_buf )
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "inp.fh"
#include "eaf.fh"
      integer munit                       ! I/O unit
      integer oseglo, oseghi              ! Segment range
      integer vlo                         ! Lowest virtual index
      integer pvlo, pvhi                  ! Virtual index range for this processor
      integer nnbf, rlen            
      double precision iocnt              ! file header offset
      integer g_buf
c
      double precision faddr
      integer o, vv, pv, nocc, npv, ooffset
      integer myid, rlo, rhi, clo, chi, k_local, ld
#ifdef  BAD_GACCESS
      integer l_local,howmy
#endif
      integer ierr
      character*80 errmsg
c
      nocc = oseghi - oseglo + 1
      npv = pvhi - pvlo + 1
      myid = ga_nodeid()
      call ga_distribution(g_buf, myid, rlo, rhi, clo, chi )
c
      do o=1,nocc
        ooffset = (o-1)*rlen
        do pv=pvlo,pvhi	
	  vv = pv - vlo + 1
          faddr = iocnt + 
     $            ma_sizeof(MT_DBL,((o-1)*npv+pv-pvlo)*nnbf,MT_BYTE)
          if ((vv.ge.clo).and.(vv.le.chi)) then
#ifdef  BAD_GACCESS
      ld=rhi-rlo+1
      howmy=max(ld,ooffset+rlen)
      if(.not.ma_push_get(MT_DBL,howmy,
     $  'scratch buff', l_local, k_local)) call
     $  errquit('mointsemi: pushget failed',howmy,0)
            call ga_get(g_buf,rlo,rhi,vv,vv,dbl_mb(k_local),ld)
#else
            call ga_access(g_buf, rlo, rhi, vv, vv, k_local, ld )
#endif
            ierr=eaf_write(munit, faddr, dbl_mb(k_local+ooffset), 
     $                    ma_sizeof(MT_DBL,rlen,MT_BYTE))
            if(ierr.ne.0) then
               call eaf_errmsg(ierr, errmsg)
               write(6,*) ' IO offset ', faddr
               write(6,*) ' IO error message ',
     ,              errmsg(1:inp_strlen(errmsg))
               call errquit('moints_wrbuf: failed on write',0, DISK_ERR)
            endif
#ifdef  BAD_GACCESS
      if(.not.ma_pop_stack(l_local)) call 
     $  errquit('mointsemi: popstack failed',0,0)
#else
            call ga_release(g_buf, rlo, rhi, vv, vv )
#endif
          else
            call errquit('moints_semi: wrong distrib for GA buffer',0,
     &       GA_ERR)
          endif
        enddo
      enddo
      return
      end
          








      subroutine moints_readintK( nbf, oseglo, oseghi, olo, ohi,
     $                            vlo, vhi, c, orbe, g_epair )
      implicit none
#include "errquit.fh"
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
      double precision c1,c2       ! constants for moints_epair_eval (because of SCS-MP2)

      tmp2 = util_cpusec()
c
      c1=4.0d0
      c2=2.0d0
c
      call util_file_name('kh',.true.,.true.,fname)
#ifdef NOIO
      if (eaf_open(fname, 200, kunit).ne.0)
#else
      if (eaf_open(fname, EAF_RW, kunit).ne.0)
#endif
     $  call errquit('moints_readintK: cannot open half int file',0,
     &       DISK_ERR)
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
*ga:1:0
      if (.not.ga_create( MT_DBL, (nvir*nvir), nocc, 'exch',
     $                    (nvir*nvir), 0, g_exch))
     $    call errquit('moints_semimp2: cannot allocate exch',0, GA_ERR)
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

c
c     The following call is suspect since the definition of this
c     function has 8 arguments (4 integers at the start, not 3.
c     Added dummy arg to get it to link under WIN32, but someone
c     should fix this properly.
c     BGJ (9/99)
c
c            xx = moints_epair_eval( nvir, 0, nvir, dbl_mb(k_local),
c     $                              orbe(ohi+1), denom )
            call errquit('suspect call to moints_epair_eval',0,
     &       UNKNOWN_ERR)
            xx = moints_epair_eval( nvir, 0, nvir, 0, dbl_mb(k_local),
     $                              orbe(ohi+1), denom, c1, c2 )
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
     $  call errquit('moints_semimp2: failed to destroy',g_exch, GA_ERR)
      if (.not. ma_pop_stack(l_x))
     $  call errquit('moints_readint: failed to pop', l_x, MA_ERR)
      if (.not. ma_pop_stack(l_t))
     $  call errquit('moints_readint: failed to pop', l_t, MA_ERR)
      if (.not. ma_pop_stack(l_v))
     $  call errquit('moints_readint: failed to pop', l_v, MA_ERR)
c
      if (eaf_close(kunit).ne.0)
     $  call errquit('moints_readintK: failed to close file',kunit,
     &       DISK_ERR)
c
      call ga_sync()
      
      return
      end





c
c  Read half-integrals from local file
c
      subroutine moints_rdhfint( munit, plo, phi, o, p, nnbf, ioff, t )
      implicit none
#include "errquit.fh"
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
     $  call errquit('moints_rdhfint:cannot read integral record',0,
     &       DISK_ERR)

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
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "eaf.fh"
      integer munit
      integer plo, phi, nnbf, ioff
      integer vtmp(10)
c$$$      INTEGER I

      if (eaf_read(munit, 0.d0, vtmp, 
     $             ma_sizeof(MT_INT,10,MT_BYTE)).ne.0)
     $  call errquit('moints_vrange: cannot read header info',0,
     &       DISK_ERR)
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
#include "errquit.fh"
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
     $  call errquit('moints_getscattv: cannot read scatter v',0,
     &       DISK_ERR)

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

       call intb_mem_2e4c(imax2e, imem2)
       imax2e = max(imax2e,min(50*maxbfsh**4,21**4)) ! Enuf room for 1 cartesian H shell
       
*      call int_mem_2e4c(imax2e, imem2) For non-blocking integrals
       memi = nsh + 4*ngrp + nbf*nbf + nsh*(nsh+1) +
     $        2*ga_nnodes()
       memd = nbf**2 + blen*nbf*maxbfsh**2 + 
     $      max(maxbfsh**2*blen**2+3*imax2e+imem2,
     $          2*nvir*maxbfsh**2+nbf**2)
c
c     Add on the ubiquitous 10% 
c
       moints_lmem = (memd + ma_sizeof(MT_INT, memi, MT_DBL))*1.1
c
       if (oprint_mem) then
          write(6,*)
          write(6,1) ' Integer Lmemory ',
     $         nsh , 4*ngrp , nbf*nbf , nsh*(nsh+1) ,
     $        2*ga_nnodes()
 1        format(a,1x,10i8)
          write(6,1) ' Real Lmemory    ',
     $         nbf**2, blen*nbf*maxbfsh**2,
     $         max(maxbfsh**2*blen**2+3*imax2e+imem2,
     $         2*nvir*maxbfsh**2+nbf**2)
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

#ifdef NEW_SPARSE
      subroutine moints_trf2Kynew( nbf, qlo, qhi, oseglo, oseghi, 
     $                             vlo, vhi, ilo, ihi, jlo, jhi, 
     $                             rlen, gloc, ssni, h, h2, ct,
     $                             g_buf )
      implicit none
#include "global.fh"
      integer nbf, qlo, qhi, oseglo, oseghi, vlo, vhi
      integer ilo, ihi, jlo, jhi
      integer rlen
      integer gloc(nbf,nbf)
      double precision ssni(nbf,jlo:jhi,ilo:ihi,oseglo:oseghi)
      double precision h(vlo:vhi,jlo:jhi,ilo:ihi),h2(*)
      double precision ct(qlo:qhi,nbf)
      integer g_buf
c
      double precision s
      integer ijolo, ijohi, nij
      integer ilen, jlen, vlen
      integer glo, ghi, i, j, ij, k, jtop, klo, khi
      integer o, v, vtop, vbot
      integer kchunk, vchunk
      parameter (kchunk=32, vchunk=96) ! For cache reuse
c
      ilen = ihi - ilo + 1
      jlen = jhi - jlo + 1
      vlen = vhi - vlo + 1
      glo  = gloc(ilo,jlo)
      ghi  = gloc(ihi,jhi)
      nij = ilen*jlen
      if (ihi.eq.jhi) nij = (ilen*(ilen+1))/2

      do o = oseglo, oseghi
        ijolo = (o - oseglo)*rlen + glo
        ijohi = (o - oseglo)*rlen + ghi
        call dfill((vhi-vlo+1)*ilen*jlen, 0.0d0, h, 1 )
        do klo = 1, nbf, kchunk                       ! Blocking for cache
          khi = min(nbf,klo+kchunk-1)
          do vbot = vlo, vhi, vchunk                  ! Blocking for cache
            vtop = min(vhi,vbot+vchunk-1)
            do i = ilo, ihi
              jtop = jhi
              if (ilo .eq. jlo) jtop = i
              do j = jlo, jtop
                do k = klo, khi
                  s = ssni(k,j,i,o)
                  if (abs(s) .gt. 1d-16) then
                    do v = vbot,vtop
                      h(v,j,i) = h(v,j,i) + s*ct(v,k)
                    enddo
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
        ij = 0
        s = 0.d0
        do v=vlo,vhi
          do i=ilo,ihi
            if (ihi.eq.jhi) jtop = i
            do j=jlo,jtop
              ij = ij + 1
              h2(ij) = h(v,j,i)
              s = s + h(v,j,i)*h(v,j,i)
            enddo
          enddo
        enddo
        if (s.gt.1.d-20) call ga_put(g_buf,ijolo,ijohi,1,vlen,h2,nij)
      enddo
c
      return
      end
#endif
