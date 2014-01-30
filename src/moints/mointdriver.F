       subroutine moint_driver( rtdb, geom, basis, movecs_in )
C$Id$
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "schwarz.fh"
#include "sym.fh"
#include "util.fh"
#include "moints.fh"
c
c
c
       integer rtdb                        ! [input] database handle
       integer geom, basis
       character*255 movecs_in
c
c
       character*255 title_vecs, basis_vecs
       character*20 scftype_vecs
       double precision nuclear_charge, charge
       integer nbf_vecs, nsets
       integer k_occ, l_occ, k_evals, l_evals, g_vecs
*skip       integer k_epair, l_epair
       integer k_x, l_x
*skip       integer npairs
       integer nbf, nsh, maxbfsh, nelec
       integer nact, nocc, nopen, nocc_correl
       integer occ_frozen, occ_lo, occ_hi
       integer vir_frozen, vir_lo, vir_hi, nmo, nmo2
       integer qlo, qhi
       integer noper, oseg_lo, oseg_hi
*skip       integer npass 
       integer blen, max_noper, nocct
       integer noper_core, seg_incr
       integer g_coul, g_exch, g_epair
       integer i, j, ij, ipass
*skip       integer ii
       integer lmemreq, gmem, mxrlen, nvir
       logical oexch, ocoul, osym, oprint, status, ohalf, omemprint
       logical oblk, odisk, oaoreuse, oaokeep
       double precision emp2, e2pair
       double precision tpass, tmp2
       double precision tol2e
       character*7 algo
c
       logical movecs_read_header, movecs_read
       external movecs_read_header, movecs_read
       integer ga_create_atom_blocked
       external ga_create_atom_blocked
       integer ga_create_JKblocked
       external ga_create_JKblocked
       logical int_normalize
       external int_normalize
c
       data g_coul,g_exch/-1,-1/
       data oexch,ocoul/.true.,.false./
       data ohalf/.true./
       data charge/0.d0/
c
c Print settings
c
       call util_print_push()
       call util_print_rtdb_load(rtdb,'mointdriver')
       oprint = util_print('information',print_default)
       omemprint = util_print('memory',print_high)
c
c Memory settings ... leave this to the user input
c
**       status = ma_set_hard_fail(.true.)
**       status = ma_set_auto_verify(.true.)
c
c
       status = geom_nuc_charge(geom, nuclear_charge)
       status = status.and.bas_numbf(basis, nbf)
       status = status.and.bas_numcont(basis,nsh)
       status = status.and.bas_nbf_cn_max(basis,maxbfsh)
       if (.not. status)
     $   call errquit('moint_driver: failed to get info', 0,
     &       INPUT_ERR)
c
c Get algorithm choice
c
       if (.not. rtdb_cget(rtdb, 'fourindex:method ', 1, algo))
     $      algo = 'twofold'
      if (.not. rtdb_get(rtdb, 'mp2:aotol2e', mt_dbl, 1, tol2e))
     $     tol2e = 1d-9  ! For screening integrals
c
c Check disk caching enabled
c
       if (.not.rtdb_get(rtdb,'fourindex:aodisk', MT_LOG,
     $                   1, odisk)) odisk = .false.
       if (.not.rtdb_get(rtdb,'fourindex:aoreuse', MT_LOG,
     $                   1, oaoreuse)) oaoreuse = .false.
       if (.not.rtdb_get(rtdb,'fourindex:aokeep', MT_LOG,
     $                   1, oaokeep)) oaokeep = .false.
c
c Reconcile number of electrons and occupied/virtual orbitals
c
       nelec = int(nuclear_charge - charge)
       status = MA_push_get(MT_DBL, nbf, 'MO occ', l_occ, k_occ)
       call dfill(nbf, 0.0d0, dbl_mb(k_occ), 1)
       if (.not. rtdb_get(rtdb, 'scfoccupation', MT_DBL, nbf,
     $      dbl_mb(k_occ))) then
         nocc = nelec/2
         nopen = mod(nelec,2)
         nact = nopen
         do i = 1, nocc
            dbl_mb(i-1+k_occ) = 2.0d0
         enddo
         do i = 1, nact
           dbl_mb(i-1+nocc+k_occ) = 1.d0
         enddo
       else
         nocc = 0
         nact = 0
         do i=1,nbf
           if (dbl_mb(k_occ+i-1).eq.2.d0) then
             nocc = nocc + 1
           else if (dbl_mb(k_occ+i-1).ne.0.d0) then
             nact = nact + 1
           endif
         enddo
       endif
C       if (nact.le.0) return
c
c Find MO vector filename in the database.
c
       if (.not.rtdb_cget(rtdb,'moints:input vectors',1,movecs_in))
     $      then
         if (.not.rtdb_cget(rtdb,'scf:output vectors',1,movecs_in))
     $        call util_file_name('movecs',.false.,.false.,movecs_in)
       endif
c
c Read in MO coefficients
c
       g_vecs = ga_create_atom_blocked(geom,basis,'Molecular orbital')
       if (.not. movecs_read_header(movecs_in, title_vecs,
     $    basis_vecs, scftype_vecs, nbf_vecs, nsets, nmo, 2))
     $    call errquit('mointdriver: cannot read MO vectors header',0,
     &       DISK_ERR)
       if (nbf.ne.nbf_vecs)
     $    call errquit(
     $    'mointdriver: number of basis fun does not match input',0,
     &       INPUT_ERR)
       status = MA_push_get(MT_DBL, nbf, 'Evals', l_evals, k_evals)
       if (.not. movecs_read(movecs_in, 1, dbl_mb(k_occ),
     $                       dbl_mb(k_evals), g_vecs)) call
     $      errquit('mointdriver: cannot read MO vectors',0, DISK_ERR)
c
c Get frozen core & virtual offsets
c
c
c Get frozen core and virtual offset
c
       if(.not. geom_freeze(rtdb, geom, 'mp2', occ_frozen)) then
          call errquit('direct_mp2: unable to get frozen core',0,
     &       GEOM_ERR)
       endif

       if (.not.rtdb_get(rtdb, 'mp2:frozen virtual', MT_INT, 1,
     $                   vir_frozen)) vir_frozen = 0

c$$$       if (.not.rtdb_get(rtdb,'fourindex:occ_frozen',MT_INT,1,
c$$$     $     occ_frozen)) occ_frozen = 0
c$$$       if (.not.rtdb_get(rtdb,'fourindex:vir_frozen',MT_INT,1,
c$$$     $     vir_frozen)) vir_frozen = 0
c
c Set index ranges
c
       occ_lo = occ_frozen + 1
       occ_hi = nocc
       vir_lo = nocc + 1
       vir_hi = nmo - vir_frozen 
       qlo = occ_lo
       qhi = vir_hi
       nmo2 = vir_hi - vir_lo + 1
       nocc_correl = occ_hi - occ_lo + 1
c
c Get Block-length
c
       if (.not.rtdb_get(rtdb, 'fourindex:block', MT_INT,
     $      1, blen)) blen = 16
       blen = max(blen,maxbfsh)
       if (.not.rtdb_get(rtdb,'fourindex:aoblock', MT_LOG,
     $      1, oblk)) oblk = .true.
c
c Initialise integrals and screening
c
       if (.not.int_normalize(rtdb,basis))
     &     call errquit('moint_driver: int_normalize failed',911,
     &       INT_ERR)
       call int_init(rtdb,1, basis)
       call schwarz_init(geom,basis)
c
c Maximum number of operator matrices we can allocate
c
       nocct = (nocc_correl*(nocc_correl+1))/2
       if ((ocoul).and.(oexch)) nocct = 2*nocct
       max_noper = moints_numoper( basis, algo, nbf,
     $                             nocc, nmo2, blen, omemprint )
C       MAX_NOPER = 10
       if (max_noper.gt.nocct) then
         noper =  (nocc_correl*(nocc_correl+1))/2
       else
         if ((ocoul).and.(oexch)) then
           noper = max_noper/2
         else
           noper = max_noper
         endif
       endif
       noper_core = noper
       if ((ocoul).and.(oexch)) noper_core = noper*2
c
c Allocate operator matrices
c
       if (ocoul)
     $      g_coul = ga_create_JKblocked(noper,nbf,nbf,'Coulomb oper')
       if (oexch) then
         if ((algo.eq.'twofold').and.(.not.(ohalf))) then
           g_exch = ga_create_JKblocked(noper,nbf,nmo2,'X oper')
         else
           g_exch = ga_create_JKblocked(noper,nbf,nbf,'X oper')
         endif
       endif
c
c Allocate pair energy array
c
*ga:1:0
       if (.not.ga_create(MT_DBL,1,((nocc_correl*(nocc_correl+1))/2),
     $      'epair',0,0,g_epair))
     $      call errquit('mointdriver: cannot create epair array',0,
     &       GA_ERR)
c
c Symmetry toggle
c
       osym = (sym_number_ops(geom) .gt. 0) 
c
c Print info
c
       if ((ga_nodeid().eq.0).and.(oprint)) then
         write(6,*)
         write(6,*)
         call util_print_centered(6,'Four-Index Transformation',
     $                            30,.true.)
         write(6,966) nbf,nsh,nocc,nocc_correl,blen,occ_lo,occ_hi,
     $                vir_lo,vir_hi,movecs_in, noper_core
 966     format(10x,'Number of basis functions:',10x,i5,/,
     $        10x,'Number of shells:',19x,i5,/,
     $        10x,'Number of occupied orbitals:',8x,i5,/,
     $        10x,'Number of occ. correlated orbitals:',1x,i5,/,
     $        10x,'Block length:',23x,i5,/,
     $        10x,'Superscript MO index range:',6x,i3,1x,'-',i3,/,
     $        10x,'Subscript MO index range:',8x,i3,1x,'-',i3,/,
     $        10x,'MO coefficients read from:',2x,a25,/,
     $        10x,'Number of operator matrices in core:',i5)
         if (algo.eq.'sixfold') then
           write(6,975) 
 975       format(10x,'Sixfold algorithm used')
         else
           write(6,976) 
 976       format(10x,'Twofold algorithm used')
         endif
         if (ohalf) then
           write(6,957) 
 957       format(10x,'Half-transformed integrals produced')
         else
           write(6,958)
 958       format(10x,'Half-transformed integrals not produced')
         endif
         if (oblk) then
           write(6,961)
 961       format(10x,'Use blocked AO integrals')
         else
           write(6,962)
 962       format(10x,'Use single quartet AO integrals')
         endif
         if (odisk) then
           write(6,963)
 963       format(10x,'AO integrals cached to disk')
         endif
           
         write(6,*)
         call util_flush(6)
       endif
*skip:       goto 333    ! skip direct MP2 section
*skip:c
*skip:c  Initialize aodisk
*skip:c
*skip:       if (.not.moints_aodisk_init( odisk, oaoreuse ))
*skip:     $    call errquit('direct_mp2: cannot initiliaze ao cache',0)
*skip:c
*skip:c   Clear 4-index statistics
*skip:c
*skip:       call moints_stats_clear()
*skip:
*skip:c
*skip:c   Loop over segments of occupied orbitals
*skip:c
*skip:       npass = 0
*skip:       oseg_hi = occ_lo - 1
*skip: 22    continue
*skip:       tpass = tcgtime()
*skip:       npass = npass + 1
*skip:       oseg_lo = oseg_hi + 1
*skip:       oseg_hi = moints_occrange( noper, occ_lo, oseg_lo, occ_hi )
*skip:c
*skip:c  J & K build routines
*skip:c
*skip:       if (algo.eq.'sixfold') then
*skip:         call moints_build_6x( basis, osym, occ_lo,
*skip:     $                         oseg_lo, oseg_hi,
*skip:     $                         vir_lo, vir_hi, g_vecs,
*skip:     $                         g_coul, ocoul, g_exch, oexch,
*skip:     $                         blen, oblk )
*skip:       else
*skip:         call moints_build_2x( basis,  ohalf, osym, occ_lo,
*skip:     $                         oseg_lo, oseg_hi,
*skip:     $                         vir_lo, vir_hi, g_vecs,
*skip:     $                         g_coul, ocoul, g_exch, oexch,
*skip:     $                         blen, oblk )
*skip:         call moints_build_2x( basis,  ohalf, osym, occ_lo,
*skip:     $                         oseg_lo, oseg_hi,
*skip:     $                         vir_lo, vir_hi, g_vecs,
*skip:     $                         g_coul, ocoul, g_exch, oexch,
*skip:     $                         blen, oblk )
*skip:       endif
*skip:c
*skip:c Evaluate and store pair energies
*skip:c
*skip:c$$$       call moints_e2pair( nocc, occ_lo, oseg_lo, oseg_hi, nmo2,
*skip:c$$$     $                     nmo2, nbf, g_exch, dbl_mb(k_evals),
*skip:c$$$     $                     g_epair )
*skip:       call moints_e2pair( nocc, occ_lo, oseg_lo, oseg_hi, 
*skip:     $                     vir_lo, vir_hi, vir_lo, vir_hi,
*skip:     $                     g_exch, dbl_mb(k_evals), g_epair )
*skip:c
*skip:c
*skip:c
*skip:       tpass = tcgtime() - tpass
*skip:       if ((ga_nodeid().eq.0).and.
*skip:     $     (util_print('progress reports',print_default)))
*skip:     $     write(6,772) npass, oseg_lo, oseg_hi, tpass
*skip: 772   format(5x,'Pass:',i5,5x,'Index range:',i5,2x,'-',i5,
*skip:     $        5x,'Time:',f10.2)
*skip:       if (oseg_hi.lt.occ_hi) goto 22
*skip:c
*skip:c End multipass loop
*skip:c
*skip:       call moints_aodisk_close( oaokeep )
*skip:       if (util_print('statistics',print_high))
*skip:     $   call moints_stats_print(algo)
*skip:c
*skip:c Evaluate mp2 correction from pair energies
*skip:c
*skip:       if ((oexch).and.(ga_nodeid().eq.0).
*skip:     $    and.(util_print('mp2',print_default))) then
*skip:         npairs = (nocc_correl*(nocc_correl+1))/2
*skip:         status = ma_push_get(MT_DBL, npairs, 'Pair energies',
*skip:     $                        l_epair, k_epair)
*skip:         ii = 0
*skip:         emp2 = 0.0d0
*skip:         do i=occ_lo,nocc
*skip:           do j=occ_lo,i
*skip:             ij = ((i-occ_frozen)*(i-occ_frozen-1))/2 + j - occ_frozen
*skip:             call ga_get(g_epair,1,1,ij,ij,e2pair,1)
*skip:             dbl_mb(k_epair+ii) = e2pair
*skip:             emp2 = emp2 + e2pair
*skip:             ii = ii + 1
*skip:           enddo
*skip:         enddo
*skip:c
*skip:c Print out pair energies
*skip:c
*skip:         if (util_print('pair energies',print_debug)) then
*skip:           write(6,997)
*skip: 997       format(//,15x,'Pair Energies',/)
*skip:           write(6,998) (dbl_mb(k_epair+ii),ii=0,npairs-1)
*skip: 998       format(5x,f12.6)
*skip:         endif
*skip:c
*skip:c Print out total correlation correction
*skip:c
*skip:         write(6,999) emp2
 999     format(//,15x,31('-'),
     $          /,15x,'MP2 Correction:',f16.10,
     $          /,15x,31('-'),/)
*skip:         status = ma_pop_stack(l_epair)
*skip:       endif
*skip:       call ga_sync()
*skip:c
*skip:c
*skip: 333   continue
       if (ocoul) status = ga_destroy(g_coul)
       if (oexch) status = ga_destroy(g_exch)
C       goto 444          !skip semi-direct code
c
c ---------------------------
c  Test semi-direct routine
c ---------------------------
c
       tmp2 = tcgtime()
       call ga_zero(g_epair)
       oseg_hi = occ_lo - 1
       nvir = vir_hi - vir_lo + 1
       ipass = 0
c
c  Memory arithmetic to determine occ segment length
c
       lmemreq = moints_lmem(basis, nocc, nvir, blen)
       if ((.not.ga_uses_ma()).and.(.not.(ga_memory_limited())))
     $   call errquit('cannot determine memory limit',0, GA_ERR)
       gmem = (ga_memory_avail() - lmemreq)*ga_nnodes()
       mxrlen=maxbfsh*maxbfsh*ga_nnodes()*10
       seg_incr = gmem/(nvir*mxrlen)
       seg_incr = min(max(seg_incr,1),nocc_correl)

       if (ga_nodeid().eq.0) then
         write(6,738) lmemreq, gmem, mxrlen, seg_incr
 738     format(' Local memory required:',2x,i20,/,
     $          ' Global memory remaining:',i20,/,
     $          ' Record length:',10x,i20,/,
     $        ' Occupied segment length:',i20)
       endif

       if (ga_nodeid().eq.0) then
         write(6,978) seg_incr
 978     format(//,10x,25('='),
     $           /,10x,'Test semi-direct MP2 code'
     $           /,10x,25('='),
     $        //, 'Segment length = ',I5)
         call util_flush(6)
       endif
 222   continue

       tpass = tcgtime()
       ipass = ipass + 1
       oseg_lo = oseg_hi + 1
       oseg_hi = min((oseg_lo + seg_incr - 1),occ_hi)
c
cedo  open a 128mb file enough?
       call moints_semi( basis, tol2e, oseg_lo, oseg_hi,
     $                  occ_lo, occ_hi, vir_lo, vir_hi, g_vecs, oblk ,
     c      128*1024*1024)

       status = ma_push_get(MT_DBL, (nbf*nmo), 'mos', l_x, k_x )
       call ga_get(g_vecs,1,nbf,1,nmo,dbl_mb(k_x),nbf)
       call moints_readintK( nbf, oseg_lo, oseg_hi, occ_lo, occ_hi,
     $                       vir_lo, vir_hi, dbl_mb(k_x),
     $                       dbl_mb(k_evals), g_epair )
       status = ma_pop_stack(l_x)
       tpass = tcgtime() - tpass
       if (ga_nodeid().eq.0) then
         write(6,332) ipass,oseg_lo,oseg_hi,tpass
 332     format('Pass:',i3,5x,'Segment:',i3,'-',i3,5x,f12.2)
         call util_flush(6)
       endif
       if (oseg_hi.lt.occ_hi) goto 222
c
c  Sum pair energies       
c
       emp2 = 0.d0
       do i=occ_lo,occ_hi
         do j=occ_lo,i
           ij = ((i-occ_lo+1)*(i-occ_lo))/2 + j - occ_lo + 1
           call ga_get(g_epair,1,1,ij,ij,e2pair,1)
           emp2 = emp2 + e2pair
         enddo
       enddo
       tmp2 = tcgtime() - tmp2
       if (ga_nodeid().eq.0) then
         write(6,999) emp2
         write(6,881) tmp2
 881     format(/,'Time for semi-direct energy:',f12.3)
         call util_flush(6)
       endif
c
c
c
 444   continue
       status = ga_destroy(g_epair)
       status = ma_pop_stack(l_evals)
       call schwarz_tidy()
       status = status.and.ma_pop_stack(l_occ)
       if (.not.status) call errquit('moints: cannot pop stack?',0,
     &       MA_ERR)
       status = ga_destroy(g_vecs)
c
c
c
       call int_terminate()
       call util_print_pop()
c
       return
       end











       integer function moints_occrange( noper, ostart, ilo, ohi )
       implicit none
#include "errquit.fh"
       integer noper, ostart, ilo, ohi
       integer trioff, i, ihi, ii
       integer ofroz

       ofroz = ostart - 1
       ihi = 0
       trioff = ((ilo-ofroz)*(ilo-ofroz-1))/2
       do i=ilo,ohi
         ii = ((i-ofroz)*(i-ofroz+1))/2
         if ((ii-trioff).le.noper) ihi = i
       enddo
       if (ihi.eq.0)
     $   call errquit('moints_occrange: not enough memory to do MP2',0,
     &       MEM_ERR)
       moints_occrange = ihi
       return
       end







