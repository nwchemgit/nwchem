      subroutine ccsd_driver( rtdb, geom, basis, movecs_in,
     $     oconverged)
C     $Id$
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "schwarz.fh"
#include "inp.fh"
#include "msgids.fh"
#include "ccsd_len.fh"
#include "stdio.fh"
#include "util.fh"
#include "ccsdps.fh"
c
c
c
       integer rtdb                        ! [input] database handle
       integer geom, basis
       character*255 movecs_in
       logical oconverged
c
       character*255 title_vecs, basis_vecs
       integer nbf_vecs, nsets
       integer k_occ, l_occ, k_evals, l_evals, g_vecs
       integer nbf, nsh, nelec, nocc, nopen, nvir, nocc_correl,
     &         nvir_correl, lnof
       integer occ_frozen, occ_lo, occ_hi, nmo_lo
       integer vir_frozen, vir_lo, vir_hi, nmo_hi, nmo, nmo2
       integer noper, oseg_lo, oseg_hi, npass, blen, max_noper, nocct
       integer noper_core
       integer g_coul, g_exch, g_epair, g_ncoul, g_nexch
       logical oexch,ocoul,status,ohalf,do_trip
       double precision emp2, RefEner, TotEner
       double precision tpass
       character*7 algo
       Integer Me
       Integer ALo, AHi, IJLo, IJHi
       Integer MaxVec
C
       Logical OUseSym, OBlk, ODisk, occd
       Integer MaxBfSh
       Character*80 CC_Theory
       character*20 scftype_vecs
c
       integer maxit,iprt,ncor,ndel,l_mo,k_mo, l_scr,k_scr
       Double Precision Convi, Tol2E
c
       logical movecs_read_header, movecs_read
       external movecs_read_header, movecs_read
       integer ga_create_atom_blocked
       external ga_create_atom_blocked
       integer ga_create_JKblocked
       external ga_create_JKblocked
       logical ccsd_memory_fit
       external ccsd_memory_fit
       integer moints_occrange
       external moints_occrange
       logical int_normalize
       external int_normalize,ccsd_drivacc
       logical ccsd_drivacc,oprint, oprintdef
c
       double precision zip
       data zip/0.0d00/
c
       data g_coul,g_exch/-1,-1/
       data oexch,ocoul/.true.,.true./
       data ohalf/.true./
       Data OUseSym /.FALSE./
capr
       data algo/'twofold'/
capr
       Me = GA_NodeID()
       oprint = util_print('meminfo',print_debug)
       oprintdef = util_print('information',print_low)
c
       status = bas_numbf(basis, nbf)
       status = status.and.bas_numcont(basis,nsh)
       if (.not. status)
     $   call errquit('ccsd_driver: failed to get info', 0, BASIS_ERR)
c
c Reconcile number of electrons and occupied/virtual orbitals
c
       if (.not. rtdb_get(rtdb, 'scf:nclosed', MT_INT, 1, nocc))
     $      call errquit('direct_mp2: get of nocc failed',0, RTDB_ERR)
       if (.not. rtdb_get(rtdb, 'scf:nopen', MT_INT, 1, nopen))
     $      call errquit('ccsd: get of nopen failed',0, RTDB_ERR)
       if (nopen.ne.0) 
     $      call errquit('ccsd: nopen is not zero',nopen, INPUT_ERR)
       nelec = 2*nocc
       status = ma_push_get(MT_DBL, nbf, 'MO occ', l_occ, k_occ)
c
c Find MO vector filename in the database.
c
       if (.not.rtdb_cget(rtdb,'moints:input vectors',1,movecs_in))
     $      then
         if (.not.rtdb_cget(rtdb,'scf:output vectors',1,movecs_in))
     $        movecs_in = 'movecs'
       endif
c
c Read in MO coefficients
c
       g_vecs = ga_create_atom_blocked(geom,basis,'Molecular orbital')
       if (.not. movecs_read_header(movecs_in, title_vecs,
     $    basis_vecs, scftype_vecs, nbf_vecs, nsets, nmo, 2))
     $    call errquit('ccsddriver: cannot read MO vectors header',0,
     &       INPUT_ERR)
       if (nbf.ne.nbf_vecs)
     $    call errquit(
     $    'ccsddriver: number of basis fun does not match input',0,
     &       INPUT_ERR)
       status = MA_push_get(MT_DBL, nbf, 'Evals', l_evals, k_evals)
       if (.not. movecs_read(movecs_in, 1, dbl_mb(k_occ),
     $                       dbl_mb(k_evals), g_vecs)) call
     $      errquit('ccsddriver: cannot read MO vectors',0, INPUT_ERR)

       If ( Util_Print('reference', Print_High) )
     $    Call GA_Print(G_Vecs)
c
c Get frozen core & virtual offsets
c
      if(.not. geom_freeze(rtdb, geom, 'ccsd', occ_frozen)) then
        call errquit('ccsd: unable to freeze core orbitals', 0,
     &       INPUT_ERR)
      endif

      if (.not. rtdb_put(rtdb,'ccsd:number frozen core',MT_INT,1,
     $     occ_frozen))
     $   call errquit('ccsddriver: rtdb failure for occ_frozen',555,
     &       RTDB_ERR)
      if (.not.rtdb_get(rtdb, 'ccsd:frozen virtual', MT_INT, 1,
     $     vir_frozen)) then
         if (.not.rtdb_get(rtdb,'vir_frozen',MT_INT,1,vir_frozen))
     $        vir_frozen = 0
      endif
c
c Set index ranges
c
       occ_lo = occ_frozen + 1
       occ_hi = nocc
       nocc_correl = occ_hi - occ_lo + 1
       vir_lo = nocc + 1
       vir_hi = nmo - vir_frozen
       nvir_correl = vir_hi - vir_lo + 1
       nvir = nmo - nocc
       nmo_lo = occ_frozen + 1
       nmo_hi = nmo - vir_frozen
       nmo2 = nmo_hi - nmo_lo + 1
       noper = occ_hi - occ_lo + 1
c
c Get Block-length
c
       if (.not. bas_nbf_cn_max(basis,maxbfsh))
     $    call errquit('ccsd_driver: cannot get basis info',0,
     &       BASIS_ERR)
       if (.not.rtdb_get(rtdb,'fourindex:block',MT_INT,
     $    1,blen)) blen = 16
       blen = max(blen,maxbfsh)
*       if (.not.rtdb_get(rtdb,'fourindex:aoblock',MT_LOG,
*     $    1,oblk)) oblk = .false.
       oblk = .true.            ! non-blocked 4-index currently broken
c
c Disk caching disabled
c
       ODisk = .false.
c     
c Initialise integrals and screening
c
       if (.not.int_normalize(rtdb,basis))
     &     call errquit('ccsd_driver: int_normalize failed',911,
     &       UNKNOWN_ERR)
       call int_init(rtdb,1, basis)
       call schwarz_init(geom,basis)
c
c Maximum number of operator matrices we can allocate
c
       if(oprint) then
          call ga_summarize(.true.)
          call ma_summarize_allocated_blocks()
       endif
       nocct = (nocc_correl*(nocc_correl+1))/2
       if ((ocoul).and.(oexch)) nocct = 2*nocct
c
      lnoo=nocc_correl*nocc_correl
      lnov=nocc_correl*nvir_correl
      lnvv=nvir_correl*nvir_correl
      lnoov=lnoo*nvir_correl
      lnooo=lnoo*nocc_correl
      lnovv=lnvv*nocc_correl
c
      CC_Theory = " "
      call ccsd_get_param(rtdb,maxit,convi,iprt, CC_Theory, MaxVec,
     $     Tol2E, occd)
      do_trip =  Inp_Compare(.FALSE., 'ccsd(t)', CC_Theory)
     $   .OR. Inp_Compare(.FALSE., 'ccsd+t(ccsd)', CC_Theory)
c
      if (.not. ccsd_memory_fit(basis,algo,nmo,nocc,nmo2,blen,maxvec,
     $   do_trip)) call errquit('ccsd_driver:memory error',0, GA_ERR)
c
      noper =  (nocc_correl*(nocc_correl+1))/2
      noper_core = noper
      if ((ocoul).and.(oexch)) noper_core = noper*2
c
      lnof=nmo2*nocc
      if (.not.ga_create(MT_DBL,lnof,lnof,'ncoul',
     &                   nmo2,nmo2,g_ncoul))
     &     call errquit('ga_create g_ncoul failed',0, GA_ERR)
      if (.not.ga_create(MT_DBL,lnof,lnof,'nexch',
     &                   nmo2,nmo2,g_nexch))
     &     call errquit('ga_create g_nexch failed',0, GA_ERR)
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
     $      call errquit('cannot create epair array',0, GA_ERR)
c
c Print info
c
       if (Me.eq.0.and.oprintdef) then
         write(6,*)
         write(6,*)
         call util_print_centered(6,'Four-Index Transformation',
     $                            30,.true.)
         write(6,966) nbf,nsh,nocc,nocc_correl,blen,occ_lo,occ_hi,
     $                nmo_lo,nmo_hi,movecs_in, noper_core
 966     format(10x,'Number of basis functions:',10x,i5,/,
     $        10x,'Number of shells:',19x,i5,/,
     $        10x,'Number of occupied orbitals:',8x,i5,/,
     $        10x,'Number of occ. correlated orbitals:',1x,i5,/,
     $        10x,'Block length:',23x,i5,/,
     $        10x,'Superscript MO index range:',2x,i5,1x,'-',i5,/,
     $        10x,'Subscript MO index range:',4x,i5,1x,'-',i5,/,
     $        10x,'MO coefficients read from:',2x,a25,/,
     $        10x,'Number of operator matrices in core:',i5)
         call util_flush(6)
         if (ohalf) then
           write(6,957) 
 957       format(10x,'Half-transformed integrals produced')
         else
           write(6,958)
 958       format(10x,'Half-transformed integrals not produced')
         endif
       endif
c
c Clear 4-index statistics
c
       call moints_stats_clear()
c
c Loop over segments of occupied orbitals
c
       if (Me.eq.0.and.oprintdef) write(6,*)
       npass = 0
       oseg_hi = occ_lo - 1
 22    continue
       tpass = tcgtime()
       npass = npass + 1
       oseg_lo = oseg_hi + 1
       oseg_hi = moints_occrange( noper, occ_lo, oseg_lo, occ_hi )
c
c  J & K build routines
c
      if (occsdps) call pstat_on(ps_aomo)
       if (algo.eq.'sixfold') then
         call moints_build_6x(basis,  OUseSym, occ_lo, oseg_lo, oseg_hi,
     $                        nmo_lo, nmo_hi, g_vecs,
     $                        g_coul, ocoul, g_exch, oexch,
     $                        blen, OBlk )
       else
         call moints_build_2x(basis,  ohalf, OUseSym,
     $       occ_lo, oseg_lo, oseg_hi, nmo_lo, nmo_hi, g_vecs,
     $                        g_coul, ocoul, g_exch, oexch,
     $                        blen, OBlk )
       endif
      if (occsdps) call pstat_off(ps_aomo)
c
c
c Evaluate and store pair energies
c
capr insert oexch
       if (oexch)then
        call moints_e2pair(rtdb, nocc, occ_lo, oseg_lo, oseg_hi, 
     $                      nmo_lo, nmo_hi, vir_lo, vir_hi,
     $                      g_exch, dbl_mb(k_evals), g_epair )
       endif

c
c
c
       tpass = tcgtime() - tpass
       if (Me.eq.0.and.oprintdef) write(6,772) 
     $    npass, oseg_lo, oseg_hi, tpass
 772   format(5x,'Pass:',i5,5x,'Index range:',i5,2x,'-',i5,
     $        5x,'Time:',f10.2)
       if (oseg_hi.lt.occ_hi) goto 22
c
c End multipass loop
c
C       call moints_stats_print(algo)
c
c      Evaluate mp2 correction from pair energies
C
c      Eventually this should be replaced with a derivative of
C      rimp2/sum_pairs (currently understands only rectangular
C      g_epair arrays, not triangular).
C
C      Note that these are not strictly pair energies -- the diagonals
C      are off by a factor of two.
C
       call ga_sync
       If ( OExch ) then
          If ( Util_Print('guess pair energies', Print_Debug) ) then
             If ( Me .eq. 0 ) Call Banner( LuOut,
     $          'MP2 Pair Correlation Energies', '-',
     $          .TRUE., .TRUE., .FALSE.)
             Call GA_Sync
             Call GA_Print( G_EPair )
          EndIf
C
          EMP2 = zip
          Call GA_Distribution(G_EPair, GA_NodeID(),
     $       Alo, AHi, IJLo, IJHi)
*          If ( ALo .ne. 1 .OR. AHi .ne. 1) Call ErrQuit(
*     $       'CCSD: Internal error: G_EPair has wrong shape', AHi-ALo+1)
C
C         GA_Access gets upset if there is nothing local
C
          If ( IJHi - IJLo + 1 .gt. 0 ) then
             call ga_access_callback_release(
     $        g_EPair, Alo, Ahi, Ijlo, Ijhi,ccsd_drivacc, EMP2,
     $        0, 0, 0, 0, 0,0) 
          EndIf
          Call GA_DGOp(msg_cc_emp2, EMP2, 1, '+')
          Call GA_Sync
       EndIF
       call ga_sync
C
C     Also get the SCF energy so we can report the total energy
C
      If ( .NOT. RTDB_Get(RTDB, 'scf:energy', MT_Dbl, 1, RefEner)
     $   ) then
         RefEner = Zip
         If ( Me .eq. 0) Write ( LuOut, *)
     $      'CCSD: Warning! Unable to get scf:energy from RTDB.'
      EndIf
C
      TotEner = RefEner + EMP2
c
c     Print out total correlation correction
c
      If ( Util_Print('byproduct energies', Print_Default)
     $   .AND. Me .eq. 0) then
         Call Banner(LuOut,
     $      'MP2 Energy (coupled cluster initial guess)', '-',
     $      .TRUE., .TRUE., .FALSE.)
         Write (LuOut, '(1X, A, F30.15)') 'Reference energy: ',
     $      RefEner
         Write (LuOut, '(1X, A, F30.15)') 'MP2 Corr. energy: ',
     $      EMP2
         Write (LuOut, '(1X, A, F30.15)') 'Total MP2 energy: ',
     $      TotEner
      EndIf
C
C     Write our energy out to the database
C
      If ( .NOT. RTDB_Put( RTDB, 'ccsd:mp2 correlation energy',
     $   MT_Dbl, 1, EMP2) )
     $   Call ErrQuit('CCSD: Can''t write to RTDB', 1,0)
      If ( .NOT. RTDB_Put( RTDB, 'ccsd:mp2 energy', MT_Dbl, 1,
     $   TotEner) ) 
     $   Call ErrQuit('CCSD: Can''t write to RTDB', 2,0)
c
      call ecce_print1('mp2 energy', mt_dbl, totener, 1)
      call ga_sync()
c
      status = ga_destroy(g_epair)
capr
c
c Start CCSD part here
c
c
      ncor=occ_frozen
      ndel=vir_frozen
c
      lnobj=nmo2*nmo2+nmo2*nmo2
      status=status.and.ma_push_get(MT_DBL,lnobj,'scr',l_scr,k_scr)
      if (.not.status)call errquit('error before ccsd_sortmo',0, MA_ERR)
      call ccsd_sortmo(ncor,nocc_correl,nvir_correl,nmo2,nbf,
     $                 dbl_mb(k_scr),dbl_mb(k_scr+nmo2*nmo2),g_coul,
     $                 g_exch,g_ncoul,g_nexch,iprt)
          status=status.and.ga_destroy(g_coul)
          status=status.and.ga_destroy(g_exch)
      status=status.and.ma_pop_stack(l_scr)
      if (.not.status) call errquit('error after ccsd_sortmo',0, GA_ERR)
      call rddcom()
      status=ma_push_get(MT_DBL,(nbf*nbf),'movecs',l_mo,k_mo)
      if (.not.status)call errquit('movecs: failed to alloc', l_mo,
     &       MA_ERR)
      call ga_get(g_vecs,1,nbf,1,nbf,dbl_mb(k_mo),nbf)
      call aoccsd(basis,ncor,nocc_correl,nvir_correl,ndel,nmo2,nbf,
     &            maxit,convi,iprt,
     &            dbl_mb(k_mo),dbl_mb(k_evals),
     &            blen,g_ncoul,g_nexch, RefEner,
     $   CC_Theory, RTDB, MaxVec, Geom, Tol2E, occd, oconverged)
      status = ma_pop_stack(l_mo)
      if (.not.ga_destroy(g_nexch))
     &    call errquit('ga_dest g_nexch fail',0, GA_ERR)
      if (.not.ga_destroy(g_ncoul))
     &    call errquit('ga_dest g_ncoul fail',0, GA_ERR)
      if (.not.status)call errquit('deallocate after aoccsd',0, MA_ERR)
c
      if (.not.status) call errquit('error after aoccsd',0, MA_ERR)
      status = ma_pop_stack(l_evals)
c
c Cleanup
c
      call schwarz_tidy()
      status = status.and.ma_pop_stack(l_occ)
      if (.not.status) call errquit('moints: cannot pop stack?',0,
     &       MA_ERR)
      status = ga_destroy(g_vecs)
      call int_terminate()
      call ga_sync()
c
      return
      end
      integer function ccsd_drivacc(
     $        g_EPair, Alo, Ahi, Ijlo, Ijhi, a,ld, EMP2, 
     $     idum2, idum3, idum4, idum5, idum6,idum7)
      implicit none                                                             
      integer g_EPair, Alo, Ahi, Ijlo, Ijhi, ld,
     $     idum2, idum3, idum4, idum5, idum6,idum7
      double precision a(Alo:Alo+ld-1,ijlo:ijhi), EMP2
      integer j 
      do j=ijlo,ijhi
         emp2 = emp2 + a(alo,J)
      enddo
      ccsd_drivacc = 0    ! no update necessary
c
      end     
