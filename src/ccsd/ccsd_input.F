      subroutine ccsd_input(rtdb)
C     $Id$
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
C
C     We need the LenDbg parameter only, not the arrays themselves
C     (the data goes into the RTDB!)
C
#include "ccsd_debug.fh"
c
      integer rtdb
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 19)
      character*8 dirs(num_dirs)
      character*255 test
      character*80 key
      integer maxit, iprt, MaxDIIS
      Double Precision convi, tol2e, fss, fos
C
C     Length of debugging arrays (LenDbg) is defined in ccsd_debug.fh
C
      Integer IList(LenDbg), NGiven
      data dirs /'maxiter','thresh','tol2e','iprt',
     $   'print','noprint',
     $   'doa', 'dob', 'dog', 'doh', 'dojk', 'dos', 'dod',
     $   'diisbas', 'freeze', 'nodisk', 'fss', 'fos', 'end'/
c
c     Remove previous entries so empty directive forces freeze nothing.
c     To be consistent with the rest of the code, a new ccsd block
c     should reset no the default of no freezing
c
      key = ' '
      write(key,'(a,'':'',a)') 'ccsd','frozen core'
      call rtdb_clean(rtdb,key)
      key = ' '
      write(key,'(a,'':'',a,'':'')') 'ccsd','frozen virtual'
      call rtdb_clean(rtdb,key)
c
 10   if (.not. inp_read()) call errquit
     $     ('ccsd_input: inp_read failed', 0, INPUT_ERR)
      if (.not. inp_a(test))
     $     call errquit('ccsd_input: failed to read keyword', 0,
     &       INPUT_ERR)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     $     call errquit('ccsd_input: unknown directive', 0,
     &       INPUT_ERR)
c     
      goto (100, 200, 300, 400, 500, 600, 700,
     $   800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600,
     $  1700, 1800, 9999) ind
      call errquit('ccsd_input: unimplemented directive', ind,
     &       INPUT_ERR)
c     
c     Maximum iterations
c
 100  if (.not. inp_i(maxit)) call errquit 
     $     ('ccsd_input: failed reading maxiter', 0, INPUT_ERR)
      if (.not. rtdb_put(rtdb, 'ccsd:maxiter', mt_int, 1, maxit))
     $     call errquit('ccsd-input: rtdb_put failed', 0, RTDB_ERR)
      goto 10
c     
c     thresh <convergence threshold>
c     
 200  if (.not. inp_f(convi)) call errquit
     $     ('ccsd_input: failed reading convergence', 0, INPUT_ERR)
      if (.not. rtdb_put(rtdb, 'ccsd:thresh', mt_dbl, 1, convi))
     $     call errquit('ccsd_input: rtdb_put failed', 0, RTDB_ERR)
      goto 10
c     
c     tol2e <two-electron integral tolerance>
c     
 300  if (.not. inp_f(tol2e)) call errquit
     $     ('ccsd_input: failed reading tolerance', 0, INPUT_ERR)
      if (.not. rtdb_put(rtdb, 'ccsd:tol2e', mt_dbl, 1, tol2e))
     $     call errquit('ccsd_input: rtdb_put failed', 0, RTDB_ERR)
      goto 10
c
c     iprt is a carry over from the old titan stuff.
c
 400  if (.not. inp_i(iprt)) call errquit
     $     ('ccsd_input: failed reading convergence', 0, INPUT_ERR)
      if (.not. rtdb_put(rtdb, 'ccsd:iprt', mt_int, 1, iprt))
     $     call errquit('ccsd-input: rtdb_put failed', 0, RTDB_ERR)
      goto 10
C
C     Use the NWChem printing control mechanisms
C
 500  Continue
 600  Continue
      Call Util_Print_Input(RTDB, 'ccsd')
      Goto 10
C
C     Debugging term switches
C
 700  Continue                  ! doa
      if ( .not. Inp_IList(LenDbg, IList, Ngiven) ) Call ErrQuit(
     $   'ccsd: bad format for debugging term switch', 0, INPUT_ERR)
      If ( NGiven .ne. LenDbg ) Call ErrQuit(
     $   'ccsd: number of elements is fixed', LenDbg, INPUT_ERR)
      If ( .NOT. RTDB_Put( RTDB, 'ccsd:doa', MT_Int, NGiven,
     $   IList)) Call ErrQuit('ccsd: rtdb_put failed', 0, RTDB_ERR)
      Goto 10
 800  Continue                  ! dob
      if ( .not. Inp_IList(LenDbg, IList, Ngiven) ) Call ErrQuit(
     $   'ccsd: bad format for debugging term switch', 0, INPUT_ERR)
      If ( NGiven .ne. LenDbg ) Call ErrQuit(
     $   'ccsd: number of elements is fixed', LenDbg, INPUT_ERR)
      If ( .NOT. RTDB_Put( RTDB, 'ccsd:dob', MT_Int, NGiven,
     $   IList) ) Call ErrQuit('ccsd: rtdb_put failed', 0, RTDB_ERR)
      Goto 10
 900  Continue                  ! dog
      if ( .not. Inp_IList(LenDbg, IList, Ngiven) ) Call ErrQuit(
     $   'ccsd: bad format for debugging term switch', 0, INPUT_ERR)
      If ( NGiven .ne. LenDbg ) Call ErrQuit(
     $   'ccsd: number of elements is fixed', LenDbg, INPUT_ERR)
      If ( .NOT. RTDB_Put( RTDB, 'ccsd:dog', MT_Int, NGiven,
     $   IList) ) Call ErrQuit('ccsd: rtdb_put failed', 0, RTDB_ERR)
      Goto 10
 1000  Continue                  ! doh
      if ( .not. Inp_IList(LenDbg, IList, Ngiven) ) Call ErrQuit(
     $   'ccsd: bad format for debugging term switch', 0, INPUT_ERR)
      If ( NGiven .ne. LenDbg ) Call ErrQuit(
     $   'ccsd: number of elements is fixed', LenDbg, INPUT_ERR)
      If ( .NOT. RTDB_Put( RTDB, 'ccsd:doh', MT_Int, NGiven,
     $   IList) ) Call ErrQuit('ccsd: rtdb_put failed', 0, RTDB_ERR)
      Goto 10
 1100  Continue                  ! dojk
      if ( .not. Inp_IList(LenDbg, IList, Ngiven) ) Call ErrQuit(
     $   'ccsd: bad format for debugging term switch', 0, INPUT_ERR)
      If ( NGiven .ne. LenDbg ) Call ErrQuit(
     $   'ccsd: number of elements is fixed', LenDbg, INPUT_ERR)
      If ( .NOT. RTDB_Put( RTDB, 'ccsd:dojk', MT_Int, NGiven,
     $   IList) ) Call ErrQuit('ccsd: rtdb_put failed', 0, RTDB_ERR)
      Goto 10
 1200  Continue                  ! dos
      if ( .not. Inp_IList(LenDbg, IList, Ngiven) ) Call ErrQuit(
     $   'ccsd: bad format for debugging term switch', 0, INPUT_ERR)
      If ( NGiven .ne. LenDbg ) Call ErrQuit(
     $   'ccsd: number of elements is fixed', LenDbg, INPUT_ERR)
      If ( .NOT. RTDB_Put( RTDB, 'ccsd:dos', MT_Int, NGiven,
     $   IList) ) Call ErrQuit('ccsd: rtdb_put failed', 0, RTDB_ERR)
      Goto 10
 1300  Continue                  ! dod
      if ( .not. Inp_IList(LenDbg, IList, Ngiven) ) Call ErrQuit(
     $   'ccsd: bad format for debugging term switch', 0, INPUT_ERR)
      If ( NGiven .ne. LenDbg ) Call ErrQuit(
     $   'ccsd: number of elements is fixed', LenDbg, INPUT_ERR)
      If ( .NOT. RTDB_Put( RTDB, 'ccsd:dod', MT_Int, NGiven,
     $   IList) ) Call ErrQuit('ccsd: rtdb_put failed', 0, RTDB_ERR)
      Goto 10
C
C     Maximum size of DIIS basis
C
 1400  if (.not. inp_i(maxdiis)) call errquit
     $     ('ccsd_input: failed reading maxdiis', 0, INPUT_ERR)
      if (.not. rtdb_put(rtdb, 'ccsd:maxdiis', mt_int, 1, maxdiis))
     $     call errquit('ccsd-input: rtdb_put failed', 0, RTDB_ERR)
      goto 10
c
c     Freeze
c
 1500 call freeze_input(rtdb,'ccsd')
      goto 10
C
C     nodisk means do not use scratch disk to store integrals
C
 1600 print*,'usedisk called input'
      if (.not. rtdb_put(rtdb, 'ccsd:usedisk', mt_log, 1, .false.))
     $     call errquit('ccsd-input: rtdb_put failed', 0, RTDB_ERR)
      goto 10
c     
c     fss <same spin scaling factor for SCS-CCSD>
c     
 1700 if (.not. inp_f(fss)) call errquit
     $     ('ccsd_input: failed reading fss', 0, INPUT_ERR)
      if (.not. rtdb_put(rtdb, 'ccsd:fss', mt_dbl, 1, fss))
     $     call errquit('ccsd_input: rtdb_put failed', 0, RTDB_ERR)
      goto 10
c     
c     fos <opposite spin scaling factor for SCS-CCSD>
c     
 1800 if (.not. inp_f(fos)) call errquit
     $     ('ccsd_input: failed reading fos', 0, INPUT_ERR)
      if (.not. rtdb_put(rtdb, 'ccsd:fos', mt_dbl, 1, fos))
     $     call errquit('ccsd_input: rtdb_put failed', 0, RTDB_ERR)
      goto 10
c
c     end
c     
 9999 continue
c     
      end

