      subroutine ccsd_input(rtdb)
      implicit none
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      integer rtdb
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 15)
      character*8 dirs(num_dirs)
      character*255 test
      integer maxit, convi, iprt, MaxDIIS
      Integer ListLen
      Parameter (ListLen = 15) ! Length of the arrays in the /debug/ common
      Integer IList(ListLen), NGiven
      data dirs /'maxiter','thresh','        ','iprt',
     $   'print','noprint',
     $   'doa', 'dob', 'dog', 'doh', 'dojk', 'dos', 'dod',
     $   'dissbas', 'end'/
c
 10   if (.not. inp_read()) call errquit
     $     ('ccsd_input: inp_read failed', 0)
      if (.not. inp_a(test))
     $     call errquit('ccsd_input: failed to read keyword', 0)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     $     call errquit('ccsd_input: unknown directive', 0)
c     
      goto (100, 200, 300, 400, 500, 600,
     $   700, 705, 710, 715, 720, 725, 730, 800, 9999) ind
      call errquit('ccsd_input: unimplemented directive', ind)
c     
c     Maximum iterations
c
 100  if (.not. inp_i(maxit)) call errquit
     $     ('ccsd_input: failed reading maxiter', 0)
      if (.not. rtdb_put(rtdb, 'ccsd:maxiter', mt_int, 1, maxit))
     $     call errquit('ccsd-input: rtdb_put failed', 0)
      goto 10
c     
c     thresh <convergence threshold>
c     
 200  if (.not. inp_i(convi)) call errquit
     $     ('ccsd_input: failed reading convergence', 0)
      if (.not. rtdb_put(rtdb, 'ccsd:thresh', mt_int, 1, convi))
     $     call errquit('ccsd-input: rtdb_put failed', 0)
      goto 10
c
c     room for more
c
 300  goto 10
c
c     iprt is a carry over from the old titan stuff.
c
 400  if (.not. inp_i(iprt)) call errquit
     $     ('ccsd_input: failed reading convergence', 0)
      if (.not. rtdb_put(rtdb, 'ccsd:iprt', mt_int, 1, iprt))
     $     call errquit('ccsd-input: rtdb_put failed', 0)
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
      if ( .not. Inp_IList(ListLen, IList, Ngiven) ) Call ErrQuit(
     $   'ccsd: bad format for debugging term switch', 0)
      If ( NGiven .ne. ListLen ) Call ErrQuit(
     $   'ccsd: number of elements is fixed', ListLen)
      If ( .NOT. RTDB_Put( RTDB, 'ccsd:doa', MT_Int, NGiven,
     $   IList)) Call ErrQuit('ccsd: rtdb_put failed', 0)
      Goto 10
 705  Continue                  ! dob
      if ( .not. Inp_IList(ListLen, IList, Ngiven) ) Call ErrQuit(
     $   'ccsd: bad format for debugging term switch', 0)
      If ( NGiven .ne. ListLen ) Call ErrQuit(
     $   'ccsd: number of elements is fixed', ListLen)
      If ( .NOT. RTDB_Put( RTDB, 'ccsd:dob', MT_Int, NGiven,
     $   IList) ) Call ErrQuit('ccsd: rtdb_put failed', 0)
      Goto 10
 710  Continue                  ! dog
      if ( .not. Inp_IList(ListLen, IList, Ngiven) ) Call ErrQuit(
     $   'ccsd: bad format for debugging term switch', 0)
      If ( NGiven .ne. ListLen ) Call ErrQuit(
     $   'ccsd: number of elements is fixed', ListLen)
      If ( .NOT. RTDB_Put( RTDB, 'ccsd:dog', MT_Int, NGiven,
     $   IList) ) Call ErrQuit('ccsd: rtdb_put failed', 0)
      Goto 10
 715  Continue                  ! doh
      if ( .not. Inp_IList(ListLen, IList, Ngiven) ) Call ErrQuit(
     $   'ccsd: bad format for debugging term switch', 0)
      If ( NGiven .ne. ListLen ) Call ErrQuit(
     $   'ccsd: number of elements is fixed', ListLen)
      If ( .NOT. RTDB_Put( RTDB, 'ccsd:doh', MT_Int, NGiven,
     $   IList) ) Call ErrQuit('ccsd: rtdb_put failed', 0)
      Goto 10
 720  Continue                  ! dojk
      if ( .not. Inp_IList(ListLen, IList, Ngiven) ) Call ErrQuit(
     $   'ccsd: bad format for debugging term switch', 0)
      If ( NGiven .ne. ListLen ) Call ErrQuit(
     $   'ccsd: number of elements is fixed', ListLen)
      If ( .NOT. RTDB_Put( RTDB, 'ccsd:dojk', MT_Int, NGiven,
     $   IList) ) Call ErrQuit('ccsd: rtdb_put failed', 0)
      Goto 10
 725  Continue                  ! dos
      if ( .not. Inp_IList(ListLen, IList, Ngiven) ) Call ErrQuit(
     $   'ccsd: bad format for debugging term switch', 0)
      If ( NGiven .ne. ListLen ) Call ErrQuit(
     $   'ccsd: number of elements is fixed', ListLen)
      If ( .NOT. RTDB_Put( RTDB, 'ccsd:dos', MT_Int, NGiven,
     $   IList) ) Call ErrQuit('ccsd: rtdb_put failed', 0)
      Goto 10
 730  Continue                  ! dod
      if ( .not. Inp_IList(ListLen, IList, Ngiven) ) Call ErrQuit(
     $   'ccsd: bad format for debugging term switch', 0)
      If ( NGiven .ne. ListLen ) Call ErrQuit(
     $   'ccsd: number of elements is fixed', ListLen)
      If ( .NOT. RTDB_Put( RTDB, 'ccsd:dod', MT_Int, NGiven,
     $   IList) ) Call ErrQuit('ccsd: rtdb_put failed', 0)
      Goto 10
C
C     Maximum size of DIIS basis
C
 800  if (.not. inp_i(maxdiis)) call errquit
     $     ('ccsd_input: failed reading maxdiis', 0)
      if (.not. rtdb_put(rtdb, 'ccsd:maxdiis', mt_int, 1, maxdiis))
     $     call errquit('ccsd-input: rtdb_put failed', 0)
      goto 10
c
c     end
c     
 9999 continue
c     
      end

