      subroutine mcscf_input(rtdb)
c
c $Id$
c
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      integer rtdb
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      integer num_profopt       ! No. of known profiling options
      parameter (num_dirs = 17)
      parameter (num_profopt = 1)
      character*12 dirs(num_dirs)
      character*255 test, keyword
      double precision conv, levelshift, tol2e
      logical status
      integer maxiter, nact, multip, nel, pflg, statesym
      data dirs / 'thresh',
     $            'maxiter',
     $            'level',
     $            'tol2e',
     $            'inactive',
     $            'active',
     $            'actelec',
     $            'multiplicity',
     $            'vectors',
     $            'profile',
     $            'symmetry',
     $            'hessian',
     $            'print',
     $            'noprint',
     $            'lock',
     $            'state',
     $            'end'/
c
c     Any input causes convergence flag to be deleted
c
      status = rtdb_delete(rtdb, 'mcscf:converged')
c
c
 10   if (.not. inp_read()) call errquit
     $     ('mcscf_input: inp_read failed', 0, INPUT_ERR)
      if (.not. inp_a(test))
     $     call errquit('mcscf_input: failed to read keyword', 0,
     &       INPUT_ERR)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     $     call errquit('mcscf_input: unknown directive', 0,
     &       INPUT_ERR)
c
c
c
      goto ( 100,  200, 300, 400, 500, 600, 700, 800, 900, 1000,
     $      1100, 1200, 1300, 1400, 1500, 1600,
     $      9999) ind
      call errquit('mcscf_input: unimplemented directive', ind,
     &       INPUT_ERR)
c     
c  thresh <convergence threshold>
c     
 100  if (.not. inp_f(conv)) call errquit
     $     ('mcscf_input: failed reading convergence', 0, INPUT_ERR)
      if (.not. rtdb_put(rtdb, 'mcscf:thresh', mt_dbl, 1, conv))
     $     call errquit('mcscf_input: rtdb_put failed', 0, RTDB_ERR)
      goto 10
c
c  maximum macro-iterations
c
 200  if (.not. inp_i(maxiter)) call errquit
     $     ('mcscf_input: failed reading maxiter', 0, INPUT_ERR)
      if (.not. rtdb_put(rtdb, 'mcscf:maxiter', mt_int, 1, maxiter))
     $     call errquit('mcscf_input: rtdb_put failed', 0,
     &       RTDB_ERR)
      goto 10
c
c  level shift 
c 
 300  if (.not. inp_f(levelshift)) call errquit
     $     ('mcscf_input: level <levelshift>',0, INPUT_ERR)
      if (.not. rtdb_put(rtdb,'mcscf:lshift', MT_DBL, 1, levelshift))
     $     call errquit('mcscf_input: rtdb_put of level failed',0,
     &       RTDB_ERR)
      goto 10
c
c  tol2e <real tol2e>
c
 400  if (.not. inp_f(tol2e)) call errquit
     $     ('mcscf_input: tol2e <tolerance>',0, INPUT_ERR)
      if (.not. rtdb_put(rtdb,'mcscf:tol2e', MT_DBL, 1, tol2e))
     $     call errquit('mcscf_input: rtdb_put of tol2e failed',0,
     &       RTDB_ERR)
      goto 10
c
c  inactive <integer nclosed>
c
 500  call errquit('mcscf_input: inactive directive not used',0,
     &       INPUT_ERR)
c
c  active <integer nact>
c
 600  if (.not. inp_i(nact)) call
     $  errquit('mcscf_input: active <active>',0, INPUT_ERR)
      if (.not. rtdb_put(rtdb,'mcscf:nact', MT_INT, 1, nact))
     $     call errquit('mcscf_input: rtdb_put of nact failed',0,
     &       RTDB_ERR)
      goto 10
c
c  active electrons <integer actel>
c
 700  if (.not. inp_i(nel)) call
     $  errquit('mcscf_input: actelec <activeelec>',0, INPUT_ERR)
      if (.not. rtdb_put(rtdb,'mcscf:nactelec', MT_INT, 1, nel))
     $     call errquit('mcscf_input: rtdb_put of nactel failed',0,
     &       RTDB_ERR)
      goto 10
c
c  spin multiplicity
c
 800  if (.not. inp_i(multip)) call
     $  errquit('mcscf_input: multiplicity <multi>',0, INPUT_ERR)
      if (.not. rtdb_put(rtdb,'mcscf:multiplicity', MT_INT, 1, multip))
     $     call errquit('mcscf_input: rtdb_put of multip failed',0,
     &       RTDB_ERR)
      goto 10
c
c  vectors
c
 900  call vectors_input(rtdb,'mcscf')
      goto 10
c
c  profiling options
c
 1000 if (.not. rtdb_put(rtdb,'mcscf:prof', MT_LOG, 1, .true.))
     $  call errquit('mcscf_input: rtdb_put of profile failed',0,
     &       RTDB_ERR)
 1001 continue
      if (.not.inp_a(keyword)) goto 10
      if (keyword.eq.'fock') then
        if (.not.rtdb_put(rtdb,'mcscf:fockprof',MT_LOG,1,.true.))
     $     call errquit('mcscf_input: rtdb_put of fockprof failed',0,
     &       RTDB_ERR)
      else if (keyword.eq.'ci') then
        if (.not.rtdb_put(rtdb,'mcscf:ciprof',MT_LOG,1,.true.))
     $     call errquit('mcscf_input: rtdb_put of ciprof failed',0,
     &       RTDB_ERR)
      else
        call errquit('mcscf_input: unrecognized profiling option',0,
     &       INPUT_ERR)
      endif
      goto 1001
c
c  state symmetry
c
 1100 if (.not.(inp_i(statesym)))
     $  call errquit('mcscf_input: symmetry <irrep>',0, INPUT_ERR)
      if (.not.rtdb_put(rtdb,'mcscf:statesymmetry',MT_INT,1,statesym))
     $  call errquit('mcscf_input: rtdb_put of statesym failed',0,
     &       RTDB_ERR)
      goto 1001
c
c  orbital Hessian exact/approx.
c
 1200 if (.not. inp_a(keyword))
     $  call errquit('mcscf_input: hessian (type>',0, INPUT_ERR)
      if (keyword.eq.'exact') then
        pflg = 2
      else if (keyword.eq.'onel') then
        pflg = 1
      else
        call errquit('mcscf_input: hessian: unknown type',0, INPUT_ERR)
      endif
      if (.not. rtdb_put(rtdb,'mcscf:hessiantype', MT_INT, 1, pflg))
     $     call errquit('mcscf_input: rtdb_put of pflg failed',0,
     &       RTDB_ERR)
      goto 1001
c
c     print/noprint
c
 1300 continue
 1400 continue
      call util_print_input(rtdb,'mcscf')
      goto 1001
c
c     lock
c
 1500 if (.not. rtdb_put(rtdb,'mcscf:movecs lock',mt_log, 1, .true.))
     $     call errquit('mcscf_input: failed putting lock',0, RTDB_ERR)
      goto 1001
c
c     state multiplictyxsym
c
 1600 if (.not. inp_a(keyword)) call errquit
     $     ('mcscf_input: failed reading the state',0, INPUT_ERR)
      if (.not. rtdb_cput(rtdb,'mcscf:state',1,keyword)) 
     $     call errquit('mcscf_input: failed storing the state',0,
     &       RTDB_ERR)
      status = rtdb_delete(rtdb,'mcscf:statesymmetry') ! Delete alternate form
      status = rtdb_delete(rtdb,'mcscf:multiplicity')
      goto 1001
c
c
 9999 continue
c     
      end
