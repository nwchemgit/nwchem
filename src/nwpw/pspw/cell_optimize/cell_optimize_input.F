      subroutine cell_optimize_input(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "bafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
c
      integer rtdb
      logical value
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 8)
 
      character*30 dirs(num_dirs)
      character*255 test

      data dirs / 'cell_name',
     >            'optimize',
     >            'optimize_lattice',
     >            'optimize_lattice_vectors',
     >            'lattice_tolerances',
     >            'lattice_vectors_tolerances',
     >            'cycles',
     >            'end'/
c
      integer ntokens
      parameter (ntokens=6)
      character*10 tokens(ntokens)
      data tokens /'lat_a', 'lat_b', 'lat_c',
     >     'alpha', 'beta', 'gamma' /


      integer ntokens2
      parameter (ntokens2=9)
      character*10 tokens2(ntokens2)
      data tokens2 /'a11', 'a21', 'a31',
     >              'a12', 'a22', 'a32',
     >              'a13', 'a23', 'a33' /

      character*50     cell_name
      integer          optimize_strategy,cycles
      integer          optimize_lattice(6)
      integer          optimize_lattice_vectors(9)
      double precision lattice_tolerances(6),fe
       


*     ***** initializations ****
      cell_name             = 'cell_default'

      optimize_strategy = 0  ! 0-complete lattice, 1-specified lattice, 2-specified lattice_vectors
      cycles            = 10
      optimize_lattice(1) = 1
      optimize_lattice(2) = 1
      optimize_lattice(3) = 1
      optimize_lattice(4) = 1
      optimize_lattice(5) = 1
      optimize_lattice(6) = 1
      optimize_lattice_vectors(1) = 1
      optimize_lattice_vectors(2) = 1
      optimize_lattice_vectors(3) = 1
      optimize_lattice_vectors(4) = 1
      optimize_lattice_vectors(5) = 1
      optimize_lattice_vectors(6) = 1
      optimize_lattice_vectors(7) = 1
      optimize_lattice_vectors(8) = 1
      optimize_lattice_vectors(9) = 1
   

      lattice_tolerances(1) = 1.0d-3
      lattice_tolerances(2) = 1.0d-3
      lattice_tolerances(3) = 1.0d-3
      lattice_tolerances(4) = 1.0d-3
      lattice_tolerances(5) = 1.0d-3
      lattice_tolerances(6) = 1.0d-3


 10   if (.not. inp_read()) 
     >     call errquit(
     >           'cell_optimize_input: inp_read failed', 0, INPUT_ERR)
      if (.not. inp_a(test))
     >     call errquit(
     >           'cell_optimize_input: failed to read keyword', 0,
     &       INPUT_ERR)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit(
     >           'cell_optimize_input: unknown directive', 0, INPUT_ERR)


      goto ( 100, 200, 300, 400, 500, 600,700,
     >      9999) ind
      call errquit(
     >      'cell_optimize_input: unimplemented directive', ind,
     &       INPUT_ERR)


c     
c  cell_name
c     
 100  if (.not. inp_a(cell_name))
     >     call errquit(
     >          'cell_optimize_input: failed to read cell_name', 0,
     &       INPUT_ERR)
      goto 10


c
c  optimize
c
 200  if (.not. inp_a(test))
     >     call errquit(
     >          'cell_optimize_input: failed to read optimize:', 0,
     &       INPUT_ERR)
      if (test.eq.'all')             optimize_strategy = 0
      if (test.eq.'lattice')         optimize_strategy = 1
      if (test.eq.'lattice_vectors') optimize_strategy = 2
      goto 10

c
c  optimize_lattice
c
 300  optimize_lattice(1) = 0
      optimize_lattice(2) = 0
      optimize_lattice(3) = 0
      optimize_lattice(4) = 0
      optimize_lattice(5) = 0
      optimize_lattice(6) = 0
 301  if (.not. inp_a(test)) goto 320

      if (.not. inp_match(ntokens, .false., test, tokens, ind))
     >     call errquit(
     >           'optimize_lattice: unknown directive', 0, INPUT_ERR)

      goto (311,312,313,314,315,316) ind
c
 311      optimize_lattice(1) = 1
          goto 301
 312      optimize_lattice(2) = 1
          goto 301
 313      optimize_lattice(3) = 1
          goto 301
 314      optimize_lattice(4) = 1
          goto 301
 315      optimize_lattice(5) = 1
          goto 301
 316      optimize_lattice(6) = 1
          goto 301

 320      goto 10

c
c  optimize_lattice_vectors
c
 400  optimize_lattice_vectors(1) = 0
      optimize_lattice_vectors(2) = 0
      optimize_lattice_vectors(3) = 0
      optimize_lattice_vectors(4) = 0
      optimize_lattice_vectors(5) = 0
      optimize_lattice_vectors(6) = 0
      optimize_lattice_vectors(7) = 0
      optimize_lattice_vectors(8) = 0
      optimize_lattice_vectors(9) = 0
 401  if (.not. inp_a(test)) goto 420

      if (.not. inp_match(ntokens2, .false., test, tokens2, ind))
     >     call errquit(
     >           'optimize_lattice_vectors: unknown directive', 0,
     &       INPUT_ERR)

      goto (411,412,413,414,415,416,417,418,419) ind
c
 411      optimize_lattice_vectors(1) = 1
          goto 401
 412      optimize_lattice_vectors(2) = 1
          goto 401
 413      optimize_lattice_vectors(3) = 1
          goto 401
 414      optimize_lattice_vectors(4) = 1
          goto 401
 415      optimize_lattice_vectors(5) = 1
          goto 401
 416      optimize_lattice_vectors(6) = 1
          goto 401
 417      optimize_lattice_vectors(7) = 1
          goto 401
 418      optimize_lattice_vectors(8) = 1
          goto 401
 419      optimize_lattice_vectors(9) = 1
          goto 401

 420      goto 10

c
c  lattice_tolerances
c
 500  if (.not. inp_f(fe)) goto 501
      lattice_tolerances(1) = fe
      if (.not. inp_f(fe)) goto 501
      lattice_tolerances(2) = fe
      if (.not. inp_f(fe)) goto 501
      lattice_tolerances(3) = fe
      if (.not. inp_f(fe)) goto 501
      lattice_tolerances(4) = fe
      if (.not. inp_f(fe)) goto 501
      lattice_tolerances(5) = fe
      if (.not. inp_f(fe)) goto 501
      lattice_tolerances(6) = fe
 501  goto 10


c
c  lattice_vectors_tolerances
c
 600  if (.not. inp_f(fe)) goto 601
      lattice_tolerances(1) = fe
      if (.not. inp_f(fe)) goto 601
      lattice_tolerances(2) = fe
      if (.not. inp_f(fe)) goto 601
      lattice_tolerances(3) = fe
 601  goto 10

c
c  cycles
c
 700  if (.not. inp_i(cycles)) goto 701
 701  goto 10


*     ***** cell_optimize parameters to rtdb ****
 9999 continue
      if (cell_name.eq.' ') 
     >  call errquit('cell_optimize_input: nead a cell_name', 0,
     &       INPUT_ERR)


      ind = index(cell_name,' ') - 1
      value = rtdb_cput(rtdb,'cell_optimize:cell_name',
     >                  1,cell_name(1:ind))

      value = value.and. 
     >        rtdb_put(rtdb,'cell_optimize:optimize_strategy',
     >                 mt_int,1,optimize_strategy)
      value = value.and. 
     >        rtdb_put(rtdb,'cell_optimize:cycles',
     >                 mt_int,1,cycles)
      value = value.and. 
     >        rtdb_put(rtdb,'cell_optimize:optimize_lattice',
     >                 mt_int,6,optimize_lattice)
      value = value.and. 
     >        rtdb_put(rtdb,'cell_optimize:optimize_lattice_vectors',
     >                 mt_int,9,optimize_lattice_vectors)
      value = value.and. 
     >        rtdb_put(rtdb,'cell_optimize:lattice_tolerances',
     >                 mt_dbl,6,lattice_tolerances)



      return
      end
c $Id$
