      subroutine xtb_input(rtdb)
#ifndef USE_TBLITE
      implicit none
#include "errquit.fh"
      integer rtdb
      call errquit('Recompile NWChem with TBLITE support',0,0)
#else
c
c     Input parser for TBLITE@NWCHEM
c
c     Daniel Mejia-Rodriguez
c     November 2020
c
      implicit none
c
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util_params.fh"
c
      integer rtdb
      character*255 test
c
      integer, parameter :: num_dirs = 11
      character*11 dirs(num_dirs)
      data dirs /'acc','uhf', 'method', 'etemp', 'guess', 
     $           'verbosity', 'print', 'noprint', 'nspin', 
     $           'broyden', 'end'/
c
      integer uhf, ind, verb, nspin
      double precision acc, etemp, broyden

      character(len=4) method
c
c     Parse the input file
c
 10   if (.not. inp_read()) then
        call errquit('xtb_input: failed to read input',0,INPUT_ERR)
      endif
c
      if (.not. inp_a(test)) then
        call errquit('xtb_input: failed to read keyword',0,INPUT_ERR)
      endif
c
      if (.not.inp_match(num_dirs, .false., test, dirs, ind)) then
        call errquit('xtb_input: unknown directive', 0, INPUT_ERR)
      end if
c
      goto (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 9999) ind

      call errquit('xtb_input: unknown directive',ind, INPUT_ERR)
c
c     -----------------------------------------------------------------
c
c     accuracy
c
 100  if (.not. inp_f(acc))
     $   call errquit('xtb_input: accuracy value not found',0,INPUT_ERR)
      if (.not.rtdb_put(rtdb, 'xtb:acc', mt_dbl, 1, acc))
     $  call errquit('xtb_input: rtdb_put failed',100,RTDB_ERR)
      goto 10
c
c     number of unpaired electrons
c
 200  if (.not. inp_i(uhf))
     $  call errquit('xtb_input: no. unpaired e not found',0,INPUT_ERR)
      if (.not.rtdb_put(rtdb, 'xtb:uhf', mt_int, 1, uhf))
     $  call errquit('xtb_input: rtdb_put failed', 200, RTDB_ERR)
      goto 10
c
c     method gff, gfn1 or gfn2
c
 300  if (.not.inp_a(test))
     $  call errquit('xtb_input: method not found',0,INPUT_ERR)
      if (inp_compare(.false.,'gfn1',test)) then
        method = 'gfn1'
      elseif (inp_compare(.false.,'gfn2',test)) then   
        method = 'gfn2'
      else
        call errquit('xtb_input: unknown xTB method',0,INPUT_ERR)
      endif
      if (.not.rtdb_cput(rtdb, 'xtb:method', 1, method))
     $  call errquit('xtb_input: rtdb_put failed',303, RTDB_ERR)
      goto 10
c
c     electronic temperature
c
 400  if (.not. inp_f(etemp))
     $   call errquit('xtb_input: temperature not found',0,INPUT_ERR)
      etemp = etemp*kb_au
      if (.not.rtdb_put(rtdb, 'xtb:kt', mt_dbl, 1, etemp))
     $  call errquit('xtb_input: rtdb_put failed',400,RTDB_ERR)
      goto 10
c
c     guess
c
 500  if (.not.rtdb_put(rtdb, 'xtb:guess', mt_log, 1, .true.))
     $  call errquit('xtb_input: rtdb_put failed',500,RTDB_ERR) 
      goto 10
c
c     verbosity level
c
  600 if(.not.inp_i(verb))
     $  call errquit('xtb_input: verb level not found',600,INPUT_ERR)
      if (.not.rtdb_put(rtdb, 'xtb:verb', mt_int, 1, verb))
     $  call errquit('xtb_input: rtdb_put failed',600,RTDB_ERR)
      goto 10
c
c     print
c
  700 call util_print_input(rtdb, 'xtb')
      goto 10

  800 call util_print_input(rtdb, 'xtb')
      goto 10
c
c     spin-polarization treatment
c
  900 if(.not.inp_i(nspin))
     $  call errquit('xtb_input: nspin value not found',900,INPUT_ERR)
      if (.not.(nspin.eq.1 .or. nspin.eq.2)) then
        call errquit('xtb_input: nspin can be 1 or 2',900,INPUT_ERR)
      endif
      if (.not.rtdb_put(rtdb, 'xtb:nspin', mt_int, 1, nspin))
     $  call errquit('xtb_input: rtdb_put failed',900,RTDB_ERR)
      goto 10
c
c     broyden mixer parameter
c
 1000 if(.not.inp_f(broyden))
     $  call errquit('xtb_input: broyden mixer value not found',1000,
     $                INPUT_ERR)
      if (broyden.le.0.0d0 .or. broyden .gt. 1.0d0)
     $  call errquit('xtb_input: broyden value should lie (0,1]',1001,
     $                INPUT_ERR)
      if (.not.rtdb_put(rtdb, 'xtb:broyden', mt_dbl, 1, broyden))
     $  call errquit('xtb_input: rtdb_put failed',1002,RTDB_ERR)
      goto 10 
c
c     normal termination
c
 9999 return
c
#endif
      end subroutine
