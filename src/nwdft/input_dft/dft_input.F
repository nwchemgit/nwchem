      subroutine dft_input(rtdb)
c
c     $Id: dft_input.F,v 1.26 2001-03-30 00:00:47 edo Exp $
c
c     The main purpose of this subroutine is to parse the input
c     and place any "found" parameters into the rtdb.
c
      implicit none
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "inp.fh"
#include "stdio.fh"
c     
c     read input for DFT calculation.  All input is of the form:
c     <keyword>  or   <keyword> <value>  or   <keyword> <value><units>
c     
      integer rtdb ! [input] handle to database
c
      logical status
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 23)
      character*12 dirs(num_dirs)
      character*255 test
      integer idecomp, mult, iterations, ipol, noio,
     &        idirect, imull
      integer nopen
      logical pseudo, direct, lmaxov, incore
      double precision smear_sigma
      data dirs /'convergence', 'decomp', 'dft', 'direct', 'ecp', 
     $           'grid', 'vectors', 'incore', 'iterations', 'max_ovl', 
     $           'mulliken', 'mult', 'noio', 'noprint', 'odft', 
     $           'print', 'tolerances', 'xc', 'sym', 'adapt', 'sic', 
     .           'smear',
     &           'end'/
c
c     Any DFT input causes calculation to be marked as not converged
c
      status = rtdb_delete(rtdb,'dft:converged')
c
c     Lebedev on by default
c
c
 10   if (.not. inp_read()) call errquit
     $     ('dft_input: inp_read failed', 0)
      if (.not. inp_a(test))
     $     call errquit('dft_input: failed to read keyword', 0)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     $     call errquit('dft_input: unknown directive', 0)
c     
      goto (100,  200,  300,  400,  500,  600,  700,  800,  900, 
     $     1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800,
     $     1900, 2000, 2100, 2200, 9999) ind
      call errquit('dft_input: unimplemented directive', ind)
c
c     convergence;
c             
  100 call conv_input(rtdb, 'dft')
c
      goto 10
c             
c     decomp;
c             
  200 idecomp = 1
      if (.not. rtdb_put(rtdb, 'dft:idecomp', mt_int, 1, idecomp))
     &   call errquit('dft_input: rtdb_put failed', 200)
c
      goto 10
c
c     dft;
c             
  300 ipol=1
      if (.not. rtdb_put(rtdb, 'dft:ipol', mt_int, 1, ipol))
     &   call errquit('dft_input: rtdb_put failed', 300)
c
      goto 10
c
c     direct;
c             
  400 direct = .true.
      if (.not. rtdb_put(rtdb, 'dft:direct', mt_log, 1, direct))
     &   call errquit('dft_input: rtdb_put failed', 400)
      idirect = 1
      if (.not. rtdb_put(rtdb, 'dft:idirect', mt_int, 1, idirect))
     &   call errquit('dft_input: rtdb_put failed', 400)
         if (.not. rtdb_put(rtdb,'int2e:filesize', mt_int, 1, -1))
     $        call errquit('dft_input: failed to store filesize',0)
         if (.not. rtdb_put(rtdb,'int2e:memsize', mt_int, 1, -1))
     $        call errquit('dft_input: failed to store filesize',0)
c
      goto 10
c
c     ecp;
c             
  500 pseudo = .true.
c
      goto 10
c
c     grid;
c             
  600 call grid_input(rtdb, 'dft')
c
      goto 10
c             
c     vectors; get input vectors, output vectors, initial guess, etc.
c             
  700 call vectors_input(rtdb,'dft')
c
      goto 10
c
c     incore; try to do an incore coulomb construction
c             
  800 incore = .true.
      if (.not.rtdb_put(rtdb,'dft:incore', mt_log, 1, incore)) 
     &   call errquit('dft_input: rtdb_put failed', 800)
c
      goto 10
c
c     iterations;
c             
  900 if (.not. inp_i(iterations)) then
         write(LuOut,*)' dft_input: iterations value not found;',
     &             ' defaulting to 30'
         iterations = 30
      endif 
      if (.not. rtdb_put(rtdb, 'dft:iterations', mt_int, 1, iterations))
     &   call errquit('dft_input: rtdb_put failed', 900) 
c
      goto 10
c
c     max_ovl; Move orbitals around in movecs to maintain max 
c              overlap between iterations
c             
 1000 lmaxov = .true.
      if (.not. rtdb_put(rtdb, 'dft:max_ovl', mt_log, 1, lmaxov))
     &   call errquit('dft_input: rtdb_put failed', 1000)
c
      goto 10
c
c     mulliken; do Mulliken pop analysis
c             
 1100 imull=1
      if (.not. rtdb_put(rtdb, 'dft:imull', mt_int, 1, imull))
     &   call errquit('dft_input: rtdb_put failed', 1100)
c
      goto 10
c             
c     mult;
c             
 1200 if (.not. inp_i(mult)) then
         write(LuOut,*)' dft_input: spin multiplicity value not found;',
     &             ' ignoring option using default of singlet'
         mult = 1
      endif 
      if (.not. rtdb_put(rtdb, 'dft:mult', mt_int, 1, mult))
     &   call errquit('dft_input: rtdb_put failed', 1200)
      nopen = mult - 1
      if (.not.rtdb_put(rtdb, 'dft:nopen',mt_int,1,nopen))
     &   call errquit('dft_input: rtdb_put failed', 1201)
c
      goto 10
c             
c     noio;
c             
 1300 noio = 1
      if (.not. rtdb_put(rtdb, 'dft:noio', mt_int, 1, noio))
     &   call errquit('dft_input: rtdb_put failed', 1300)
c
      goto 10
c
c     noprint;
c
 1400 call util_print_input(rtdb, 'dft')
c
      goto 10
c
c     odft;
c             
 1500 ipol=2
      if (.not. rtdb_put(rtdb, 'dft:ipol', mt_int, 1, ipol))
     &   call errquit('dft_input: rtdb_put failed', 1500)
c
      goto 10
c
c     print;
c
 1600 call util_print_input(rtdb, 'dft')
c
      goto 10
c
c     tolerances;
c             
 1700 call tol_input(rtdb, 'dft')
c
      goto 10
c
c     xc;
c             
 1800 call xc_input(rtdb, 'dft')
c
      goto 10
c
c     sym [off|on] (defaults to on)
c
 1900 if (.not. inp_a(test)) test = 'on'
      status = inp_compare(.false.,test,'on')
      if (.not. rtdb_put(rtdb,'dft:skeleton', MT_LOG, 1, status))
     &   call errquit('dft_input: rtdb_put failed', 1900)
      goto 10
c
c     adapt [off|on] (defaults to on)
c
 2000 if (.not. inp_a(test)) test = 'on'
      status = inp_compare(.false.,test,'on')
      if (.not. rtdb_put(rtdb,'dft:adapt', MT_LOG, 1, status))
     &   call errquit('dft_input: rtdb_put failed', 2000)
      goto 10
c
 2100 call sic_input(rtdb, 'dft')
c
      goto 10
c
c     smearing
c
 2200 if (.not. inp_f(smear_sigma)) then
         write(LuOut,*)'  smearing value not found;',
     &             ' defaulting to 1d-3 Hartree.'
         smear_sigma = 1.d-3
      endif
      if (.not.rtdb_put(rtdb, 'dft:smear_sigma',mt_dbl,1,smear_sigma))
     &   call errquit('dft_inp: rtdb_put failed', 2200)
      goto 10
c
c             
c     end;
c             
 2300 continue
 9999 continue
      return 
c       
      end
