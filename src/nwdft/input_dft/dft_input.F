      subroutine dft_input(rtdb)
c
c     $Id: dft_input.F,v 1.16 1997-02-17 17:58:32 d3g681 Exp $
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
c     
c     read input for DFT calculation.  All input is of the form:
c     <keyword>  or   <keyword> <value>  or   <keyword> <value><units>
c     
      integer rtdb ! [input] handle to database
c
      logical status
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 19)
      character*12 dirs(num_dirs)
      character*255 test
      data dirs /'convergence', 'decomp', 'dft', 'direct', 'ecp', 
     $           'grid', 'guess', 'incore', 'itrscf', 'max_ovl', 
     $           'mulliken', 'mult', 'noio', 'noprint', 'odft', 
     $           'print', 'tolerances', 'xc', 'end'/
c      integer nopt,  nopt_quad, nopt_munit, nopt_grid,
c     &     nopt_system
c      parameter (nopt=57,  nopt_quad=2, nopt_munit=6,
c     &     nopt_grid=5, nopt_system=4)
      integer  n,
     &     idecomp, mult, 
     &     itrscf, ipol, noio,
     &     iguess,
     &     idirect, nfock, ifield,
     &     imull
      integer nfield, iline
      logical pseudo
      logical direct
      logical lmaxov, incore
c
c     Any DFT input causes calculation to be marked as not converged
c
      status = rtdb_delete(rtdb,'dft:converged')
c
 10   if (.not. inp_read()) call errquit
     $     ('dft_input: inp_read failed', 0)
      if (.not. inp_a(test))
     $     call errquit('dft_input: failed to read keyword', 0)
      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     $     call errquit('dft_input: unknown directive', 0)
c     
*      write(6,*)' test = ',test     
*      write(6,*)' ind = ',ind     
      goto (100,  200,  300,  400,  500,  600,  700,  800,  900, 
     $     1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800,
     $     1900, 9999) ind
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
c     guess; next keyword gives type of initial guess
c             
  700 if (inp_a(test)) then
         if (inp_compare(.false., 'atomic', test)) then
            iguess = 0
         elseif (inp_compare(.false., 'hcore', test)) then
            iguess = 1
         else 
            write(6,*)' dft_input: guess option not',
     &                ' recognized'
         endif 
      else 
         write (6,*)' dft_input:  WARNING:',
     &              ' guess option not found;',
     &              ' defaulting to 0'
         iguess = 0
      endif 
      if (.not. rtdb_put(rtdb, 'dft:iguess', mt_int, 1, iguess))
     &   call errquit('dft_input: rtdb_put failed', 700)
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
c     itrscf;
c             
  900 if (.not. inp_i(itrscf)) then
         write(6,*)' dft_input: itrscf value not found;',
     &             ' defaulting to 30'
         itrscf = 30
      endif 
      if (.not. rtdb_put(rtdb, 'dft:itrscf', mt_int, 1, itrscf))
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
         write(6,*)' dft_input: spin multiplicity value not found;',
     &             ' ignoring option'
         mult = 0
      endif 
      if (.not. rtdb_put(rtdb, 'dft:mult', mt_int, 1, mult))
     &   call errquit('dft_input: rtdb_put failed', 1200)
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
c     end;
c             
 1900 continue
 9999 continue
      return 
c       
      end
