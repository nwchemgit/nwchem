      subroutine tol_input(rtdb, module)
C$Id$
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      integer rtdb
      character*(*) module
c     
c     Parse the tolerance directive which customizes the screening
c     parameters for the DFT procedure.  Insert options into the 
c     database prefixed with module:...
c     
c     Possible variables are:
c
c     tight
c     tol_rho
c     accAOfunc
c     accCDfunc
c     accXCfunc
c     accCoul
c     accQrad
c     accXCfunc
c
c     Any default assumptions are based upon "no input" defaults.
c     
      integer num_dirs, ind, mlen
      parameter (num_dirs = 8)
      character*9 dirs(num_dirs)
      character*255 test
      logical ltight
      integer iAOacc, iCDacc, itol2e, irqmax, iXCacc
      double precision tol_rho,ddumm
      data dirs /'tight', 'tol_rho', 'accAOfunc', 'accCDfunc', 
     $           'accCoul', 'accQrad', 'accXCfunc', 'radius'/
c     
      mlen = inp_strlen(module)
c     
 10   if (.not. inp_a(test)) goto 999
c     
      if (.not. inp_match(num_dirs, .false., test, dirs, ind)) then
c     
c     Does not match a keyword ... 
c     
         goto 10000
      endif
c     
      goto (100, 200, 300, 400, 500, 600, 700, 800, 999) ind
      call errquit('tol_inp: unimplemented directive', ind, INPUT_ERR)
c     
c     tight; set screening threshholds tight
c             
  100 ltight = .true.
      if (.not. rtdb_put(rtdb, 'dft:ltight', mt_log, 1, ltight))
     &   call errquit('tol_inp: rtdb_put failed', 100, RTDB_ERR)
c
      goto 10
c     
c     tol_rho; a density screening parameter used to avoid computation 
c              of xc functional on a grid point and than be subsequently 
c              multiplied by a density element less than tol_rho.
c     
  200 if (.not. inp_f(tol_rho)) then
         write(LuOut,*)' tol_inp: tol_rho value not found; defaulting ',
     &             ' to 1.d-10 '
         tol_rho=1.D-10
      endif            
      if (.not. rtdb_put(rtdb, 'dft:tol_rho', mt_dbl, 1, tol_rho))
     &   call errquit('tol_inp: rtdb_put failed', 200, RTDB_ERR)
c
      goto 10
c     
c     accAOfunc; for use in int_eval_gbsets (called by quadvxcs),
c                accuracy for evaluating AO gaussians
c             
  300 if (.not. inp_i(iAOacc)) then
         write(LuOut,*)' tol_inp: accAOfunc value not found; ',
     &                 'defaulting to 14.'
         iAOacc = 14
      endif            
      if (.not. rtdb_put(rtdb, 'dft:iAOacc', mt_int, 1, iAOacc))
     &   call errquit('tol_inp: rtdb_put failed', 300, RTDB_ERR)
c
      goto 10
c     
c     accCDfunc; for use in screening evaluation of CD gaussians
c                on grid
c             
  400 if (.not. inp_i(iCDacc)) then
         write(LuOut,*)' tol_inp: accCDfunc value not found; ',
     &                 'defaulting to 20.'
         iCDacc = 20
      endif            
      if (.not. rtdb_put(rtdb, 'dft:iCDacc', mt_int, 1, iCDacc))
     &   call errquit('tol_inp: rtdb_put failed', 400, RTDB_ERR)
c
      goto 10
c     
c     accCoul; for use in fitcd & fitvcoul for schwarz screening 
c             
  500 if (.not. inp_i(itol2e)) then
      call errquit(' cannot found value for acccoul ',0, INPUT_ERR)
      endif
      if (.not. rtdb_put(rtdb, 'dft:itol2e', mt_int, 1, itol2e))
     &   call errquit('tol_inp: rtdb_put failed', 500, RTDB_ERR)
c
      goto 10
c     
c     accQrad; for use in setquad - r_max for quadrature points (distance in bohr)
c             
  600 if (.not. inp_i(irqmax)) then
         write(LuOut,*)' tol_inp: accQrad value not found; defaulting ',
     &             ' to 12.'
         irqmax = 12
      endif            
      if (.not. rtdb_put(rtdb, 'dft:irqmax', mt_int, 1, irqmax))
     &   call errquit('tol_inp: rtdb_put failed', 600, RTDB_ERR)
c
      goto 10
c     
c     accXCfunc; for use in int_eval_gbsets (called by quadvxc1),
c                accuracy for evaluating XC gaussians
c             
  700 if (.not. inp_i(iXCacc)) then
         write(LuOut,*)' tol_inp: accXCfunc value not found; ',
     &                 'defaulting to 20.'
         iXCacc = 20
      endif            
      if (.not. rtdb_put(rtdb, 'dft:iXCacc', mt_int, 1, iXCacc))
     &   call errquit('tol_inp: rtdb_put failed', 700, RTDB_ERR)
c
      goto 10
c     
c     radius; screening radius used in computation of spatial weights (distance in bohr)
c     
  800 if (.not. inp_f(ddumm)) then
cedo         radius=25.d0
      endif            
      write(LuOut,*)' WARNING: radius keyword NOT recognized '
c
      goto 10
c     
 999  continue
c     
      return
c     
10000 write(LuOut,10001)
10001 format(/,' tolerances [[tight], [tol_rho dtol], ',
     &       /,' [accAOfunc itol], [accCDfunc itol], ',
     &       /,' [accCoul itol], [accQrad itol], ',
     &       /,' [accXCfunc itol] ')
      call util_flush(LuOut)
      call errquit('tol_input: invalid format', 0, INPUT_ERR)
c     
      end

