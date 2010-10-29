c===============================================================c
c                                                               c
c       Socrates - an NWChem module for teaching                c
c       Copyright © 2009 Jeff Hammond                           c
c                                                               c
c       Developed by:                                           c
c                                                               c
c               Jeff R. Hammond                                 c
c               Leadership Computing Facility                   c
c               Argonne National Laboratory                     c
c               jhammond@mcs.anl.gov                            c
c                                                               c
c===============================================================c
      logical function socrates_input(rtdb)
c
c $Id$
c
      implicit none
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
c
      integer rtdb
      integer itemp
      integer scf_alg
      integer nopen
c
      double precision dtemp
      double precision tol2e
      double precision int_thresh
c
      character*20 test
      character*20 ctemp
c
      logical debug
      logical nodezero
      logical ltemp
c
#ifdef DEBUG_PRINT
      debug = (GA_NodeId().eq.0) ! debug print on nodezero only
c      debug = .true. ! debug print everywhere
#else
      parameter debug = .false.
#endif
      if (debug) write(LuOut,*) 'top of socrates_input'
c
      nodezero = (GA_NodeId().eq.0) 
c
      socrates_input = .false.
c
c ----------
c Read input
c ----------
c
c
 10   if (.not. inp_read()) then
        call errquit(__FILE__,__LINE__,INPUT_ERR)
      endif
      if (.not. inp_a(test)) then
        call errquit(__FILE__,__LINE__,INPUT_ERR)
      endif
c
c     PRINT
c
      if (inp_compare(.false.,test,'print')) then
        if (debug) write(LuOut,100) 'print'
        call util_print_input(rtdb,'socrates')
c
c     SCF_ALG
c
      elseif (inp_compare(.false.,test,'scf_alg')) then
        if (debug) write(LuOut,100) 'scf_alg'
        if (.not.inp_i(scf_alg)) then
           if (nodezero) write(LuOut,*) 'usage: scf_alg <int>'
           call errquit(__FILE__,__LINE__,INPUT_ERR)
        endif
        if (.not.rtdb_put(rtdb,'socrates:scf_alg',mt_int,1,
     1                    scf_alg)) then
           call errquit(__FILE__,__LINE__,RTDB_ERR)
        endif
c
c     NOPEN
c
      elseif (inp_compare(.false.,test,'nopen')) then
        if (debug) write(LuOut,100) 'nopen'
        if (.not.inp_i(nopen)) then
           if (nodezero) write(LuOut,*) 'usage: nopen <int>'
           call errquit(__FILE__,__LINE__,INPUT_ERR)
        endif
        if (.not.rtdb_put(rtdb,'socrates:nopen',mt_int,1,
     1                    nopen)) then
           call errquit(__FILE__,__LINE__,RTDB_ERR)
        endif
c
c     INT_THRESH
c
      elseif (inp_compare(.false.,test,'int_thresh')) then
        if (debug) write(LuOut,100) 'int_thresh'
        if (.not.inp_f(int_thresh)) then
           if (nodezero) write(LuOut,*) 'usage: int_thresh <float>'
           call errquit(__FILE__,__LINE__,INPUT_ERR)
        endif
        if (.not.rtdb_put(rtdb,'socrates:int_thresh',mt_dbl,1,
     1                    int_thresh)) then
           call errquit(__FILE__,__LINE__,RTDB_ERR)
        endif
c
c     TOL2E
c
      elseif (inp_compare(.false.,test,'tol2e')) then
        if (debug) write(LuOut,100) 'tol2e'
        if (.not.inp_f(tol2e)) then
           if (nodezero) write(LuOut,*) 'usage: tol2e <float>'
           call errquit(__FILE__,__LINE__,INPUT_ERR)
        endif
        if (.not.rtdb_put(rtdb,'socrates:tol2e',mt_dbl,1,tol2e)) then
           call errquit(__FILE__,__LINE__,RTDB_ERR)
        endif
c
c     END
c
      else if (inp_compare(.false.,test,'end')) then
        goto 20
      else
        call errquit(__FILE__,__LINE__,INPUT_ERR)
      endif
      goto 10
c
c ------
c Return
c ------
c
 20   continue
c
      socrates_input = .true.
c
      if (debug) write(LuOut,*) 'end of socrates_input'
c
      return
c
 100  format(8x,'keyword found:         ',a20)
 150  format(8x,'using value:           ',i20)
 160  format(8x,'using value:           ',f20.8)
 170  format(8x,'using value:           ',a20)
 200  format(8x,'value not found for:   ',a20,/,
     &          'using default value of ',i20)
c
      end
