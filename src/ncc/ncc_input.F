c===============================================================c
c                                                               c
c       NCC - a New Coupled-Cluster code for NWChem             c
c                                                               c
c       Developed by:                                           c
c                                                               c
c               Jeff R. Hammond                                 c
c               Leadership Computing Facility                   c
c               Argonne National Laboratory                     c
c               jhammond@mcs.anl.gov                            c
c                                                               c
c               Karol Kowalski                                  c
c               Environmental Molecular Sciences Laboratory     c
c               Pacific Northwest National Laboratory           c
c               karol.kowalski@pnl.gov                          c
c                                                               c
c               Marta Włoch                                     c
c               Department of Chemistry                         c
c               Michigan Technological University               c
c               wloch@mtu.edu                                   c
c                                                               c
c===============================================================c
      logical function ncc_input(rtdb)
c
c $Id$
c
      implicit none
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"
c===============================================================c
c                                                               c
#include "ncc_input.fh"                                         
c                                                               c
c      ...defines these common variables:                       c
c                                                               c
c      integer debugmode                                        c
c      double precision thresh                                  c
c      character*20 method                                      c
c                                                               c
c      I am using common variables to pass from this routine    c
c      into the NCC module so that I don't have to use the      c
c      RTDB since that would involve disk access.               c
c                                                               c
c===============================================================c
c
      integer rtdb
      integer itemp
c
      double precision dtemp
c
      character*20 test
      character*20 ctemp
c
      logical debug
      logical ltemp
c
#ifdef DEBUG_PRINT
      print*,'top of ncc_input'
      debug = .true.
#else
      debug = .false.
#endif
c
      ncc_input = .false.
c
c ----------
c Read input
c ----------
c
c
 10   if (.not. inp_read()) then
        call errquit('ncc_input: failed reading input',0,RTDB_ERR)
      endif
      if (.not. inp_a(test)) then
        call errquit('ncc_input: failed reading keyword',0,RTDB_ERR)
      endif
c
c     DEBUGMODE
c
      if (inp_compare(.false.,test,'debugmode')) then
        if (debug) write(LuOut,100) 'debugmode'
        if (.not.inp_i(itemp)) then
          write(LuOut,200) 'debugmode',100
          debugmode = 100
        else
          debugmode = itemp
        endif
#ifndef NO_RTDB_CALLS
        if (.not.rtdb_put(rtdb,'ncc:debugmode',mt_int,1,thresh)) then
          call errquit('ncc_input: rtdb_put failed',0,RTDB_ERR)
        endif
#endif
        if (debug) write(LuOut,150) debugmode
c
c     FREEZE (frozen cores/virtuals)
c
      elseif (inp_compare(.false.,test,'freeze')) then
        if (debug) write(LuOut,100) 'freeze'
        call freeze_input(rtdb,'ncc')
c
c     METHOD (ccsd, etc.)
c
      elseif (inp_compare(.false.,test,'method')) then
        if (debug) write(LuOut,100) 'method'
        if (.not.inp_a(ctemp)) then
          call errquit('ncc_input: no method given',0,RTDB_ERR)
        endif
        method = ctemp
#ifndef NO_RTDB_CALLS
        if (.not.rtdb_cput(rtdb,'ncc:method',1,ctemp)) then
          call errquit('ncc_input: rtdb_cput failed',0,RTDB_ERR)
        endif
#endif
        if (debug) write(LuOut,170) method
c
c     THRESH (convergence threshold for Davidson iteration)
c
      elseif (inp_compare(.false.,test,'thresh')) then
        if (debug) write(LuOut,100) 'thresh'
        if (.not.inp_f(dtemp)) then
          call errquit('ncc_input: no thresh given',0,RTDB_ERR)
        endif
        thresh = dtemp
#ifndef NO_RTDB_CALLS
        if (.not.rtdb_put(rtdb,'ncc:thresh',mt_dbl,1,thresh)) then
          call errquit('ncc_input: rtdb_put failed',0,RTDB_ERR)
        endif
#endif
        if (debug) write(LuOut,160) thresh
c
c     TILING
c
      elseif (inp_compare(.false.,test,'occtile')) then
        if (debug) write(LuOut,100) 'occtile'
        if (.not.inp_i(itemp)) then
          write(LuOut,200) 'occtile',10
          occtile = 10
        else
          occtile = itemp
        endif
#ifndef NO_RTDB_CALLS
        if (.not.rtdb_put(rtdb,'ncc:occtile',mt_int,1,thresh)) then
          call errquit('ncc_input: rtdb_put failed',0,RTDB_ERR)
        endif
#endif
        if (debug) write(LuOut,150) occtile
c
      elseif (inp_compare(.false.,test,'virtile')) then
        if (debug) write(LuOut,100) 'virtile'
        if (.not.inp_i(itemp)) then
          write(LuOut,200) 'virtile',25
          virtile = 25
        else
          virtile = itemp
        endif
#ifndef NO_RTDB_CALLS
        if (.not.rtdb_put(rtdb,'ncc:virtile',mt_int,1,thresh)) then
          call errquit('ncc_input: rtdb_put failed',0,RTDB_ERR)
        endif
#endif
        if (debug) write(LuOut,150) virtile
c
c     PRINT
c
      elseif (inp_compare(.false.,test,'print')) then
        if (debug) write(LuOut,100) 'print'
        call util_print_input(rtdb,'ncc')
c
c     END
c
      elseif (inp_compare(.false.,test,'end')) then
        if (debug) write(LuOut,100) 'end'
        goto 20
      else
        call errquit('ncc_input: unknown directive',0,INPUT_ERR)
      endif
      goto 10
c
c ------
c Return
c ------
c
 20   continue
c
      ncc_input = .true.
c
#ifdef DEBUG_PRINT
      print*,'end of ncc_input'
#endif
c
      return
c
 100  format(8x,'keyword found:         ',a20)
 150  format(8x,'using value:           ',i20)
 160  format(8x,'using value:           ',f20.8)
 170  format(8x,'using value:           ',a20)
 200  format(8x,'value not found for:   ',a20,/,
     &       8x,'using default value of ',i20)
c
      end
