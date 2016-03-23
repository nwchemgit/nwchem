c
c     QMD NAMD input handler
c
      subroutine qmd_namd_input(rtdb)
c
c     NAMD
c     END
c     
      implicit none
c
#include "errquit.fh"
#include "global.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "inp.fh"
#include "stdio.fh"
#include "case.fh"
c
      integer rtdb
c
      integer nstates
      integer state
      integer ks_spin
      double precision dt_elec
      logical do_tdks
      logical deco
c
      character*32 test
      character*30 pname
c
      pname = "qmd_namd_input: "
c
c     ------------------
c     Set default values
c     ------------------
c
c     Electronic dynamics time step: dt_elec
      dt_elec = 0.01d0
      if (.not.rtdb_put(rtdb,'qmd:dt_elec',mt_dbl,1,dt_elec))
     & call errquit(pname//'failed to write dt_elec',0,RTDB_ERR)
c
c     Number of states
      nstates = 2
      if (.not.rtdb_put(rtdb,'qmd:nstates',mt_int,1,nstates))
     & call errquit(pname//'failed to write nstates',0,RTDB_ERR)
c
c     Initial state
c      corresponds to first excited state
      state = 2
      if (.not.rtdb_put(rtdb,'qmd:init_state',mt_int,1,state))
     & call errquit(pname//'failed to write state',0,RTDB_ERR)
c
c     decoherence flag
      deco = .false.
      if (.not.rtdb_put(rtdb,'qmd:deco',mt_log,1,deco))
     & call errquit(pname//'failed to write deco',0,RTDB_ERR)
c
c     TDKS flag
      do_tdks = .false.
      if (.not.rtdb_put(rtdb,'qmd:tdks',mt_log,1,do_tdks))
     & call errquit(pname//'failed to write do_tdks',0,RTDB_ERR)
c
c     TDKS spin channel
      ks_spin = 1
      if (.not.rtdb_put(rtdb,'qmd:ks_spin',mt_int,1,ks_spin))
     & call errquit(pname//'failed to write ks_spin',0,RTDB_ERR)
c
c     ----------
c     Read input
c     ----------
c
 10   if (.not. inp_read()) 
     1  call errquit(pname//'failed reading input',0,INPUT_ERR)
      if (.not. inp_a(test)) 
     1  call errquit(pname//'failed reading keyword',0,INPUT_ERR)
c
c     Electronic dynamics time step: dt_elec
      if (inp_compare(.false.,test,'dt_elec')) then
       if (.not.inp_f(dt_elec)) then
          write(luOut,*) 'dt_elec: value not found; ',
     1      'default value of 0.01d0 will be used'
          dt_elec = 0.01d0
       endif
       if (.not.rtdb_put(rtdb,'qmd:dt_elec',mt_dbl,1,dt_elec))
     &  call errquit(pname//'failed to write dt_elec',0,RTDB_ERR)
c
c     Number of states
      else if (inp_compare(.false.,test,'nstates')) then
       if (.not.inp_i(nstates)) then
         write(luOut,*) 'nstates: value not found; ',
     1     'default value of 2 will be used'
         nstates=2
       end if
       if (.not.rtdb_put(rtdb,'qmd:nstates',mt_int,1,nstates))
     &  call errquit(pname//'failed to write nstates',0,
     &     RTDB_ERR)
c
c     Initial state
      else if (inp_compare(.false.,test,'init_state')) then
       if (.not.inp_i(state)) then
         write(luOut,*) 'init_state: value not found; ',
     1     'default value of 1 will be used'
         state=1
       end if
c      For input/output state runs from 0 to nstates-1
c      inside the code state runs from 1 to nstates
       state=state+1
       if (.not.rtdb_put(rtdb,'qmd:init_state',mt_int,1,state))
     &  call errquit(pname//'failed to write state',0,
     &     RTDB_ERR)
c
c     Decoherence?
      else if (inp_compare(.false.,test,'deco')) then
          deco=.true.
       if (.not.rtdb_put(rtdb,'qmd:deco',mt_log,1,deco))
     &   call errquit(pname//'failed to write deco',0,RTDB_ERR)
c
c     TDKS?
      else if (inp_compare(.false.,test,'tdks')) then
          do_tdks=.true.
       if (.not.rtdb_put(rtdb,'qmd:tdks',mt_log,1,do_tdks))
     &   call errquit(pname//'failed to write do_tdks',0,RTDB_ERR)
       if (.not.inp_i(ks_spin)) then
         write(luOut,*) 'tdks: value not found; ',
     1     'default value of 1 will be used'
         ks_spin=1
       end if
       if (.not.rtdb_put(rtdb,'qmd:ks_spin',mt_int,1,ks_spin))
     &  call errquit(pname//'failed to write ks_spin',0,
     &     RTDB_ERR)
c
c     END
      else if (inp_compare(.false.,test,'end')) then
        goto 20
c
      else
        call errquit(pname//'unknown directive',0, INPUT_ERR)
      endif
      goto 10
c
c ------
c Return
c ------
c
 20   return
      end
