      subroutine rt_tddft_input(rtdb)
c
c Input parser for RT-TDDFT called by input/input_parse.F.
c Also sets default values for input parameters.
c
c     RT_TDDFT
c     END
c     
      implicit none
c
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"
c
      integer rtdb
      character*20 test
c
      integer ntime
      double precision delta
      logical sodft
      logical dft
      logical rk4
      logical euler
      logical magnus
      double precision field_Emax
      double precision field_w
c
c     == Set default values ==
c
c     == NTIME (number of time steps) ==
c
      ntime=100
      if (.not.rtdb_put(rtdb,'rt_tddft:ntime',mt_int,1,ntime))
     & call errquit('rt_tddft_input: write failed to rtdb',0, RTDB_ERR)
c
c     == DELTA (time step in atomic units) ==
c
      delta=0.02
      if (.not.rtdb_put(rtdb,'rt_tddft:delta',mt_dbl,1,delta))
     & call errquit('rt_tddft_input: write failed to rtdb',0, RTDB_ERR)
c
c     == SODFT (dynamics with spin-orbit interactions) ==
c
      sodft=.true.
      if (.not.rtdb_put(rtdb,'rt_tddft:sodft',mt_log,1,sodft))
     & call errquit('rt_tddft_input: write failed to rtdb',0, RTDB_ERR)
c
c     == DFT (dynamics without spin-orbit interactions) ==
c
      dft=.false.
      if (.not.rtdb_put(rtdb,'rt_tddft:dft',mt_log,1,dft))
     & call errquit('rt_tddft_input: write failed to rtdb',0, RTDB_ERR)
c
c     == RK4 (4th order Runge-Kutta propagator) ==
c
      rk4= .true.
      if (.not.rtdb_put(rtdb,'rt_tddft:rk4',mt_log,1,rk4))
     & call errquit('rt_tddft_input: write failed to rtdb',0, RTDB_ERR)
c
c     == Euler (Euler propagator) ==
c
      euler = .false.
      if (.not.rtdb_put(rtdb,'rt_tddft:euler',mt_log,1,euler))
     & call errquit('rt_tddft_input: write failed to rtdb',0, RTDB_ERR)
c
c     == Magnus (Magnus propagator) ==
c
      magnus = .false.
      if (.not.rtdb_put(rtdb,'rt_tddft:magnus',mt_log,1,magnus))
     & call errquit('rt_tddft_input: write failed to rtdb',0, RTDB_ERR)
c
c     == Electric field (time-dependent electric field influencing the dynamics) ==
c
      field_Emax = 0.05d0    ! field maximum
      field_w = 0.4764d0     ! field omega
      if (.not.rtdb_put(rtdb,'rt_tddft:field_Emax',mt_dbl,1,field_Emax))
     & call errquit('rt_tddft_input: write failed to rtdb',0, RTDB_ERR)
      if (.not.rtdb_put(rtdb,'rt_tddft:field_w',mt_dbl,1,field_w))
     & call errquit('rt_tddft_input: write failed to rtdb',0, RTDB_ERR)
c
c     == Parse the input ==
c
 10   if (.not. inp_read()) 
     &  call errquit('rt_tddft_input: read failed input',0, INPUT_ERR)
      if (.not. inp_a(test)) 
     &  call errquit('rt_tddft_input: read failed keyword',0, INPUT_ERR)
c
c     == NTIME (number of time steps) ==
c
      if (inp_compare(.false.,test,'ntime')) then
        if (.not.inp_i(ntime)) then
          write(LuOut,*) 'rt_tddft_input: ntime value not found; ',
     &      'default value of 100 will be used'
          ntime=100
        endif
        if (.not.rtdb_put(rtdb,'rt_tddft:ntime',mt_int,1,ntime))
     &   call errquit('rt_tddft_input: write failed to rtdb',0,RTDB_ERR)
c
c     == DELTA ==
c
      else if (inp_compare(.false.,test,'delta')) then
        if (.not.inp_f(delta)) then
          write(LuOut,*) 'rt_tddft_input: delta value not found; ',
     &      'default value of 0.02 atomic units will be used'
          delta=0.02
        endif
        if (.not.rtdb_put(rtdb,'rt_tddft:delta',mt_dbl,1,delta))
     &   call errquit('rt_tddft_input: write failed to rtdb',0,RTDB_ERR)
c
c     == SODFT ==
c
      else if (inp_compare(.false.,test,'sodft')) then
        sodft=.true.
        if (.not.rtdb_put(rtdb,'rt_tddft:sodft',mt_log,1,sodft))
     &   call errquit('rt_tddft_input: write failed to rtdb',0,RTDB_ERR)
c
c     == DFT ==
c
      else if (inp_compare(.false.,test,'dft')) then
        dft=.true.
        if (.not.rtdb_put(rtdb,'rt_tddft:dft',mt_log,1,dft))
     &   call errquit('rt_tddft_input: write failed to rtdb',0,RTDB_ERR)
c
c     == RK4 ==
c
      else if (inp_compare(.false.,test,'rk4')) then
        rk4=.true.
        if (.not.rtdb_put(rtdb,'rt_tddft:rk4',mt_log,1,rk4))
     &   call errquit('rt_tddft_input: write failed to rtdb',0,RTDB_ERR)
c
c     == Euler ==
c
      else if (inp_compare(.false.,test,'euler')) then
        euler=.true.
        if (.not.rtdb_put(rtdb,'rt_tddft:euler',mt_log,1,euler))
     &   call errquit('rt_tddft_input: write failed to rtdb',0,RTDB_ERR)
c
c     == Magnus ==
c
      else if (inp_compare(.false.,test,'magnus')) then
        magnus=.true.
        if (.not.rtdb_put(rtdb,'rt_tddft:magnus',mt_log,1,magnus))
     &   call errquit('rt_tddft_input: write failed to rtdb',0,RTDB_ERR)
c
c     == Field Description (maximum value) ==
c
      else if (inp_compare(.false.,test,'field_Emax')) then
      if (.not.inp_f(field_Emax)) then
          write(LuOut,*) 'rt_tddft_input: delta value not found; ',
     &      'default value of 0.05 atomic units will be used'
          field_Emax=0.05
      endif
      if (.not.rtdb_put(rtdb,'rt_tddft:field_Emax',mt_dbl,1,field_Emax))
     & call errquit('rt_tddft_input: write failed to rtdb',0,RTDB_ERR)
c
c     == Field Description (omega) ==
c
      else if (inp_compare(.false.,test,'field_w')) then
      if (.not.inp_f(field_w)) then
          write(LuOut,*) 'rt_tddft_input: delta value not found; ',
     &      'default value of 0.5 atomic units will be used'
          field_w=0.5
      endif
      if (.not.rtdb_put(rtdb,'rt_tddft:field_Emx',mt_dbl,1,field_w))
     & call errquit('rt_tddft_input: write failed to rtdb',0,RTDB_ERR)
c
c     == Print ==
c
      else if (inp_compare(.false.,test,'print')) then
        call util_print_input(rtdb,'rt_tddft')
c
c     == End ==
c
      else if (inp_compare(.false.,test,'end')) then
        goto 20
      else
        call errquit('rt_tddft_input: unknown directive',0, INPUT_ERR)
      endif
      goto 10
c
c     == Return ==
c
 20   return
      end
