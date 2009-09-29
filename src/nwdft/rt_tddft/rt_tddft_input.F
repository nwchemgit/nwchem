c Input parser for RT_TDDFT called by input/input_parse.F.
c Sets default values for input parameters.
c
c     RT_TDDFT
c     END
c
      subroutine rt_tddft_input(rtdb)
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
      double precision tmax
      double precision dt
      logical sodft
      logical dft
      logical rk4
      logical euler
      logical magnus
      double precision field_Emax
      double precision field_w
      logical nuclear
      character*32 pname
c      
      pname = "rt_tddft_input:"
c
c     == Set default values ==
c
c     == TMAX (maximum time) ==
      tmax=100
      if (.not.rtdb_put(rtdb,'rt_tddft:tmax',mt_dbl,1,tmax))
     & call errquit(pname//'Write failed to rtdb',0, RTDB_ERR)
c
c     == DELTA (time step in atomic units) ==
      dt=0.02
      if (.not.rtdb_put(rtdb,'rt_tddft:dt',mt_dbl,1,dt))
     & call errquit(pname//'Write failed to rtdb',0, RTDB_ERR)
c
c     == SODFT (dynamics with spin-orbit interactions) ==
      sodft=.true.
      if (.not.rtdb_put(rtdb,'rt_tddft:sodft',mt_log,1,sodft))
     & call errquit(pname//'Write failed to rtdb',0, RTDB_ERR)
c
c     == DFT (dynamics without spin-orbit interactions) ==
      dft=.false.
      if (.not.rtdb_put(rtdb,'rt_tddft:dft',mt_log,1,dft))
     & call errquit(pname//'Write failed to rtdb',0, RTDB_ERR)
c
c     == RK4 (4th order Runge-Kutta propagator) ==
      rk4= .true.
      if (.not.rtdb_put(rtdb,'rt_tddft:rk4',mt_log,1,rk4))
     & call errquit(pname//'Write failed to rtdb',0, RTDB_ERR)
c
c     == EULER (Euler propagator) ==
      euler = .false.
      if (.not.rtdb_put(rtdb,'rt_tddft:euler',mt_log,1,euler))
     & call errquit(pname//'Write failed to rtdb',0, RTDB_ERR)
c
c     == MAGNUS (Magnus propagator) ==
      magnus = .false.
      if (.not.rtdb_put(rtdb,'rt_tddft:magnus',mt_log,1,magnus))
     & call errquit(pname//'Write failed to rtdb',0, RTDB_ERR)
c
c     == ELECTRIC FIELD ==
      field_Emax = 0.05d0    ! maximum value
      field_w = 0.4764d0     ! omega
      if (.not.rtdb_put(rtdb,'rt_tddft:field_Emax',mt_dbl,1,field_Emax))
     & call errquit(pname//'Write failed to rtdb',0, RTDB_ERR)
      if (.not.rtdb_put(rtdb,'rt_tddft:field_w',mt_dbl,1,field_w))
     & call errquit(pname//'Write failed to rtdb',0, RTDB_ERR)
c
c     == NUCLEAR (include nuclear motion) ==
      nuclear = .false.
      if (.not.rtdb_put(rtdb,'rt_tddft:nuclear',mt_log,1,nuclear))
     & call errquit(pname//'Write failed to rtdb',0, RTDB_ERR)
c
c     == Parse the input ==
c
 10   if (.not. inp_read()) 
     &  call errquit(pname//'Read failed input',0, INPUT_ERR)
      if (.not. inp_a(test)) 
     &  call errquit(pname//'Read failed keyword',0, INPUT_ERR)
c
c     == TMAX ==
      if (inp_compare(.false.,test,'tmax')) then
        if (.not.inp_f(tmax)) then
          write(LuOut,*) pname//'tmax value not found; ',
     &      'default value of 100 au will be used'
          tmax=100.d0
        endif
        if (.not.rtdb_put(rtdb,'rt_tddft:tmax',mt_dbl,1,tmax))
     &   call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
c
c     == DELTA ==
      else if (inp_compare(.false.,test,'dt')) then
        if (.not.inp_f(dt)) then
          write(LuOut,*) pname//'dt value not found; ',
     &      'default value of 0.02 atomic units will be used'
          dt=0.02
        endif
        if (.not.rtdb_put(rtdb,'rt_tddft:dt',mt_dbl,1,dt))
     &   call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
c
c     == SODFT ==
      else if (inp_compare(.false.,test,'sodft')) then
        sodft=.true.
        if (.not.rtdb_put(rtdb,'rt_tddft:sodft',mt_log,1,sodft))
     &   call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
c
c     == DFT ==
      else if (inp_compare(.false.,test,'dft')) then
        dft=.true.
        if (.not.rtdb_put(rtdb,'rt_tddft:dft',mt_log,1,dft))
     &   call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
c
c     == RK4 ==
      else if (inp_compare(.false.,test,'rk4')) then
        rk4=.true.
        if (.not.rtdb_put(rtdb,'rt_tddft:rk4',mt_log,1,rk4))
     &   call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
c
c     == EULER ==
      else if (inp_compare(.false.,test,'euler')) then
        euler=.true.
        if (.not.rtdb_put(rtdb,'rt_tddft:euler',mt_log,1,euler))
     &   call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
c
c     == MAGNUS ==
      else if (inp_compare(.false.,test,'magnus')) then
        magnus=.true.
        if (.not.rtdb_put(rtdb,'rt_tddft:magnus',mt_log,1,magnus))
     &   call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
c
c     == ELECTRIC FIELD DESCRIPTION (maximum value) ==
      else if (inp_compare(.false.,test,'field_Emax')) then
      if (.not.inp_f(field_Emax)) then
          write(LuOut,*) pname//'field_Emax value not found; ',
     &      'default value of 0.05 atomic units will be used'
          field_Emax=0.05
      endif
      if (.not.rtdb_put(rtdb,'rt_tddft:field_Emax',mt_dbl,1,field_Emax))
     & call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
c
c     == ELECTRIC FIELD DESCRIPTION (omega) ==
      else if (inp_compare(.false.,test,'field_w')) then
      if (.not.inp_f(field_w)) then
          write(LuOut,*) pname//'field_w value not found; ',
     &      'default value of 0.5 atomic units will be used'
          field_w=0.5
      endif
      if (.not.rtdb_put(rtdb,'rt_tddft:field_Emx',mt_dbl,1,field_w))
     & call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
c
c     == NUCLEAR ==
      else if (inp_compare(.false.,test,'nuclear')) then
        nuclear=.true.
        if (.not.rtdb_put(rtdb,'rt_tddft:magnus',mt_log,1,nuclear))
     &   call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
c
c     == PRINT ==
      else if (inp_compare(.false.,test,'print')) then
        call util_print_input(rtdb,'rt_tddft')
c
c     == END ==
      else if (inp_compare(.false.,test,'end')) then
        goto 20
      else
        call errquit(pname//'Unknown directive',0, INPUT_ERR)
      endif
      goto 10
c
c     == RETURN ==
c
 20   return
      end
