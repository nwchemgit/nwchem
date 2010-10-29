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
      logical function socrates_driver(rtdb)
c
c $Id$
c
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "bas.fh"
#include "geom.fh"
#include "sym.fh"
#include "schwarz.fh"
c
c     object handles
c
      integer rtdb             ! RTDB handle
      integer geom             ! GEOM handle
c
c     variables
c
      integer scf_alg
      integer nuclear_charge
      integer total_charge
      integer nopen
      integer num_elec
      integer num_alfa
      integer num_beta
c
      double precision dtemp
c
      logical debug
      logical nodezero
      logical stat
      logical use_symm
c
      character*8 group
c
c     function declarations
c
      logical socrates_scf_serial
      external socrates_scf_serial
      logical socrates_scf_ga
      external socrates_scf_ga
c
#ifdef DEBUG_PRINT
      debug = (GA_NodeId().eq.0) ! debug print on nodezero only
c      debug = .true. ! debug print everywhere
#else
      parameter debug = .false.
#endif
c
c===============================================================c
c                                                               c
c                        INITIALIZATION                         c
c                                                               c
c===============================================================c
c
      nodezero=(ga_nodeid().eq.0)
c
      if (debug) write(LuOut,*) 'top of socrates_driver'
c
      socrates_driver = .false.
c
      if (.not.rtdb_get(rtdb,'socrates:scf_alg',mt_int,1,scf_alg)) then
        scf_alg = 2
      endif
c
c===============================================================c
c                                                               c
c     Setup molecular parameters                                c
c                                                               c
c===============================================================c
c
c     ---------------------
c     Geometry and symmetry
c     ---------------------
c
      if (.not.geom_create(geom,'geometry')) then
        call errquit(__FILE__,__LINE__,GEOM_ERR)
      endif
      if (.not.geom_rtdb_load(rtdb,geom,'geometry')) then
        call errquit(__FILE__,__LINE__,GEOM_ERR)
      endif
c
      call sym_group_name(geom,group)
      use_symm = (group(1:2) .ne. 'C1')
c
      if (use_symm) then
#ifndef SYMMETRY_ENABLED
        write(LuOut,*) 'Please implement symmetry'
        call errquit(__FILE__,__LINE__,GEOM_ERR)
#else
        if(nodezero) then
          write(LuOut,*) ' BE CAREFUL!  Using symmetry: ',group
        endif
#endif
        call util_flush(LuOut)
      endif
c
c     ------
c     Charge
c     ------
c
c     total charge
c
      if (.not.rtdb_get(rtdb, 'charge',mt_dbl,1,dtemp)) then
        dtemp = 0.0d0
      endif
      total_charge = int(dtemp)
      if ((dble(total_charge)-dtemp).gt.1e-7) then
        call errquit(__FILE__,__LINE__,INPUT_ERR)
      endif
c
c     nuclear charge
c
      if (.not. geom_nuc_charge(geom, dtemp)) then
        call errquit(__FILE__,__LINE__,GEOM_ERR)
      endif
      nuclear_charge = int(dtemp)
      if ((dble(nuclear_charge)-dtemp).gt.1e-7) then
        call errquit(__FILE__,__LINE__,INPUT_ERR)
      endif
c
c     electron count
c
      num_elec = nuclear_charge - total_charge
      if (num_elec.le.0) then
        if (nodezero) then
          write(LuOut,*) 'Electrons required!'
        endif
        call util_flush(LuOut)
        call errquit(__FILE__,__LINE__,INPUT_ERR)
      endif
c
c     open-shell
c
      if (.not.rtdb_get(rtdb,'socrates:nopen',mt_int,1,nopen)) then
        nopen = 0
      endif
      if (mod(num_elec-nopen,2).eq.0) then
        num_beta = (num_elec-nopen)/2
        num_alfa = num_elec-num_beta
      else
        if (nodezero) write(LuOut,*) 'nopen is wrong'
        call util_flush(LuOut)
        call errquit(__FILE__,__LINE__,INPUT_ERR)
        num_beta = 0 ! I hate compiler warnings 
        num_alfa = 0 ! I hate compiler warnings
      endif
      if (.not.rtdb_put(rtdb,'socrates:num_alfa',mt_int,1,
     1                  num_alfa)) then
           call errquit(__FILE__,__LINE__,RTDB_ERR)
      endif
      if (.not.rtdb_put(rtdb,'socrates:num_beta',mt_int,1,
     1                  num_beta)) then
           call errquit(__FILE__,__LINE__,RTDB_ERR)
      endif
c
c===============================================================c
c                                                               c
c                           RUN SCF                             c
c                                                               c
c===============================================================c
c
      if (scf_alg.eq.1) then
        if (nodezero) write(LuOut,102)
  102   format(/,'SCF with GA algorithm',/)
        stat = socrates_scf_ga(rtdb,geom,num_alfa,num_beta)
        if (.not.stat) then
          write(LuOut,*) 'socrates_scf_ga failed'
          call util_flush(LuOut)
          call errquit(__FILE__,__LINE__,GEOM_ERR)
        endif
      else
        write(LuOut,*) 'invalid scf_alg'
        call util_flush(LuOut)
        call errquit(__FILE__,__LINE__,GEOM_ERR)
      endif
c
c===============================================================c
c                                                               c
c                         TERMINATION                           c
c                                                               c
c===============================================================c
c
      if (.not.geom_destroy(geom)) then
        call errquit(__FILE__,__LINE__,GEOM_ERR)
      endif
c
      socrates_driver = .true.
c
      if (debug) then
        write(LuOut,*) 'end of socrates_driver'
        call util_flush(LuOut)
      endif
c
      return
c
      end
