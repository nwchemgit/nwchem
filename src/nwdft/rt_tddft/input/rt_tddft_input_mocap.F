C
C     Input sub-block parser for MO absorbing boundary condition.
C

      subroutine rt_tddft_input_mocap (rtdb)
      implicit none

#include "rt_tddft.fh"
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"

      
C     == Inputs ==
      integer, intent(in)      :: rtdb

C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_input_mocap: "

C     == Variables ==
      logical done
      character*255 test
      double precision emin, expconst, maxval, prefac
      logical lprint, lchecks, lactive


C     (default values)
      maxval = 100d0    ! exponential max value (clamped to this)
      emin = 0d0        ! vacuum energy level
      prefac = 1d0      ! exp prefactor
      expconst = 1d0    ! exp constant
      lchecks = .true.  ! print potential at each step?  This is slow.
      lprint = .false.  ! check projectors at each step?
      lactive = .false.  ! check projectors at each step?
      

C     XXX FIX LOGICALS:  read in true/false
C     XXX CHANGE NAME: MOCAP

      done = .false.
      do while (.not. done)
         
         if (.not. inp_read()) 
     $        call errquit(pname//'Read failed input',0, INPUT_ERR)
         if (.not. inp_a(test)) 
     $        call errquit(pname//'Read failed keyword',0, INPUT_ERR)
         
         if (inp_compare(.false.,test,'emin')) then
            if (.not.inp_f(emin)) call errquit (pname//
     $           "emin takes a float", 0, 0)

         elseif (inp_compare(.false.,test,'maxval')) then
            if (.not.inp_f(maxval)) call errquit (pname//
     $           "maxval takes a float", 0, 0)

         elseif (inp_compare(.false.,test,'prefac')) then
            if (.not.inp_f(prefac)) call errquit (pname//
     $           "prefac takes a float", 0, 0)
            
         elseif (inp_compare(.false.,test,'expconst')) then
            if (.not.inp_f(expconst)) call errquit (pname//
     $           "expconst takes a float", 0, 0)

         elseif (inp_compare(.false.,test,'print')) then
            lprint = .true.

         elseif (inp_compare(.false.,test,'noprint')) then
            lprint = .false.

         elseif (inp_compare(.false.,test,'checks')) then
            lchecks = .true.

         elseif (inp_compare(.false.,test,'nochecks')) then
            lchecks = .false.

         elseif (inp_compare(.false.,test,'on')) then
            lactive = .true.

         elseif (inp_compare(.false.,test,'off')) then
            lactive = .false.
            
            
C     (end of parse)
         elseif (inp_compare(.false.,test,'end')) then
            done = .true.
         else
            call errquit(pname//'Unknown directive: '//trim(test),
     $           0, INPUT_ERR)
         endif
      enddo


C
C     check that we have required parameters
C
C     not needed, as we have defaults for all params

C
C     save in rtdb
C
      if (.not.rtdb_put(rtdb,'rt_tddft:mocap_maxval',mt_dbl,1,maxval))
     $     call errquit(pname//
     $     'Failed to write tstart to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_tddft:mocap_emin',mt_dbl,1,emin))
     $     call errquit(pname//
     $     'Failed to write tstart to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_tddft:mocap_prefac',mt_dbl,1,prefac))
     $     call errquit(pname//
     $     'Failed to write tstart to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_tddft:mocap_expconst',mt_dbl,
     $     1,expconst))
     $     call errquit(pname//
     $     'Failed to write tstart to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_tddft:mocap_print',mt_log,1,lprint))
     $     call errquit(pname//
     $     'Failed to write dplot to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_tddft:mocap_checks',mt_log,1,lchecks))
     $     call errquit(pname//
     $     'Failed to write dplot to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_tddft:mocap_active',mt_log,1,lactive))
     $     call errquit(pname//
     $     'Failed to write dplot to rtdb',0,RTDB_ERR)
         
      end subroutine
