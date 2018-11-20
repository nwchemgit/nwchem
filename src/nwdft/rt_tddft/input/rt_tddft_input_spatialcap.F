C
C     Input sub-block parser for spatial complex absorbing potential
C     (CAP).
C
C     TODO:
C     - print CAP(x,y,z)
C     - light polarization angle
C     - test polyatomic 
C      
      subroutine rt_tddft_input_spatialcap (rtdb)
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
      character(*), parameter :: pname = "rt_tddft_input_spatialcap: "

C     == Variables ==
      logical done
      logical rangedone
      character*255 test, rangetest, units, atomname, tag
      double precision potmax  ! max of CAP potential in Hartree
      logical lprint, lchecks, lactive
      logical lhave_potmax, lhave_qmin, lhave_qmax, lhave_nq
      logical lhave_units
      double precision :: qmin(3), qmax(3)
      integer :: nq(3), nfound
      double precision :: posconv
      double precision :: rstart, rend
      

C     (default values)
      lactive = .true.
C      units = "au"    ! 0 = au, 1 = angstroms, default: user specfies geom in au
C      posconv = 1d0

      lhave_potmax = .false.
      lhave_qmin = .false.
      lhave_qmax = .false.
      lhave_nq = .false.
      lhave_units = .false.

      
      done = .false.
      do while (.not. done)
         
         if (.not. inp_read()) 
     $        call errquit(pname//'Read failed input',0, INPUT_ERR)
         if (.not. inp_a(test)) 
     $        call errquit(pname//'Read failed keyword',0, INPUT_ERR)
         
         if (inp_compare(.false.,test,'potmax')) then
            if (.not.inp_f(potmax)) call errquit (pname//
     $           "potmax takes a float", 0, INPUT_ERR)
            lhave_potmax = .true.

         elseif (inp_compare(.false.,test,'nq')) then
            if (.not. inp_ilist (3, nq, nfound))
     $           call errquit (pname//"nq takes three ints",
     $           0, INPUT_ERR)
            lhave_nq = .true.

         elseif (inp_compare(.false.,test,'qmin')) then
            if (.not. lhave_units) call errquit (pname//
     $           "units must be specified before qmin",0,INPUT_ERR)

            if (.not. inp_f (qmin(1))) call errquit (pname//
     $           "qmin takes three floats",0,0)
            if (.not. inp_f (qmin(2))) call errquit (pname//
     $           "qmin takes three floats",0,0)
            if (.not. inp_f (qmin(3))) call errquit (pname//
     $           "qmin takes three floats",0,0)
            lhave_qmin = .true.

         elseif (inp_compare(.false.,test,'qmax')) then
            if (.not. lhave_units) call errquit (pname//
     $           "units must be specified before qmax",0,INPUT_ERR)
            
            if (.not. inp_f (qmax(1))) call errquit (pname//
     $           "qmax takes three floats",0,0)
            if (.not. inp_f (qmax(2))) call errquit (pname//
     $           "qmax takes three floats",0,0)
            if (.not. inp_f (qmax(3))) call errquit (pname//
     $           "qmax takes three floats",0,0)
            lhave_qmax = .true.
            
            
         elseif (inp_compare(.false.,test,'on')) then
            lactive = .true.

         elseif (inp_compare(.false.,test,'units')) then
            if (.not. inp_a (units))
     $           call errquit (pname//"failed to read unit type",
     $           0, INPUT_ERR)

C     input is in au, internal in au, no need to convert            
            if (inp_compare(.false.,units,'au')) then
               posconv = 1d0
               lhave_units = .true.
            elseif (inp_compare(.false.,units,'bohr')) then
               posconv = 1d0
               lhave_units = .true.               
            elseif (inp_compare(.false.,units,'bohrs')) then
               posconv = 1d0
               lhave_units = .true.               

C     input is in Angstroms, internal in au, need to convert
            elseif (inp_compare(.false.,units,'angstroms')) then
               posconv = ang2au
               lhave_units = .true.               
            elseif (inp_compare(.false.,units,'angstrom')) then
               posconv = ang2au
               lhave_units = .true.               
            elseif (inp_compare(.false.,units,'ang')) then
               posconv = ang2au
               lhave_units = .true.

            else
               call errquit (pname//"invalid units: "//trim(units),
     $              0, INPUT_ERR)
            endif
            
         elseif (inp_compare(.false.,test,'off')) then
            lactive = .false.


C     read in CAP range, i.e., atom-centered Rstart and Rend
         elseif (inp_compare(.false.,test,"range")) then
            if (.not. lhave_units) call errquit (pname//
     $           "units must be specified before range sub-block",
     $           0, INPUT_ERR)
            
            rangedone = .false.
            do while (.not. rangedone)

               if (.not. inp_read()) 
     $              call errquit(pname//'Read failed input',
     $              0, INPUT_ERR)
               if (.not. inp_a(rangetest)) 
     $              call errquit(pname//'Read failed keyword',
     $              0, INPUT_ERR)

               if (inp_compare(.false.,rangetest,'end')) then
                  rangedone = .true.
               else
                  atomname = rangetest
                  if (.not. inp_f (rstart)) call errquit (pname//
     $                 "range line should be: atom, Rstart, Rend",0,0)
                  if (.not. inp_f (rend)) call errquit (pname//
     $                 "range line should be: atom, Rstart, Rend",0,0)

                  rstart = rstart * posconv
                  rend = rend * posconv

                  tag = "rt_tddft:cap_rstart_"//trim(atomname) ! e.g. "rt_tddft:cap_rstart_O"
                  if (.not. rtdb_put (rtdb, trim(tag),mt_dbl,1,rstart))
     $                 call errquit (pname//
     $                 "failed to write rstart "//trim(tag), 0,RTDB_ERR)

                  tag = "rt_tddft:cap_rend_"//trim(atomname) ! e.g. "rt_tddft:cap_rend_O"
                  if (.not. rtdb_put (rtdb, trim(tag),mt_dbl,1,rend))
     $                 call errquit (pname//
     $                 "failed to write rend "//trim(tag), 0,RTDB_ERR)
               endif
               
            enddo  ! done with "range" sub-block



C     (end of parse)
         elseif (inp_compare(.false.,test,'end')) then
            done = .true.
         else
            call errquit(pname//'Unknown directive: '//trim(test),
     $           0, INPUT_ERR)
         endif
      enddo


C
C     Check that we have all required parameters
C      
         if (.not. lhave_potmax) call errquit(pname//
     $     "must specify potential max for CAP", 0, INPUT_ERR)

         if (potmax .lt. 0d0)
     $        call errquit(pname//
     $        "by potmax should be positive for CAP", 0, INPUT_ERR)

         if (.not. lhave_nq) call errquit(pname//
     $     "must specify nq for CAP", 0, INPUT_ERR)

         if (.not. lhave_qmin) call errquit(pname//
     $     "must specify qmin for CAP", 0, INPUT_ERR)

         if (.not. lhave_qmax) call errquit(pname//
     $     "must specify qmax for CAP", 0, INPUT_ERR)


C
C     Convert all units to au before saving to rtdb
C         
         qmin(1) = qmin(1)*posconv
         qmin(2) = qmin(2)*posconv
         qmin(3) = qmin(3)*posconv
         qmax(1) = qmax(1)*posconv
         qmax(2) = qmax(2)*posconv
         qmax(3) = qmax(3)*posconv
         

C
C     Save in rtdb
C
      if (.not.rtdb_put(rtdb,'rt_tddft:cap_potmax',mt_dbl,1,potmax))
     $     call errquit(pname//
     $     'Failed to write cap_potmax to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_tddft:cap_active',mt_log,1,lactive))
     $     call errquit(pname//
     $     'Failed to write cap_active to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_tddft:cap_nq',mt_int,3,nq))
     $     call errquit(pname//
     $     'Failed to write cap_nq to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_tddft:cap_qmin',mt_dbl,3,qmin))
     $     call errquit(pname//
     $     'Failed to write cap_qmin to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_tddft:cap_qmax',mt_dbl,3,qmax))
     $     call errquit(pname//
     $     'Failed to write cap_qmax to rtdb',0,RTDB_ERR)
      
      end subroutine
