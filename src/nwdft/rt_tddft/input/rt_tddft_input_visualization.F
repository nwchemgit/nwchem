      subroutine rt_tddft_input_visualization (rtdb)
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
      character(*), parameter :: pname="rt_tddft_input_visualization: "

C     == Variables ==
      logical done
      character*255 test
      double precision tstart, tend, treference
      logical lhave_tstart, lhave_tend, lhave_treference
      logical lsubgs
      logical ldplot
      

      lhave_tstart = .false.
      lhave_tend = .false.
      lhave_treference = .false.
      ldplot = .false.


      done = .false.
      do while (.not. done)

         if (.not. inp_read()) 
     $        call errquit(pname//'Read failed input',0, INPUT_ERR)
         if (.not. inp_a(test)) 
     $        call errquit(pname//'Read failed keyword',0, INPUT_ERR)


C
C     nsnapshots XXXTODO
C


C
C     tstart and tend
C
         if (inp_compare(.false.,test,'tstart')) then
            if (.not.inp_f(tstart)) call errquit (pname//
     $           "tstart takes a float", 0, 0)
            lhave_tstart = .true.

         elseif (inp_compare(.false.,test,'tend')) then
            if (.not.inp_f(tend)) call errquit (pname//
     $           "tend takes a float", 0, 0)
            lhave_tend = .true.


C
C     time value for reference density matrix (to be subtracted), e.g, 0.0 for ground state
C
         elseif (inp_compare(.false.,test,'treference')) then
            if (.not.inp_f(treference)) call errquit (pname//
     $           "treference takes a float", 0, 0)
            lhave_treference = .true.


C
C     subtract ground state (XXX only used in dplot call, dens mats NOT subtracted)
C            
C     xxx no longer used
C
C         elseif (inp_compare(.false.,test,'subgs')) then
C            lsubgs = .true.


C
C     do dplot post-processing
C            
         elseif (inp_compare(.false.,test,'dplot')) then
            ldplot = .true.


C     
C     end of parse
C     
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
      if (.not. lhave_tstart)
     $     call errquit (pname//"must supply tstart",0,0)

      if (.not. lhave_tend)
     $     call errquit (pname//"must supply tend",0,0)



      if (lhave_treference) then
         lsubgs = .true.
      else
         lsubgs = .false.
      endif
         
      
C
C     save in rtdb
C
      if (.not.rtdb_put(rtdb,'rt_tddft:viz_tstart',mt_dbl,1,tstart))
     $     call errquit(pname//
     $     'Failed to write tstart to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_tddft:viz_tend',mt_dbl,1,tend))
     $     call errquit(pname//
     $     'Failed to write tend to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_tddft:viz_treference',mt_dbl,1,
     $     treference)) call errquit(pname//
     $     'Failed to write tend to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_tddft:viz_subgs',mt_log,1,lsubgs))
     $     call errquit(pname//
     $     'Failed to write subgs to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_tddft:viz_dplot',mt_log,1,ldplot))
     $     call errquit(pname//
     $     'Failed to write dplot to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_tddft:viz_active',mt_log,1,.true.))
     $     call errquit(pname//
     $     'Failed to write viz_active to rtdb',0,RTDB_ERR)

         
      end subroutine
c $Id$
