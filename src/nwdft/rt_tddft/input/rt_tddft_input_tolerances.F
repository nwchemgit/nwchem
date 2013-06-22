C
C     rt_tddft_input_tolerances.F
C
      subroutine rt_tddft_input_tolerances (rtdb)
      implicit none

#include "rt_tddft.fh"
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"

      
C     == Inputs ==
      integer, intent(in) :: rtdb


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_input_tolerances: "


C     == Variables ==
      logical got_opt
      character*255 opt
      double precision val
      


C
C     Scan over all words found after "tolerances" directive, make sure
C     it is a valid tolerance keyword, then read in the value and store
C     in rtdb.
C
      got_opt = .true.
      
      do while (got_opt)
         got_opt = inp_a (opt)
         
         if ( trim(opt) .eq. "") exit

         if (opt .eq. "zero") then
            if (.not.inp_f(val)) call errquit (pname//
     $           "zero tolerance should be a positive float", 0, 0)
            if (.not.rtdb_put(rtdb,'rt_tddft:tol_zero',mt_dbl,1,val))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

         elseif (opt .eq. "series") then
            if (.not.inp_f(val)) call errquit (pname//
     $           "series tolerance should be a positive float", 0, 0)
            if (.not.rtdb_put(rtdb,'rt_tddft:tol_series',mt_dbl,1,val))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

         elseif (opt .eq. "interpol") then
            if (.not.inp_f(val)) call errquit (pname//
     $           "interpol tolerance should be a positive float", 0, 0)
            if (.not.rtdb_put(rtdb,'rt_tddft:tol_interpol',
     $           mt_dbl,1,val))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

         else
            call errquit(pname//'unknown tolerance keyword: '//
     $           trim(opt),0,RTDB_ERR)

         endif
      enddo
      
      end subroutine

      
c $Id$
