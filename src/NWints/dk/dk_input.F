      subroutine dk_input(rtdb)
*
* $Id: dk_input.F,v 1.1 1999-11-11 21:30:03 bert Exp $
*
      implicit none
#include "rtdb.fh"
#include "context.fh"
#include "mafdecls.fh"
#include "inp.fh"
c     
c     Read input for calculations with Douglas-Kroll
c     
c     All input is of the form:
c     <keyword>  or   <keyword> <on>  or   <keyword> <off>
c     
      integer rtdb ! [input] handle to database
      
      character*255 test
      logical dk_run, dk_status, status
      integer dk_type
c
c     Default settings
c
      dk_status = .false.
      dk_run    = .true.
      dk_type   = 2
c
c     Check if there has been a calculation before. If so, check if it 
c     was a Douglas-Kroll calculation
c
c     status = rtdb_get(rtdb, 'doug_kroll', MT_LOG, 1,dk_status)
c
c     Read input data beyond douglas-kroll; store in rtdb
c     
  10  if (inp_a(test)) then
         if (inp_compare(.false., 'off', test)) then
            dk_run  = .false.
            dk_type = 1
         endif
         goto 10
      endif
c
c     If the type of calculation has changed compared to previous one
c     the calculation cannot be converged
c
c     if (.not.(dk_status.and.dk_run))
c    $   status = rtdb_delete(rtdb,'scf:converged')
c
c     Put parameters in rtdb
c       
      if (.not. rtdb_put(rtdb, 'doug_kroll', mt_log, 1, dk_run))
     $   call errquit('dk_input: rtdb_put failed', 0)
      if (.not. rtdb_put(rtdb, 'doug_kroll:type', mt_int, 1, dk_type))
     $   call errquit('dk_input: rtdb_put failed', 0)
        return
c       
        end
