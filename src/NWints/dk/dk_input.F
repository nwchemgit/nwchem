      subroutine dk_input(rtdb)
*
* $Id: dk_input.F,v 1.2 1999-11-19 19:57:29 bert Exp $
*
      implicit none
#include "rtdb.fh"
#include "context.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"
c     
c     Read input for calculations with Douglas-Kroll
c     
c     All input is of the form:
c     <keyword>  or   <keyword> <on>  or   <keyword> <off>
c     
      integer rtdb ! [input] handle to database
      
      character*255 test
      logical dk_run, status
      integer dk_type
c
c     Default settings
c
      dk_run    = .true.
      dk_type   = 2
c
c     Check if there has been a calculation before. If so, check if it 
c     was a Douglas-Kroll calculation
c
      status = rtdb_get(rtdb, 'doug_kroll', MT_LOG, 1,dk_run)
      status = rtdb_get(rtdb, 'doug_kroll:type', MT_LOG, 1,dk_type)
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
c     If we are going to do Douglas-Kroll, write out a line informing user
c
      if ((ga_nodeid() .eq. 0) .and .dk_run) then
         write(luout,20)
  20     format(/12x, 'Douglas-Kroll integrals are switched on !'/)
         if (dk_type .eq. 1) then
            write(luout,21)
  21        format(/12x, 'Integrals for the no-pair equation with ',
     $             'free-particle projection operators'/
     $             12x,'will be computed.'//)
         elseif (dk_type .eq. 2) then
            write(luout,22)
  22        format(/12x, 'Integrals for the no-pair equation with ',
     $             'external-field projection operators'/
     $             12x,'will be computed (Douglas-Kroll-Hess).'//)
         else  
            write(luout,23)
  23        format(/12x, 'Integrals for the no-pair equation with ',
     $             'external-field projection operators'/12x,
     $             'including cross-product terms will be computed.'//)
         endif
         call util_flush(luout)
      endif
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
