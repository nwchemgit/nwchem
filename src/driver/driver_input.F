      subroutine driver_input(rtdb)
      implicit none
#include "inp.fh"
#include "coptopt.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
c
      integer rtdb
      character*255 field
c
c     So that task knows that driver is in control
c
      if (.not. rtdb_put(rtdb,'opt:driver',mt_log,1,.true.))
     $     call errquit('driver_input: rtdb_put failed',0)
c
c     ----- read values from input -----
c
 10   if(.not.inp_read())
     $     call errquit('driver_input: unexpected eof',911)
      if (.not. inp_a(field))
     $     call errquit('driver_input: failed to read field',911)
c
      if(inp_compare(.false.,'end',field)) then
         return
c
      else if(inp_compare(.false.,'cvgopt',field) .or. 
     $        inp_compare(.false.,'thresh',field)) then
         if(.not. inp_f(cvgopt))
     $        call errquit('driver_input: cvgopt/thresh',911)
         if (.not. rtdb_put(rtdb,'driver:cvgopt',mt_dbl,1,cvgopt))
     $        call errquit('driver_input: rtdb put failed',0)
c
      else if(inp_compare(.false.,'sadstp',field)) then
         if(.not. inp_f(sadstp)) 
     $        call errquit('driver_input: sadstp',911)
         if (.not. rtdb_put(rtdb,'driver:sadstp',mt_dbl,1,sadstp))
     $        call errquit('driver_input: rtdb put failed',0)
c
      else if(inp_compare(.false.,'trust',field)) then
         if(.not. inp_f(sadstp)) 
     $        call errquit('driver_input: trust',911)
         if (.not. rtdb_put(rtdb,'driver:trust',mt_dbl,1,trust))
     $        call errquit('driver_input: rtdb put failed',0)
c
      else if(inp_compare(.false.,'nptopt',field) .or. 
     $        inp_compare(.false.,'maxiter',field)) then
         if(.not. inp_i(nptopt))
     $        call errquit('driver_input: nptopt/maxiter',911)
         if (.not. rtdb_put(rtdb,'driver:nptopt',mt_int,1,nptopt))
     $        call errquit('driver_input: rtdb put failed',0)
c
      else if(inp_compare(.false.,'linopt',field)) then
         if(.not. inp_i(linopt))
     $        call errquit('driver_input: linopt',911)
         if (.not. rtdb_put(rtdb,'driver:linopt',mt_int,1,linopt))
     $        call errquit('driver_input: rtdb put failed',0)
c
      else if (inp_compare(.false.,'clear',field)) then
         call driver_opt_del_hess()
c
      else if (inp_compare(.false.,'redoautoz',field)) then
         if (.not. rtdb_put(rtdb,'driver:redoautoz',mt_log,1,.true.))
     $        call errquit('driver_input: rtdb put failed',0)
c
      else if(inp_compare(.false.,'inhess',field)) then
         if(.not. inp_i(inhess)) 
     $        call errquit('driver_input: inhess',911)
         if (.not. rtdb_put(rtdb,'driver:inhess',mt_int,1,inhess))
     $        call errquit('driver_input: rtdb put failed',0)
c
      else if(inp_compare(.false.,'modupd',field)) then
         if(.not. inp_i(modupd))
     $        call errquit('driver_input: modupd',911)
         if (.not. rtdb_put(rtdb,'driver:modupd',mt_int,1,modupd))
     $        call errquit('driver_input: rtdb put failed',0)
c
      else if(inp_compare(.false.,'modsad',field)) then
         if(.not. inp_i(modsad))
     $        call errquit('driver_input: modsad',911)
         if (.not. rtdb_put(rtdb,'driver:modsad',mt_int,1,modsad))
     $        call errquit('driver_input: rtdb put failed',0)
c
      else if(inp_compare(.false.,'moddir',field)) then
         if(.not. inp_i(moddir)) 
     $        call errquit('driver_input: moddir',911)
         if (.not. rtdb_put(rtdb,'driver:moddir',mt_int,1,moddir))
     $        call errquit('driver_input: rtdb put failed',0)
c
      else if(inp_compare(.false.,'print',field)) then
         call util_print_input(rtdb,'driver')
c
      else if(inp_compare(.false.,'noprint',field)) then
         call util_print_input(rtdb,'driver')
c
      else if (inp_compare(.false.,'qstep',field)) then
         oqstep = .true.
         if (.not. rtdb_put(rtdb,'driver:oqstep',mt_log,1,oqstep))
     $        call errquit('driver_input: rtdb put failed',0)
c
      else if (inp_compare(.false.,'noqstep',field)) then
         oqstep = .false.
         if (.not. rtdb_put(rtdb,'driver:oqstep',mt_log,1,oqstep))
     $        call errquit('driver_input: rtdb put failed',0)
      endif
c
      goto 10                   ! While loop
c
      end
