      subroutine tropt_input(rtdb)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "nwc_const.fh"
#include "ctropt.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
c
      integer rtdb
      character*255 field
      double precision ascale, bscale, tscale, hscale
      logical ignore
      integer mh
c
c     So that task knows that tropt is in control
c
      if (.not. rtdb_put(rtdb,'opt:tropt',mt_log,1,.true.))
     $     call errquit('tropt_input: rtdb_put failed',0,
     &       RTDB_ERR)
c
c     Turn driver off
c
      if (.not.rtdb_put(rtdb, 'opt:driver', mt_log, 1, .false.))
     &    call errquit
     &    ('tropt_input: unable to turn opt:driver off',911, RTDB_ERR)
c
c     ----- read values from input -----
c
 10   if(.not.inp_read())
     $     call errquit('tropt_input: unexpected eof',911, INPUT_ERR)
      if (.not. inp_a(field))
     $     call errquit('tropt_input: failed to read field',911,
     &       INPUT_ERR)
c
      if(inp_compare(.false.,'end',field)) then
         return
c
      else if(inp_compare(.false.,'cvgopt',field) .or. 
     $        inp_compare(.false.,'gmax',field) .or.
     $        inp_compare(.false.,'thresh',field)) then
         if(.not. inp_f(gmax_tol))
     $        call errquit('tropt_input: cvgopt/thresh',911,
     &       INPUT_ERR)
         if (.not. rtdb_put(rtdb,'tropt:gmax_tol',mt_dbl,1,gmax_tol))
     $        call errquit('tropt_input: rtdb put failed',0,
     &       RTDB_ERR)
c
      else if(inp_compare(.false.,'grms',field)) then
         if(.not. inp_f(grms_tol))
     $        call errquit('tropt_input: grms',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'tropt:grms_tol',mt_dbl,1,grms_tol))
     $        call errquit('tropt_input: rtdb put failed',0,
     &       RTDB_ERR)
c
      else if(inp_compare(.false.,'xrms',field)) then
         if(.not. inp_f(xrms_tol))
     $        call errquit('tropt_input: xrms',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'tropt:xrms_tol',mt_dbl,1,xrms_tol))
     $        call errquit('tropt_input: rtdb put failed',0,
     &       RTDB_ERR)
c
      else if(inp_compare(.false.,'xmax',field)) then
         if(.not. inp_f(xmax_tol))
     $        call errquit('tropt_input: xmax',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'tropt:xmax_tol',mt_dbl,1,xmax_tol))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
c
      else if (inp_compare(.false.,'mh', field)) then
         if(.not. inp_i(mh))
     $        call errquit('tropt_input: mh',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'tropt:mh',mt_int,1,mh))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
c
      else if (inp_compare(.false.,'lbfgs', field)) then
         if (.not. rtdb_put(rtdb,'tropt:lbfgs',mt_log,1,.true.))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
c
      else if (inp_compare(.false.,'loose', field)) then
         gmax_tol = 0.0045d0    ! default*10
         grms_tol = 0.0030d0    ! default*10
         xmax_tol = 0.0180d0    ! default*10
         xrms_tol = 0.0120d0    ! default*10
         if (.not. rtdb_put(rtdb,'tropt:gmax_tol',mt_dbl,1,gmax_tol))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'tropt:grms_tol',mt_dbl,1,grms_tol))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'tropt:xrms_tol',mt_dbl,1,xrms_tol))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'tropt:xmax_tol',mt_dbl,1,xmax_tol))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
c
      else if (inp_compare(.false.,'default', field)) then
         gmax_tol = 0.00045d0
         grms_tol = 0.00030d0
         xmax_tol = 0.00180d0
         xrms_tol = 0.00120d0
         if (.not. rtdb_put(rtdb,'tropt:gmax_tol',mt_dbl,1,gmax_tol))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'tropt:grms_tol',mt_dbl,1,grms_tol))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'tropt:xrms_tol',mt_dbl,1,xrms_tol))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'tropt:xmax_tol',mt_dbl,1,xmax_tol))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
c
      else if (inp_compare(.false.,'tight', field)) then
         gmax_tol = 0.000015d0  ! default / 30
         grms_tol = 0.00001d0   ! default / 30
         xmax_tol = 0.00006d0   ! default / 30
         xrms_tol = 0.00004d0   ! default / 30
         if (.not. rtdb_put(rtdb,'tropt:gmax_tol',mt_dbl,1,gmax_tol))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'tropt:grms_tol',mt_dbl,1,grms_tol))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'tropt:xrms_tol',mt_dbl,1,xrms_tol))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'tropt:xmax_tol',mt_dbl,1,xmax_tol))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
CJMC
      else if(inp_compare(.false.,'opttol',field)) then
         if(.not. inp_f(opt_tol))
     $        call errquit('tropt_input: opttol',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'tropt:opt_tol',mt_dbl,1,opt_tol))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
         grms_tol = 1.0d0*opt_tol
         gmax_tol = 1.5d0*opt_tol
         xrms_tol = 4.0d0*opt_tol
         xmax_tol = 6.0d0*opt_tol
         if (.not. rtdb_put(rtdb,'tropt:gmax_tol',mt_dbl,1,gmax_tol))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'tropt:grms_tol',mt_dbl,1,grms_tol))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'tropt:xrms_tol',mt_dbl,1,xrms_tol))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'tropt:xmax_tol',mt_dbl,1,xmax_tol))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
CJMC
c
      else if(inp_compare(.false.,'eprec',field)) then
         if(.not. inp_f(eprec))
     $        call errquit('tropt_input: eprec',911, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'tropt:eprec',mt_dbl,1,eprec))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'trust',field)) then
         if(.not. inp_f(trust)) 
     $        call errquit('tropt_input: trust',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'tropt:trust',mt_dbl,1,trust))
     $        call errquit('tropt_input: rtdb put failed',0, INPUT_ERR)
CJMC
         if (.not. rtdb_put(rtdb,'tropt:maxdr',mt_dbl,1,trust))
     $        call errquit('tropt_input: rtdb put failed',0, INPUT_ERR)
CJMC
c
      else if(inp_compare(.false.,'sadstp',field)) then
         if(.not. inp_f(sadstp)) 
     $        call errquit('tropt_input: sadstp',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'tropt:sadstp',mt_dbl,1,sadstp))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'nptopt',field) .or. 
     $        inp_compare(.false.,'maxiter',field)) then
         if(.not. inp_i(nptopt))
     $        call errquit('tropt_input: nptopt/maxiter',911,
     &       INPUT_ERR)
         if (.not. rtdb_put(rtdb,'tropt:nptopt',mt_int,1,nptopt))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'linopt',field)) then
         if(.not. inp_i(linopt))
     $        call errquit('tropt_input: linopt',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'tropt:linopt',mt_int,1,linopt))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
c
      else if (inp_compare(.false.,'clear',field)) then
         call tropt_del_hess()
         ignore = rtdb_delete(rtdb,'tropt:moddir')
c
      else if (inp_compare(.false.,'redoautoz',field)) then
         if (.not. rtdb_put(rtdb,'tropt:redoautoz',mt_log,1,.true.))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'inhess',field)) then
         if(.not. inp_i(inhess)) 
     $        call errquit('tropt_input: inhess',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'tropt:inhess',mt_int,1,inhess))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'modupd',field)) then
         if(.not. inp_i(modupd))
     $        call errquit('tropt_input: modupd',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'tropt:modupd',mt_int,1,modupd))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'modsad',field)) then
         if(.not. inp_i(modsad))
     $        call errquit('tropt_input: modsad',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'tropt:modsad',mt_int,1,modsad))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'moddir',field)) then
         if(.not. inp_i(moddir)) 
     $        call errquit('tropt_input: moddir',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'tropt:moddir',mt_int,1,moddir))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'print',field)) then
         call util_print_input(rtdb,'tropt')
c
      else if(inp_compare(.false.,'noprint',field)) then
         call util_print_input(rtdb,'tropt')
c
      else if (inp_compare(.false.,'qstep',field)) then
         oqstep = .true.
         if (.not. rtdb_put(rtdb,'tropt:oqstep',mt_log,1,oqstep))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
c
      else if (inp_compare(.false.,'noqstep',field)) then
         oqstep = .false.
         if (.not. rtdb_put(rtdb,'tropt:oqstep',mt_log,1,oqstep))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'hscale',field)) then
         if(.not. inp_f(hscale)) 
     $        call errquit('tropt_input: hscale',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'tropt:hscale',mt_dbl,1,hscale))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'ascale',field)) then
         if(.not. inp_f(ascale)) 
     $        call errquit('tropt_input: ascale',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'tropt:ascale',mt_dbl,1,ascale))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'bscale',field)) then
         if(.not. inp_f(bscale)) 
     $        call errquit('tropt_input: bscale',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'tropt:bscale',mt_dbl,1,bscale))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'tscale',field)) then
         if(.not. inp_f(tscale)) 
     $        call errquit('tropt_input: tscale',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'tropt:tscale',mt_dbl,1,tscale))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'xyz',field)) then
         if (.not. inp_a(xyz)) then
            call util_file_prefix('x',xyz)
            xyz(inp_strlen(xyz)-1:) = ' '
         endif
         if (.not. rtdb_cput(rtdb,'tropt:xyz',1,xyz))
     $        call errquit('tropt_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'noxyz',field)) then
         ignore = rtdb_delete(rtdb,'tropt:xyz')
c
      else
         call errquit('tropt_input: unknown directive',0, INPUT_ERR)
      endif
c
      goto 10                   ! While loop
c
      end
