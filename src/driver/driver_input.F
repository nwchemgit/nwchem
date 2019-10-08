C>
C> \ingroup driver
C> @{
C>
C> \file driver_input.F
C> The Driver input reader
C>
C> @}
C>
C> \ingroup driver
C> @{
C>
      subroutine driver_input(rtdb)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "nwc_const.fh"
#include "coptopt.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
c
      integer rtdb
      character*255 field
      character*50  f2
      character*128 ipname
      double precision ascale, bscale, tscale, hscale
      logical ignore,diagh
      integer mh,ind

*     **** external functions ****
      logical  driver_parse_boolean
      external driver_parse_boolean
c
c     So that task knows that driver is in control
c
      if (.not. rtdb_put(rtdb,'opt:driver',mt_log,1,.true.))
     $     call errquit('driver_input: rtdb_put failed',0,
     &       RTDB_ERR)
c
c     ----- read values from input -----
c
 10   if(.not.inp_read())
     $     call errquit('driver_input: unexpected eof',911, INPUT_ERR)
      if (.not. inp_a(field))
     $     call errquit('driver_input: failed to read field',911,
     &       INPUT_ERR)
c
      if(inp_compare(.false.,'end',field)) then
         return
c
      else if(inp_compare(.false.,'cvgopt',field) .or. 
     $        inp_compare(.false.,'gmax',field) .or.
     $        inp_compare(.false.,'thresh',field)) then
         if(.not. inp_f(gmax_tol))
     $        call errquit('driver_input: cvgopt/thresh',911,
     &       INPUT_ERR)
         if (.not. rtdb_put(rtdb,'driver:gmax_tol',mt_dbl,1,gmax_tol))
     $        call errquit('driver_input: rtdb put failed',0,
     &       RTDB_ERR)
c
      else if(inp_compare(.false.,'grms',field)) then
         if(.not. inp_f(grms_tol))
     $        call errquit('driver_input: grms',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'driver:grms_tol',mt_dbl,1,grms_tol))
     $        call errquit('driver_input: rtdb put failed',0,
     &       RTDB_ERR)
c
      else if(inp_compare(.false.,'xrms',field)) then
         if(.not. inp_f(xrms_tol))
     $        call errquit('driver_input: xrms',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'driver:xrms_tol',mt_dbl,1,xrms_tol))
     $        call errquit('driver_input: rtdb put failed',0,
     &       RTDB_ERR)
c
      else if(inp_compare(.false.,'xmax',field)) then
         if(.not. inp_f(xmax_tol))
     $        call errquit('driver_input: xmax',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'driver:xmax_tol',mt_dbl,1,xmax_tol))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
c
      else if (inp_compare(.false.,'mh', field)) then
         if(.not. inp_i(mh))
     $        call errquit('driver_input: mh',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'driver:mh',mt_int,1,mh))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
c
      else if (inp_compare(.false.,'lbfgs', field)) then
         if (.not. rtdb_put(rtdb,'driver:lbfgs',mt_log,1,.true.))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)

      else if (inp_compare(.false.,'diagonal_hessian', field)) then
         diagh = .true.
         if (inp_a(f2)) diagh = driver_parse_boolean(f2,.true.)
         if (.not. rtdb_put(rtdb,'driver:diagonal_hessian',
     $                      mt_log,1,diagh))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
c
c socket ipi_client ip:port
c
      else if (inp_compare(.false.,'socket', field)) then
         if (inp_a(f2)) then
            if (driver_parse_boolean(f2,.true.)) then
               ind = index(f2,' ') - 1
               if (.not.rtdb_cput(rtdb,'driver:socket_type',1,
     >                      f2(1:ind)))
     >          call errquit('driver_input: rtdb put failed',7,RTDB_ERR)
      
               if (inp_a(ipname)) then
                  ind = index(ipname,' ') - 1
                  if (.not.rtdb_cput(rtdb,'driver:socket_ip',1,
     >                           ipname(1:ind)))
     >          call errquit('driver_input: rtdb put failed',8,RTDB_ERR)
               end if
               diagh = .true.
               if (.not. rtdb_put(rtdb,'driver:socket',
     $                      mt_log,1,diagh))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
c
            else
               if (.not.rtdb_delete(rtdb,'driver:socket_type'))
     >         call errquit('driver_input:rtdbdelete failed',9,RTDB_ERR)
               diagh = .false.
               if (.not. rtdb_put(rtdb,'driver:socket',
     $                      mt_log,1,diagh))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
            end if
         end if
c
      else if (inp_compare(.false.,'deloc', field)) then
         if (.not. rtdb_put(rtdb,'driver:deloc',mt_log,1,.true.))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
c
      else if (inp_compare(.false.,'loose', field)) then
         gmax_tol = 0.0045d0    ! default*10
         grms_tol = 0.0030d0    ! default*10
         xmax_tol = 0.0180d0    ! default*10
         xrms_tol = 0.0120d0    ! default*10
         if (.not. rtdb_put(rtdb,'driver:gmax_tol',mt_dbl,1,gmax_tol))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'driver:grms_tol',mt_dbl,1,grms_tol))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'driver:xrms_tol',mt_dbl,1,xrms_tol))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'driver:xmax_tol',mt_dbl,1,xmax_tol))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
c
      else if (inp_compare(.false.,'default', field)) then
         gmax_tol = 0.00045d0
         grms_tol = 0.00030d0
         xmax_tol = 0.00180d0
         xrms_tol = 0.00120d0
         if (.not. rtdb_put(rtdb,'driver:gmax_tol',mt_dbl,1,gmax_tol))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'driver:grms_tol',mt_dbl,1,grms_tol))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'driver:xrms_tol',mt_dbl,1,xrms_tol))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'driver:xmax_tol',mt_dbl,1,xmax_tol))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
c
      else if (inp_compare(.false.,'tight', field)) then
         gmax_tol = 0.000015d0  ! default / 30
         grms_tol = 0.00001d0   ! default / 30
         xmax_tol = 0.00006d0   ! default / 30
         xrms_tol = 0.00004d0   ! default / 30
         if (.not. rtdb_put(rtdb,'driver:gmax_tol',mt_dbl,1,gmax_tol))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'driver:grms_tol',mt_dbl,1,grms_tol))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'driver:xrms_tol',mt_dbl,1,xrms_tol))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'driver:xmax_tol',mt_dbl,1,xmax_tol))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'eprec',field)) then
         if(.not. inp_f(eprec))
     $        call errquit('driver_input: eprec',911, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'driver:eprec',mt_dbl,1,eprec))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'trust',field)) then
         if(.not. inp_f(trust)) 
     $        call errquit('driver_input: trust',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'driver:trust',mt_dbl,1,trust))
     $        call errquit('driver_input: rtdb put failed',0, INPUT_ERR)
c
      else if(inp_compare(.false.,'sadstp',field)) then
         if(.not. inp_f(sadstp)) 
     $        call errquit('driver_input: sadstp',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'driver:sadstp',mt_dbl,1,sadstp))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'nptopt',field) .or. 
     $        inp_compare(.false.,'maxiter',field)) then
         if(.not. inp_i(nptopt))
     $        call errquit('driver_input: nptopt/maxiter',911,
     &       INPUT_ERR)
         if (.not. rtdb_put(rtdb,'driver:nptopt',mt_int,1,nptopt))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'linopt',field)) then
         if(.not. inp_i(linopt))
     $        call errquit('driver_input: linopt',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'driver:linopt',mt_int,1,linopt))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
c
      else if (inp_compare(.false.,'clear',field)) then
         call driver_opt_del_hess()
         ignore = rtdb_delete(rtdb,'driver:moddir')
c
      else if (inp_compare(.false.,'redoautoz',field)) then
         if (.not. rtdb_put(rtdb,'driver:redoautoz',mt_log,1,.true.))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'inhess',field)) then
         if(.not. inp_i(inhess)) 
     $        call errquit('driver_input: inhess',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'driver:inhess',mt_int,1,inhess))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'modupd',field)) then
         if(.not. inp_i(modupd))
     $        call errquit('driver_input: modupd',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'driver:modupd',mt_int,1,modupd))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'modsad',field)) then
         if(.not. inp_i(modsad))
     $        call errquit('driver_input: modsad',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'driver:modsad',mt_int,1,modsad))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'moddir',field)) then
         if(.not. inp_i(moddir)) 
     $        call errquit('driver_input: moddir',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'driver:moddir',mt_int,1,moddir))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'vardir',field)) then
         if(.not. inp_i(vardir)) 
     $        call errquit('driver_input: vardir',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'driver:vardir',mt_int,1,vardir))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
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
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
c
      else if (inp_compare(.false.,'firstneg',field)) then
         ofirstneg = .true.
         if (.not. rtdb_put(rtdb,'driver:ofirstneg',mt_log,1,ofirstneg))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
c
      else if (inp_compare(.false.,'nofirstneg',field)) then
         ofirstneg = .false.
         if (.not. rtdb_put(rtdb,'driver:ofirstneg',mt_log,1,ofirstneg))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
c
      else if (inp_compare(.false.,'noqstep',field)) then
         oqstep = .false.
         if (.not. rtdb_put(rtdb,'driver:oqstep',mt_log,1,oqstep))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'hscale',field)) then
         if(.not. inp_f(hscale)) 
     $        call errquit('driver_input: hscale',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'driver:hscale',mt_dbl,1,hscale))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'ascale',field)) then
         if(.not. inp_f(ascale)) 
     $        call errquit('driver_input: ascale',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'driver:ascale',mt_dbl,1,ascale))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'bscale',field)) then
         if(.not. inp_f(bscale)) 
     $        call errquit('driver_input: bscale',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'driver:bscale',mt_dbl,1,bscale))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'tscale',field)) then
         if(.not. inp_f(tscale)) 
     $        call errquit('driver_input: tscale',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'driver:tscale',mt_dbl,1,tscale))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'xyz',field)) then
         if (.not. inp_a(xyz)) then
            call util_file_prefix('x',xyz)
            xyz(inp_strlen(xyz)-1:) = ' '
         endif
         if (.not. rtdb_cput(rtdb,'driver:xyz',1,xyz))
     $        call errquit('driver_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'noxyz',field)) then
         ignore = rtdb_delete(rtdb,'driver:xyz')
c
      else
         call errquit('driver_input: unknown directive',0, INPUT_ERR)
      endif
c
      goto 10                   ! While loop
c
      end
C>
C> @}

*     *********************************************
*     *                                           *
*     *         driver_parse_boolean              *
*     *                                           *
*     *********************************************
      logical function driver_parse_boolean(word,default)
      implicit none
      character*(*) word
      logical default

#include "inp.fh"
#include "errquit.fh"

      logical value
      value = default
      if (inp_compare(.false.,word,'on'))      value = .true.
      if (inp_compare(.false.,word,'yes'))     value = .true.
      if (inp_compare(.false.,word,'.true.'))  value = .true.
      if (inp_compare(.false.,word,'True'))    value = .true.
      if (inp_compare(.false.,word,'T'))       value = .true.
      if (inp_compare(.false.,word,'off'))     value = .false.
      if (inp_compare(.false.,word,'no'))      value = .false.
      if (inp_compare(.false.,word,'.false.')) value = .false.
      if (inp_compare(.false.,word,'False'))   value = .false.
      if (inp_compare(.false.,word,'F'))       value = .false.
      if (inp_compare(.false.,word,'fix'))     value = .false.

      driver_parse_boolean = value
      return
      end

