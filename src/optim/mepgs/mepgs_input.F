      subroutine mepgs_input(rtdb)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "nwc_const.fh"
#include "cgsopt.fh"
#include "cmepgs.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
c
      integer rtdb
      character*255 field
      double precision ascale, bscale, tscale, hscale
      logical ignore
      integer mh
c
c     === set defaults ===
c
      mswg = .false.
c
c     ----- read values from input -----
c
 10   if(.not.inp_read())
     $     call errquit('mepgs_input: unexpected eof',911, INPUT_ERR)
      if (.not. inp_a(field))
     $     call errquit('mepgs_input: failed to read field',911,
     &       INPUT_ERR)
c
      if(inp_compare(.false.,'end',field)) then
         return
c
      else if(inp_compare(.false.,'evib',field)) then
         if(.not. inp_f(evib))
     $        call errquit('mepgs_input: evib',911, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'mepgs:evib',mt_dbl,1,evib))
     $        call errquit('mepgs_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'stride',field)) then
         if(.not. inp_f(stride)) 
     $        call errquit('mepgs_input: stride',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'mepgs:stride',mt_dbl,1,stride))
     $        call errquit('mepgs_input: rtdb put failed',0, INPUT_ERR)
      else if(inp_compare(.false.,'maxst',field)) then
         if(.not. inp_f(maxst))
     $        call errquit('mepgs_input: maxst',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'mepgs:maxst',mt_dbl,1,maxst))
     $        call errquit('mepgs_input: rtdb put failed',0, INPUT_ERR)
c
      else if(inp_compare(.false.,'nircopt',field) .or. 
     $        inp_compare(.false.,'maxmep',field)) then
         if(.not. inp_i(nircopt))
     $        call errquit('mepgs_input: nircopt/maxmep',911,
     &       INPUT_ERR)
         if (.not. rtdb_put(rtdb,'mepgs:nircopt',mt_int,1,nircopt))
     $        call errquit('mepgs_input: rtdb put failed',0, RTDB_ERR)
c
      else if (inp_compare(.false.,'mswg',field)) then
         mswg = .true.
         if (.not. rtdb_put(rtdb,'ircgs:mswg',mt_log,1,mswg))
     $        call errquit('ircgs_input: rtdb put failed',0, RTDB_ERR)
c
      else if (inp_compare(.false.,'nomswg',field)) then
         mswg = .false.
         if (.not. rtdb_put(rtdb,'ircgs:mswg',mt_log,1,mswg))
     $        call errquit('ircgs_input: rtdb put failed',0, RTDB_ERR)
c
      else if (inp_compare(.false.,'backward',field)) then
         backward = .true.
         if (.not. rtdb_put(rtdb,'mepgs:backward',mt_log,1,backward))
     $        call errquit('mepgs_input: rtdb put failed',0, RTDB_ERR)
      else if (inp_compare(.false.,'forward',field)) then
         forward = .true.
         if (.not. rtdb_put(rtdb,'mepgs:forward',mt_log,1,forward))
     $        call errquit('mepgs_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'print',field)) then
         call util_print_input(rtdb,'mepgs')
c
      else if(inp_compare(.false.,'noprint',field)) then
         call util_print_input(rtdb,'mepgs')
c
      else if(inp_compare(.false.,'xyz',field)) then
         if (.not. inp_a(xyz)) then
            call util_file_prefix('x',xyz)
            xyz(inp_strlen(xyz)-1:) = ' '
         endif
         if (.not. rtdb_cput(rtdb,'mepgs:xyz',1,xyz))
     $        call errquit('mepgs_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'noxyz',field)) then
         ignore = rtdb_delete(rtdb,'mepgs:xyz')
c
c     Constrained optimization control 
c
      else if(inp_compare(.false.,'nptopt',field) .or.
     $        inp_compare(.false.,'maxiter',field)) then
         if(.not. inp_i(nptopt))
     $        call errquit('mepgs_input: nptopt/maxiter',911,
     &       INPUT_ERR)
         if (.not. rtdb_put(rtdb,'gsopt:nptopt',mt_int,1,nptopt))
     $        call errquit('mepgs_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'inhess',field)) then
         if(.not. inp_i(inhess)) 
     $        call errquit('mepgs_input: inhess',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'mepgs:inhess',mt_int,1,inhess))
     $        call errquit('mepgs_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'opttol',field)) then
         if(.not. inp_f(opt_tol))
     $        call errquit('mepgs_input: opttol',911, INPUT_ERR)
         if (.not. rtdb_put(rtdb,'gsopt:opt_tol',mt_dbl,1,opt_tol))
     $        call errquit('mepgs_input: rtdb put failed',0, RTDB_ERR)
         grms_tol = 1.0d0*opt_tol
         gmax_tol = 1.5d0*opt_tol
         xrms_tol = 4.0d0*opt_tol
         xmax_tol = 6.0d0*opt_tol
         if (.not. rtdb_put(rtdb,'gsopt:gmax_tol',mt_dbl,1,gmax_tol))
     $        call errquit('mepgs_input: rtdb put failed',0, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'gsopt:grms_tol',mt_dbl,1,grms_tol))
     $        call errquit('mepgs_input: rtdb put failed',0, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'gsopt:xrms_tol',mt_dbl,1,xrms_tol))
     $        call errquit('mepgs_input: rtdb put failed',0, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'gsopt:xmax_tol',mt_dbl,1,xmax_tol))
     $        call errquit('mepgs_input: rtdb put failed',0, RTDB_ERR)
c
      else if(inp_compare(.false.,'eprec',field)) then
         if(.not. inp_f(eprec))
     $        call errquit('mepgs_input: eprec',911, RTDB_ERR)
         if (.not. rtdb_put(rtdb,'gsopt:eprec',mt_dbl,1,eprec))
     $        call errquit('mepgs_input: rtdb put failed',0, RTDB_ERR)

c
      else
         call errquit('mepgs_input: unknown directive',0, INPUT_ERR)
      endif
c
      goto 10                   ! While loop
c
      end
