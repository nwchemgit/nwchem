c $Id$

*     *************************************
*     *                                   *
*     *           neb_input               *
*     *                                   *
*     *************************************

      subroutine neb_input(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
c     
      logical impose
      integer rtdb
      character*255 test,test2
      integer nbeads
      double precision timestep,kbeads

*     **** external functions ****
      logical  neb_parse_boolean
      external neb_parse_boolean

c
 100  if (.not. inp_read()) call errquit('neb_input: inp_read failed',0,
     &       0)
c
      if (.not. inp_a(test))
     $     call errquit('neb_input: failed to read keyword', 0,
     &       0)
c     
c     print level is read
c     -------------------
      if (inp_compare(.false.,'print', test)) then
         call util_print_input(rtdb, "neb")
c
      goto 100

c
c     algorihm - 
c     ----------------------------
      else if (inp_compare(.false.,'algorithm', test)) then
       if(.not. inp_i(nbeads))
     &    call errquit('neb_input: algorithm',911,
     &    INPUT_ERR)
       if (.not. rtdb_put(rtdb,'neb:algorithm',mt_int,1,nbeads))
     $     call errquit('neb_input: rtdb put failed',0,
     &     RTDB_ERR)
      goto 100

c
c     nhist = m - number histories
c     ----------------------------
      else if (inp_compare(.false.,'nhist', test)) then
       if(.not. inp_i(nbeads))
     &    call errquit('neb_input: algorithm',911,
     &    INPUT_ERR)
       if (.not. rtdb_put(rtdb,'neb:m',mt_int,1,nbeads))
     $     call errquit('neb_input: rtdb put failed',0,
     &     RTDB_ERR)
      goto 100


c
c     stepsize
c     ----------------------------
      else if (inp_compare(.false.,'stepsize', test)) then
       if(.not. inp_f(timestep))
     &    call errquit('neb_input: timestep',911,
     &    INPUT_ERR)
       if (.not. rtdb_put(rtdb,'neb:stepsize',mt_dbl,1,timestep))
     $     call errquit('neb_input: rtdb put failed',0,
     &     RTDB_ERR)

      goto 100

c
c     maxiter
c     ----------------------------
      else if (inp_compare(.false.,'maxiter', test)) then
       if(.not. inp_i(nbeads))
     &    call errquit('neb_input: maxiter',911,
     &    INPUT_ERR)
       if (.not. rtdb_put(rtdb,'neb:steps',mt_int,1,nbeads))
     $     call errquit('neb_input: rtdb put failed',0,
     &     RTDB_ERR)

      goto 100

c
c     nbeads - number of beads
c     ----------------------------
      else if (inp_compare(.false.,'nbeads', test)) then
       if(.not. inp_i(nbeads))
     &    call errquit('neb_input: nbeads',911,
     &    INPUT_ERR)
       if (.not. rtdb_put(rtdb,'neb:nbeads',mt_int,1,nbeads))
     $     call errquit('neb_input: rtdb put failed',0,
     &     RTDB_ERR)

      goto 100

c
c     kbeads - neb spring constant
c     ----------------------------
      else if (inp_compare(.false.,'kbeads', test)) then
       if(.not. inp_f(kbeads))
     &    call errquit('neb_input: kbeads',911,
     &    INPUT_ERR)
       if (.not. rtdb_put(rtdb,'neb:kbeads',mt_dbl,1,kbeads))
     $     call errquit('neb_input: rtdb put failed',0,
     &     RTDB_ERR)

      goto 100

c
c     reset  
c     ----------------------------
      else if (inp_compare(.false.,'reset', test)) then
       if (.not. rtdb_put(rtdb,'neb:nebnew',mt_log,1,.true.))
     $     call errquit('neb_input: rtdb put failed',0,
     &     RTDB_ERR)

      goto 100

c
c     gmax - neb spring constant
c     ----------------------------
      else if (inp_compare(.false.,'gmax', test)) then
       if(.not. inp_f(kbeads))
     &    call errquit('neb_input: gmax',911,
     &    INPUT_ERR)
       if (.not. rtdb_put(rtdb,'neb:gmax',mt_dbl,1,kbeads))
     $     call errquit('neb_input: rtdb put failed',0,
     &     RTDB_ERR)
      goto 100
c
c     grms - neb spring constant
c     ----------------------------
      else if (inp_compare(.false.,'grms', test)) then
       if(.not. inp_f(kbeads))
     &    call errquit('neb_input: gmax',911,
     &    INPUT_ERR)
       if (.not. rtdb_put(rtdb,'neb:grms',mt_dbl,1,kbeads))
     $     call errquit('neb_input: rtdb put failed',0,
     &     RTDB_ERR)
      goto 100
c
c     xmax - neb spring constant
c     ----------------------------
      else if (inp_compare(.false.,'xmax', test)) then
       if(.not. inp_f(kbeads))
     &    call errquit('neb_input: xmax',911,
     &    INPUT_ERR)
       if (.not. rtdb_put(rtdb,'neb:xmax',mt_dbl,1,kbeads))
     $     call errquit('neb_input: rtdb put failed',0,
     &     RTDB_ERR)
      goto 100
c
c     xrms - neb spring constant
c     ----------------------------
      else if (inp_compare(.false.,'xrms', test)) then
       if(.not. inp_f(kbeads))
     &    call errquit('neb_input: xrms',911,
     &    INPUT_ERR)
       if (.not. rtdb_put(rtdb,'neb:xrms',mt_dbl,1,kbeads))
     $     call errquit('neb_input: rtdb put failed',0,
     &     RTDB_ERR)
      goto 100

c
c    loose - 
c     ----------------------------
      else if (inp_compare(.false.,'loose', test)) then
       if (.not. rtdb_put(rtdb,'neb:gmax',mt_dbl,1,0.0045d0))
     >     call errquit('neb_input: rtdb put failed',0,
     >     RTDB_ERR)
       if (.not. rtdb_put(rtdb,'neb:grms',mt_dbl,1,0.0030d0))
     >     call errquit('neb_input: rtdb put failed',0,
     >     RTDB_ERR)
       if (.not. rtdb_put(rtdb,'neb:xmax',mt_dbl,1,0.0054d0))
     >     call errquit('neb_input: rtdb put failed',0,
     >     RTDB_ERR)
       if (.not. rtdb_put(rtdb,'neb:xrms',mt_dbl,1,0.0036d0))
     >     call errquit('neb_input: rtdb put failed',0,
     >     RTDB_ERR)
      goto 100
c
c    default -
c     ----------------------------
      else if (inp_compare(.false.,'default', test)) then
       if (.not. rtdb_put(rtdb,'neb:gmax',mt_dbl,1,0.00045d0))
     >     call errquit('neb_input: rtdb put failed',0,
     >     RTDB_ERR)
       if (.not. rtdb_put(rtdb,'neb:grms',mt_dbl,1,0.00030d0))
     >     call errquit('neb_input: rtdb put failed',0,
     >     RTDB_ERR)
       if (.not. rtdb_put(rtdb,'neb:xmax',mt_dbl,1,0.0018d0))
     >     call errquit('neb_input: rtdb put failed',0,
     >     RTDB_ERR)
       if (.not. rtdb_put(rtdb,'neb:xrms',mt_dbl,1,0.0012d0))
     >     call errquit('neb_input: rtdb put failed',0,
     >     RTDB_ERR)
      goto 100
c
c    tight -
c     ----------------------------
      else if (inp_compare(.false.,'tight', test)) then
       if (.not. rtdb_put(rtdb,'neb:gmax',mt_dbl,1,0.00015d0))
     >     call errquit('neb_input: rtdb put failed',0,
     >     RTDB_ERR)
       if (.not. rtdb_put(rtdb,'neb:grms',mt_dbl,1,0.00001d0))
     >     call errquit('neb_input: rtdb put failed',0,
     >     RTDB_ERR)
       if (.not. rtdb_put(rtdb,'neb:xmax',mt_dbl,1,0.00006d0))
     >     call errquit('neb_input: rtdb put failed',0,
     >     RTDB_ERR)
       if (.not. rtdb_put(rtdb,'neb:xrms',mt_dbl,1,0.00004d0))
     >     call errquit('neb_input: rtdb put failed',0,
     >     RTDB_ERR)
      goto 100

c
c    impose 
c     ----------------------------
      else if (inp_compare(.false.,'impose', test)) then
         impose = .true.
         if(inp_a(test2)) impose = neb_parse_boolean(test2,.true.)
         if (.not. rtdb_put(rtdb,'neb:impose',mt_log,1,impose))
     >     call errquit('neb_input: rtdb put failed',0,
     >     RTDB_ERR)
      
      goto 100

c
c     end of input block
c     -------------------
      else if (inp_compare(.false.,'end', test)) then
       goto 101
      endif
c     
c
101   continue
      
      return
      end
c



*     *********************************************
*     *                                           *
*     *          neb_parse_boolean                *
*     *                                           *  
*     *********************************************  
      logical function neb_parse_boolean(word,default)
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
      if (inp_compare(.false.,word,'T'))       value = .true.
      if (inp_compare(.false.,word,'off'))     value = .false.
      if (inp_compare(.false.,word,'no'))      value = .false.
      if (inp_compare(.false.,word,'.false.')) value = .false.
      if (inp_compare(.false.,word,'F'))       value = .false.
      if (inp_compare(.false.,word,'fix'))     value = .false.

      neb_parse_boolean = value
      return
      end 

