*
* $Id$
*
      subroutine util_set_rtdb_state(rtdb,os,oc,or)
      implicit none
#include "errquit.fh"
*
* sets the rtdb state 
*
c:: includes
#include "stdio.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c:: passed
      integer rtdb ! [input] rtdb handle
      logical os   ! [input] T -> nwchem had start    directive
      logical oc   ! [input] T -> nwchem had continue directive
      logical or   ! [input] T -> nwchem had restart  directive
c:: local
      integer iset 
      integer ncond
      parameter (ncond = 3)
      logical ostate(ncond)
c
c sanity checking
c
      iset = 0
      if (os) iset = iset + 1
      if (oc) iset = iset + 1
      if (or) iset = iset + 1
      if (iset.ne.1) then
        write(luout,*)' os[tart]    = ',os
        write(luout,*)' oc[ontinue] = ',oc
        write(luout,*)' or[estart]  = ',or
        call errquit(' util_set_rtdb_state: fatal error',911, RTDB_ERR)
      endif
      ostate(1) = os
      ostate(2) = oc
      ostate(3) = or
      if(.not.
     &    rtdb_put(rtdb,'rtdb:stored:state',mt_log,ncond,ostate)) then
        call errquit(' util_set_rtdb_state: rtdb_put failed',911,
     &       RTDB_ERR)
      endif
      end
      subroutine util_get_rtdb_state(rtdb,os,oc,or)
      implicit none
#include "errquit.fh"
*
* gets the rtdb state 
*
c:: includes
#include "stdio.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c:: passed
      integer rtdb ! [input] rtdb handle
      logical os   ! [output] T -> nwchem had start    directive
      logical oc   ! [output] T -> nwchem had continue directive
      logical or   ! [output] T -> nwchem had restart  directive
c:: local
      integer iget 
      integer ncond
      parameter (ncond = 3)
      logical ostate(ncond)
c
      if(.not.
     &    rtdb_get(rtdb,'rtdb:stored:state',mt_log,ncond,ostate)) then
        call errquit(' util_get_rtdb_state: rtdb_put failed',911,
     &       RTDB_ERR)
      endif
c
      os = ostate(1) 
      oc = ostate(2) 
      or = ostate(3) 
c
c sanity checking
c
      iget = 0
      if (os) iget = iget + 1
      if (oc) iget = iget + 1
      if (or) iget = iget + 1
      if (iget.ne.1) then
        write(luout,*)' os[tart]    = ',os
        write(luout,*)' oc[ontinue] = ',oc
        write(luout,*)' or[estart]  = ',or
        call errquit(' util_get_rtdb_state: fatal error',911, RTDB_ERR)
      endif
      end
