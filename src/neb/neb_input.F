c $Id: neb_input.F,v 1.1 2004-02-11 02:37:09 marat Exp $
      subroutine neb_input(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
c     
      integer rtdb
      character*255 test
      integer nbeads
      double precision timestep

      call banner(6, 
     & "Processing input for NEB module", 
     & "*",.true.,.true.,.true.)
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
c     time step
c     ----------------------------
      else if (inp_compare(.false.,'timestep', test)) then
       if(.not. inp_f(timestep))
     &    call errquit('neb_input: timestep',911,
     &    INPUT_ERR)
       if (.not. rtdb_put(rtdb,'neb:timestep',mt_dbl,1,timestep))
     $     call errquit('driver_input: rtdb put failed',0,
     &     RTDB_ERR)
       write(6,*) 'neb:timestep',timestep

      goto 100

c
c     number of beads
c     ----------------------------
      else if (inp_compare(.false.,'nbeads', test)) then
       if(.not. inp_i(nbeads))
     &    call errquit('neb_input: nbeads',911,
     &    INPUT_ERR)
       if (.not. rtdb_put(rtdb,'neb:nbeads',mt_int,1,nbeads))
     $     call errquit('driver_input: rtdb put failed',0,
     &     RTDB_ERR)
       write(6,*) 'neb:nbeads',nbeads

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

