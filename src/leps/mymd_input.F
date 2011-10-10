      subroutine mymd_input(rtdb)
c
      implicit none
#include "errquit.fh"
c
#include "stdio.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "rtdb.fh"
c
      integer rtdb
c
      character*32 tag
      character*32 pname
      character*255 token
      double precision temp,timestep 
      integer nsteps
c
      pname = "mymd_input: "
c
c      write(luout,*) "in ",pname
c
      call inp_set_field(0)
c
c     start parsing input
c     ------------------
      if (.not.inp_a(token))
     + call errquit('input: no input available',0, INPUT_ERR)
C      write(*,*) 'token is: ',token
c      write(*,*) 'token is: again' 
      if (.not.inp_compare(.false.,token,'mymd'))
     + call errquit('mymd_input: no input available',0, INPUT_ERR)
      goto 2
    1 continue
      if (.not.inp_read()) call errquit('my_input: premature EOF',0,
     &       INPUT_ERR)
    2 continue
      if(.not.inp_a(token)) goto 1
c
c     temperature 
c     -----------------
      if (inp_compare(.false.,"temp", token)) then
        if (.not.inp_f(temp))
     >    call errquit('qmmm_input: failed to read temp',0, INPUT_ERR)
        tag = "qmmm:temp"
        if (.not.rtdb_put(rtdb,"mymd:temp",mt_dbl,1,temp)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
        goto 2
      endif
c
c     number of  steps 
c     -----------
      if (inp_compare(.false.,"nsteps", token)) then
        if (.not.inp_i(nsteps))
     >  call errquit(pname//'failed to read num. of steps',0, INPUT_ERR)
        tag = "mymd:nsteps"
        if (.not.rtdb_put(rtdb,tag(1:inp_strlen(tag)),mt_int,1,nsteps)) 
     >    call errquit(pname//'failed to store '//tag,0, RTDB_ERR)
        goto 2
      endif
c
c     size of timestep
c     ----------------
      if (inp_compare(.false.,"timestep", token) ) then
        if (.not.inp_f(timestep))
     >    call errquit('failed to read timestep',0, INPUT_ERR)
        tag = "mymd:timestep"
        if (.not.rtdb_put(rtdb,tag,mt_dbl,1,timestep)) 
     >      call errquit(pname//'failed to store'//tag,0,RTDB_ERR)
        goto 2
      endif
c      write(*,*) 'End of mymd_input.'
      end
c $Id$
