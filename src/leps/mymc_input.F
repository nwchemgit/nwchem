      subroutine mymc_input(rtdb)
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
      double precision temp,stepsize,ratio,seed
      integer mcsteps,mode,gseed
c
      pname = "mymc_input: "
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
      if (.not.inp_compare(.false.,token,'mymc'))
     + call errquit('mymc_input: no input available',0, INPUT_ERR)
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
     >    call errquit('mc_input: failed to read temp',0, INPUT_ERR)
        tag = "mc:temp"
        if (.not.rtdb_put(rtdb,"mymc:temp",mt_dbl,1,temp)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
        goto 2
      endif
c
c     number of  Monte carlo steps 
c     -----------
      if (inp_compare(.false.,"mcsteps", token)) then
        if (.not.inp_i(mcsteps))
     >  call errquit(pname//'failed to read MC steps',0, INPUT_ERR)
        tag = "mymc:mcsteps"
        if (.not.rtdb_put(rtdb,tag(1:inp_strlen(tag)),mt_int,1,mcsteps))
     >    call errquit(pname//'failed to store '//tag,0, RTDB_ERR)
        goto 2
      endif

c
c     number of  Monte carlo steps 
c     -----------
      if (inp_compare(.false.,"mode", token)) then
        if (.not.inp_i(mode))
     >  call errquit(pname//'failed to read mode',0, INPUT_ERR)
        tag = "mymc:mode"
        if (.not.rtdb_put(rtdb,tag(1:inp_strlen(tag)),mt_int,1,mode))
     >    call errquit(pname//'failed to store '//tag,0, RTDB_ERR)
        goto 2
      endif
c
c     size of stepsize
c     ----------------
      if (inp_compare(.false.,"stepsize", token) ) then
        if (.not.inp_f(stepsize))
     >    call errquit('failed to read stepsize',0, INPUT_ERR)
        tag = "mymc:stepsize"
        if (.not.rtdb_put(rtdb,tag,mt_dbl,1,stepsize)) 
     >      call errquit(pname//'failed to store'//tag,0,RTDB_ERR)
        goto 2
      endif

c
c     size of stepsize
c     ----------------
      if (inp_compare(.false.,"ratio", token) ) then
        if (.not.inp_f(ratio))
     >    call errquit('failed to read ratio',0, INPUT_ERR)
        tag = "mymc:ratio"
        if (.not.rtdb_put(rtdb,tag,mt_dbl,1,ratio))
     >      call errquit(pname//'failed to store'//tag,0,RTDB_ERR)
        goto 2
      endif

c
c     seed for dynrnd
c     ----------------
      if (inp_compare(.false.,"seed", token) ) then
        if (.not.inp_f(seed))
     >    call errquit('failed to read seed',0, INPUT_ERR)
        tag = "mymc:seed"
        if (.not.rtdb_put(rtdb,tag,mt_dbl,1,seed))
     >      call errquit(pname//'failed to store'//tag,0,RTDB_ERR)
        goto 2
      endif

c
c     seed for  gasdev
c     ----------------
      if (inp_compare(.false.,"gseed", token)) then
        if (.not.inp_i(gseed))
     >  call errquit(pname//'failed to read gseed',0, INPUT_ERR)
        tag = "mymc:gseed"
        if (.not.rtdb_put(rtdb,tag(1:inp_strlen(tag)),mt_int,1,gseed))
     >    call errquit(pname//'failed to store '//tag,0, RTDB_ERR)
        goto 2
      endif

c      write(*,*) 'End of mymc_input.'
      end


c $Id$
