      subroutine qmmm_input(rtdb)
c
      implicit none
c
#include "stdio.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "rtdb.fh"
c
      integer rtdb
c
      character*255 token
      double precision Eatoms
c
      call inp_set_field(0)
c
      if (.not.inp_a(token))
     + call errquit('qmmm_input: no input available',911)
      if (.not.inp_compare(.false.,token,'qmmm')) goto 4
      goto 2
    1 continue
      if (.not.inp_read()) call errquit('qmmm_input: premature EOF',911)
    2 continue
      if(.not.inp_a(token)) goto 1
      if (token.eq.'end') then
      goto 3
      else if (token.eq.'eatoms') then
      if (inp_f(Eatoms)) then
      if (.not. rtdb_put(rtdb,'qmmm:uqmatm',mt_dbl,1,Eatoms))
     + call errquit
     + ('qmmm_input: failed to store qmmm:uqmatm',911)
      write(luout,*)' QM/MM separated atom energy:',Eatoms
      call util_flush(luout)
      else
      call errquit
     + ('qmmm_input: failed to read eatoms value as a real number',911)
      endif
      goto 2
      else
      write(luout,*)' unrecognized token in qmmm input:',
     + token(1:inp_strlen(token))
      goto 4
      endif
    3 continue
      return
    4 continue
      write(luout,1000)
 1000 format(/,/,'QM/MM input format:',/,'qmmm',/,'eatoms <real>',/,
     + 'end',///)
      call errquit ('qmmm_input: called with bad input',911)
      return
      end
