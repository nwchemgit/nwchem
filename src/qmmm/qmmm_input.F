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
      real*8 eatoms,qatoms
      integer linkatm,nobq
c
      call inp_set_field(0)
c
      linkatm=0
      nobq=1
c
      if (.not.inp_a(token))
     + call errquit('qmmm_input: no input available',0)
      if (.not.inp_compare(.false.,token,'qmmm'))
     + call errquit('qmmm_input: no input available',0)
      goto 2
    1 continue
      if (.not.inp_read()) call errquit('qmmm_input: premature EOF',0)
    2 continue
      if(.not.inp_a(token)) goto 1
c
      if (token.eq.'eatoms') then
      if (.not.inp_f(eatoms))
     +  call errquit('qmmm_input: failed to read eatoms',0)
      write(luout,*)' QM/MM separated atom energy:',eatoms
      call util_flush(luout)
      goto 2
      endif
c
      if (token.eq.'qatoms') then
      if (.not.inp_f(qatoms))
     +  call errquit('qmmm_input: failed to read qatoms',0)
      write(luout,*)' QM/MM atom charge:',qatoms
      call util_flush(luout)
      if (.not.rtdb_put(rtdb,'qmmm:qqmatm',mt_dbl,1,qatoms))
     + call errquit('qmmm_input: failed to store qmmm:uqmatm',0)
      goto 2
      endif
c
      if(token.eq.'link') then
      if(.not.inp_a(token)) call errquit('qmmm_input: no link atom',0)
      if(token.eq.'hydrogen') then
      linkatm=1
      goto 2
      endif
      call errquit('qmmm_input: invalid link atom',0)
      endif
c
      if(token.eq.'nobq') then
      if(.not.inp_a(token)) call errquit('qmmm_input: no nobq atom',0)
      if(token.eq.'hydrogen') then
      nobq=1
      goto 2
      endif
      if(token.eq.'all') then
      nobq=2
      goto 2
      endif
      if(token.eq.'none') then
      nobq=0
      goto 2
      endif
      call errquit('qmmm_input: invalid nobq atom',0)
      endif
c
      if (token.eq.'end') then
      if (.not.rtdb_put(rtdb,'qmmm:uqmatm',mt_dbl,1,eatoms))
     + call errquit('qmmm_input: failed to store qmmm:uqmatm',0)
      if (.not.rtdb_put(rtdb,'qmmm:linkatm',mt_int,1,linkatm))
     + call errquit('qmmm_input: failed to store qmmm:linkatm',linkatm)
      if (.not.rtdb_put(rtdb,'qmmm:nobq',mt_int,1,nobq))
     + call errquit('qmmm_input: failed to store qmmm:nobq',nobq)
      return
      endif
c
      write(luout,*)' unrecognized token in qmmm input:',
     + token(1:inp_strlen(token))
      write(luout,1000)
 1000 format(/,/,'QM/MM input format:',/,'qmmm',/,'eatoms <real>',/,
     + 'end',///)
      call errquit ('qmmm_input: called with bad input',911)
      return
      end
