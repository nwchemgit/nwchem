      subroutine errquit(string, icode)
      implicit none
C$Id: errquit.F,v 1.13 2000-04-10 12:05:34 hvd Exp $
#include "global.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      character*(*) string
      integer icode
c     
c     error termination 
c     
c     string = error message printed to stdout
c     icode  = informative value printed to stdout
c
c     try to minimize the amount of output if all nodes 
c     detect an error
c     
      call ecce_print_module_entry('errquit')
      call ecce_print1_char('errquit message',string,1)
      call ecce_print1('errquit code', mt_int,icode,1)
      call ecce_print_module_exit('errquit', 'failed')
c
      if (ga_nodeid() .eq. 0) write(LuOut,1)
 1    format(1x,72('-'))
      write(LuOut, 2) string, icode
 2    format(1x,a,i9)
      if (ga_nodeid() .eq. 0) then
        write(LuOut,1)
        call inp_errout
        write(LuOut,1)
        write(LuOut,*) ' current input line : '
        call inp_outrec
        write(LuOut,1)
      endif
      call util_flush(LuOut)
      call ga_error(string, icode)
c
      end
