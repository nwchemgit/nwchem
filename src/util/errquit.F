      subroutine errquit(string, icode)
C$Id: errquit.F,v 1.9 1996-08-19 15:28:56 d3g681 Exp $
#include "global.fh"
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
      if (ga_nodeid() .eq. 0) write(6,1)
 1    format(1x,72('-'))
      write(6, 2) string, icode
 2    format(1x,a,i4)
      if (ga_nodeid() .eq. 0) then
        write(6,1)
        call inp_errout
        write(6,1)
        write(6,*) ' current input line : '
        call inp_outrec
        write(6,1)
      endif
      call util_flush(6)
      call ga_error(string, icode)
c
      end
