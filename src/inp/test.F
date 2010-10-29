      program test
C$Id$
      implicit none
#include "inp.fh"
      character*30 aval
      integer i, ival, ival2, ival3, line
      double precision dval
c
      line = 0
      call inp_init(5,6)
 10   if (inp_read()) then
         line = line + 1
         call inp_outrec()
         do i = 1, inp_n_field()
            if (inp_i(ival)) then
               write(6,1) line, i, ival
 1             format(' line=',i2,', field=',i2,', integer=',i8)
            Elseif (inp_irange(ival, ival2, ival3)) then
               write(6,5) line, i, ival, ival2, ival3
 5             format(' line=',i2,', field=',i2,', integers=',3i8)
            else if (inp_f(dval)) then
               write(6,2) line, i, dval
 2             format(' line=',i2,', field=',i2,', double=',1pd9.2)
            else if (inp_a(aval)) then
               write(6,3) line, i, aval
 3             format(' line=',i2,', field=',i2,', string=',a)
            else
               write(6,4) line, i
 4             format(' line=',i2,', field=',i2,', error!')
               call inp_errout()
            endif
         enddo
         goto 10
      else 
         if (inp_eof()) then
            write(6,*) ' EOF detected at line ', line
         else
            write(6,*) ' input failed at line ', line
         endif
         call inp_errout()
      endif
c
      end
