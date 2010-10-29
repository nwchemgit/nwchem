      program test
C$Id$
      implicit none
#include "inp.fh"
      character*30 aval
      integer i, ival, ival2, ival3, line, N, maxlist, j
      Parameter (MaxList = 100)
      Integer List(MaxList)
      double precision dval
c
      line = 0
      call inp_init(5,6)
 10   if (inp_read()) then
         line = line + 1
         call inp_outrec()
c$$$         do i = 1, inp_n_field()
            if (inp_ilist(MaxList, List, N)) then
               write(6,6) line, i, N, (List(j), j=1, N)
 6             format(' line=',i2,', field=',i2,', elements=',i8,
     $            /' data = ', 7i8/(8X,7i8))
            Else
               Write (6, 6) line, i, n
               Call Inp_ErrOut()
            EndIf
c$$$            Elseif (inp_i(ival)) then
c$$$               write(6,1) line, i, ival
c$$$ 1             format(' line=',i2,', field=',i2,', integer=',i8)
c$$$            Elseif (inp_irange(ival, ival2, ival3)) then
c$$$               write(6,5) line, i, ival, ival2, ival3
c$$$ 5             format(' line=',i2,', field=',i2,', integers=',3i8)
c$$$            else if (inp_f(dval)) then
c$$$               write(6,2) line, i, dval
c$$$ 2             format(' line=',i2,', field=',i2,', double=',1pd9.2)
c$$$            else if (inp_a(aval)) then
c$$$               write(6,3) line, i, aval
c$$$ 3             format(' line=',i2,', field=',i2,', string=',a)
c$$$            else
c$$$               write(6,4) line, i
c$$$ 4             format(' line=',i2,', field=',i2,', error!')
c$$$               call inp_errout()
c$$$            endif
c$$$         enddo
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
