      subroutine util_pname(name, pname)
      implicit none
#include "global.fh"
#include "inp.fh"
      character*(*) name, pname
c
c     construct a unique parallel name by appending the process
c     number after the stub name
c     i.e. <pname> = <name>.<mynode>
c
c     e.g., fred.0001, fred.0002, ...
c
c     The number of leading zeroes are adjusted so that there are
c     none in front of the highest numbered processor
c
      integer i, p
c
      i = inp_strlen(name)
      if (i+5.gt.len(pname))
     &  call errquit
     $     ('util_pname: pname too short for name.id',len(pname))
      pname = name
      p = ga_nnodes()
c
      if (p .lt. 10) then
         write(pname(i+1:i+2),1) ga_nodeid()
 1       format('.',i1)
      else if (p .lt. 100) then
         write(pname(i+1:i+3),2) ga_nodeid()
 2       format('.',i2.2)
      else if (p .lt. 1000) then
         write(pname(i+1:i+4),3) ga_nodeid()
 3       format('.',i3.3)
      else 
         write(pname(i+1:i+5),4) ga_nodeid()
 4       format('.',i4.4)
      endif
c
      end
