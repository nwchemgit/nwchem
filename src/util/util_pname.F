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
      integer i
c
      i = inp_strlen(name)
      if (i+5.gt.len(pname))
     &  call errquit
     $     ('util_pname: pname too short for name.id',len(pname))
      pname = name
      write(pname(i+1:i+5),1) ga_nodeid()
1     format('.',i4.4)
c
      end
