      subroutine util_pname(name, pname)
*
* $Id$
*
      implicit none
#include "errquit.fh"
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
      integer nnumber,howmany
      nnumber=ga_nodeid()
      howmany=ga_nnodes()
      call util_pname0(name, pname,nnumber,howmany)
      
      end
      subroutine util_pname0(name, pname,nnumber,howmany)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "inp.fh"
      integer nnumber,howmany
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
      if (i+5.gt.len(pname)) then
         write(6,*) ' util_pname: |',name,'|'
         call util_flush(6)
         call errquit
     $     ('util_pname: pname too short for name.id',len(pname),
     &       INPUT_ERR)
      endif
      pname = name
      p = howmany
c
      if (p .lt. 10) then
         write(pname(i+1:i+2),1) nnumber
 1       format('.',i1)
      else if (p .lt. 100) then
         write(pname(i+1:i+3),2) nnumber
 2       format('.',i2.2)
      else if (p .lt. 1000) then
         write(pname(i+1:i+4),3) nnumber
 3       format('.',i3.3)
      else if (p .lt. 10000) then
         write(pname(i+1:i+5),4) nnumber
 4       format('.',i4.4)
      else if (p .lt. 100000) then
         write(pname(i+1:i+6),5) nnumber
 5       format('.',i5.5)
      else 
         write(pname(i+1:i+7),6) nnumber
 6       format('.',i6.6)
      endif
c
      end
