      logical function sym_mol(zname, ipntmol)
c     
c     implicit double precision (a-h,p-w)
      implicit none
#include "errquit.fh"
#include "util.fh"
c     
      character*8 zname
      character*8 zgroup
      integer maxireps
      parameter (maxireps=20)
      integer ipntmol(maxireps)
c
      zgroup = zname
      call inp_lcase(zgroup)
c     
      if(zgroup.eq.'c1') then
         ipntmol(1)=1            
         goto 10
      endif
c     
      if(zgroup.eq.'c2'.or.zgroup.eq.'ci'.or.zgroup.eq.'cs') then
         ipntmol(1)=1
         ipntmol(1)=2
         goto 10
      endif
c     
      if(zgroup.eq.'c2h'.or.zgroup.eq.'c2v') then
         ipntmol(1)=1
         ipntmol(2)=3
         ipntmol(3)=2
         ipntmol(4)=4
         goto 10
      endif
c     
      if(zgroup.eq.'d2') then
         ipntmol(1)=1
         ipntmol(2)=2
         ipntmol(3)=3
         ipntmol(4)=4
         goto 10
      endif
c
      if(zgroup.eq.'d2h') then
         ipntmol(1)=1
         ipntmol(2)=8
         ipntmol(3)=4
         ipntmol(4)=5
         ipntmol(5)=6
         ipntmol(6)=3
         ipntmol(7)=7
         ipntmol(8)=2
         goto 10
      endif
c     
c     here it will come if group was not recognized
c     
      sym_mol = .false.
      return
c
c     here it comes on success
c
 10   continue
      sym_mol = .true.
c
      end
c $Id$
