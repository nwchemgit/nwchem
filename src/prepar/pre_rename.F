      logical function pre_rename(lfnout,ffield,
     + lseq,cseq,mseq,nseq,latm,catm,matm,natm)
c
c $Id: pre_rename.F,v 1.4 1998-06-08 23:29:37 d3j191 Exp $
c
c     function to comply with ffield segment names
c
      implicit none
c
#include "util.fh"
c
      integer lfnout
      character*80 ffield
      integer mseq,nseq,matm,natm
      integer lseq(4,mseq),latm(5,matm)
      character*5 cseq(2,mseq)
      character*6 catm(3,matm)
c
      integer i,j,k,l
c
      do 1 i=1,nseq
      cseq(2,i)=cseq(1,i)
    1 continue
c
      if(ffield(1:5).eq.'amber') then
      do 2 i=1,natm
      if(cseq(2,latm(1,i))(1:3).eq.'HIS'.and.catm(2,i).eq.' HD1  ')
     + cseq(2,latm(1,i))(1:3)='HID'
      if(cseq(2,latm(1,i))(1:3).eq.'HIS'.and.catm(2,i).eq.' HE2  ')
     + cseq(2,latm(1,i))(1:3)='HIE'
      if(cseq(2,latm(1,i))(1:3).eq.'HID'.and.catm(2,i).eq.' HE2  ')
     + cseq(2,latm(1,i))(1:3)='HIP'
      if(cseq(2,latm(1,i))(1:3).eq.'HIE'.and.catm(2,i).eq.' HD1  ')
     + cseq(2,latm(1,i))(1:3)='HIP'
      if(cseq(2,latm(1,i))(1:3).eq.'CYS'.and.catm(2,i).eq.' SG   '.and.
     + latm(5,i).eq.3) cseq(2,latm(1,i))(1:3)='CYX'
    2 continue
      endif
c
      if(ffield(1:5).eq.'charm') then
      do 3 i=1,natm
      if(cseq(2,latm(1,i))(1:3).eq.'HIS'.and.catm(2,i).eq.' HD1  ')
     + cseq(2,latm(1,i))(1:3)='HSD'
      if(cseq(2,latm(1,i))(1:3).eq.'HIS'.and.catm(2,i).eq.' HE2  ')
     + cseq(2,latm(1,i))(1:3)='HSE'
      if(cseq(2,latm(1,i))(1:3).eq.'HSD'.and.catm(2,i).eq.' HE2  ')
     + cseq(2,latm(1,i))(1:3)='HSP'
      if(cseq(2,latm(1,i))(1:3).eq.'HSE'.and.catm(2,i).eq.' HD1  ')
     + cseq(2,latm(1,i))(1:3)='HSP'
      if(cseq(2,latm(1,i))(1:3).eq.'CYS'.and.catm(2,i).eq.' SG   '.and.
     + latm(5,i).eq.3) cseq(2,latm(1,i))(1:3)='CYX'
    3 continue
      endif
c
      if(util_print('connectivity',print_high)) then
      write(lfnout,1000)
 1000 format(
     + '    Atom Sgmnt Sgmnt Sgmnt   Atom   Atom Link   ',
     + '    Atom Sgmnt Sgmnt Sgmnt   Atom   Atom Link',/,
     + '  number number name  name   name   name        ',
     + '  number number name  name   name   name     ',/,
     + '                 PDB   top    PDB    top        ',
     + '                 PDB   top    PDB    top     ',/)
      k=natm/2
      if(natm.gt.2*k) k=k+1
      do 4 i=1,k
      l=i+k
      if(l.le.natm) then
      write(lfnout,1001) i,latm(1,i),(cseq(j,latm(1,i)),j=1,2),
     + (catm(j,i),j=1,2),latm(5,i),
     + l,latm(1,l),(cseq(j,latm(1,l)),j=1,2),
     + (catm(j,l),j=1,2),latm(5,l)
      else
      write(lfnout,1001) i,latm(1,i),(cseq(j,latm(1,i)),j=1,2),
     + (catm(j,i),j=1,2),latm(5,i)
      endif
 1001 format(2(i8,i6,2(1x,a5),2(1x,a6),i5,3x))
    4 continue
      write(lfnout,1002)
 1002 format(' ')
      endif
c
      pre_rename=.true.
      return
c
 9999 continue
      pre_rename=.false.
      return
      end
