      logical function pre_rename(lfnout,lprint,ffield,
     + lseq,cseq,mseq,nseq,latm,catm,matm,natm)
c
c $Id: pre_rename.F,v 1.2 1998-06-04 22:28:38 d3j191 Exp $
c
c     function to comply with ffield segment names
c
      implicit none
c
      integer lfnout
      integer lprint
      character*80 ffield
      integer mseq,nseq,matm,natm
      integer lseq(4,mseq),latm(5,matm)
      character*5 cseq(2,mseq)
      character*6 catm(3,matm)
c
      integer i
c
      do 1 i=1,nseq
      cseq(2,i)=cseq(1,i)
    1 continue
c
      if(ffield(1:5).eq.'amber') then
      do 2 i=1,natm
      if(cseq(1,latm(1,i))(1:3).eq.'HIS'.and.catm(2,i).eq.' HD1  ')
     + cseq(2,i)(1:3)='HID'
      if(cseq(1,latm(1,i))(1:3).eq.'HIS'.and.catm(2,i).eq.' HE2  ')
     + cseq(2,i)(1:3)='HIE'
      if(cseq(1,latm(1,i))(1:3).eq.'HID'.and.catm(2,i).eq.' HE2  ')
     + cseq(2,i)(1:3)='HIP'
      if(cseq(1,latm(1,i))(1:3).eq.'HIE'.and.catm(2,i).eq.' HD1  ')
     + cseq(2,i)(1:3)='HIP'
      if(cseq(1,latm(1,i))(1:3).eq.'CYS'.and.lseq(4,latm(1,i)).ne.0)
     + cseq(2,i)(1:3)='CYX'
    2 continue
      endif
c
      if(ffield(1:5).eq.'charm') then
      do 3 i=1,natm
      if(cseq(1,latm(1,i))(1:3).eq.'HIS'.and.catm(2,i).eq.' HD1  ')
     + cseq(2,i)(1:3)='HSD'
      if(cseq(1,latm(1,i))(1:3).eq.'HIS'.and.catm(2,i).eq.' HE2  ')
     + cseq(2,i)(1:3)='HSE'
      if(cseq(1,latm(1,i))(1:3).eq.'HID'.and.catm(2,i).eq.' HE2  ')
     + cseq(2,i)(1:3)='HSP'
      if(cseq(1,latm(1,i))(1:3).eq.'HIE'.and.catm(2,i).eq.' HD1  ')
     + cseq(2,i)(1:3)='HSP'
      if(cseq(1,latm(1,i))(1:3).eq.'CYS'.and.lseq(4,latm(1,i)).ne.0)
     + cseq(2,i)(1:3)='CYX'
    3 continue
      endif
c
      pre_rename=.true.
      return
c
 9999 continue
      pre_rename=.false.
      return
      end
