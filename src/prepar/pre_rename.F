      logical function pre_rename(lfnout,lprint,ffield,
     + lseq,cseq,mseq,nseq,latm,catm,matm,natm)
c
c $Id: pre_rename.F,v 1.3 1998-06-05 22:44:00 d3j191 Exp $
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
      integer i,j
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
      if(lprint.ge.5) write(lfnout,'(i5,2a5,3a6,i5)')
     + i,(cseq(j,latm(1,i)),j=1,2),(catm(j,i),j=1,3),latm(5,i)
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
      if(lprint.ge.5) write(lfnout,'(i5,2a5,3a6,i5)')
     + i,(cseq(j,latm(1,i)),j=1,2),(catm(j,i),j=1,3),latm(5,i)
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
