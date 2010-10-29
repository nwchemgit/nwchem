      logical function pre_rename(lfnout,lfnmod,filmod,
     + lseq,cseq,mseq,nseq,latm,catm,matm,natm,
     + llnk,clnk,mlnk,nlnk,icyren,his)
c
c $Id$
c
c     function to comply with ffield segment names
c
      implicit none
c
#include "util.fh"
#include "pre_common.fh"
c
      integer lfnout,lfnmod
      character*255 filmod
      character*80 card
      integer mseq,nseq,matm,natm,mlnk,nlnk,icyren
      integer lseq(6,mseq),latm(5,matm)
      character*5 cseq(2,mseq),string
      character*6 catm(3,matm)
      integer llnk(4,mlnk)
      character*4 clnk(2,mlnk)
      character*1 dorr,his
      character*2 ctest
c
      integer i,j,k,l
c
      dorr='D'
c
      do 1 i=1,nseq
      cseq(2,i)=cseq(1,i)
      if(cseq(1,i)(1:3).eq.'U  ') dorr='R'
      if(cseq(1,i)(1:2).eq.'U_') dorr='R'
    1 continue
c
      if(ffield(1:5).eq.'amber') then
      do 2 i=1,natm
c
      if(cseq(2,latm(1,i))(1:3).eq.'HIS'.and.catm(2,i).eq.' HD1  ') then
      cseq(2,latm(1,i))(1:3)='HID'
      if(util_print('connectivity',print_high)) then
      write(lfnout,1004) cseq(2,latm(1,i)),lseq(1,latm(1,i))
 1004 format(' Renamed ',a,i3.3,' (HD1 found)')
      endif
      endif
      if(cseq(2,latm(1,i))(1:3).eq.'HIS'.and.catm(2,i).eq.' HE2  ') then
      cseq(2,latm(1,i))(1:3)='HIE'
      if(util_print('connectivity',print_high)) then
      write(lfnout,1005) cseq(2,latm(1,i)),lseq(1,latm(1,i))
 1005 format(' Renamed ',a,i3.3,' (HE2 found)')
      endif
      endif
      if(cseq(2,latm(1,i))(1:3).eq.'HID'.and.catm(2,i).eq.' HE2  ') then
      cseq(2,latm(1,i))(1:3)='HIP'
      if(util_print('connectivity',print_high)) then
      write(lfnout,1006) cseq(2,latm(1,i)),lseq(1,latm(1,i))
 1006 format(' Renamed ',a,i3.3,' (HE2 found)')
      endif
      endif
      if(cseq(2,latm(1,i))(1:3).eq.'HIE'.and.catm(2,i).eq.' HD1  ') then
      cseq(2,latm(1,i))(1:3)='HIP'
      if(util_print('connectivity',print_high)) then
      write(lfnout,1007) cseq(2,latm(1,i)),lseq(1,latm(1,i))
 1007 format(' Renamed ',a,i3.3,' (HD1 found)')
      endif
      endif
      if(icyren.eq.0) then
      if(cseq(2,latm(1,i))(1:3).eq.'CYS'.and.catm(2,i).eq.' SG   '.and.
     + latm(5,i).eq.3) then
      cseq(2,latm(1,i))(1:3)='CYX'
      if(util_print('connectivity',print_high)) then
      write(lfnout,1008) cseq(2,latm(1,i)),lseq(1,latm(1,i))
 1008 format(' Renamed ',a,i3.3,' (SG found with link id 3)')
      endif
      endif
      else
      if(cseq(2,latm(1,i))(1:3).eq.'CYS'.and.catm(2,i).eq.' SG   ') then
      if(latm(5,i).ge.3) then
      if((cseq(2,llnk(1,latm(5,i)-2))(1:3).eq.'CYS'.or.
     + cseq(2,llnk(1,latm(5,i)-2))(1:3).eq.'CYX').and.
     + (cseq(2,llnk(1,latm(5,i)-2))(1:3).eq.'CYS'.or.
     + cseq(2,llnk(1,latm(5,i)-2))(1:3).eq.'CYX').and.
     + catm(2,llnk(3,latm(5,i)-2)).eq.' SG   '.and.
     + catm(2,llnk(4,latm(5,i)-2)).eq.' SG   ') then
      cseq(2,latm(1,i))(1:3)='CYX'
      if(util_print('connectivity',print_high)) then
      write(lfnout,1009) cseq(2,latm(1,i)),lseq(1,latm(1,i))
 1009 format(' Renamed ',a,i3.3,' (SG-SG found with link id 3)')
      endif
      endif
      endif
      endif
      endif
c
      ctest=cseq(2,latm(1,i))(1:2)
c
      if(ctest.eq.'A '.or.ctest.eq.'A_'.or.
     + ctest.eq.'C '.or.ctest.eq.'C_'.or.
     + ctest.eq.'G '.or.ctest.eq.'G_')
     + cseq(2,latm(1,i))='D'//cseq(2,latm(1,i))(1:4)
      if(ctest.eq.'T '.or.ctest.eq.'T_')
     + cseq(2,latm(1,i))='D'//cseq(2,latm(1,i))(1:4)
      if(ctest.eq.'U '.or.ctest.eq.'U_')
     + cseq(2,latm(1,i))='R'//cseq(2,latm(1,i))(1:4)
c
      if(catm(2,i).eq.' O2*  '.and.
     + (cseq(2,latm(1,i))(1:3).eq.'DA '.or.
     +  cseq(2,latm(1,i))(1:3).eq.'DA_'.or.
     +  cseq(2,latm(1,i))(1:3).eq.'DC '.or.
     +  cseq(2,latm(1,i))(1:3).eq.'DC_'.or.
     +  cseq(2,latm(1,i))(1:3).eq.'DG '.or.
     +  cseq(2,latm(1,i))(1:3).eq.'DG_'))
     + cseq(2,latm(1,i))(1:1)='R'
c
    2 continue
      do 3 i=1,nseq
      if(cseq(2,i)(1:3).eq.'HIS'.and.his.eq.'D') cseq(2,i)(1:3)='HID'
      if(cseq(2,i)(1:3).eq.'HIS'.and.his.eq.'E') cseq(2,i)(1:3)='HIE'
      if(cseq(2,i)(1:3).eq.'HIS'.and.his.eq.'P') cseq(2,i)(1:3)='HIP'
    3 continue
      endif
c
      if(ffield(1:5).eq.'charm') then
      do 4 i=1,natm
      if(cseq(2,latm(1,i))(1:3).eq.'HIS'.and.catm(2,i).eq.' HD1  ')
     + cseq(2,latm(1,i))(1:3)='HSD'
      if(cseq(2,latm(1,i))(1:3).eq.'HIS'.and.catm(2,i).eq.' HE2  ')
     + cseq(2,latm(1,i))(1:3)='HSE'
      if(cseq(2,latm(1,i))(1:3).eq.'HSD'.and.catm(2,i).eq.' HE2  ')
     + cseq(2,latm(1,i))(1:3)='HSP'
      if(cseq(2,latm(1,i))(1:3).eq.'HSE'.and.catm(2,i).eq.' HD1  ')
     + cseq(2,latm(1,i))(1:3)='HSP'
      if(icyren.eq.0) then
      if(cseq(2,latm(1,i))(1:3).eq.'CYS'.and.catm(2,i).eq.' SG   '.and.
     + latm(5,i).eq.3) cseq(2,latm(1,i))(1:3)='CYX'
      else
      if(cseq(2,latm(1,i))(1:3).eq.'CYS'.and.catm(2,i).eq.' SG   ') then
      if(latm(5,i).ge.3) then
      if((cseq(2,llnk(1,latm(5,i)-2))(1:3).eq.'CYS'.or.
     + cseq(2,llnk(1,latm(5,i)-2))(1:3).eq.'CYX').and.
     + (cseq(2,llnk(1,latm(5,i)-2))(1:3).eq.'CYS'.or.
     + cseq(2,llnk(1,latm(5,i)-2))(1:3).eq.'CYX').and.
     + catm(2,llnk(3,latm(5,i)-2)).eq.' SG   '.and.
     + catm(2,llnk(4,latm(5,i)-2)).eq.' SG   ') then
      cseq(2,latm(1,i))(1:3)='CYX'
      endif
      endif
      endif
      endif
    4 continue
      endif
c
c     get modifications
c
      open(unit=lfnmod,file=filmod(1:index(filmod,' ')-1),
     + form='formatted',status='old',err=999)
      rewind(unit=lfnmod)
    5 continue
      read(lfnmod,1000,end=99) card
 1000 format(a)
c
      if(card(1:8).eq.'sequence') then
      card=card(9:80)
      i=index(card,':')
      if(i.le.1) call md_abort('modify sequence error',0)
      read(card(1:i-1),*) j
      do 6 k=1,nseq
      if(lseq(1,k).eq.j) then
      string=cseq(2,k)
      if(index(cseq(2,k),'_').eq.4) then
      cseq(2,k)=card(i+1:i+3)//cseq(2,k)(4:5)
      else
      cseq(2,k)=card(i+1:i+5)
      endif
      write(lfnout,1010) j,string,cseq(2,k)
 1010 format(' Sequence modification: ',i6,':',a5,' to ',a5)
      endif
    6 continue
      endif
      goto 5
   99 continue
      close(unit=lfnmod)
  999 continue
c
      if(util_print('connectivity',print_debug)) then
      if(natm.gt.0) write(lfnout,1001)
 1001 format(/,
     + '    Atom Sgmnt Sgmnt Sgmnt   Atom   Atom Link   ',
     + '    Atom Sgmnt Sgmnt Sgmnt   Atom   Atom Link',/,
     + '  number number name  name   name   name        ',
     + '  number number name  name   name   name     ',/,
     + '                 PDB   top    PDB    top        ',
     + '                 PDB   top    PDB    top     ',/)
      k=natm/2
      if(natm.gt.2*k) k=k+1
      do 7 i=1,k
      l=i+k
      if(l.le.natm) then
      write(lfnout,1002) i,lseq(1,latm(1,i)),(cseq(j,latm(1,i)),j=1,2),
     + (catm(j,i),j=1,2),latm(5,i),
     + l,lseq(1,latm(1,l)),(cseq(j,latm(1,l)),j=1,2),
     + (catm(j,l),j=1,2),latm(5,l)
      else
      write(lfnout,1002) i,lseq(1,latm(1,i)),(cseq(j,latm(1,i)),j=1,2),
     + (catm(j,i),j=1,2),latm(5,i)
      endif
 1002 format(2(i8,i6,2(1x,a5),2(1x,a6),i5,3x))
    7 continue
      if(util_print('topology',print_low)) then
      write(lfnout,1003)
 1003 format(' ')
      endif
      endif
c
      pre_rename=.true.
      return
c
 9999 continue
      pre_rename=.false.
      return
      end
