      subroutine dia_wrtxyz(lfnxxx,filxxx,sgmnam,isel,xs,qs,
     + iwsel,xw,logs,logw,ndxw)
c
c $Id$
c
      implicit none
c
#include "dia_common.fh"
c
      character*16 sgmnam(msa)
      integer isel(msa),iwsel(mwm),ndxw(nwm)
      integer lfnxxx
      real*8 xs(msa,3),qs(msa),xw(mwm,mwa,3)
      character*255 filxxx
      logical logs,logw
      character*2 atom
c
      integer i,j,k,length
c
      if(me.eq.0) then
c
      length=index(filxxx,' ')-1
      open(unit=lfnxxx,file=filxxx(1:length),form='formatted',
     + status='unknown')
c
      rewind(lfnxxx)
c
      write(lfnxxx,1000) nsa
 1000 format('* Coordinates in XYZ format',/,i5)
      if(logs) then
      do 1 i=1,nsa
      if(isel(i).gt.0) then
      atom=sgmnam(i)(6:7)
      if(atom(1:1).eq.'1') atom(1:1)=' '
      if(atom(1:1).eq.'2') atom(1:1)=' '
      if(atom(1:1).eq.'3') atom(1:1)=' '
      if(atom(1:1).eq.'4') atom(1:1)=' '
      if(atom(1:1).eq.'5') atom(1:1)=' '
      if(atom(1:1).eq.'6') atom(1:1)=' '
      if(atom(1:1).eq.'7') atom(1:1)=' '
      if(atom(1:1).eq.'8') atom(1:1)=' '
      if(atom(1:1).eq.'9') atom(1:1)=' '
      if(atom(1:1).eq.'0') atom(1:1)=' '
      if(atom(1:1).eq.' ') atom=atom(2:2)//' '
      write(lfnxxx,1001) atom,(1.0d1*xs(i,k),k=1,3)
 1001 format(a2,3f12.6)
      endif
    1 continue
      else
      do 2 i=1,nsa
      atom=sgmnam(i)(6:7)
      if(atom(1:1).eq.'1') atom(1:1)=' '
      if(atom(1:1).eq.'2') atom(1:1)=' '
      if(atom(1:1).eq.'3') atom(1:1)=' '
      if(atom(1:1).eq.'4') atom(1:1)=' '
      if(atom(1:1).eq.'5') atom(1:1)=' '
      if(atom(1:1).eq.'6') atom(1:1)=' '
      if(atom(1:1).eq.'7') atom(1:1)=' '
      if(atom(1:1).eq.'8') atom(1:1)=' '
      if(atom(1:1).eq.'9') atom(1:1)=' '
      if(atom(1:1).eq.'0') atom(1:1)=' '
      if(atom(1:1).eq.' ') atom=atom(2:2)//' '
      write(lfnxxx,1001) atom,(1.0d1*xs(i,k),k=1,3)
    2 continue
      endif
c
      if(logw) then
      do 3 i=1,nwm
      if(iwsel(i).ne.0) then
      do 4 j=1,nwa
      atom=wnam(j)(6:7)
      if(atom(1:1).eq.'1') atom(1:1)=' '
      if(atom(1:1).eq.'2') atom(1:1)=' '
      if(atom(1:1).eq.'3') atom(1:1)=' '
      if(atom(1:1).eq.'4') atom(1:1)=' '
      if(atom(1:1).eq.'5') atom(1:1)=' '
      if(atom(1:1).eq.'6') atom(1:1)=' '
      if(atom(1:1).eq.'7') atom(1:1)=' '
      if(atom(1:1).eq.'8') atom(1:1)=' '
      if(atom(1:1).eq.'9') atom(1:1)=' '
      if(atom(1:1).eq.'0') atom(1:1)=' '
      if(atom(1:1).eq.' ') atom=atom(2:2)//' '
      write(lfnxxx,1001) atom,(1.0d1*xw(i,j,k),k=1,3)
    4 continue
      endif
    3 continue
      endif
c
      close(unit=lfnxxx)
c
      length=index(filxxx,' ')-1
      open(unit=lfnxxx,file=filxxx(1:length)//'q',form='formatted',
     + status='unknown')
c
      rewind(lfnxxx)
c
      write(lfnxxx,1002) nsa,1,box
 1002 format(i7,i3,3f12.6)
      if(logs) then
      do 5 i=1,nsa
      if(isel(i).gt.0) then
      atom=sgmnam(i)(6:7)
      if(atom(1:1).eq.'1') atom(1:1)=' '
      if(atom(1:1).eq.'2') atom(1:1)=' '
      if(atom(1:1).eq.'3') atom(1:1)=' '
      if(atom(1:1).eq.'4') atom(1:1)=' '
      if(atom(1:1).eq.'5') atom(1:1)=' '
      if(atom(1:1).eq.'6') atom(1:1)=' '
      if(atom(1:1).eq.'7') atom(1:1)=' '
      if(atom(1:1).eq.'8') atom(1:1)=' '
      if(atom(1:1).eq.'9') atom(1:1)=' '
      if(atom(1:1).eq.'0') atom(1:1)=' '
      if(atom(1:1).eq.' ') atom=atom(2:2)//' '
      write(lfnxxx,1003) atom,(1.0d1*xs(i,k),k=1,3),qs(i)
 1003 format(a2,4f12.6)
      endif
    5 continue
      else
      do 6 i=1,nsa
      atom=sgmnam(i)(6:7)
      if(atom(1:1).eq.'1') atom(1:1)=' '
      if(atom(1:1).eq.'2') atom(1:1)=' '
      if(atom(1:1).eq.'3') atom(1:1)=' '
      if(atom(1:1).eq.'4') atom(1:1)=' '
      if(atom(1:1).eq.'5') atom(1:1)=' '
      if(atom(1:1).eq.'6') atom(1:1)=' '
      if(atom(1:1).eq.'7') atom(1:1)=' '
      if(atom(1:1).eq.'8') atom(1:1)=' '
      if(atom(1:1).eq.'9') atom(1:1)=' '
      if(atom(1:1).eq.'0') atom(1:1)=' '
      if(atom(1:1).eq.' ') atom=atom(2:2)//' '
      write(lfnxxx,1003) atom,(1.0d1*xs(i,k),k=1,3),qs(i)
    6 continue
      endif
c
      close(unit=lfnxxx)
c
      endif
c
      return
      end
