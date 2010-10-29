      subroutine ana_wrtpdb(lfnpdb,filpdb,sgmnam,isel,xs,
     + iwsel,xw,logs,logw,ndxw)
c
c $Id$
c
      implicit none
c
#include "ana_common.fh"
c
      character*16 sgmnam(msa)
      integer isel(msa),iwsel(mwm),ndxw(nwm)
      integer lfnpdb
      real*8 xs(msa,3),xw(mwm,mwa,3)
      character*255 filpdb
      logical logs,logw
      character*3 string
c
      integer i,j,k,length,ia,isoff
c
      if(me.eq.0) then
c
      length=index(filpdb,' ')-1
      open(unit=lfnpdb,file=filpdb(1:length),form='formatted',
     + status='unknown')
c
      rewind(lfnpdb)
c
      write(lfnpdb,1000) (1.0d1*box(i),i=1,3),90.0,90.0,90.0
 1000 format(
     + 'HEADER',/,
     + 'TITLE     ',/,
     + 'TITLE    2',/,
     + 'TITLE    3',/,
     + 'REMARK   4 XXXX COMPLIES WITH FORMAT V. 2.1, 25-OCT-1996',/,
     + 'CRYST1',3f9.3,3f7.2)
c
      if(logs) then
      do 1 i=1,nsa
      if(isel(i).gt.0) then
      write(lfnpdb,1001) i,sgmnam(i)(6:9),sgmnam(i)(1:3),
     + sgmnam(i)(13:16),(1.0d1*xs(i,k),k=1,3)
 1001 format('ATOM',i7,1x,a4,1x,a3,2x,a4,4x,3f8.3)
      endif
    1 continue
      else
      do 2 i=1,nsa
      write(lfnpdb,1001) i,sgmnam(i)(6:9),sgmnam(i)(1:3),
     + sgmnam(i)(13:16),(1.0d1*xs(i,k),k=1,3)
    2 continue
      endif
c
      if(logw) then
      write(lfnpdb,1002)
 1002 format('TER')
      read(sgmnam(nsa)(13:16),'(i4)') isoff
      ia=nsa
      string=wnam(1)(1:3)
      if(string.eq.'spc') string='HOH'
      if(string.eq.'SPC') string='HOH'
      do 3 i=1,nwm
      if(iwsel(i).ne.0) then
      do 4 j=1,nwa
      ia=ia+1
      write(lfnpdb,1003) ia,wnam(j)(6:9),string,isoff+i,
     + (1.0d1*xw(i,j,k),k=1,3)
 1003 format('ATOM',i7,1x,a4,1x,a3,i6,4x,3f8.3)
    4 continue
      endif
    3 continue
      endif
c
      write(lfnpdb,1004)
 1004 format('END')
c
      close(unit=lfnpdb)
c
      endif
c
      return
      end
