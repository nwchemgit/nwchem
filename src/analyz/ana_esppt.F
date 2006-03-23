      subroutine ana_esppt(isel,xs,xw,qs,qw)
c
      implicit none
c
#include "ana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      integer isel(msa)
      real*8 xs(msa,3),xw(mwm,mwa,3),qs(msa),qw(mwa)
      integer i,j,lfnpdb,ndx
      character*255 filpdb,string
      real*8 espval,point(3)
c
      open(unit=lfnpin,file=filpin(1:index(filpin,' ')-1),
     + form='formatted',status='old')
c
      open(unit=lfnhol,file=filhol(1:index(filhol,' ')-1),
     + form='formatted',status='old')
c
      open(unit=lfnpin,file=filpou(1:index(filpou,' ')-1),
     + form='formatted',status='unknown')
c
      rewind(lfnpin)
      rewind(lfnpou)
c
    1 continue
      read(lfnpin,1000,end=9) string
 1000 format(a)
      if(string(1:12).ne.' From file: ') goto 1
      string=string(13:255)
      filpdb=string(1:index(string,' ')-1)
c
      write(*,'(a,a)') 'PDB FILE TO READ IS ',filpdb
c
      lfnpdb=44
c
      open(unit=lfnpdb,file=filpdb,form='formatted',status='old')
c
      do 2 i=1,nsa
    3 continue
      read(lfnpdb,1001,end=9) string
 1001 format(a)
      if(string(1:6).ne.'ATOM  ') goto 3
      read(string,1002) (xs(i,j),j=1,3)
 1002 format(30x,3f8.3)
      xs(i,1)=1.0d-1*xs(i,1)
      xs(i,2)=1.0d-1*xs(i,2)
      xs(i,3)=1.0d-1*xs(i,3)
    2 continue
c
      close(unit=lfnpdb)
c
    4 continue
c
      read(lfnhol,1000,end=9) string
      if(string(1:12).eq.'LAST-REC-END') then
      read(lfnhol,1000,end=9) string
      goto 4
      endif
      if(string(1:6).ne.'ATOM  ') goto 4
      read(string,1002) point
      point(1)=1.0d-01*point(1)
      point(2)=1.0d-01*point(2)
      point(3)=1.0d-01*point(3)
c
      espval=0.0d0
      do 10 j=1,nsa
      if(isel(j).gt.0) then
      espval=espval+138.9354*qs(i)/sqrt((point(1)-xs(i,1))**2+
     + (point(2)-xs(i,2))**2+(point(3)-xs(i,3))**2)
      endif
   10 continue
      write(*,'(3f12.6,e12.5)') point,espval
c
      goto 1
c
    9 continue
c
      close(unit=lfnpin,status='keep')
      close(unit=lfnpou,status='keep')
      close(unit=lfnhol,status='keep')
c
      return
      end
