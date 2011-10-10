      subroutine argos_diana_esppt(isel,xs,xw,qs,qw)
c
      implicit none
c
#include "argos_diana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      integer isel(msa)
      real*8 xs(msa,3),xw(mwm,mwa,3),qs(msa),qw(mwa)
      real*8 val(1000),zval(1000)
      integer nval(1000)
      integer i,j,lfnpdb,last,number
      character*255 filpdb,string
      real*8 espval,point(3)
c
      write(*,'(a,a)') 'filpin=',filpin(1:index(filpin,' ')-1)
      open(unit=lfnpin,file=filpin(1:index(filpin,' ')-1),
     + form='formatted',status='old')
c
      write(*,'(a,a)') 'filhol=',filhol(1:index(filhol,' ')-1)
      open(unit=lfnhol,file=filhol(1:index(filhol,' ')-1),
     + form='formatted',status='old')
c
      write(*,'(a,a)') 'filpou=',filpou(1:index(filpou,' ')-1)
      open(unit=lfnpou,file=filpou(1:index(filpou,' ')-1),
     + form='formatted',status='unknown')
c
      write(*,'(a,a)') 'filavg=',filavg(1:index(filavg,' ')-1)
      open(unit=lfnavg,file=filavg(1:index(filavg,' ')-1),
     + form='formatted',status='unknown')
c
      print*,'lfnpin=',lfnpin
      print*,'lfnhol=',lfnhol
      print*,'lfnpou=',lfnpou
      print*,'lfnavg=',lfnavg
c
      rewind(lfnpin)
      rewind(lfnhol)
      rewind(lfnpou)
      rewind(lfnavg)
c
      last=0
      number=0
c
    1 continue
      last=0
      read(lfnpin,1000,end=9) string
 1000 format(a)
      if(string(1:12).ne.' From file: ') goto 1
      string=string(13:255)
      if(index(string,'pdb').le.0) goto 1
      filpdb=string(1:index(string,' ')-1)
c
      write(*,'(a,a)') 'Opening PDB file ',filpdb
c
      lfnpdb=44
c
      open(unit=lfnpdb,file=filpdb,form='formatted',status='old')
c
      do 2 i=1,nsa
    3 continue
      read(lfnpdb,1001,end=999) string
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
      write(*,'(a,a)') 'Read and closed PDB file ',filpdb
c
    4 continue
c
      read(lfnhol,1000,end=998) string
      if(string(23:26).eq.'-888') then
      read(lfnhol,1000,end=998) string
      if(last.lt.0) goto 1
      goto 4
      endif
      if(string(1:6).ne.'ATOM  ') goto 4
      read(string,2002) last,point
 2002 format(22x,i4,4x,3f8.3)
      point(1)=1.0d-01*point(1)
      point(2)=1.0d-01*point(2)
      point(3)=1.0d-01*point(3)
c
      espval=0.0d0
      do 10 i=1,nsa
      if(isel(i).gt.0) then
      espval=espval+138.9354*qs(i)/sqrt((point(1)-xs(i,1))**2+
     + (point(2)-xs(i,2))**2+(point(3)-xs(i,3))**2)
      endif
   10 continue
      write(*,'(3f12.6,1x,e12.5,i5)') point,espval,last
      write(lfnpou,3001) point(3),espval
 3001 format(2f12.6)
      do 5 i=1,number
      if(zval(i).eq.point(3)) then
      nval(i)=nval(i)+1
      val(i)=val(i)+espval
      goto 4
      endif
    5 continue
      number=number+1
      nval(number)=1
      zval(number)=point(3)
      val(i)=espval
c
      goto 4
c
    9 continue
c
      close(unit=lfnpin,status='keep')
      close(unit=lfnpou,status='keep')
      close(unit=lfnhol,status='keep')
c
      do 6 i=1,number
      val(i)=val(i)/dble(nval(i))
    6 continue
      do 7 i=1,number
      nval(i)=i
    7 continue
      do 8 i=1,number-1
      do 12 j=i+1,number
      if(zval(nval(i)).gt.zval(nval(j))) then
      last=nval(i)
      nval(i)=nval(j)
      nval(j)=last
      endif
   12 continue
    8 continue
c
      do 11 i=1,number
      write(lfnavg,4000) zval(nval(i)),val(nval(i))
 4000 format(2f12.6)
   11 continue
c
      close(unit=lfnavg)
c
      return
c
  998 continue
      write(*,'(a)') 'Error reading hole output file ',filhol
      return
  999 continue
      write(*,'(a)') 'Error reading pdb file ',filpdb
      return
      end
c $Id$
