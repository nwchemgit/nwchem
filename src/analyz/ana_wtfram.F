      subroutine ana_wtfram(iunit,fmt,x,sgmnam,w,isel)
c
c $Id: ana_wtfram.F,v 1.8 2000-04-27 15:22:00 d3j191 Exp $
c
      implicit none
c
#include "ana_common.fh"
#include "mafdecls.fh"
c
      real*8 x(msa,3),w(mwm,mwa,3)
      integer isel(msa)
      integer iunit
      character*3 fmt
      character*16 sgmnam(msa)
c
      character*10 today,now
      integer i,j,k
      character*2 elemnt
c
      character*1 cdummy
      integer nsb,ib,jb
c
      if(me.eq.0) then
c
      if(iunit.eq.lfncop) then
      if(ncop.gt.mcopf.and.mcopf.gt.0) then
      ncop=1
      if(fmt.eq.'eci') write(iunit,'(a)') 'End'
      close(unit=iunit)
      icopf=icopf+1
      call ana_wthdr(lfncop,fmtcop,byte_mb(i_snam),isel)
      endif
      endif
c
      if(iunit.eq.lfnsup) then
      if(nsup.gt.msupf.and.msupf.gt.0) then
      nsup=1
      if(fmt.eq.'eci') write(iunit,'(a)') 'End'
      close(unit=iunit)
      isupf=isupf+1
      call ana_wthdr(lfnsup,fmtsup,byte_mb(i_snam),isel)
      endif
      endif
c
      if(fmt.eq.'trj') then
      write(iunit,1000)
 1000 format('frame')
      write(iunit,1001) time,temp,pres
 1001 format(2f12.6,e12.5)
c      write(*,'(a,2f12.6,e12.5)') 'frame ',time,temp,pres
      write(iunit,1002) box
 1002 format(f12.6,36x,f12.6,36x,f12.6)
      write(iunit,1003) lxw,.false.,lsx,.false.,0,0,nsel
 1003 format(4l1,3i10)
      if(lxw) then
      do 1 i=1,nwm
      do 2 j=1,nwa
      write(iunit,1004) (w(i,j,k),k=1,3)
 1004 format(3f8.3)
    2 continue
    1 continue
      endif
      if(lsx) then
      do 3 i=1,nsa
      if(isel(i).ne.0) write(iunit,1005) (x(i,j),j=1,3)
 1005 format(3f8.3)
    3 continue
      endif
      endif
c
      if(fmt.eq.'arc') then
      call swatch(today,now)
      write(iunit,2000) time,today,now
 2000 format('title ',t67,f10.3,/,'!DATE ',2a10)
      if(lsx) then
      do 4 i=1,nsa
      if(sgmnam(i)(6:6).eq.' '.or.
     + sgmnam(i)(6:6).eq.'2'.or.sgmnam(i)(6:6).eq.'3'.or.
     + sgmnam(i)(6:6).eq.'4'.or.sgmnam(i)(6:6).eq.'1') then
      write(iunit,2001) sgmnam(i)(6:9),(1.0d1*x(i,j),j=1,3),
     + sgmnam(i)(1:4),sgmnam(i)(13:16),'?',sgmnam(i)(7:7)//' ',0.0
      else
      write(iunit,2001) sgmnam(i)(6:9),(1.0d1*x(i,j),j=1,3),
     + sgmnam(i)(1:4),sgmnam(i)(13:16),'?',sgmnam(i)(6:7),0.0
      endif
 2001 format(a4,1x,3f15.9,1x,a4,a4,4x,a1,7x,a2,f7.3)
    4 continue
      endif
      if(lxw) then
      do 44 i=1,nwm
      do 45 j=1,nwa
      if(wnam(j)(6:6).eq.' '.or.
     + wnam(j)(6:6).eq.'2'.or.wnam(j)(6:6).eq.'3'.or.
     + wnam(j)(6:6).eq.'4'.or.wnam(j)(6:6).eq.'1') then
      write(iunit,2001) wnam(j)(6:9),(1.0d1*w(i,j,k),k=1,3),
     + wnam(j)(1:4),wnam(j)(13:16),'?',wnam(j)(7:7)//' ',0.0
      else
      write(iunit,2001) wnam(j)(6:9),(1.0d1*w(i,j,k),k=1,3),
     + wnam(j)(1:4),wnam(j)(13:16),'?',wnam(j)(6:7),0.0
      endif
   45 continue
   44 continue
      endif
      write(iunit,2002)
 2002 format('end',/,'end')
      endif
c
      if(fmt.eq.'amb') then
      write(iunit,3000) ((1.0d1*(x(i,j)+0.5d0*box(j)),j=1,3),i=1,nsa)
      write(iunit,3000) (1.0d1*box(j),j=1,3)
 3000 format(10f8.3)
      endif
c
      if(fmt.eq.'bam') then
      write(iunit) ((real(1.0d1*(x(i,j)+0.5d0*box(j))),j=1,3),i=1,nsa),
     + (real(1.0d1*box(j)),j=1,3)
      endif
c
      if(fmt.eq.'eci') then
      write(iunit,5000) time
 5000 format('Time ',f12.4)
      write(iunit,5001) nsa
 5001 format('Geometry{',/,i7,/,' ')
      do 5 i=1,nsa
      elemnt=sgmnam(i)(6:7)
      if(elemnt(1:1).eq.'1'.or.elemnt(1:1).eq.'2'.or.
     + elemnt(1:1).eq.'3'.or.elemnt(1:1).eq.'4') elemnt(1:1)=' '
      write(iunit,5002) elemnt,(1.0d1*x(i,j),j=1,3)
 5002 format(a2,3f12.4)
    5 continue
      write(iunit,5003)
 5003 format('}')
      endif
c
      if(fmt.eq.'mvm') then
      write(iunit,6001) nsa
 6001 format('Geometry{',/,'num_atoms:',i7,/,'atom_list:')
      do 6 i=1,nsa
      elemnt=sgmnam(i)(6:7)
      if(elemnt(1:1).eq.'1'.or.elemnt(1:1).eq.'2'.or.
     + elemnt(1:1).eq.'3'.or.elemnt(1:1).eq.'4') elemnt(1:1)=' '
      write(iunit,6002) elemnt,(1.0d1*x(i,j),j=1,3)
 6002 format(a2,3f12.4)
    6 continue
      open(unit=44,file='.header',form='formatted',
     + status='unknown')
      rewind(44)
      read(44,4401) nwa,nsa,nsb
 4401 format(3i10)
      do 7 i=1,nwa+nsa
      read(44,4402) cdummy
 4402 format(a1)
    7 continue
      write(iunit,6003) nsb
 6003 format('num_bonds:      ',i7,/,'bond_list:')
      do 8 i=1,nsb
      read(44,6004) ib,jb
      write(iunit,6004) ib,jb
 6004 format(2i8)
    8 continue
      close(unit=44)
      write(iunit,6005)
 6005 format('}')
      endif
c
      if(fmt.eq.'frm') then
      write(iunit,7001)
 7001 format(/,'Geometry {')
      do 9 i=1,nsa
      if(isel(i).ne.0) write(iunit,7002) (1.0d1*x(i,j),j=1,3)
 7002 format(3f12.4)
    9 continue
      write(iunit,7003)
 7003 format('}')
      endif
c
      endif
c
      return
c
 9999 call errquit('Oops',0)
      return
      end
