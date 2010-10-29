      subroutine ana_wthdr(iunit,fmt,sgmnam,tag,isel,logw)
c
c $Id$
c
      implicit none
c
#include "ana_common.fh"
c
      character*16 cwa,sgmnam(msa)
      integer isel(msa)
      integer iunit
      character*3 fmt
      character*4 cdnum
      character*255 fname
      character*24 tag(msa,2)
      logical logw
c
      character*80 label
c
      integer i,ib,j,jb,nsb,lq,nwb
      logical binary
      character*4 dcdhdr
      character*80 dcdtit(10)
      integer*4 icntrl(20),numtit
c
      character*1 cdummy
      character*2 elemnt
      character*18 name
c
      if(me.eq.0) then
c
      if(fmt.eq.'pov') then
      if(iunit.eq.lfncop) filpov=filcop
      if(iunit.eq.lfnsup) filpov=filsup
      npov=0
      return
      endif
c
      binary=fmt(1:1).eq.'b'.or.fmt.eq.'dcd'
c
      if(iunit.eq.lfncop) then
      fname=filcop
      if(mcopf.gt.0) then
      lq=index(filcop,'.')
      write(cdnum,'(a1,i3.3)') '_',icopf
      fname=filcop(1:lq-1)//cdnum//filcop(lq:index(filcop,' ')-1)
      endif
      if(binary) then
      open(unit=lfncop,file=fname(1:index(fname,' ')-1),
     + form='unformatted',status='unknown',err=9999)
      else
      open(unit=lfncop,file=fname(1:index(fname,' ')-1),
     + form='formatted',status='unknown',err=9999)
      endif
      write(*,3333) fname(1:index(fname,' ')-1)
 3333 format(/,' Opening copy file ',a,/)
      endif
c
      if(iunit.eq.lfnsup) then
      fname=filsup
      if(msupf.gt.0) then
      lq=index(filsup,'.')
      write(cdnum,'(a1,i3.3)') '_',isupf
      fname=filsup(1:lq-1)//cdnum//filsup(lq:index(filsup,' ')-1)
      endif
      if(binary) then
      open(unit=lfnsup,file=fname(1:index(fname,' ')-1),
     + form='unformatted',status='unknown',err=9999)
      else
      open(unit=lfnsup,file=fname(1:index(fname,' ')-1),
     + form='formatted',status='unknown',err=9999)
      endif
      write(*,3333) fname(1:index(fname,' ')-1)
 3334 format(/,' Opening super file ',a,/)
      endif
c
      if(fmt.eq.'trj') then
      write(iunit,1000)
 1000 format('header')
      open(unit=44,file='.header',form='formatted',
     + status='unknown')
      rewind(44)
      read(44,4401) nwa,nsa,nsb,nwb
c      if(nsa.eq.nsel) then
      if(logw) then
      write(iunit,4401) nwa,nsa,nsb,nwb,nwsel,mcopf
      else
      write(iunit,4401) nwa,nsa,nsb,nwb,0,mcopf
      endif
c      else
c      write(iunit,4401) nwa,nsel,0,0,0
c      endif
 4401 format(6i10)
      do 1 i=1,nwa
      read(44,4402) cwa
 4402 format(a16,i10,i5)
      write(iunit,4402) cwa,i
    1 continue
      do 2 i=1,nsa
      read(44,4402) cwa
      write(iunit,4402) cwa,i,isel(i)
    2 continue
c      if(nsa.eq.nsel) then
      do 3 i=1,nsb+nwb
      read(44,4403) ib,jb
      write(iunit,4403) ib,jb
 4403 format(2i8)
    3 continue
c      endif
      close(unit=44)
      endif
c
      if(fmt.eq.'arc') then
      write(iunit,2000)
 2000 format('!BIOSYM archive 3',/,'PBC=OFF')
      endif
c
      if(fmt.eq.'amb') then
      label='AMBER trajectory file'
      write(iunit,3000) label
 3000 format(a80)
      endif
c
      if(fmt.eq.'crd') then
      label='AMBER trajectory file'
      write(iunit,3000) label
      endif
c
      if(fmt.eq.'bam') then
      label='AMBER binary trajectory file'
      write(iunit) label
      endif
c
      if(fmt.eq.'mvm') then
      label='ecce mvm trajectory file'
      write(iunit,4000) label
 4000 format('# ',a)
      write(iunit,4001)
 4001 format('type:           molecule')
      write(iunit,4002)
 4002 format('centering:      0')
      write(iunit,4003) nsa
 4003 format('num_atoms:      ',i7)
      write(iunit,4004)
 4004 format('atom_info:      symbol cart')
      open(unit=44,file='.header',form='formatted',
     + status='unknown')
      rewind(44)
      read(44,4401) nwa,nsa,nsb,nwb
      do 4 i=1,nwa+nsa
      read(44,4402) cdummy
    4 continue
      write(iunit,6003) nsb
 6003 format('num_bonds:      ',i7,/,'bond_list:')
      do 5 i=1,nsb
      read(44,4403) ib,jb
      write(iunit,4403) ib,jb
    5 continue
      close(unit=44)
      endif
c
      if(fmt.eq.'eci') then
      label='EcceImport 1.1'
      write(iunit,5000) label
 5000 format(a,/,' ')
      endif
c
      if(fmt.eq.'frm') then
      write(iunit,6001) nsel
 6001 format('Atoms {',/,i8)
      do 6 i=1,nsa
      if(isel(i).ne.0) then
      elemnt=sgmnam(i)(6:7)
      name=sgmnam(i)(11:16)//':'//sgmnam(i)(1:5)//':'//sgmnam(i)(6:10)
      do 7 j=8,14
      if(name(j:j).eq.' ') name(j:j)='_'
    7 continue
      if(elemnt(1:1).eq.'1'.or.elemnt(1:1).eq.'2'.or.
     + elemnt(1:1).eq.'3'.or.elemnt(1:1).eq.'4') elemnt(1:1)=' '
      write(iunit,6002) elemnt,name,tag(i,1)(1:index(tag(i,1),' ')-1),
     + tag(i,2)(1:index(tag(i,2),' ')-1)
 6002 format(a2,1x,a18,1x,a,1x,a)
      endif
    6 continue
      open(unit=44,file='.header',form='formatted',
     + status='unknown')
      rewind(44)
      read(44,4401) nwa,nsa,nsb,nwb
      do 8 i=1,nwa+nsa
      read(44,4402) cdummy
    8 continue
      if(nsel.eq.nsa) then
      write(iunit,6005) nsb
      else
      write(iunit,6005) 0
      endif
 6005 format('}',/,'Bonds {',/,i8)
      if(nsel.eq.nsa) then
      do 9 i=1,nsb
      read(44,6006) ib,jb
      write(iunit,6006) ib,jb
 6006 format(2i8)
    9 continue
      endif
      close(unit=44)
      write(iunit,6007)
 6007 format('}')
      endif
c
      if(fmt.eq.'dcd') then
      dcdhdr='CORD'
      do 21 i=1,20
      icntrl(i)=0
   21 continue
c     number of frames
      icntrl(1)=(ifrto-ifrfr+1)/ifrsk
      icntrl(2)=ifrfr
      icntrl(3)=1
      icntrl(4)=ifrto-ifrfr+1
      icntrl(8)=1
      icntrl(10)=1026003170
      icntrl(20)=22
      write(iunit) dcdhdr,icntrl
      dcdtit(1)='charmm dcd formatted trajectory file'
      numtit=1
      write(iunit) numtit,dcdtit(1)
      ldcd=.false.
      endif
c
      endif
c
      return
c
 9999 continue
      call md_abort('Failed to open file',iunit)
      return
      end
