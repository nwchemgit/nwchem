      subroutine ana_wthdr(iunit,fmt,sgmnam,tag)
c
c $Id: ana_wthdr.F,v 1.2 1999-03-15 21:27:53 d3j191 Exp $
c
      implicit none
c
#include "ana_common.fh"
c
      character*16 sgmnam(msa)
      integer iunit
      character*3 fmt
      character*4 cnum
      character*255 fname
      character*24 tag(msa,2)
c
      character*80 label
c
      integer i,ib,j,jb,nsb,lq
      logical binary
c
      character*1 cdummy
      character*2 elemnt
      character*18 name
c
      if(me.eq.0) then
c
      binary=fmt(1:1).eq.'b' 
c
      if(iunit.eq.lfncop) then
      fname=filcop
      if(mcopf.gt.0) then
      lq=index(filcop,'.')
      write(cnum,'(a1,i3.3)') '_',icopf
      fname=filcop(1:lq-1)//cnum//filcop(lq:index(filcop,' ')-1)
      endif
      if(binary) then
      open(unit=lfncop,file=fname(1:index(fname,' ')-1),
     + form='unformatted',status='unknown',err=9999)
      else
      open(unit=lfncop,file=fname(1:index(fname,' ')-1),
     + form='formatted',status='unknown',err=9999)
      endif
      write(*,'(a,a)') 'Opened copy file ',fname(1:index(fname,' ')-1)
      endif
c
      if(iunit.eq.lfnsup) then
      fname=filsup
      if(msupf.gt.0) then
      lq=index(filsup,'.')
      write(cnum,'(a1,i3.3)') '_',isupf
      fname=filsup(1:lq-1)//cnum//filsup(lq:index(filsup,' ')-1)
      endif
      if(binary) then
      open(unit=lfnsup,file=fname(1:index(fname,' ')-1),
     + form='unformatted',status='unknown',err=9999)
      else
      open(unit=lfnsup,file=fname(1:index(fname,' ')-1),
     + form='formatted',status='unknown',err=9999)
      endif
      write(*,'(a,a)') 'Opened super file ',fname(1:index(fname,' ')-1)
      endif
c
      if(fmt.eq.'sco') then
      write(iunit,1000)
 1000 format('title: NWChem:analysis')
      write(iunit,1001)
 1001 format('type: molecule')
      write(iunit,1002)
 1002 format('cell_units: angstroms')
      write(iunit,1003)
 1003 format('atom_info: symbol cart')
      write(iunit,1004) 0
 1004 format('num_solvent_atoms:',i8)
      write(iunit,1005) nsa
 1005 format('num_solute_atoms:',i8)
      write(iunit,1006) nsa
 1006 format('num_atoms:',i8)
      write(iunit,1007)
 1007 format('atom_names:')
      do 1 i=1,nsa
      write(iunit,1008) sgmnam(i)(6:9),sgmnam(i)(1:5),sgmnam(i)(11:16)
 1008 format(a4,a5,a6)
    1 continue
c      open(unit=44,file='bond.list',form='formatted',
c     + status='unknown')
c      rewind(44)
c      read(44,4401) nsb
c 4401 format(11x,i8)
c      write(iunit,4005) nsb
c 4005 format('num_bonds:      ',i7,/,'bond_list')
c      do 2 i=1,nsb
c      read(44,4402) ib,jb
c 4402 format(2i8)
c      write(iunit,4006) ib,jb
c 4006 format(2i8)
c    2 continue
c      close(unit=44)
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
      if(fmt.eq.'bam') then
      label='AMBER binary trajectory file'
      write(iunit) label
      endif
c
      if(fmt.eq.'mvm') then
      label='ecce mvm trajectory file'
      write(iunit,4000) label
 4000 format('# ',a)
c      write(iunit,4001)
c 4001 format('type:           molecule')
c      write(iunit,4002)
c 4002 format('centering:      0')
c      write(iunit,4003) nsa
c 4003 format('num_atoms:      ',i7)
c      write(iunit,4004)
c 4004 format('atom_info:      symbol cart')
c      open(unit=44,file='bond.list',form='formatted',
c     + status='unknown')
c      rewind(44)
c      read(44,4401) nsb
c      write(iunit,4005) nsb
c      do 3 i=1,nsb
c      read(44,4402) ib,jb
c      write(iunit,4006) ib,jb
c    3 continue
c      close(unit=44)
      endif
c
      if(fmt.eq.'eci') then
      label='EcceImport 1.1'
      write(iunit,5000) label
 5000 format(a,/,' ')
      endif
c
      if(fmt.eq.'frm') then
      write(iunit,6001) nsa
 6001 format('Atoms {',/,i8)
      do 2 i=1,nsa
      elemnt=sgmnam(i)(6:7)
      name=sgmnam(i)(11:16)//':'//sgmnam(i)(1:5)//':'//sgmnam(i)(6:10)
      do 3 j=8,14
      if(name(j:j).eq.' ') name(j:j)='_'
    3 continue
      if(elemnt(1:1).eq.'1'.or.elemnt(1:1).eq.'2'.or.
     + elemnt(1:1).eq.'3'.or.elemnt(1:1).eq.'4') elemnt(1:1)=' '
      write(iunit,6002) elemnt,name,tag(i,1)(1:index(tag(i,1),' ')-1),
     + tag(i,2)(1:index(tag(i,2),' ')-1)
 6002 format(a2,1x,a18,1x,a,1x,a)
    2 continue
      open(unit=44,file='bond.list',form='formatted',
     + status='unknown')
      rewind(44)
      read(44,4401) nsb
 4401 format(11x,i8)
      read(44,4402) cdummy
 4402 format(a1)
      write(iunit,6005) nsb
 6005 format('}',/,'Bonds {',/,i8)
      do 4 i=1,nsb
      read(44,6006) ib,jb
      write(iunit,6006) ib,jb
 6006 format(2i8)
    4 continue
      close(unit=44)
      write(iunit,6007)
 6007 format('}')
      endif
c
      endif
c
      return
c
 9999 continue
      call errquit('Failed to open file',iunit)
      return
      end
