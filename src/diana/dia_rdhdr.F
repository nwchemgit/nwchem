      subroutine dia_rdhdr(sgmnam)
c
c $Id$
c
      implicit none
c
#include "dia_common.fh"
#include "mafdecls.fh"
c
      integer ftell
      external ftell
c
      character*16 sgmnam(msa)
      character*255 fname
      integer i,lq
c
      integer ib,jb,nsb,nwb
c
      if(me.eq.0) then
c
      iscof=ifrst
      timoff=0.0d0
      time=0.0d0
      timr=0.0d0
      write(cnum,'(i3.3)') iscof
      ioffhd=0
      iofffr=0
c
      fname=filtrj
      lq=index(filtrj,'?')
      if(lq.gt.0) then
      fname=filtrj(1:lq-1)//cnum//filtrj(lq+1:index(filtrj,' ')-1)
      endif
c
      if(fmttrj.eq.'trj') then
      write(*,3333) fname(1:index(fname,' ')-1)
 3333 format(/,' Opening trj file ',a)
      open(unit=lfntrj,file=fname(1:index(fname,' ')-1),
     + status='old',err=9999)
      rewind(lfntrj)
      if(lrdf) call dia_rdfhdr(int_mb(i_rdf))
    1 continue
      read(lfntrj,1000,err=9998,end=9997) card
 1000 format(a)
      if(card(1:6).ne.'header') goto 1
      open(unit=44,file='.header',form='formatted',
     + status='unknown')
      rewind(44)
      read(lfntrj,1003) nwa,nsa,nsb,nwb,ib,nfilf
      if(nfilf.eq.0) nfilf=nframi
      write(44,1003) nwa,nsa,nsb,nwb,ib,nfilf
 1003 format(6i10)
      if(nsa.gt.msa) call md_abort('Error in trj',0)
      if(nwa.gt.0) then
      do 2 i=1,nwa
      read(lfntrj,1004) wnam(i)
      write(44,1004) wnam(i)
 1004 format(a16)
    2 continue
      endif
      if(nsa.gt.0) then
      do 3 i=1,nsa
      read(lfntrj,1005) sgmnam(i)
      write(44,1005) sgmnam(i)
 1005 format(a16)
    3 continue
      endif
      if(nsb+nwb.gt.0) then
      do 4 i=1,nsb+nwb
      read(lfntrj,1001) ib,jb
      write(44,1001) ib,jb
 1001 format(2i8)
    4 continue
      endif
      ioffhd=ftell(lfntrj)
      close(unit=44)
      elseif(fmttrj.eq.'sco'.or.fmttrj.eq.'coo') then
      open(unit=lfntrj,file=fname(1:index(fname,' ')-1),
     + status='old',err=9999)
      write(*,3333) fname(1:index(fname,' ')-1)
   11 continue
      read(lfntrj,1000,err=9998,end=9997) card
      if(card(1:18).eq.'num_solvent_atoms:')then
      read(card(19:26),'(i8)') nwa
      goto 11
      endif
      if(card(1:17).eq.'num_solute_atoms:') then
      read(card(18:25),'(i8)') nsa
      goto 11
      endif
      if(card(1:10).eq.'num_bonds:') then
      read(card(11:19),'(i9)') nsb
      else
      if(card(1:5).ne.'time:') goto 11
      if(.not.ltop) call md_abort('Format error trajectory file',0)
      endif
      rewind(lfntrj)
      nwa=0
      if(.not.ltop) then
      open(unit=44,file='.header',form='formatted',
     + status='unknown')
      rewind(44)
      write(44,1003) nwa,nsa,nsb,0
      endif
      if(nsa.gt.msa) call md_abort('Error in trj',0)
      if(fmttrj.ne.'coo') then
   12 continue
      read(lfntrj,1000,err=9998,end=9997) card
      if(card(1:11).ne.'atom_names:') goto 12
      do 13 i=1,nsa
      if(ltop) then
      read(lfntrj,1000,err=9998,end=9997) card
      else
      read(lfntrj,1005) sgmnam(i)
      write(44,1005) sgmnam(i)
      endif
   13 continue
      if(nsb.gt.0) then
   14 continue
      read(lfntrj,1000,err=9998,end=9997) card
      if(card(1:10).ne.'bond_list:') goto 14
      do 15 i=1,nsb
      read(lfntrj,1001) ib,jb
      if(.not.ltop) write(44,1001) ib,jb
   15 continue
      endif
      endif
      if(.not.ltop) close(unit=44)
      else
      call md_abort('Trajectory file format error',0)
      endif
c
      endif
c
      return
c
 9997 continue
      call md_abort('End of file reading trajectory file',0)
 9998 continue
      call md_abort('Error reading trajectory file',0)
 9999 continue
      call md_abort('Failed to open trajectory file',0)
      return
      end
      subroutine dia_getofffr()
c
      implicit none
c
#include "dia_common.fh"
c
      integer ftell
      external ftell
c
      character*5 string
c
      if(me.eq.0) then
      read(lfntrj,1000) string
    1 continue
      iofffr=ftell(lfntrj)
      read(lfntrj,1000) string
 1000 format(a5)
      if(string.ne.'frame') goto 1
      close(unit=lfntrj)
      endif
c
      return
      end

