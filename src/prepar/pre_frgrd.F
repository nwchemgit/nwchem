      logical function pre_frgrd(lfnfrg,filfrg,lfnout,
     + latm,catm,qatm,matm,natm,lbnd,mbnd,nbnd,izmat,zmat,mzmat,nzmat)
c
c $Id: pre_frgrd.F,v 1.6 2000-09-15 22:21:53 d3j191 Exp $
c
c     function to read dimensions from a frg file
c
c     in  : lfnfrg = frg file logical file number
c           filfrg = frg file name
c           lfnout = output file logical file number
c
c     out : natm = number of atoms
c           nbnd = number of bonds
c
      implicit none
c
#include "util.fh"
c
      integer lfnfrg,lfnout,matm,natm,mbnd,nbnd,mzmat,nzmat
      integer latm(5,matm),lbnd(2,mbnd),izmat(mzmat,4)
      character*6 catm(2,matm)
      real*8 qatm(2,matm),zmat(mzmat,3)
      character*255 filfrg
c
      character*80 card
      integer i,icon(16),j,l,length
c
      natm=0
      nbnd=0
c
      length=index(filfrg,' ')-1
c
      if(util_print('where',print_debug)) then
      write(lfnout,2000) filfrg(1:length)
 2000 format(' Attempt reading fragment',t40,a,/)
      endif
c
      open(unit=lfnfrg,file=filfrg(1:length),form='formatted',
     + status='old',err=9999)
c
      if(util_print('sequence',print_high)) then
      write(lfnout,2001) filfrg(1:length)
 2001 format(' Reading fragment',t40,a,/)
      endif
c
    1 continue
c
      read(lfnfrg,1000,end=9999,err=9999) card
 1000 format(a)
c
      if(card(1:1).eq.'#'.or.card(1:1).eq.'$') goto 1
c
      read(card,1001) natm
 1001 format(i5)
c
      do 2 i=1,natm
      read(lfnfrg,1002,end=9999,err=9999) catm(1,i),catm(2,i),
     + (latm(j,i),j=1,4),qatm(1,i),qatm(2,i)
 1002 format(5x,2a6,4i5,2f12.6)
    2 continue
c
    3 continue
      read(lfnfrg,1003,end=999,err=9999) icon
 1003 format(16i5)
c
      if(icon(1).ne.0) then
      do 4 i=2,16
      if(icon(i).ne.0) then
      nbnd=nbnd+1
      lbnd(1,nbnd)=icon(i-1)
      lbnd(2,nbnd)=icon(i)
      endif
    4 continue
      goto 3
      endif
c
      nzmat=0
      read(lfnfrg,1000,end=999,err=9999) card
      if(card(1:8).eq.'z-matrix') then
   10 continue
      read(lfnfrg,1004) (izmat(nzmat+1,i),i=1,4),(zmat(nzmat+1,i),i=1,3)
 1004 format(4i5,3f12.6)
      if(izmat(nzmat+1,1).gt.0) then
      nzmat=nzmat+1
      goto 10
      endif
      endif
c
  999 continue
c
      close(lfnfrg)
c
      if(nbnd.gt.0) then
      do 5 i=1,nbnd
      if(lbnd(1,i).gt.lbnd(2,i)) then
      l=lbnd(1,i)
      lbnd(1,i)=lbnd(2,i)
      lbnd(2,i)=l
      endif
    5 continue
      l=nbnd
      nbnd=1
      do 6 i=2,l
      do 7 j=1,nbnd
      if(lbnd(1,i).eq.lbnd(1,j).and.lbnd(2,i).eq.lbnd(2,j)) goto 6
    7 continue
      nbnd=nbnd+1
    6 continue
      do 8 i=1,nbnd-1
      do 9 j=i+1,nbnd
      if(lbnd(1,i).gt.lbnd(1,j).or.
     + (lbnd(1,i).eq.lbnd(1,j).and.lbnd(2,i).gt.lbnd(2,j))) then
      l=lbnd(1,i)
      lbnd(1,i)=lbnd(1,j)
      lbnd(1,j)=l
      l=lbnd(2,i)
      lbnd(2,i)=lbnd(2,j)
      lbnd(2,j)=l
      endif
    9 continue
    8 continue
      endif
c
      if(util_print('where',print_high)) then
      write(lfnout,2002) natm,nbnd
 2002 format(' Number of atoms',t40,i8,/,
     + ' Number of bonds',t40,i8,/)
      endif
c
      pre_frgrd=.true.
      return
c
 9999 continue
c
      if(util_print('sequence',print_none)) then
      write(lfnout,2003) filfrg(1:length)
 2003 format(' Error reading fragment',t40,a,/)
      endif
c
      pre_frgrd=.false.
      return
      end
