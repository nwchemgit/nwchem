      logical function pre_frgrd(lfnfrg,filfrg,lfnout,lprint,
     + latm,catm,qatm,matm,natm,lbnd,mbnd,nbnd)
c
c $Id: pre_frgrd.F,v 1.4 1998-04-29 18:35:15 d3j191 Exp $
c
c     function to read dimensions from a frg file
c
c     in  : lfnfrg = frg file logical file number
c           filfrg = frg file name
c           lfnout = output file logical file number
c           lprint = integer for print to output
c
c     out : natm = number of atoms
c           nbnd = number of bonds
c
      implicit none
c
      integer lfnfrg,lfnout,matm,natm,mbnd,nbnd
      integer latm(5,matm),lbnd(2,mbnd)
      character*6 catm(2,matm)
      real*8 qatm(2,matm)
      character*255 filfrg
      integer lprint
c
      character*80 card
      integer i,icon(16),j,l,length
c
      natm=0
      nbnd=0
c
      length=index(filfrg,' ')-1
c
      if(lprint.ge.5) then
      write(lfnout,2000) filfrg(1:length)
 2000 format(/,'pre_frgrd  ',a)
      endif
c
      open(unit=lfnfrg,file=filfrg(1:length),form='formatted',
     + status='old',err=9999)
c
      if(lprint.ge.1) write(lfnout,2001) filfrg(1:length)
 2001 format('FRAGMENT FILE ',a)
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
      if(lprint.ge.5) then
      write(lfnout,2002) natm,nbnd,filfrg(1:length)
 2002 format(i20,' atoms',/,i20,' bonds',//,
     + 'FRG CLOSE ',a)
      endif
c
      pre_frgrd=.true.
      return
c
 9999 continue
c
      if(lprint.ge.1) then
      write(lfnout,2003) filfrg(1:length)
 2003 format('FRG ERROR ',a)
      endif
c
      pre_frgrd=.false.
      return
      end
