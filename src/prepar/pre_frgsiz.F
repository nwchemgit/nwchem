      logical function pre_frgsiz(lfnfrg,filfrg,lfnout,lprint,
     + natm,nbnd)
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
      integer lfnfrg,lfnout,natm,nbnd
      character*80 filfrg
      integer lprint
c
      character*80 card
      integer i,icon(16),length
c
      natm=0
      nbnd=0
c
      length=index(filfrg,' ')-1
c
      if(lprint.ge.5) then
      write(lfnout,2000) filfrg(1:length)
 2000 format(/,'pre_frgsiz  ',a)
      endif
c
      open(unit=lfnfrg,file=filfrg(1:length),form='formatted',
     + status='old',err=9999)
c
      if(lprint.ge.5) write(lfnout,2001) filfrg(1:length)
 2001 format('FRG OPEN  ',a)
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
      read(lfnfrg,1000,end=9999,err=9999) card
    2 continue
c
    3 continue
      read(lfnfrg,1002,end=9,err=9999) icon
 1002 format(16i5)
c
      if(icon(1).ne.0) then
      do 4 i=2,16
      if(icon(i).ne.0) nbnd=nbnd+1
    4 continue
      goto 3
      endif
c
    9 continue
c
      close(lfnfrg)
c
      if(lprint.ge.5) then
      write(lfnout,2002) natm,nbnd,filfrg(1:length)
 2002 format(i20,' atoms',/,i20,' bonds',//,
     + 'FRG CLOSE ',a)
      endif
c
      pre_frgsiz=.true.
      return
c
 9999 continue
c
      if(lprint.ge.1) then
      write(lfnout,2003) filfrg(1:length)
 2003 format('FRG ERROR ',a)
      endif
c
      pre_frgsiz=.false.
      return
      end
