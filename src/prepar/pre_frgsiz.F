      logical function pre_frgsiz(lfn,fil,lfnout,
     + natm,nbnd,nzmat,nparms,ipardef)
c
c $Id$
c
c     function to read dimensions from a frg file
c
c     in  : lfn = frg file logical file number
c           fil = frg file name
c           lfnout = output file logical file number
c
c     out : natm = number of atoms
c           nbnd = number of bonds
c
      implicit none
c
#include "util.fh"
c
      integer lfn,lfnout,natm,nbnd,nzmat
      character*255 fil
c
      character*80 card
      integer i,icon(16),length
c
      integer nparms,ipardef
c
      natm=0
      nbnd=0
      nzmat=0
c
      length=index(fil,' ')-1
c
      if(util_print('where',print_debug)) then
      write(lfnout,2000) fil(1:length)
 2000 format(/,'pre_frgsiz  ',a)
      endif
c
      open(unit=lfn,file=fil(1:length),form='formatted',
     + status='old',err=9999)
c
      if(util_print('where',print_debug)) then
      write(lfnout,2001) fil(1:length)
 2001 format('FRG OPEN  ',a)
      endif
c
    1 continue
c
      read(lfn,1000,end=9999,err=9999) card
 1000 format(a)
c
      if(card(1:1).eq.'#'.or.card(1:1).eq.'$') goto 1
c
      read(card,1001) natm,nparms,ipardef,nzmat
 1001 format(4i5)
c
      do 2 i=1,nparms*(natm+1)
      read(lfn,1000,end=9999,err=9999) card
    2 continue
c
    3 continue
      read(lfn,1002,end=9,err=9999) icon
 1002 format(16i5)
c
      if(icon(1).ne.0) then
      do 4 i=2,16
      if(icon(i).ne.0) nbnd=nbnd+1
    4 continue
      goto 3
      endif
c
      read(lfn,1000,end=9,err=9999) card
      if(card(1:8).eq.'z-matrix') then
    5 continue
      read(lfn,1003) i
 1003 format(i5)
      if(i.gt.0) then
      nzmat=nzmat+1
      goto 5
      endif
      endif
c
    9 continue
c
      close(lfn)
c
      if(util_print('where',print_debug)) then
      write(lfnout,2002) natm,nbnd,fil(1:length)
 2002 format(i20,' atoms',/,i20,' bonds',//,
     + 'FRG CLOSE ',a)
      endif
c
      pre_frgsiz=.true.
      return
c
 9999 continue
c
      if(util_print('where',print_none)) then
      write(lfnout,2003) fil(1:length)
 2003 format('FRG ERROR ',a)
      endif
c
      pre_frgsiz=.false.
      return
      end
