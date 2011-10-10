      subroutine argos_prepare_quick_charge(lfnout,lfnpar,
     + ilo,ihi,
     + latm,catm,qatm,matm,natm,lbnd,qbnd,mbnd,nbnd)
c
      implicit none
c
#include "util.fh"
#include "argos_prepare_common.fh"
c
      integer lfnout,lfnpar
      integer matm,natm,mbnd,nbnd,ilo,ihi
      character*255 filnam
      integer latm(5,matm)
      real*8 qatm(matm)
      character*6 catm(3,matm)
      integer lbnd(2,mbnd)
      real*8 qbnd(mbnd)
      character*100 card
c
      integer i,j,k,l,lt,ld,len
      real*8 q
c
      do 1 i=1,nbnd
      qbnd(i)=0.0d0
    1 continue
c
      do 7 j=1,mdirpar
      do 77 k=1,nfilpar(j)
      write(filnam,'(a,a)') dirpar(j)(1:index(dirpar(j),' ')-1),
     + filpar(j,k)(1:index(filpar(j,k),' '))
      len=index(filnam,' ')-1
      open(unit=lfnpar,file=filnam(1:len),status='old',form='formatted',
     + err=2)
    3 continue
      read(lfnpar,1000,end=4,err=9999) card
      if(card(1:5).ne.'Bonds') goto 3
      if(util_print('files',print_medium)) then
      write(lfnout,2000) filnam(1:len)
 2000 format(' Bond charge definitions',t40,a)
      endif
    5 continue
      read(lfnpar,1000,end=4,err=9999) card
 1000 format(a)
      if(card(1:1).eq.'#'.or.card(1:4).eq.'    ') goto 5
      if(card(1:6).eq.'Angles') goto 4
      read(card(34:45),'(f12.6)') q
      do 6 i=1,nbnd
      if(catm(3,lbnd(1,i))(1:4).eq.card(1:4).and.
     + catm(3,lbnd(2,i))(1:4).eq.card(7:10)) qbnd(i)=q
      if(catm(3,lbnd(2,i))(1:4).eq.card(1:4).and.
     + catm(3,lbnd(1,i))(1:4).eq.card(7:10)) qbnd(i)=-q
    6 continue
      goto 5
    4 continue
      close(unit=lfnpar)
      if(util_print('where',print_debug)) then
      write(lfnout,'(a,a)') filnam(1:len),' closed'
      endif
    2 continue
   77 continue
    7 continue
c
      do 8 i=1,nbnd
      qatm(lbnd(1,i))=qatm(lbnd(1,i))+qbnd(i)
      qatm(lbnd(2,i))=qatm(lbnd(2,i))-qbnd(i)
    8 continue
c
      return
 9999 continue
      call md_abort('Error in argos_prepare_quick_charge',0)
      return
      end
c $Id$
