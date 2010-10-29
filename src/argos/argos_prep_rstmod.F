      logical function argos_prep_rstmod(lseq,mseq,nseq,
     + lfnmod,filmod,lfnout)
c
c $Id$
c
      implicit none
c
      integer mseq,nseq,lfnmod,lfnout
      integer lseq(mseq)
      character*80 card
      character*255 filmod
c
      integer mset,isgm
c
      open(unit=lfnmod,file=filmod(1:index(filmod,' ')-1),
     + form='formatted',status='old',err=999)
      rewind(unit=lfnmod)
    1 continue
      read(lfnmod,1000,end=99) card
 1000 format(a)
c
      if(card(1:7).eq.'segment') then
      read(card(9:13),'(i5)') isgm
      read(card(18:19),'(i2)') mset
c
      if(card(21:23).eq.'hop') then
cc      if(util_print('where',print_default)) then
      write(lfnout,2019) isgm,mset
 2019 format(' modify segment ',i5,' protonation state ',i5)
cc      endif
      lseq(isgm)=mset
      endif
      endif
      goto 1
 99   continue
c
      close(unit=lfnmod,status='keep')
c
 999  continue
c
      argos_prep_rstmod=.true.
c
      return
      end
