      logical function argos_prep_rdseq(lfnout,lfnseq,filseq,
     + lseq,cseq,mseq,nseq,llnk,clnk,mlnk,nlnk)
c
c $Id$
c
c     function to read the sequence file
c
c     out : lseq(1,mseq) = segment number
c                2       = terminal id
c
      implicit none
c
#include "util.fh"
c
      integer lfnout,lfnseq,mseq,nseq,mlnk,nlnk
      integer lseq(6,mseq),llnk(3,mlnk)
      character*10 cseq(mseq)
      character*4 clnk(2,mlnk)
      character*255 filseq
      character*80 card
      integer i,j,iseq,ilnk,length,imol
c
      length=index(filseq,' ')-1
      open(unit=lfnseq,file=filseq(1:length),form='formatted',
     + status='old',err=9999)
c
      if(util_print('where',print_debug)) then
      write(lfnout,2000) filseq(1:length)
 2000 format('OPENED SEQUENCE FILE ',a)
      endif
c
      iseq=0
      ilnk=0
      imol=1
c
    1 continue
      read(lfnseq,1000,end=9999,err=9999) card
 1000 format(a)
      if(card(1:1).eq.'#') goto 1
c
      if(util_print('where',print_debug)) then
      write(lfnout,2001) card
 2001 format(a)
      endif
c
      if(card(1:5).eq.'break') then
      if(iseq.gt.0) lseq(2,iseq)=1
      goto 1
      endif
      if(card(1:8).eq.'fraction') then
      imol=-(iabs(imol)+1)
      if(iseq.gt.0) lseq(2,iseq)=3
      goto 1
      endif
      if(card(1:4).eq.'link') then
      ilnk=ilnk+1
      if(ilnk.gt.mlnk) call md_abort('increase dimension mlnk',9999)
      read(card,1002,err=9999) (llnk(j,ilnk),clnk(j,ilnk),j=1,2)
 1002 format(5x,i5,a4,i5,a4)
      llnk(3,ilnk)=0
      goto 1
      endif
      if(card(1:4).eq.'LINK') then
      ilnk=ilnk+1
      if(ilnk.gt.mlnk) call md_abort('increase dimension mlnk',9999)
      read(card,1002,err=9999) (llnk(j,ilnk),clnk(j,ilnk),j=1,2)
      llnk(3,ilnk)=1
      goto 1
      endif
      if(card(1:8).eq.'molecule') then
      imol=iabs(imol)+1
      if(iseq.gt.0) lseq(2,iseq)=2
      goto 1
      endif
      if(card(1:3).eq.'pmf') then
      goto 1
      endif
      if(card(1:6).eq.'repeat') then
      goto 1
      endif
      if(card(1:7).eq.'solvent') then
      read(lfnseq,1000,end=9999,err=9999) card
      if(util_print('where',print_debug)) then
      write(lfnout,2001) card
      endif
      if(iseq+1.gt.mseq) call md_abort('increase dimension mseq',9999)
      read(card,1001,err=9999) lseq(1,iseq+1),cseq(iseq+1)
      lseq(2,iseq+1)=-2
      lseq(4,iseq+1)=0
      goto 1
      endif
      if(card(1:4).eq.'stop') goto 9
c
      iseq=iseq+1
      if(iseq.gt.mseq) call md_abort('increase dimension mseq',9999)
      read(card,1001,err=9999) lseq(1,iseq),cseq(iseq)
 1001 format(i5,a10)
      lseq(2,iseq)=0
      lseq(4,iseq)=imol
      if(iseq+1.gt.mseq) call md_abort('increase dimension mseq',9999)
      cseq(iseq+1)='          '
      lseq(1,iseq+1)=0
      lseq(2,iseq+1)=0
      lseq(4,iseq+1)=0
      imol=iabs(imol)
      goto 1
c
    9 continue
c
c     determine termini based on sequence list ( not the links yet! )
c     ----------------------------------------
c
c     -1 : N-terminus
c     -2 : integral molecule
c     >0 : C-terminus
c
      if(iseq.gt.0) then
      if(lseq(2,1).eq.0) lseq(2,1)=-1
      if(lseq(2,1).eq.2) lseq(2,1)=-2
      if(iseq.eq.1) lseq(2,1)=-2
      if(lseq(2,iseq).eq.0) lseq(2,iseq)=1
      do 2 i=2,iseq
      if(lseq(2,i).eq.0.and.lseq(2,i-1).gt.0) lseq(2,i)=-1
      if(lseq(2,i).eq.-1.and.lseq(2,i-1).eq.0) lseq(2,i-1)=1
      if(lseq(2,i).eq.-1.and.lseq(2,i-1).eq.-1) lseq(2,i-1)=-2
      if(lseq(2,i).gt.0.and.lseq(2,i-1).eq.-2) lseq(2,i)=-2
      if(lseq(2,i).gt.0.and.lseq(2,i-1).gt.0) lseq(2,i)=-2
    2 continue
      endif
c
      close(unit=lfnseq)
      nseq=iseq
      nlnk=ilnk
      argos_prep_rdseq=.true.
      return
c
 9999 continue
      argos_prep_rdseq=.false.
      return
      end
