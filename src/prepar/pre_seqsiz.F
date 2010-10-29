      logical function pre_seqsiz(lfnout,lfnseq,filseq,nseq,nlnk)
c
c $Id$
c
c     function to scan the sequence file for the number of segments
c
      implicit none
c
#include "util.fh"
c
      integer lfnout,lfnseq,nseq,nlnk
      character*255 filseq
      character*80 card
      integer length,ival,ncrd,ncnt
c
      length=index(filseq,' ')-1
      open(unit=lfnseq,file=filseq(1:length),form='formatted',
     + status='old',err=9999)
c
      nseq=0
      nlnk=0
      ncnt=1
      ncrd=0
c
    1 continue
      read(lfnseq,1000,end=9999,err=9999) card
 1000 format(a)
      if(ncrd.eq.0) ncnt=1
      if(ncrd.gt.0) ncrd=ncrd-1
      if(card(1:4).eq.'stop') goto 9
      if(card(1:6).eq.'repeat') then
      read(card,1002) ncrd,ncnt
 1002 format(6x,2i5)
      goto 1
      endif
      if(card(1:4).eq.'link') then
      nlnk=nlnk+1
      goto 1
      endif
      if(card(1:4).eq.'LINK') then
      nlnk=nlnk+1
      goto 1
      endif
c
      read(card(1:5),1001,err=1) ival
 1001 format(i5)
      nseq=nseq+ncnt
      goto 1
c
    9 continue
c
      close(unit=lfnseq)
c
      pre_seqsiz=.true.
      return
c
 9999 continue
      pre_seqsiz=.false.
      return
      end
