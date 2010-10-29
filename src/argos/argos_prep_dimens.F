      logical function argos_prep_dimens(lfnout,lfnsgm,
     + lseq,cseq,mseq,nseq,natm,nbnd,nang,ndih,nimp,nparms)
c
c $Id$
c
c     function to determine total system array dimensions
c
      implicit none
c
#include "util.fh"
c
      logical argos_prep_sgmfil,argos_prep_sgmsiz
      external argos_prep_sgmfil,argos_prep_sgmsiz
c
      integer lfnout,lfnsgm,mseq,nseq,natm,nbnd,nang,ndih,nimp,nparms
      integer lseq(6,mseq)
      character*10 cseq(mseq)
c
      character*255 filsgm
      integer i,nsatm,nsbnd,nsang,nsdih,nsimp
c
      natm=0
      nbnd=0
      nang=0
      ndih=0
      nimp=0
c
      argos_prep_dimens=.true.
c
c     loop over sequence list
c     -----------------------
c
      do 1 i=1,nseq+1
      if(i.lt.nseq+1.or.cseq(i)(1:1).ne.' ') then
c
      if(util_print('where',print_debug)) then
      write(lfnout,1000) cseq(i)
 1000 format('SEARCHING FOR SEGMENT FILE FOR ',a10)
      endif
c
c     find segment file
c     -----------------
c
      if(argos_prep_sgmfil(cseq(i),lfnsgm,filsgm)) then
c
      if(util_print('where',print_debug)) then
      write(lfnout,1001) filsgm(1:index(filsgm,' ')-1)
 1001 format('FOUND SEGMENT FILE ',a)
      endif
c
c     find segment dimensions
c     -----------------------
c
      if(.not.argos_prep_sgmsiz(lfnsgm,filsgm,nsatm,nsbnd,nsang,
     + nsdih,nsimp,nparms))
     + call md_abort('argos_prep_sgmsiz failed',9999)
c
      natm=natm+nsatm
      nbnd=nbnd+nsbnd
      nang=nang+nsang
      ndih=ndih+nsdih
      nimp=nimp+nsimp
c
      else
      if(util_print('topology',print_none)) then
      write(lfnout,1002) cseq(i)
 1002 format('A SEGMENT FILE COULD NOT BE FOUND FOR ',a10)
      endif
      argos_prep_dimens=.false.
      endif
      endif
    1 continue
c
      return
      end
