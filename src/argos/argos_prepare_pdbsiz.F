      logical function argos_prepare_pdbsiz(iunit,filnam,lfnout,
     + model,natm,nseq,nssb,nlnk,namslv,numslv,altloc,chain)
c
c $Id$
c
c     function to read dimensions from a pdb formatted file
c
c     in  : iunit  = pdb file logical file number
c           filnam = pdb file name
c           lfnout = output file logical file number
c           model  = model number
c
c     out : natm = number of atoms
c           nseq = number of segments
c           nssb = number of SS bonds
c           nlnk = number of other links
c
c     return value:
c          .false. : specified file does not exist or an error
c                    occurred during reading
c          .true.  : specified file read succesfully
c
      implicit none
c
#include "util.fh"
c
      integer iunit,lfnout,model,natm,nseq,nssb,nchain,nlnk
      character*255 filnam
      integer numslv,lrgpdb
      character*3 namslv
      character*1 ,altloc,chain,pchain
c
      character*80 card
      integer length,mdl,nsgm,isgm,i,ioff
      logical lslvnt
c
      natm=0
      nseq=0
      nssb=0
      nchain=0
      mdl=0
      isgm=0
      nsgm=0
      nlnk=0
      ioff=0
      lslvnt=.false.
      numslv=0
      pchain=' '
      lrgpdb=0
c
      length=index(filnam,' ')-1
c
      if(util_print('where',print_debug)) then
      write(lfnout,2006) filnam(1:length)
 2006 format(/,'pdb_size  ',a)
      endif
c
      open(iunit,file=filnam(1:length),form='formatted',
     + status='old',err=9999)
c
      if(util_print('where',print_debug)) then
      write(lfnout,2000) filnam(1:length)
 2000 format(' PDB OPEN  ',a)
      endif
c
    1 continue
c
      read(iunit,1000,end=999,err=9999) card
 1000 format(a)
c
c     lrgpdb option
c
      if(card(1:6).eq.'LRGPDB') then
      lrgpdb=1
      goto 1
      endif
c
c     MODEL
c
      if(model.ne.mdl) then
      if(card(1:6).eq.'MODEL ') then
      read(card(7:14),1002) mdl
 1002 format(i8)
      if(model.eq.mdl.and.util_print('where',print_high)) then
      write(lfnout,2005) model
 2005 format(10x,'MODEL ',i8)
      endif
      goto 1
      endif
      endif
c
c     END and ENDMDL
c
c
      if(card(1:6).eq.'ENDMDL') then
      if(model.eq.mdl) goto 9
      goto 1
      endif
      if(card(1:6).eq.'END   ') goto 9
c
c     HEADER
c
      if(card(1:6).eq.'HEADER') then
      if(util_print('where',print_debug)) then
      write(lfnout,2001) card
 2001 format(10x,a)
      endif
      goto 1
      endif
c
c     FORMUL
c
      if(card(1:6).eq.'FORMUL'.and.card(19:19).eq.'*') then
      namslv=card(13:15)
      endif
c
c     SEQRES
c
      if(card(1:6).eq.'SEQRES') then
      read(card(14:17),1001) nseq
 1001 format(i4)
      goto 1
      endif
c
c     SSBOND
c
      if(card(1:6).eq.'SSBOND') then
      nssb=nssb+1
      goto 1
      endif
c
c     LINK
c
      if(card(1:6).eq.'LINK  ') then
      nlnk=nlnk+1
      goto 1
      endif
c
c     ATOM
c
      if(card(1:4).eq.'ATOM') then
      if(card(17:17).ne.' '.and.
     + card(17:17).ne.altloc.and.altloc.ne.'*') goto 1
      if(card(22:22).ne.' '.and.
     + card(22:22).ne.chain.and.chain.ne.'*') goto 1
      if(mdl.eq.model) then
      if(card(18:20).eq.namslv) lslvnt=.true.
      if(.not.lslvnt) natm=natm+1
      if(lrgpdb.eq.0) then
      read(card(23:26),1003) i
 1003 format(i4)
      else
      read(card(21:26),1023) i
 1023 format(i6)
      endif
      if(i.ne.isgm.or.card(22:22).ne.pchain) then
      if(i+ioff.lt.isgm) ioff=isgm+1-i
      isgm=i+ioff
      if(lslvnt) then
      numslv=numslv+1
      else
      nsgm=nsgm+1
      nchain=nchain+1
      endif
      endif
      pchain=card(22:22)
      endif
      goto 1
      endif
c
c     HETATM
c
      if(card(1:6).eq.'HETATM') then
      if(card(17:17).ne.' '.and.
     + card(17:17).ne.altloc.and.altloc.ne.'*') goto 1
      if(card(22:22).ne.' '.and.
     + card(22:22).ne.chain.and.chain.ne.'*') goto 1
      if(mdl.eq.model) then
      if(card(18:20).eq.namslv) lslvnt=.true.
      if(.not.lslvnt) natm=natm+1
      read(card(23:26),1004) i
 1004 format(i4)
      if(i.ne.isgm.or.card(22:22).ne.pchain) then
      isgm=i
      if(lslvnt) then
      numslv=numslv+1
      else
      nsgm=nsgm+1
      nchain=nchain+1
      endif
      endif
      pchain=card(22:22)
      endif
      goto 1
      endif
c
c     TER
c
      if(card(1:6).eq.'TER   ') then
      nchain=0
      endif
c
      goto 1
c
  999 continue
c
      if(util_print('where',print_high)) then
      write(lfnout,2004)
 2004 format('PDB EOF encountered')
      endif
c
    9 continue
c
      close(iunit)
c
      if(nsgm.gt.nseq) then
      nseq=nsgm
      if(util_print('where',print_debug)) then
      write(lfnout,2008)
 2008 format(10x,'number of segments read from ATOM cards')
      endif
      endif
c
      if(util_print('where',print_high)) then
      write(lfnout,2002) filnam(1:length),natm,nseq,nssb,numslv
 2002 format(' PDB file',t40,a,/,
     + ' Number of atoms',t40,i10,/,
     + ' Number of segments',t40,i10,/,
     + ' Number of SS-bonds',t40,i10,/,
     + ' Number of solvent molecules',t40,i10,/)
      endif
c
      argos_prepare_pdbsiz=.true.
      return
c
 9999 continue
c
      if(util_print('where',print_none)) then
      write(lfnout,2003) filnam(1:length)
 2003 format('PDB ERROR ',a)
      endif
c
      argos_prepare_pdbsiz=.false.
      return
      end
