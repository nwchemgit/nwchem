      logical function pre_pdbsiz(iunit,filnam,lfnout,lprint,
     + model,natm,nseq,nssb,nlnk,namslv)
c
c     function to read dimensions from a pdb formatted file
c
c     in  : iunit  = pdb file logical file number
c           filnam = pdb file name
c           lfnout = output file logical file number
c           lprint = boolean for print to output
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
      integer iunit,lfnout,model,natm,nseq,nssb,nchain,nlnk
      character*80 filnam
      logical lprint
      character*3 namslv
c
      character*80 card
      integer length,mdl,nsgm,isgm,i
c
      natm=0
      nseq=0
      nssb=0
      nchain=0
      mdl=0
      isgm=0
      nsgm=0
      nlnk=0
c
      length=index(filnam,' ')-1
c
c      if(lprint) then
c      write(lfnout,2006) filnam(1:length)
c 2006 format(/,'pdb_size  ',a)
c      endif
c
      open(iunit,file=filnam(1:length),form='formatted',
     + status='old',err=9999)
c
c      if(lprint) write(lfnout,2000) filnam(1:length)
c 2000 format('PDB OPEN  ',a)
c
    1 continue
c
      read(iunit,1000,end=999,err=9999) card
 1000 format(a)
c
c     MODEL
c
      if(model.ne.mdl) then
      if(card(1:6).eq.'MODEL ') then
      read(card(7:14),1002) mdl
 1002 format(i8)
c      if(model.eq.mdl.and.lprint) then
c      write(lfnout,2005) model
c 2005 format(10x,'MODEL ',i8)
c      endif
      goto 1
      endif
      endif
c
c     END and ENDMDL
c
      if(card(1:3).eq.'END') goto 9
      if(card(1:6).eq.'ENDMDL') goto 9
c
c     HEADER
c
      if(card(1:6).eq.'HEADER') then
c      if(lprint) write(lfnout,2001) card
c 2001 format(10x,a)
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
      read(card(11:17),1001) nseq
 1001 format(i7)
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
      if(mdl.eq.model) then
      natm=natm+1
      read(card(23:26),1003) i
 1003 format(i4)
      if(i.ne.isgm) then
      isgm=i
      nsgm=nsgm+1
      nchain=nchain+1
      endif
      endif
      goto 1
      endif
c
c     HETATM
c
      if(card(1:6).eq.'HETATM') then
      if(mdl.eq.model) then
      natm=natm+1
      read(card(23:26),1004) i
 1004 format(i4)
      if(i.ne.isgm) then
      isgm=i
      nsgm=nsgm+1
      nchain=nchain+1
      endif
      endif
      if(card(18:20).eq.namslv) goto 9
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
      if(lprint) then
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
c      if(lprint) then
c      write(lfnout,2008)
c 2008 format(10x,'number of segments read from ATOM cards')
c      endif
      endif
c
c      if(lprint) then
c      write(lfnout,2002) natm,nseq,nssb,filnam(1:length)
c 2002 format(i20,' atoms',/,i20,' segments',/,i20,' SS-bonds',//,
c     + 'PDB CLOSE ',a)
c      endif
c
      pre_pdbsiz=.true.
      return
c
 9999 continue
c
      if(lprint) then
      write(lfnout,2003) filnam(1:length)
 2003 format('PDB ERROR ',a)
      endif
c
      pre_pdbsiz=.false.
      return
      end
