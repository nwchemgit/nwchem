      logical function pre_pdbseq(iunit,filnam,lfnout,lprint,
     + model,namslv,lseq,cseq,mseq,nseq,latm,catm,xatm,matm,natm,
     + lssb,mssb,nssb,llnk,clnk,mlnk,nlnk)
c
c     function to read sequence from a pdb formatted file
c
c     in  : iunit       = pdb logical file number
c           filnam      = pdb file name
c           lfnout      = output file logical file number
c           lprint      = integer for print to output
c           model       = model number
c           mseq        = dimension of the sequence list
c           nseq        = expected length of the sequence list
c
c     out : lseq(1,*)   = segment numbers
c           cseq(1,*)   = segment names
c           nseq        = actual length of the sequence list
c           latm(1,*)   = index to sequence list
c           catm(1,*)   = atom name
c           xatm(1:3,*) = atomic coordinates
c
      implicit none
c
      integer iunit,lfnout,model,mseq,nseq,matm,natm,mssb,nssb,mlnk,nlnk
      integer lseq(4,mseq),latm(5,matm),lssb(2,mssb),llnk(2,mlnk)
      character*6 catm(3,matm)
      character*5 cseq(2,mseq)
      character*4 clnk(2,mlnk)
      character*3 namslv
      real*8 xatm(3,matm)
      character*255 filnam
      integer lprint
c
      character*80 card
      integer length,mdl,iseq,isgm,i,j,k,nchain,iatm
c
      mdl=0
      iseq=0
      isgm=0
      nchain=0
      iatm=0
      nssb=0
c
      length=index(filnam,' ')-1
c
      if(lprint.ge.5) then
      write(lfnout,2006) filnam(1:length)
 2006 format(/,'pre_pdbseq ',a)
      endif
c
c     initialize the lists
c
      do 10 i=1,mseq
      cseq(1,i)='     '
      cseq(2,i)='     '
      lseq(1,i)=0
      lseq(2,i)=0
      lseq(3,i)=0
   10 continue
      do 11 i=1,matm
      latm(1,i)=0
      latm(2,i)=0
      latm(3,i)=0
      latm(4,i)=0
      latm(5,i)=0
      catm(1,i)='      '
      catm(2,i)='      '
      catm(3,i)='      '
   11 continue
c
c     open the pdb file
c
      open(iunit,file=filnam(1:length),form='formatted',
     + status='old',err=9999)
c
      if(lprint.ge.1) write(lfnout,2000) filnam(1:length)
 2000 format('PDB FILE ',a)
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
      if(model.eq.mdl.and.lprint.ge.5) then
      write(lfnout,2005) model
 2005 format(10x,'MODEL ',i8)
      endif
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
      if(lprint.ge.5) write(lfnout,2001) card
 2001 format(10x,a)
      goto 1
      endif
c
c     SSBOND
c
      if(card(1:6).eq.'SSBOND') then
      nssb=nssb+1
      if(nssb.gt.mssb) call errquit('increase mssb',9999)
      read(card,2007) lssb(1,nssb),lssb(2,nssb)
 2007 format(16x,i4,9x,i4)
      goto 1
      endif
c
c     LINK
c
      if(card(1:6).eq.'LINK  ') then
      nlnk=nlnk+1
      if(nlnk.gt.mlnk) call errquit('increase mlnk',9999)
      read(card,2011)
     + clnk(1,nlnk),llnk(1,nlnk),clnk(2,nlnk),llnk(2,nlnk)
 2011 format(6x,a4,6x,i4,17x,a4,6x,i4)
      goto 1
      endif
c
c     ATOM
c
      if(card(1:4).eq.'ATOM') then
      if(mdl.eq.model) then
      read(card(23:26),1003) i
 1003 format(i4)
      if(i.ne.isgm) then
      isgm=i
      iseq=iseq+1
      nchain=nchain+1
      if(iseq.gt.mseq) then
      write(lfnout,2009)
 2009 format('PDB ERROR : insufficient memory available')
      goto 9999
      endif
      lseq(1,iseq)=isgm
      cseq(1,iseq)(1:3)=card(18:20)
c      if(nchain.eq.1) cseq(1,iseq)(4:5)='_N'
      endif
      iatm=iatm+1
      if(lseq(3,iseq).eq.0) lseq(3,iseq)=iatm
      latm(1,iatm)=i
      catm(1,iatm)(1:4)=card(13:16)
      read(card(31:54),'(3f8.3)') (xatm(j,iatm),j=1,3)
      endif
      goto 1
      endif
c
c     HETATM
c
      if(card(1:6).eq.'HETATM') then
      if(mdl.eq.model) then
      read(card(23:26),1004) i
 1004 format(i4)
      if(i.ne.isgm) then
      isgm=i
      iseq=iseq+1
      nchain=nchain+1
      if(iseq.gt.mseq) then
      write(lfnout,2010)
 2010 format('PDB ERROR : insufficient memory available')
      goto 9999
      endif
      lseq(1,iseq)=isgm
      cseq(1,iseq)(1:3)=card(18:20)
c      if(nchain.eq.1) cseq(1,iseq)(4:5)='_N'
      endif
      iatm=iatm+1
      if(lseq(3,iseq).eq.0) lseq(3,iseq)=iatm
      latm(1,iatm)=i
      catm(1,iatm)(1:4)=card(13:16)
      read(card(31:54),'(3f8.3)') (xatm(j,iatm),j=1,3)
      endif
      if(card(18:20).eq.namslv) then
      iseq=iseq-1
      goto 99
      endif
      goto 1
      endif
c
c     TER
c
c      if(card(1:6).eq.'TER   ') then
c      if(iseq.gt.0) then
c      if(nchain.eq.1) then
c      cseq(1,iseq)(4:5)='_M'
c      else
c      cseq(1,iseq)(4:5)='_C'
c      endif
c      endif
c      nchain=0
c      endif
c
      goto 1
c
  999 continue
c
      if(lprint.ge.1) then
      write(lfnout,2004)
 2004 format('PDB EOF encountered')
      endif
c
    9 continue
c
   99 continue
c
      close(iunit)
c
c     take white space out of segment names
c     -------------------------------------
c
      do 12 i=1,iseq
      k=0
      do 13 j=1,5
      if(cseq(1,i)(j:j).ne.' ') then
      k=k+1
      cseq(1,i)(k:k)=cseq(1,i)(j:j)
      endif
   13 continue
      k=k+1
      do 14 j=k,5
      cseq(1,i)(k:k)=' '
   14 continue
   12 continue
c
      if(iseq.ne.nseq) then
      if(lprint.ge.5) then
      write(lfnout,2008) iseq,nseq
 2008 format(10x,'number of segments read',i8,' while expecting ',i8)
      endif
      nseq=iseq
      endif
c
      if(lprint.ge.5) then
      write(lfnout,2002) nseq,filnam(1:length)
 2002 format(i20,' segments',//,
     + 'PDB CLOSE ',a)
      endif
c
c     since mseq>nseq lseq(3,nseq+1) can hold natm+1
c
      lseq(3,nseq+1)=natm+1
c
      do 23 i=1,natm
      xatm(1,i)=1.0d-01*xatm(1,i)
      xatm(2,i)=1.0d-01*xatm(2,i)
      xatm(3,i)=1.0d-01*xatm(3,i)
   23 continue
c
      pre_pdbseq=.true.
      return
c
 9999 continue
c
      if(lprint.ge.1) then
      write(lfnout,2003) filnam(1:length)
 2003 format('PDB ERROR ',a)
      endif
c
      pre_pdbseq=.false.
      return
      end
